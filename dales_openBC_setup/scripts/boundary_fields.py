# Interpolate fields to DALES domain boundary
# Creates openboundaries.inp.xxx.nc
import numpy as np
import xarray as xr
from datetime import datetime
import pandas as pd 
import glob
import re
from collections import defaultdict
import sys
import dask.array as dask_array
import os
from numba import jit #for compiling the python 
import math
from netCDF4 import Dataset
from scipy.interpolate import interp1d
import gc
    
# single-threaded Numba implementations (safe)
@jit(nopython=True, nogil=True)
def interp1d_numba(src_x, src_y, tgt_x):
    n = len(tgt_x)
    res = np.empty((n,) + src_y.shape[1:], dtype=src_y.dtype)
    for i in range(n):
        x = tgt_x[i]
        if x <= src_x[0]:
            res[i, ...] = src_y[0, ...]
        elif x >= src_x[-1]:
            res[i, ...] = src_y[-1, ...]
        else:
            j = np.searchsorted(src_x, x) - 1
            if j < 0:
                j = 0
            elif j >= len(src_x) - 1:
                j = len(src_x) - 2
            w = (x - src_x[j]) / (src_x[j+1] - src_x[j])
            res[i, ...] = src_y[j, ...]*(1 - w) + src_y[j+1, ...]*w
    return res


@jit(nopython=True, nogil=True)
def interp2d_numba(src_x, src_y, field, tgt_x, tgt_y):
    nx = len(tgt_x)
    ny = len(tgt_y)
    res = np.empty((ny, nx), dtype=field.dtype)

    for iy in range(ny):
        y = tgt_y[iy]
        jy = np.searchsorted(src_y, y) - 1
        if jy < 0:
            jy = 0
            wy = 0.0
        elif jy >= len(src_y) - 1:
            jy = len(src_y) - 2
            wy = 1.0
        else:
            wy = (y - src_y[jy]) / (src_y[jy + 1] - src_y[jy])

        for ix in range(nx):
            x = tgt_x[ix]
            jx = np.searchsorted(src_x, x) - 1
            if jx < 0:
                jx = 0
                wx = 0.0
            elif jx >= len(src_x) - 1:
                jx = len(src_x) - 2
                wx = 1.0
            else:
                wx = (x - src_x[jx]) / (src_x[jx + 1] - src_x[jx])

            f00 = field[jy, jx]
            f10 = field[jy, jx + 1]
            f01 = field[jy + 1, jx]
            f11 = field[jy + 1, jx + 1]

            res[iy, ix] = (
                f00 * (1 - wx) * (1 - wy)
                + f10 * wx * (1 - wy)
                + f01 * (1 - wx) * wy
                + f11 * wx * wy
            )
    return res
    
def interp2d_numba_block(src_x, src_y, block, tgt_x, tgt_y):
    """
    block: can be (time, ny_src, nx_src) or (other_dims..., ny_src, nx_src)
    output: same leading dims, interpolated in last two dims
    """
    leading_shape = block.shape[:-2]
    ny, nx = len(tgt_y), len(tgt_x)
    out_shape = leading_shape + (ny, nx)
    res = np.empty(out_shape, dtype=block.dtype)

    # Flatten leading dims for easier looping
    block_flat = block.reshape(-1, block.shape[-2], block.shape[-1])
    res_flat = res.reshape(block_flat.shape[0], ny, nx)

    for i in range(block_flat.shape[0]):
        res_flat[i] = interp2d_numba(src_x, src_y, block_flat[i], tgt_x, tgt_y)

    return res


def fast_interp(arr, src_x, src_y, tgt_x, tgt_y):
    """Dask-aware 2D interpolation along last two dims, keeps leading dims.
       Operates on arr.data (dask array) and returns a dask array.
    """
    # arr.data is a dask array with shape e.g. (time, ny_src, nx_src) or (other_lead..., ny_src, nx_src)
    darr = arr.data  # dask.array

    # compute new chunks: keep leading chunks, set last-2 dims to target sizes
    leading_chunks = darr.chunks[:-2]   # tuple of tuples
    # create new chunks tuples for last two dims (single-block each is fine)
    new_last_chunks = (len(tgt_y), len(tgt_x))

    # construct chunks argument for dask.map_blocks: must match number of dims
    chunks = leading_chunks + ( (len(tgt_y),), (len(tgt_x),) )

    # Use dask.array.map_blocks with a top-level function that receives numpy ndarray blocks.
    out = dask_array.map_blocks(
    lambda block: interp2d_numba_block(
        np.asarray(src_x),
        np.asarray(src_y),
        np.asarray(block),
        np.asarray(tgt_x),
        np.asarray(tgt_y)
    ),
    darr,
    dtype=arr.dtype,
    chunks=chunks,
    )
    return out
    
def postprocess_time_and_coords(ds, grid, input):
    """
    Convert time + assign standard DALES-like coordinates.
    """

    # =====================================================
    # TIME: datetime → seconds since time0
    # =====================================================
    if "time" in ds.coords:

        time_values = ds.time.values

        dt_seconds = (
            time_values.astype("datetime64[s]")
            - np.datetime64(input["time0"], "s")
        ) / np.timedelta64(1, "s")

        ds = ds.assign_coords(
            time=dt_seconds.astype("float64")
        )

    # =====================================================
    # COORDINATES
    # =====================================================
    coord_map = {
        "xt": grid.xt,
        "xm": grid.xm,
        "yt": grid.yt,
        "ym": grid.ym,
        "zt": grid.zt,
        "zm": grid.zm,
    }

    coord_map = {k: v for k, v in coord_map.items() if k in ds.dims}

    ds = ds.assign_coords(coord_map)

    return ds
    
def postprocess_metadata_and_encoding(ds, input):

    # =====================================================
    # GLOBAL ATTRIBUTES
    # =====================================================
    ds.attrs = {
        "title": f"openboundaries.inp.{input['iexpnr']:03d}.nc",
        "history": f"Created on {datetime.utcnow():%Y-%m-%d %H:%M:%S} UTC",
        "author": input["author"],
        "time0": input["time0"],
    }

    # =====================================================
    # ENCODING
    # =====================================================
    encoding = {
        v: {
            "dtype": "float32",
            "zlib": False,
        }
        for v in ds.data_vars
    }

    return ds, encoding
    
def interpolate_vertical(arr, grid, var_name=None):
    """
    Interpolate a DALES field from its source vertical grid
    onto the target fine-grid vertical coordinates.

    Scalars/u/v: zt
    w: zm

    Extrapolation is allowed because the fine grid may extend
    below the lowest coarse level.
    """

    if 'zt' in arr.dims:
        zdim = 'zt'
        target_z = np.asarray(grid.zt)

    elif 'zm' in arr.dims:
        zdim = 'zm'
        target_z = np.asarray(grid.zm)

    else:
        return arr

    source_z = np.asarray(arr[zdim].values)

    # Check monotonicity
    if not np.all(np.diff(source_z) > 0):
        raise ValueError(
            f"{var_name}: {zdim} is not strictly increasing"
        )

    if not np.all(np.diff(target_z) > 0):
        raise ValueError(
            f"Target grid {zdim} is not strictly increasing"
        )

    # Check where extrapolation is required
    n_below = np.sum(target_z < source_z.min())
    n_above = np.sum(target_z > source_z.max())

    if n_below > 0 or n_above > 0:
        print(
            f"{var_name}: vertical extrapolation: "
            f"{n_below} below, {n_above} above"
        )

    result = arr.interp(
        {zdim: target_z},
        kwargs={"fill_value": "extrapolate"}
    )

    return result.assign_coords(
        {zdim: target_z}
    )

def boundary_fields(input,grid,data):
  data = data.drop(['lat','lon'])
  # West boundary
  uwest   = data['u'].interp(z=grid.zt, y=grid.yt, x=grid.xm[0], assume_sorted=True).rename({'z': 'zt', 'y': 'yt'}).rename('uwest').drop(['x'])
  vwest   = data['v'].interp(z=grid.zt, y=grid.ym, x=grid.xm[0], assume_sorted=True).rename({'z': 'zt', 'y': 'ym'}).rename('vwest').drop(['x'])
  wwest   = data['w'].interp(z=grid.zm, y=grid.yt, x=grid.xm[0], assume_sorted=True).rename({'z': 'zm', 'y': 'yt'}).rename('wwest').drop(['x'])
  thlwest = data['thl'].interp(z=grid.zt, y=grid.yt, x=grid.xm[0], assume_sorted=True).rename({'z': 'zt', 'y': 'yt'}).rename('thlwest').drop(['x'])
  qtwest  = data['qt'].interp(z=grid.zt, y=grid.yt, x=grid.xm[0], assume_sorted=True).rename({'z': 'zt', 'y': 'yt'}).rename('qtwest').drop(['x'])
  qrwest  =  data['qr'].interp(z=grid.zt, y=grid.yt, x=grid.xm[0], assume_sorted=True).rename({'z': 'zt', 'y': 'yt'}).rename('qrwest').drop(['x'])
  e12west = (xr.ones_like(thlwest)*input['e12']).rename('e12west')
  uwest.attrs.clear(); vwest.attrs.clear(); wwest.attrs.clear(); thlwest.attrs.clear(); qtwest.attrs.clear()
  # East boundary
  ueast   = data['u'].interp(z=grid.zt, y=grid.yt, x=grid.xm[-1], assume_sorted=True).rename({'z': 'zt', 'y': 'yt'}).rename('ueast').drop(['x'])
  veast   = data['v'].interp(z=grid.zt, y=grid.ym, x=grid.xm[-1], assume_sorted=True).rename({'z': 'zt', 'y': 'ym'}).rename('veast').drop(['x'])
  weast   = data['w'].interp(z=grid.zm, y=grid.yt, x=grid.xm[-1], assume_sorted=True).rename({'z': 'zm', 'y': 'yt'}).rename('weast').drop(['x'])
  thleast = data['thl'].interp(z=grid.zt, y=grid.yt, x=grid.xm[-1], assume_sorted=True).rename({'z': 'zt', 'y': 'yt'}).rename('thleast').drop(['x'])
  qteast  = data['qt'].interp(z=grid.zt, y=grid.yt, x=grid.xm[-1], assume_sorted=True).rename({'z': 'zt', 'y': 'yt'}).rename('qteast').drop(['x'])
  e12east = (xr.ones_like(thleast)*input['e12']).rename('e12east')
  qreast  = data['qr'].interp(z=grid.zt, y=grid.yt, x=grid.xm[-1], assume_sorted=True).rename({'z': 'zt', 'y': 'yt'}).rename('qreast').drop(['x'])

  ueast.attrs.clear(); veast.attrs.clear(); weast.attrs.clear(); thleast.attrs.clear(); qteast.attrs.clear()
  # South boundary
  usouth   = data['u'].interp(z=grid.zt, y=grid.ym[0], x=grid.xm,assume_sorted=True).rename({'z': 'zt', 'x': 'xm'}).rename('usouth').drop(['y'])
  vsouth   = data['v'].interp(z=grid.zt, y=grid.ym[0], x=grid.xt,assume_sorted=True).rename({'z': 'zt', 'x': 'xt'}).rename('vsouth').drop(['y'])
  wsouth   = data['w'].interp(z=grid.zm, y=grid.ym[0], x=grid.xt,assume_sorted=True).rename({'z': 'zm', 'x': 'xt'}).rename('wsouth').drop(['y'])
  thlsouth = data['thl'].interp(z=grid.zt, y=grid.ym[0], x=grid.xt,assume_sorted=True).rename({'z': 'zt', 'x': 'xt'}).rename('thlsouth').drop(['y'])
  qtsouth  = data['qt'].interp(z=grid.zt, y=grid.ym[0], x=grid.xt,assume_sorted=True).rename({'z': 'zt', 'x': 'xt'}).rename('qtsouth').drop(['y'])
  e12south = (xr.ones_like(thlsouth)*input['e12']).rename('e12south')
  qrsouth  = data['qr'].interp(z=grid.zt, y=grid.ym[0], x=grid.xt,assume_sorted=True).rename({'z': 'zt', 'x': 'xt'}).rename('qrsouth').drop(['y'])
  usouth.attrs.clear(); vsouth.attrs.clear(); wsouth.attrs.clear(); thlsouth.attrs.clear(); qtsouth.attrs.clear()
  # North boundary
  unorth   = data['u'].interp(z=grid.zt, y=grid.ym[-1], x=grid.xm, assume_sorted=True).rename({'z': 'zt', 'x': 'xm'}).rename('unorth').drop(['y'])
  vnorth   = data['v'].interp(z=grid.zt, y=grid.ym[-1], x=grid.xt, assume_sorted=True).rename({'z': 'zt', 'x': 'xt'}).rename('vnorth').drop(['y'])
  wnorth   = data['w'].interp(z=grid.zm, y=grid.ym[-1], x=grid.xt, assume_sorted=True).rename({'z': 'zm', 'x': 'xt'}).rename('wnorth').drop(['y'])
  thlnorth = data['thl'].interp(z=grid.zt, y=grid.ym[-1], x=grid.xt, assume_sorted=True).rename({'z': 'zt', 'x': 'xt'}).rename('thlnorth').drop(['y'])
  qtnorth  = data['qt'].interp(z=grid.zt, y=grid.ym[-1], x=grid.xt, assume_sorted=True).rename({'z': 'zt', 'x': 'xt'}).rename('qtnorth').drop(['y'])
  e12north = (xr.ones_like(thlnorth)*input['e12']).rename('e12north')
  qrnorth  = data['qr'].interp(z=grid.zt, y=grid.ym[-1], x=grid.xt, assume_sorted=True).rename({'z': 'zt', 'x': 'xt'}).rename('qrnorth').drop(['y'])

  unorth.attrs.clear(); vnorth.attrs.clear(); wnorth.attrs.clear(); thlnorth.attrs.clear(); qtnorth.attrs.clear()
  # Top boundary
  utop   = data['u'].interp(z=grid.zm[-1], y=grid.yt, x=grid.xm, assume_sorted=True).rename({'y': 'yt', 'x': 'xm'}).rename('utop').drop(['z'])
  vtop   = data['v'].interp(z=grid.zm[-1], y=grid.ym, x=grid.xt, assume_sorted=True).rename({'y': 'ym', 'x': 'xt'}).rename('vtop').drop(['z'])
  wtop   = data['w'].interp(z=grid.zm[-1], y=grid.yt, x=grid.xt, assume_sorted=True).rename({'y': 'yt', 'x': 'xt'}).rename('wtop').drop(['z'])
  thltop = data['thl'].interp(z=grid.zm[-1], y=grid.yt, x=grid.xt, assume_sorted=True).rename({'y': 'yt', 'x': 'xt'}).rename('thltop').drop(['z'])
  qttop  = data['qt'].interp(z=grid.zm[-1], y=grid.yt, x=grid.xt, assume_sorted=True).rename({'y': 'yt', 'x': 'xt'}).rename('qttop').drop(['z'])
  e12top = (xr.ones_like(thltop)*input['e12']).rename('e12top')
  qrtop  = data['qr'].interp(z=grid.zm[-1], y=grid.yt, x=grid.xt, assume_sorted=True).rename({'y': 'yt', 'x': 'xt'}).rename('qrtop').drop(['z'])
  utop.attrs.clear(); vtop.attrs.clear(); wtop.attrs.clear(); thltop.attrs.clear(); qttop.attrs.clear()
  # Add fields to dataset
  openboundaries = xr.merge([uwest, vwest, wwest, thlwest, qtwest, e12west, qrwest,
                              ueast, veast, weast, thleast, qteast, e12east, qreast,
                              usouth,vsouth,wsouth,thlsouth,qtsouth,e12south, qrsouth,
                              unorth,vnorth,wnorth,thlnorth,qtnorth,e12north, qrnorth,
                              utop,  vtop,  wtop,  thltop,  qttop,  e12top, qrtop],
                              combine_attrs='drop')
  # Adjust time variable to seconds since initial field
  ts = openboundaries['time'].values.astype('datetime64[s]')
  dts = (ts-np.datetime64(input['time0'],'s'))/np.timedelta64(1, 's')
  openboundaries = openboundaries.assign_coords({'time':('time', dts)})
  openboundaries['time'].attrs.clear()
  # Add variable attributes
  openboundaries['time'] = openboundaries['time'].assign_attrs({'longname': 'Time', 'units': f"seconds since {input['time0']}"})
  openboundaries['xt'] = openboundaries['xt'].assign_attrs({'longname': 'West-East displacement of cell centers','units': 'm'})
  openboundaries['xm'] = openboundaries['xm'].assign_attrs({'longname': 'West-East displacement of cell edges','units': 'm'})
  openboundaries['yt'] = openboundaries['yt'].assign_attrs({'longname': 'South-North displacement of cell centers','units': 'm'})
  openboundaries['ym'] = openboundaries['ym'].assign_attrs({'longname': 'South-North displacement of cell edges','units': 'm'})
  openboundaries['zt'] = openboundaries['zt'].assign_attrs({'longname': 'Vertical displacement of cell centers','units': 'm'})
  openboundaries['zm'] = openboundaries['zm'].assign_attrs({'longname': 'Vertical displacement of cell edges','units': 'm'})
  variables = ['u','v','w','thl','qt','e12']
  units     = ['m/s','m/s','m/s','K','kg/kg','m/s','kg/kg']
  long_names= ['West-East velocity at ',
                'South-North velocity at ',
                'Vertical velocity at ',
                'Liquid water potential temperature at ',
                'Total water specific humidity at ',
                'Square root of turbulent kinetic energy at ',
                'Rain water mixing ratio at ']
  for ivar in range(len(variables)):
    var = variables[ivar]
    unit = units[ivar]
    long_name = long_names[ivar]
    for boundary in ['West','East','South','North','top']:
      openboundaries[var+boundary.lower()] = openboundaries[var+boundary.lower()]\
      .assign_attrs({'longname': long_name+boundary+' boundary', 'units': unit})
  # Add global attributes
  openboundaries = openboundaries.assign_attrs({'title': f"openboundaries.inp.{input['iexpnr']:03d}.nc",
                                        'history': f"Created on {datetime.utcnow().strftime('%Y-%m-%d %H:%M:%S')} UTC",
                                        'author': input['author'],
                                        'time0': input['time0']})
  openboundaries.to_netcdf(path=input['outpath']+openboundaries.attrs['title'], mode='w', format="NETCDF4")
  return openboundaries
  
def boundary_fields_fine(input,grid):
  ix_west = int(input['x_offset']/input['dx_coarse'])
  ix_east = int(ix_west+grid.xsize/input['dx_coarse'])
  iy_south = int(input['y_offset']/input['dy_coarse'])
  iy_north = int(iy_south+grid.ysize/input['dy_coarse'])
  
  ix_west=ix_west+2
  ix_east=ix_east+2
  iy_south=iy_south+2
  iy_north=iy_north+2
  print(ix_west,ix_east,iy_south,iy_north)
  
  # Get initial boundary fields from initial fields
  if(input['time0']==input['start']):
    with xr.open_mfdataset(f"{input['inpath']}initfields.inp.*.nc") as ds:
      
      exclude = ['lat', 'lon', 'time', 'transform']
      all_vars = [v for v in ds.variables if v not in ds.coords and v not in exclude]
      
      boundary_fields_west0 = {}
      boundary_fields_east0 = {}
      boundary_fields_north0 = {}
      boundary_fields_south0 = {}
      boundary_fields_top0 = {}
      
      for var in all_vars:
         
        # West boundary
        if  var == 'u0':
            boundary_fields_west0['uwest'] = ds[var].isel(xm=ix_west,drop=True).interp(yt=grid.yt+input['y_offset'], zt=grid.zt, kwargs={"fill_value": "extrapolate"}).rename('uwest').expand_dims({'time':[pd.Timestamp(input['time0'])]},axis=0).assign_coords(yt=grid.yt,zt=grid.zt)
            boundary_fields_east0['ueast'] = ds[var].isel(xm=ix_east,drop=True).interp(yt=grid.yt+input['y_offset'], zt=grid.zt, kwargs={"fill_value": "extrapolate"}).rename('ueast').expand_dims({'time':[pd.Timestamp(input['time0'])]},axis=0).assign_coords(yt=grid.yt,zt=grid.zt)
            boundary_fields_south0['usouth'] = ds[var].isel(yt=iy_south,drop=True).interp(xm=grid.xm+input['x_offset'], zt=grid.zt, kwargs={"fill_value": "extrapolate"}).rename('usouth').expand_dims({'time':[pd.Timestamp(input['time0'])]},axis=0).assign_coords(xm=grid.xm,zt=grid.zt)
            boundary_fields_north0['unorth'] = ds[var].isel(yt=iy_north,drop=True).interp(xm=grid.xm+input['x_offset'], zt=grid.zt, kwargs={"fill_value": "extrapolate"}).rename('unorth').expand_dims({'time':[pd.Timestamp(input['time0'])]},axis=0).assign_coords(xm=grid.xm,zt=grid.zt)
            boundary_fields_top0['utop'] = ds[var].interp(zt=grid.zt[-1], xm=grid.xm+input['x_offset'],yt=grid.yt+input['y_offset'], kwargs={"fill_value": "extrapolate"}).rename('utop').expand_dims({'time':[pd.Timestamp(input['time0'])]},axis=0).assign_coords(xm=grid.xm,yt=grid.yt)
        elif var == 'v0':
            boundary_fields_west0['vwest'] = ds[var].isel(xt=ix_west,drop=True).interp(ym=grid.ym+input['y_offset'], zt=grid.zt, kwargs={"fill_value": "extrapolate"}).rename('vwest').expand_dims({'time':[pd.Timestamp(input['time0'])]},axis=0).assign_coords(ym=grid.ym,zt=grid.zt)
            boundary_fields_east0['veast'] = ds[var].isel(xt=ix_east,drop=True).interp(ym=grid.ym+input['y_offset'], zt=grid.zt, kwargs={"fill_value": "extrapolate"}).rename('veast').expand_dims({'time':[pd.Timestamp(input['time0'])]},axis=0).assign_coords(ym=grid.ym,zt=grid.zt)
            boundary_fields_south0['vsouth'] = ds[var].isel(ym=iy_south,drop=True).interp(xt=grid.xt+input['x_offset'], zt=grid.zt, kwargs={"fill_value": "extrapolate"}).rename('vsouth').expand_dims({'time':[pd.Timestamp(input['time0'])]},axis=0).assign_coords(xt=grid.xt,zt=grid.zt)
            boundary_fields_north0['vnorth'] = ds[var].isel(ym=iy_north,drop=True).interp(xt=grid.xt+input['x_offset'], zt=grid.zt, kwargs={"fill_value": "extrapolate"}).rename('vnorth').expand_dims({'time':[pd.Timestamp(input['time0'])]},axis=0).assign_coords(xt=grid.xt,zt=grid.zt)
            boundary_fields_top0['vtop'] = ds[var].interp(zt=grid.zt[-1],xt=grid.xt+input['x_offset'],ym=grid.ym+input['y_offset'], kwargs={"fill_value": "extrapolate"}).rename('vtop').expand_dims({'time':[pd.Timestamp(input['time0'])]},axis=0).assign_coords(xt=grid.xt,ym=grid.ym)
        elif var == 'w0':
            boundary_fields_west0['wwest']  = ds[var].isel(xt=ix_west,drop=True).interp(yt=grid.yt+input['y_offset'],zm=grid.zm,kwargs={"fill_value": "extrapolate"}).rename('wwest').expand_dims({'time':[pd.Timestamp(input['time0'])]},axis=0).assign_coords(yt=grid.yt,zm=grid.zm)
            boundary_fields_east0['weast'] = ds[var].isel(xt=ix_east,drop=True).interp(yt=grid.yt+input['y_offset'],zm=grid.zm,kwargs={"fill_value": "extrapolate"}).rename('weast').expand_dims({'time':[pd.Timestamp(input['time0'])]},axis=0).assign_coords(yt=grid.yt,zm=grid.zm)
            boundary_fields_south0['wsouth'] = ds[var].isel(yt=iy_south,drop=True).interp(xt=grid.xt+input['x_offset'],zm=grid.zm,kwargs={"fill_value": "extrapolate"}).rename('wsouth').expand_dims({'time':[pd.Timestamp(input['time0'])]},axis=0).assign_coords(xt=grid.xt,zm=grid.zm)
            boundary_fields_north0['wnorth'] = ds[var].isel(yt=iy_north,drop=True).interp(xt=grid.xt+input['x_offset'],zm=grid.zm,kwargs={"fill_value": "extrapolate"}).rename('wnorth').expand_dims({'time':[pd.Timestamp(input['time0'])]},axis=0).assign_coords(xt=grid.xt,zm=grid.zm)
            boundary_fields_top0['wtop'] = ds[var].interp(zm=grid.zm[-1], xt=grid.xt+input['x_offset'],yt=grid.yt+input['y_offset'], kwargs={"fill_value": "extrapolate"}).rename('wtop').expand_dims({'time':[pd.Timestamp(input['time0'])]},axis=0).assign_coords(xt=grid.xt,yt=grid.yt)
        
        else:
            base = var[:-1] if var.endswith("0") else var
            boundary_fields_west0[f'{base}west']=ds[var].isel(xt=ix_west,drop=True).interp(yt=grid.yt+input['y_offset'], zt=grid.zt, kwargs={"fill_value": "extrapolate"}).rename(f'{base}west').expand_dims({'time':[pd.Timestamp(input['time0'])]},axis=0).assign_coords(yt=grid.yt,zt=grid.zt)
            boundary_fields_east0[f'{base}east']=ds[var].isel(xt=ix_east,drop=True).interp(yt=grid.yt+input['y_offset'], zt=grid.zt, kwargs={"fill_value": "extrapolate"}).rename(f'{base}east').expand_dims({'time':[pd.Timestamp(input['time0'])]},axis=0).assign_coords(yt=grid.yt,zt=grid.zt)
            boundary_fields_south0[f'{base}south']=ds[var].isel(yt=iy_south,drop=True).interp(xt=grid.xt+input['x_offset'], zt=grid.zt, kwargs={"fill_value": "extrapolate"}).rename(f'{base}south').expand_dims({'time':[pd.Timestamp(input['time0'])]},axis=0).assign_coords(xt=grid.xt,zt=grid.zt)
            boundary_fields_north0[f'{base}north']=ds[var].isel(yt=iy_north,drop=True).interp(xt=grid.xt+input['x_offset'], zt=grid.zt, kwargs={"fill_value": "extrapolate"}).rename(f'{base}north').expand_dims({'time':[pd.Timestamp(input['time0'])]},axis=0).assign_coords(xt=grid.xt,zt=grid.zt)
            boundary_fields_top0[f'{base}top']=ds[var].interp(zt=grid.zt[-1], xt=grid.xt+input['x_offset'],yt=grid.yt+input['y_offset'], kwargs={"fill_value": "extrapolate"}).rename(f'{base}top').expand_dims({'time':[pd.Timestamp(input['time0'])]},axis=0).assign_coords(xt=grid.xt,yt=grid.yt)
        
  for name, da in boundary_fields_west0.items():
    print(f"{name}: {da.dims}, shape={da.shape}")
  for name, da in boundary_fields_east0.items():
    print(f"{name}: {da.dims}, shape={da.shape}")
  for name, da in boundary_fields_south0.items():
    print(f"{name}: {da.dims}, shape={da.shape}")
  for name, da in boundary_fields_north0.items():
    print(f"{name}: {da.dims}, shape={da.shape}")
  for name, da in boundary_fields_top0.items():
    print(f"{name}: {da.dims}, shape={da.shape}")
    
    
    
  # Get later time steps from corresponding coarse simulation output
  # Exclude list
  exclude = [
    'qlyz', 'thvyz', 'buoyyz', 'qlxz', 'thvxz', 'buoyxz',
    'buoyyz', 'qlxy', 'thvxy', 'buoyxy', 'co2ags', 'co2veg'
  ]

  # Boundary definitions
  boundaries = {
    'west':  {'ix': ix_west,  'file': 'crossyz', 'suffix': 'west'},
    'east':  {'ix': ix_east,  'file': 'crossyz', 'suffix': 'east'},
    'south': {'ix': iy_south, 'file': 'crossxz', 'suffix': 'south'},
    'north': {'ix': iy_north, 'file': 'crossxz', 'suffix': 'north'},
    'top':   {'ix': None,     'file': 'crossxy', 'suffix': 'top'},
  }
  
  # Map cross-section suffixes to boundaries
  section_suffix = {
    'crossyz': 'yz',   # west, east
    'crossxz': 'xz',   # north, south
    'crossxy': 'xy'    # top
  }

  # Dictionaries for storing fields
  all_boundary_fields = {}

  boundary_dict_map = {
    'west': boundary_fields_west0,
    'east': boundary_fields_east0,
    'south': boundary_fields_south0,
    'north': boundary_fields_north0,
    'top': boundary_fields_top0
  }

  # Load the *coarse* openboundaries file to grab variables list
  coarse_fn = f"{input['inpath']}openboundaries.inp.{input['iexpnr']:03d}.nc"
  openb_coarse_file = xr.open_dataset(coarse_fn)
  openb_coarse_file = openb_coarse_file.sel(time=slice(input['start'], input['end']))
  variables = list(openb_coarse_file.data_vars)

  # Filter variables
  variables = [v for v in variables if v not in exclude]

  print("Variables to process:", variables)
  
  # Target boundaries file:
  outfn = os.path.join(input['outpath'], f"openboundaries.inp.{input['iexpnr']:03d}.nc")

  # Create empty file with coordinates only (write once)
  coord_ds = xr.Dataset(
    {
        "xt": (("xt",), grid.xt),
        "xm": (("xm",), grid.xm),
        "yt": (("yt",), grid.yt),
        "ym": (("ym",), grid.ym),
        "zt": (("zt",), grid.zt),
        "zm": (("zm",), grid.zm),
    }
  )

  coord_ds.to_netcdf(
    outfn,
    mode="w",
    format="NETCDF4",
  )
  print(f"Created output file with coords: {outfn}")

  # 2) Add coordinate attributes (open in append mode and set attrs)
  ds_tmp = xr.open_dataset(outfn, mode="a")
  ds_tmp["xt"].attrs = {"long_name": "West-East displacement of cell centers", "units": "m"}
  ds_tmp["xm"].attrs = {"long_name": "West-East displacement of cell edges",   "units": "m"}
  ds_tmp["yt"].attrs = {"long_name": "South-North displacement of cell centers","units": "m"}
  ds_tmp["ym"].attrs = {"long_name": "South-North displacement of cell edges",  "units": "m"}
  ds_tmp["zt"].attrs = {"long_name": "Vertical displacement of cell centers",   "units": "m"}
  ds_tmp["zm"].attrs = {"long_name": "Vertical displacement of cell edges",     "units": "m"}
  ds_tmp.close()
  
  #make sure that there is no outfn in the target folder:
  if os.path.exists(outfn):
    os.remove(outfn)
      
  # Process each boundary
  for bname, binfo in boundaries.items():
    print(f"\nProcessing boundary: {bname}")
    boundary_fields = {}

    # Determine cross-section folder (e.g. crossyz, crossxz, crossxy)
    level = f"{binfo['ix']:04d}" if bname != 'top' else f"{grid.kmax:04d}"
    #folder = os.path.join(input['inpath'], 'crossections', binfo['file'], level)
    folder = os.path.join(input['inpath'], 'crossections', binfo['file'], level)
    if not os.path.isdir(folder):
        print(f"⚠️  Folder not found: {folder}")
        continue

    # List .nc files available in that folder
    nc_files = sorted(glob.glob(os.path.join(folder, "*.nc")))
    available_vars = [os.path.basename(f).split('.')[0] for f in nc_files]
    print(f"  Found {len(available_vars)} files in {binfo['file']}/{level}")

    # Determine which cross-section suffix to look for (yz/xz/xy)
    sect_sfx = section_suffix.get(binfo['file'], "")

    # Build list of coarse variable *base names* (remove boundary suffix)
    coarse_base_vars = []
    for v in variables:
        for suffix in ['west', 'east', 'north', 'south', 'top']:
            if v.endswith(suffix):
                coarse_base_vars.append(v[:-len(suffix)])
                break
        else:
            coarse_base_vars.append(v)
            
    coarse_base_vars = list(dict.fromkeys(coarse_base_vars))

    # Select variables that exist as either base or base+cross-section suffix
    selected_vars = []

    for base in coarse_base_vars:
        if base in exclude:
            continue

        candidates = []

        # ONLY add "0" for e12-style variables (not for others)
        if base.startswith("e12"):  
            candidates.append(base + "0")
        else:
            candidates.append(base)

        match = next((c for c in candidates if c in available_vars), None)

        if match is not None:
            selected_vars.append(match)

    print(f"  Variables to read from {binfo['file']} ({bname}): {selected_vars}")
    
    # Now iterate over selected variables
    for var_name in selected_vars:
        
        var_file = os.path.join(folder, f"{var_name}.{level}.001.nc")
    
        if not os.path.exists(var_file):
            print(f"⚠️  File not found: {var_file}")
            continue

        # Lazy open with chunking
        ds = xr.open_dataset(var_file, chunks={"time": input['tchunk']}, engine="netcdf4")
        
        if not np.issubdtype(ds.time.dtype, np.datetime64):
            ds = ds.assign_coords(
                time=(
                    np.datetime64(input["time0"])
                    + ds.time.values.astype("timedelta64[s]")
                )
            )
            
        print("time dtype:", ds.time.dtype)
        print("first time:", ds.time.values[0])
        print("last time:", ds.time.values[-1])

        print("start:", input['start'])
        print("end:", input['end'])
        
        ds = ds.sel(time=slice(input['start'], input['end']))

        arr = ds[var_name]

        # Determine proper horizontal coordinates for interpolation
        if bname == 'top':
              if var_name.startswith('u'):
                horiz_dims = ('yt', 'xm')
              elif var_name.startswith('v'):
                horiz_dims = ('ym', 'xt')
              else:
                horiz_dims = ('yt', 'xt')
        elif bname in ['west', 'east']:
            horiz_dims = ('ym',) if var_name.startswith('v') else ('yt',)
        elif bname in ['south', 'north']:
            horiz_dims = ('xm',) if var_name.startswith('u') else ('xt',)
        else:
            raise ValueError(f"Unknown boundary: {bname}")

        print('Final var_name:', var_name)

        # Determine vertical coordinate (if any)
        vert_coord = {}
        if 'zm' in arr.dims:
            vert_coord['zm'] = grid.zm
        elif 'zt' in arr.dims:
            vert_coord['zt'] = grid.zt

        # Coordinates for interpolation
        coords_to_assign = {c: getattr(grid, c) for c in horiz_dims if hasattr(grid, c)}
        coords_to_assign.update(vert_coord)

        # Interpolate 2d
        if len(horiz_dims) == 2 and all(d in arr.dims for d in horiz_dims):
            #xname, yname = horiz_dims
            yname, xname = horiz_dims
            src_x = arr[xname].values
            src_y = arr[yname].values

            # Ensure both src axes are not empty
            if len(src_x) < 2 or len(src_y) < 2:
                print(f"⚠️ Skipping 2D interpolation for {var_name}, insufficient points: "
                      f"{xname}={len(src_x)}, {yname}={len(src_y)}")
                arr_interp = arr
            else:
                # Adjust fine grid coordinates for physical offset
                phys_tgt_x = getattr(grid, xname) + input['x_offset']
                phys_tgt_y = getattr(grid, yname) + input['y_offset']
                
                # Target coordinates  (no offset start from 0)
                tgt_x = getattr(grid, xname)
                tgt_y = getattr(grid, yname)
                
                arr = arr.chunk({
                    yname: -1,
                    xname: -1
                })
                
                print("dims:", arr.dims)
                print("shape:", arr.shape)
                print("chunks:", arr.chunks)

                arr_interp_data = fast_interp(arr, src_x, src_y, phys_tgt_x, phys_tgt_y)

                # Determine leading dims (e.g., time, zt, zm)
                leading_dims = [d for d in arr.dims if d not in horiz_dims]
                all_dims = leading_dims + [yname, xname]

                # Build coords
                coords = {d: arr.coords[d] for d in leading_dims if d in arr.coords}
                coords.update({yname: tgt_y, xname: tgt_x})

                arr_interp = xr.DataArray(
                    arr_interp_data,
                    dims=all_dims,
                    coords=coords,
                    attrs=arr.attrs,
                    name=arr.name,
                )

        else:
            # =========================================================
            # 3-D boundary field:
            #     (time, zt/zm, horizontal)
            #
            # Step 1: vertical interpolation
            #         coarse zt/zm -> fine grid.zt/grid.zm
            #
            # Step 2: horizontal interpolation
            #         coarse horizontal -> fine horizontal
            # =========================================================

            # ---------------------------------------------------------
            # 1. Vertical interpolation
            # ---------------------------------------------------------
            #
            # u, v, scalars: zt
            # w:            zm
            #
            # This also handles bottom extrapolation because the fine
            # grid may extend below the lowest coarse vertical level.
            #
            arr = interpolate_vertical(
                arr,
                grid,
                var_name=var_name
            )

            print(
                f"{var_name}: after vertical interpolation: "
                f"dims={arr.dims}, shape={arr.shape}"
            )

            # ---------------------------------------------------------
            # 2. Determine horizontal coordinate
            # ---------------------------------------------------------
            dim_h = [d for d in horiz_dims if d in arr.dims][0]

            src_x = np.asarray(arr[dim_h].values)

            if dim_h in ['xt', 'xm']:
                offset = input['x_offset']
            else:
                offset = input['y_offset']

            phys_tgt_x = np.asarray(
                getattr(grid, dim_h) + offset
            )

            tgt_x = np.asarray(
                getattr(grid, dim_h)
            )

            # ---------------------------------------------------------
            # 3. Reorder array so horizontal dimension is last
            #
            # Result:
            #     (time, z, horizontal)
            # ---------------------------------------------------------
            other_dims = [
                d for d in arr.dims
                if d != dim_h
            ]        

            arr_reordered = arr.transpose(
                *other_dims,
                dim_h
            )

            # ---------------------------------------------------------
            # 4. Horizontal interpolation
            # ---------------------------------------------------------
            f = interp1d(
                src_x,
                arr_reordered.data,
                axis=-1,
                bounds_error=False,
                fill_value="extrapolate"
            )

            arr_interp_data = f(phys_tgt_x)

            # ---------------------------------------------------------
            # 5. Build output DataArray
            # ---------------------------------------------------------
            all_dims = other_dims + [dim_h]

            coords = {}

            for d in all_dims:

                if d == dim_h:
                    # Fine horizontal coordinate
                    coords[d] = tgt_x

                elif d == 'zt':
                    # Fine scalar/u/v vertical coordinate
                    coords[d] = grid.zt

                elif d == 'zm':
                    # Fine w vertical coordinate
                    coords[d] = grid.zm

                elif d == 'time':
                    # Keep original time coordinate
                    coords[d] = arr['time']

                elif d in arr.coords:
                    coords[d] = arr.coords[d]

            arr_interp = xr.DataArray(
                arr_interp_data,
                dims=all_dims,
                coords=coords,
                attrs=arr.attrs,
                name=arr.name
            )
            
        # Rename variable
        if var_name.endswith(('xy', 'yz', 'xz')):
              base = var_name[:-2]
              if base.startswith('e12'):
                 base = base[:-1]
              new_name = f"{base}{bname}"
        else:
              if var_name.startswith('e12'):
                 var_name = var_name[:-1]
              new_name = f"{var_name}{bname}"
              
        boundary_dict_0 = boundary_dict_map[bname]
        if new_name in boundary_dict_0:
            arr_interp = xr.concat([boundary_dict_0[new_name], arr_interp], dim='time')

        # Materialize computed data to NumPy (reduce Dask overhead) 
        #    This ensures a clean write per variable and avoids keeping Dask graphs.
        arr_to_write = arr_interp #.compute()

        # Build single-variable dataset and copy attributes from coarse file if present
        #ds_single = xr.Dataset({new_name: arr_to_write}).reset_coords(drop=True)
        #f new_name in openb_coarse_file.data_vars:
            #ds_single[new_name].attrs = openb_coarse_file[new_name].attrs.copy()
        #else:
            #ds_single[new_name].attrs = {"long_name": new_name, "units": ""}

        # Append variable to final file
        #    Use mode="a" - this will add the variable and its dimensions (time, zm/zt, xt/yt/xm/ym)
        #ds_single.to_netcdf(outfn, mode="a")
        #boundary_output[new_name] = arr_interp
        #print(f"Appended variable {new_name} to {outfn}")
        
        
        # ============================================================
        # WRITE VARIABLE IMMEDIATELY (avoids storing all variables)
        # ============================================================
        ds_single = xr.Dataset({new_name: arr_interp}).reset_coords(drop=True)
        ds_single = postprocess_time_and_coords(ds_single, grid, input)

        if new_name in openb_coarse_file.data_vars:
            ds_single[new_name].attrs = openb_coarse_file[new_name].attrs.copy()
        else:
            ds_single[new_name].attrs = {
                "long_name": new_name,
                "units": ""
            }

        if not os.path.exists(outfn):
            ds_single, encoding = postprocess_metadata_and_encoding(
                ds_single, input
            )

            ds_single.to_netcdf(
                outfn,
                mode="w",
                encoding=encoding
            )
        else:
            ds_single.to_netcdf(
                outfn,
                mode="a",
                encoding={
                    new_name: {
                        "dtype": "float32",
                        "zlib": False,
                    }
                }
            )
        
        gc.collect()

        print(f"Appended variable {new_name} to {outfn}")

        ds_single.close()
        ds.close()

        for var in ['arr', 'arr_interp', 'ds_single', 'f', 'arr_interp_data']:
            if var in locals():
                del locals()[var]
        

    print(f"✅ Saved open boundaries to {outfn}")

  return outfn