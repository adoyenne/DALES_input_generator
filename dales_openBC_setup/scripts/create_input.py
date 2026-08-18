#%% Load modules
import json
from GridDales import GridDales
from prep_harmonie import prep_harmonie
import prep_harmonie_WINS50
import prep_harmonie_rd
from initial_fields import initial_fields, initial_fields_fine
from boundary_fields import boundary_fields, boundary_fields_fine
from profiles import profiles
from surface_temperature import surface_temperature, surface_temperature_fine
from synthetic_turbulence import synthetic_turbulence
from gaussian_filter import gaussian_filter

#from lsm.create_dales_input import create_lsm_input    # older LSM input format
from lsm.spatial_transforms import proj4_rd, proj4_hm

# new LSM input compatible with DEPAC
from land_surface.create_dales_input import create_lsm_input as create_lsm_input_tno
#from land_surface.spatial_transforms import proj4_rd, proj4_hm  # identical to the one in lsm
import f90nml
from datetime import datetime
import sys
import os
import shutil
import xarray as xr

def setup_DALES(dalesdir, spatial_data_dir):
  print('Creating symlinks')
  for f in ('rrtmg_lw.nc', 'rrtmg_sw.nc', 'van_genuchten_parameters.nc'):
    try:
      os.symlink(os.path.join(spatial_data_dir, f), os.path.join(dalesdir, f))
    except Exception as e:
      print(repr(e))

def patch_namelist(outdir, inp,
                   ix_west_nested, ix_east_nested,
                   iy_south_nested, iy_north_nested,
                   domain_type):

    # --- Read base namelist ---
    namelist = f90nml.read(f'setting_files/namoptions_{domain_type}.001')

    # --- Grid / domain ---
    grid = inp['grid']
    namelist['domain']['itot']  = grid['itot']
    namelist['domain']['jtot']  = grid['jtot']
    namelist['domain']['kmax']  = grid['kmax']
    namelist['domain']['xsize'] = grid['xsize']
    namelist['domain']['ysize'] = grid['ysize']

    namelist['domain']['xlat'] = round(inp['xlat'], 3)
    namelist['domain']['xlon'] = round(inp['xlon'], 3)

    # --- Time handling ---
    start_dt = datetime.fromisoformat(inp['start'])
    namelist['domain']['xtime'] = start_dt.hour
    namelist['domain']['xday']  = start_dt.timetuple().tm_yday
    namelist['domain']['xyear'] = start_dt.year

    # --- Open boundary settings ---
    namelist['openbc']['dxint']  = 4 * (grid['xsize'] / grid['itot'])
    namelist['openbc']['dyint']  = 4 * (grid['ysize'] / grid['jtot'])
    namelist['openbc']['dxturb'] = grid['xsize']
    namelist['openbc']['dyturb'] = grid['ysize']
    
    if 'namdatetime' in namelist:
        namelist['namdatetime']['startyear']  = start_dt.year
        namelist['namdatetime']['startmonth'] = start_dt.month
        namelist['namdatetime']['startday']   = start_dt.day

        # only if hour exists in namelist
        if 'starthour' in namelist['namdatetime']:
            namelist['namdatetime']['starthour'] = start_dt.hour
        
    # --- Domain-specific logic ---
    if domain_type == 'coarse':
        ps   = inp['ps']
        thls = inp['thls']

        namelist['namsurface']['ps']   = ps
        namelist['namsurface']['thls'] = thls
        namelist['openbc']['lsynturb'] = True

        if 'namchem' in namelist:
            namelist['namchem']['t_ref'] = thls
            namelist['namchem']['p_ref'] = ps
            if 'qts' in inp:
                namelist['namchem']['q_ref'] = inp['qts']

    elif domain_type == 'fine':
        profile_path = os.path.join(
            inp['inpath'],
            f"profiles.{inp['iexpnr']:03d}.nc"
        )

        ds = xr.open_dataset(profile_path)

        # Extract surface values (time=0, lowest level)
        ps   = ds['presh'].isel(time=0, zt=0).item()
        thls = ds['thl'].isel(time=0, zt=0).item()
        qts  = ds['qt'].isel(time=0, zt=0).item()

        ds.close()

        namelist['namsurface']['ps']   = ps
        namelist['namsurface']['thls'] = thls
        namelist['openbc']['lsynturb'] = False

        if 'namchem' in namelist:
            namelist['namchem']['t_ref'] = thls
            namelist['namchem']['p_ref'] = ps
            namelist['namchem']['q_ref'] = qts

    else:
        raise ValueError(f"Unknown domain_type: {domain_type}")

    # --- Parallel settings ---
    namelist['run']['nprocx'] = inp['nprocx']
    namelist['run']['nprocy'] = inp['nprocy']

    # --- Cross sections (only if provided) ---
    if iy_south_nested not in [None, ''] and iy_north_nested not in [None, '']:
        namelist['namcrosssection']['crossplane'] = [int(iy_south_nested), int(iy_north_nested)]
    else:
        namelist['namcrosssection'].pop('crossplane', None)

    if ix_west_nested not in [None, ''] and ix_east_nested not in [None, '']:
        namelist['namcrosssection']['crossortho'] = [int(ix_west_nested), int(ix_east_nested)]
    else:
        namelist['namcrosssection'].pop('crossortho', None)

    # --- Write output ---
    os.makedirs(outdir, exist_ok=True)
    namelist.write(f'{outdir}/namoptions.001', force=True)
    

    shutil.copy(f"setting_files/tracerdata.inp", f"{outdir}/tracerdata.inp")
    if 'namchem' in namelist:
        shutil.copy(f"setting_files/chem.inp.{inp['iexpnr']:03d}", f"{outdir}/chem.inp.{inp['iexpnr']:03d}")

#%% Read input file
with open(sys.argv[1]) as f: input = json.load(f)
#%% Create input for outer simulation

if 'coarse' in input:
  domain_type='coarse'
  input_coarse = input['coarse']
  os.makedirs(input_coarse['outpath'], exist_ok=True)
  input_coarse['author'] = input['author']
  #%% Create DALES grid
  grid = GridDales(input_coarse['grid'])
  #%% Transfor input data to rectilinear grid and to prognostic variables of DALES
  if(input_coarse['source'].lower() == 'harmonie'):
    data,transform = prep_harmonie(input_coarse,grid)
  elif(input_coarse['source'].lower() == 'harmonie_rd'):
    data,transform = prep_harmonie_rd.prep_harmonie(input_coarse,grid)
  elif(input_coarse['source'].lower() == 'harmonie_wins50'):
    data,transform = prep_harmonie_WINS50.prep_harmonie(input_coarse,grid)
  elif(input_coarse['source'].lower() == 'none'): # useful to test LSM alone
    pass
  else:
    print('unvalid source type')
    exit()
  #%% Apply spatial horizontal Gaussian filter to data
  if('filter' in input_coarse and input_coarse['source'].lower() != 'none'):
    data = gaussian_filter(data,input_coarse)
    data.to_netcdf(f"{input_coarse['outpath']}/input_data.nc") # experiment: save the processed data

  if 'LSM' in input_coarse:
  
    x_sw, y_sw = proj4_rd(input_coarse['lon_sw'], input_coarse['lat_sw'], inverse=False)
    print(f'LSM {x_sw}, {y_sw}.')
    dx = input_coarse['grid']['xsize'] / input_coarse['grid']['itot']
    dy = input_coarse['grid']['ysize'] / input_coarse['grid']['jtot']
    start_date = datetime.fromisoformat(input_coarse['start'])

    lsm_kind = 'old'
    if 'lsm_kind' in input_coarse['LSM']:
      lsm_kind = input_coarse['LSM']['lsm_kind']
    lwrite_ags = input_coarse['LSM'].get('lwrite_ags', False)
    
    if lsm_kind == 'TNO':
      create_lsm_input_tno(x_sw, y_sw, input_coarse['grid']['itot'], input_coarse['grid']['jtot'], dx, dy,
                       input_coarse['nprocx'], input_coarse['nprocy'], start_date,
                       input_coarse['outpath'], input_coarse['LSM']['ERA5_path'], input_coarse['LSM']['spatial_data_path'],
                       input_coarse['iexpnr'], lwrite_ags=lwrite_ags)
      print('Finished creating LSM input (TNO flavor)')
    else:
      create_lsm_input(x_sw, y_sw, input_coarse['grid']['itot'], input_coarse['grid']['jtot'], dx, dy,
                       input_coarse['nprocx'], input_coarse['nprocy'], start_date,
                       input_coarse['outpath'], input_coarse['LSM']['ERA5_path'], input_coarse['LSM']['spatial_data_path'],
                       input_coarse['iexpnr'])
      print('Finished creating LSM input')

  setup_DALES(input_coarse['outpath'], input_coarse['data_path'])
  # relies on the data_path containing rrtmg*.nc.
  
  grid_param = input_coarse['grid']

  if 'x_offset' in grid_param and 'y_offset' in grid_param:
    x_offset = grid_param['x_offset']
    y_offset = grid_param['y_offset']
    xres = grid_param['xsize'] / grid_param['itot']
    yres = grid_param['ysize'] / grid_param['jtot']

    ix_west = int(x_offset / xres)
    ix_east = int(ix_west + grid_param['xsize_fine'] / xres)
    iy_south = int(y_offset / yres)
    iy_north = int(iy_south + grid_param['ysize_fine'] / yres)

    patch_namelist(input_coarse['outpath'], input_coarse, ix_west + 2, ix_east + 2, iy_south + 2, iy_north + 2, domain_type)
  else:
    patch_namelist(input_coarse['outpath'], input_coarse, '', '', '', '', 'coarse')
  

  
  #%% Advective time interpolation of input data (optional, to be implemented)

  #%% Create initial fields > initfields.inp.xxx.nc
  if(input_coarse['start']==input_coarse['time0']): # Not required for warmstarts
    initfields = initial_fields(input_coarse,grid,data,transform)
    print('finished initial fields')
    #%% Create profiles > prof.inp.xxx, lscale.inp.xxx scalar.inp.xxx
    profiles(input_coarse,grid,initfields,data)
    print('finished profiles')
  #%% Create boundary input > openboundaries.inp.xxx.nc
  openboundaries = boundary_fields(input_coarse,grid,data)
  print('finished boundary fields')
  #%% Create synthetic turbulence for boundary input (optional) > openboundaries.inp.xxx.nc
  if('synturb' in input_coarse):
    synturb = synthetic_turbulence(input_coarse,grid,data,transform)
    print('finished synthetic turbulence')
  #%% Create heterogeneous and time dependend skin temperature > tskin.inp.xxx.nc (if ltskin==true)
  if('tskin' in input_coarse):
    tskin = surface_temperature(input_coarse,grid,data,transform)
    print('finished surface temperature')
  

#%% Write data to input files
if('fine' in input):
  domain_type='fine'
  input_fine = input['fine']
  input_fine['author'] = input['author']
  #%% Create DALES grid
  grid = GridDales(input_fine['grid'])
  #%% Create initial fields > initfields.inp.xxx.nc
  if(input_fine['start']==input_fine['time0']): # Not required for warmstarts
    initfields_fine = initial_fields_fine(input_fine,grid)
    print('finished initial fields')
    #%% Create profiles > prof.inp.xxx, lscale.inp.xxx scalar.inp.xxx
    profiles(input_fine,grid,initfields_fine)
    print('finished profiles')
  #%% Create boundary input > openboundaries.inp.xxx.nc
  openboundaries_fine = boundary_fields_fine(input_fine,grid)
  print('finished boundary fields')
  
  if 'LSM' in input_fine:
  
    x_sw, y_sw = input_fine['x_fine'], input_fine['y_fine']
    print(f'LSM {x_sw}, {y_sw}.')
    dx = input_fine['grid']['xsize'] / input_fine['grid']['itot']
    dy = input_fine['grid']['ysize'] / input_fine['grid']['jtot']
    start_date = datetime.fromisoformat(input_fine['start'])

    lsm_kind = 'old'
    
    if 'lsm_kind' in input_fine['LSM']:
      lsm_kind = input_fine['LSM']['lsm_kind']

    if lsm_kind == 'TNO':
      create_lsm_input_tno(x_sw, y_sw, input_fine['grid']['itot'], input_fine['grid']['jtot'], dx, dy,
                       input_fine['nprocx'], input_fine['nprocy'], start_date,
                       input_fine['outpath'], input_fine['LSM']['ERA5_path'], input_fine['LSM']['spatial_data_path'],
                       input_fine['iexpnr'])
      print('Finished creating LSM input (TNO flavor)')
    else:
      create_lsm_input(x_sw, y_sw, input_fine['grid']['itot'], input_fine['grid']['jtot'], dx, dy,
                       input_fine['nprocx'], input_fine['nprocy'], start_date,
                       input_fine['outpath'], input_fine['LSM']['ERA5_path'], input_fine['LSM']['spatial_data_path'],
                       input_fine['iexpnr'])
      print('Finished creating LSM input')

  setup_DALES(input_fine['outpath'], input_fine['data_path'])
  
  #%% Apply spatial horizontal Gaussian filter to data
  #if('filter' in input_fine and input_fine['source'].lower() != 'none'):
    #data = gaussian_filter(data,input_fine)
    #data.to_netcdf(f"{input_fine['outpath']}/input_data.nc") # experiment: save the processed data

  #%% Create heterogeneous and time dependend skin temperature > tskin.inp.xxx.nc
  if('tskin' in input_fine):
    tskin_fine = surface_temperature_fine(input_fine,grid)
    print('finished surface temperature')
    
  patch_namelist(input_fine['outpath'], input_fine, '', '', '', '', domain_type)
  
  files_to_copy = [
        "backrad.inp.001.nc",
        "exnr.inp.001",
        "tracerdata.inp"
    ]

  for fname in files_to_copy:
    src = os.path.join(input_fine['inpath'], fname)
    dst = os.path.join(input_fine['outpath'], fname)

    if os.path.exists(src):
        shutil.copy2(src, dst)  # preserves metadata
# %%
