import xarray as xr
import numpy as np
import pyproj
from pyproj import Transformer, CRS, Proj, transform
from scipy.interpolate import griddata
from scipy.ndimage import gaussian_filter
import os
from scipy.spatial import Delaunay
from numba import jit
from pathlib import Path

# Define transformation class
class Transform:
    def __init__(self, parameters):
        self.parameters = parameters
        self.crs_latlon = 'epsg:4326'
        self.crs_rd = 'epsg:28992'
        # Construct transformation objects
        self.latlon_to_xy_transform = Transformer.from_crs(self.crs_latlon, self.parameters['proj4'])
        self.xy_to_latlon_transform = Transformer.from_crs(self.parameters['proj4'], self.crs_latlon)
        self.rd_to_latlon_transform = Transformer.from_crs(self.crs_rd, self.crs_latlon)
        self.latlon_to_rd_transform = Transformer.from_crs(self.crs_latlon, self.crs_rd)
        self.rd_to_lcc_transform = Transformer.from_crs(self.crs_rd, self.parameters['proj4'])

    def latlon_to_xy(self, lat, lon):
        return self.latlon_to_xy_transform.transform(lat, lon)

    def xy_to_latlon(self, x, y):
        return self.xy_to_latlon_transform.transform(x, y)
    
    def rd_to_latlon(self, x, y):
        return self.rd_to_latlon_transform.transform(x, y)
    
    def latlon_to_rd(self, lat, lon):
        return self.latlon_to_rd_transform.transform(lat, lon)

    def rd_to_lcc(self, x, y):
        return self.rd_to_lcc_transform.transform(x, y)
        
proj_params = {
    'proj4': '+proj=sterea +lat_0=52.15616055555555 +lon_0=5.38763888888889 +k=0.9999079 +x_0=155000 +y_0=463000 +ellps=bessel +towgs84=565.417,50.3319,465.552,-0.398957,0.343988,-1.8774,4.0725 +units=m +no_defs'
}

proj4_str_rd = '+proj=sterea +lat_0=52.15616055555555 +lon_0=5.38763888888889 +k=0.9999079 +x_0=155000 +y_0=463000 +ellps=bessel +towgs84=565.417,50.3319,465.552,-0.398957,0.343988,-1.8774,4.0725 +units=m +no_defs'

def precompute_barycentric_weights(lon2d, lat2d, x_target, y_target):
    points = np.column_stack((lon2d.ravel(), lat2d.ravel()))
    tri = Delaunay(points)

    target_points = np.column_stack((x_target.ravel(), y_target.ravel()))
    simplex = tri.find_simplex(target_points)  # triangle indices for each target point

    X = tri.transform[simplex, :2]
    Y = target_points - tri.transform[simplex, 2]
    bary = np.einsum('ijk,ik->ij', X, Y)  # compute barycentric coords
    bary_coords = np.c_[bary, 1 - bary.sum(axis=1)]  # w0, w1, w2

    vertices = tri.simplices[simplex]  # triangle vertex indices

    mask = simplex != -1  # mask for valid triangles

    return vertices, bary_coords, mask



@jit(nopython=True, nogil=True)
def interpolate_curvilinear_numba(values, vertices, bary_coords, mask, out):
    for i in range(out.size):
        if mask[i]:
            v0, v1, v2 = vertices[i]
            w0, w1, w2 = bary_coords[i]
            out[i] = w0 * values[v0] + w1 * values[v1] + w2 * values[v2]
        else:
            out[i] = np.nan

transform = Transform(proj_params)

# Projection parameters dictionary
proj_params_output = {
    "grid_mapping_name": "oblique_stereographic",  # EPSG 28992 projection name
    "latitude_of_projection_origin": float(proj4_str_rd.split("+lat_0=")[1].split()[0]),
    "longitude_of_central_meridian": float(proj4_str_rd.split("+lon_0=")[1].split()[0]),
    "scale_factor_at_projection_origin": float(proj4_str_rd.split("+k=")[1].split()[0]),
    "false_easting": float(proj4_str_rd.split("+x_0=")[1].split()[0]),
    "false_northing": float(proj4_str_rd.split("+y_0=")[1].split()[0]),
    "ellipsoid": proj4_str_rd.split("+ellps=")[1].split()[0],
    "towgs84": proj4_str_rd.split("+towgs84=")[1].split()[0],
    "units": proj4_str_rd.split("+units=")[1].split()[0],
    "proj4": proj4_str_rd
}

# Metadata dictionary for variable attributes
variable_metadata = {
    'clw': {'long_name': 'Specific cloud liquid water content', 'units': 'kg kg-1'},
    'tke': {'long_name': 'Turbulent Kinetic Energy', 'units': 'J kg-1'},
    'ta': {'long_name': 'Temperature', 'units': 'K'},
    'ua': {'long_name': 'u-component of wind', 'units': 'm s-1'},
    'va': {'long_name': 'v-component of wind', 'units': 'm s-1'},
    'wa': {'long_name': 'Geometrical vertical velocity', 'units': 'm s-1'},
    'hus': {'long_name': 'Specific humidity', 'units': 'kg kg-1'},
    'hfss': {'long_name': 'Accumulated Surface Upward Sensible Heat Flux', 'units': 'J m-2'},
    'huss': {'long_name': 'Specific humidity at 2m', 'units': 'kg kg-1'},
    'ps': {'long_name': 'Surface Air Pressure', 'units': 'Pa'},
    'tas': {'long_name': 'Temperature at 2m', 'units': 'K'},
    'tauu': {'long_name': 'Accumulated Momentum flux, u-component', 'units': 'N m-2'},
    'tauv': {'long_name': 'Accumulated Momentum flux, v-component', 'units': 'N m-2'},
    'cb': {'long_name': 'Cloud base', 'units': 'm'},
    'qr': {'long_name': 'Rain water mixing ratio', 'units': 'kg kg-1'}
}

#Load ROOTDIR from bash script (where data is stored), default = ~ (home directory)
rootdir = os.environ.get('PATH_NC_COMP','~')

#Load OUTDIR from bash script (where output is stored), default = ROOTDIR/merger
outdir = os.environ.get('OUTDIR', f'{rootdir}/RD/')

outdir = Path(outdir)

HA_config = os.environ.get('HA_CONFIG', 'N20')
HA_cycle= os.environ.get('HA_CYCLE', 'HA43')

# Create folder if it doesn't exist
outdir.mkdir(parents=True, exist_ok=True)

#Load RUNNR from bash script, default = 001
startdate = os.environ.get('STARTDATE', f'20000101')

# Define the file path template for each variable
input_file_template = f"{rootdir}/{{}}_{HA_cycle}_{HA_config}_NETHERLANDS_start_at_{startdate}00.1hr.nc"
output_folder = f"{outdir}/"

# Project lat/lon to x/y coordinates in meters
example_ds = xr.open_dataset(input_file_template.format('ps'))
lat_orig = example_ds["lat"].values
lon_orig = example_ds["lon"].values

lat_min=lat_orig.min() + 1
lon_min=lon_orig.min() + 1

x_proj, y_proj = transform.latlon_to_rd(lat_orig, lon_orig )

# Compute approximate spacing from data
dx_orig = np.median(np.diff(x_proj, axis=1))  # spacing along x
dy_orig = np.median(np.diff(y_proj, axis=0))  # spacing along y

# Round spacing to nearest 1000 m
dx = np.round(dx_orig / 1000) * 1000
dy = np.round(dy_orig / 1000) * 1000

x_proj_min, y_proj_min = transform.latlon_to_rd(lat_min, lon_min )
x_proj_min = np.round(x_proj_min / dx) * dx
y_proj_min = np.round(y_proj_min / dy) * dy

x_proj = np.round(x_proj / 5) * 5
y_proj = np.round(y_proj / 5) * 5

x_grid = np.arange(x_proj_min, x_proj_min+(250*dx), dx)
y_grid = np.arange(y_proj_min, y_proj_min+(250*dy), dy)

x_grid_2d, y_grid_2d = np.meshgrid(x_grid, y_grid)
#
lat_rd,lon_rd= transform.rd_to_latlon(x_grid_2d, y_grid_2d)

target_mesh_x = x_proj[0,:]
target_mesh_y = y_proj[:,0]

# Precompute barycentric weights
vertices, bary_coords, mask = precompute_barycentric_weights(
    lon_orig, lat_orig, lon_rd, lat_rd
)

# Loop over each variable in variable_metadata
for var_name, var_meta in variable_metadata.items():
    input_file_path = input_file_template.format(var_name)
     
    # Check if the input file exists
    if not os.path.exists(input_file_path):
        print(f"File for variable '{var_name}' not found: {input_file_path}, skipping.")
        continue
        
    print(f"variable '{var_name}' is processing.....")

    # Load the data for the current variable
    ds = xr.open_dataset(input_file_path)
    var_data = ds[var_name].values
    dimensions = var_data.shape
    
    # For 3D: (time, lat, lon)
    if len(var_data.shape) == 3:
        t_steps = var_data.shape[0]
        interpolated_var_new = np.zeros((t_steps, y_grid.shape[0], x_grid.shape[0]))
        for t in range(t_steps):
            
            # Prepare field
            field = var_data[t, :, :].ravel()  # shape (nx * ny)
            
            field=np.where(np.isnan(field), 0.0, field)
            
            out = np.empty(lon_rd.size)

            
            interpolate_curvilinear_numba(field, vertices, bary_coords, mask, out)

            # Reshape back
            interpolated = out.reshape(lon_rd.shape)

            interpolated_var_new[t,:,:] = interpolated
            
    # For 4D: (time, level, lat, lon)
    elif len(var_data.shape) == 4:
        t_steps, levels = var_data.shape[:2]
        interpolated_var_new = np.zeros((t_steps, levels, y_grid.shape[0], x_grid.shape[0]))
        
        for t in range(t_steps):
            for lev in range(levels):
                # Prepare field
                field = var_data[t, lev, :, :].ravel()  # shape (nx * ny)
                out = np.empty(lon_rd.size)

            
                interpolate_curvilinear_numba(field, vertices, bary_coords, mask, out)

                # Reshape back
                interpolated = out.reshape(lon_rd.shape)

                interpolated_var_new[t,lev,:,:] = interpolated



    # Create a new xarray dataset with interpolated data
    coords = {
        "time": (["time"], ds['time'].values),
        "y": (["y"], y_grid, {
            "units": "meters", "axis": "Y", "standard_name": "projection_y_coordinate", "long_name": "Y Coordinate Of Projection"
        }),
        "x": (["x"], x_grid, {
            "units": "meters", "axis": "X", "standard_name": "projection_x_coordinate", "long_name": "X Coordinate Of Projection"
        }),
        "lat": (["y", "x"], lat_rd, {"long_name": "Latitude", "units": "degrees_north"}),
        "lon": (["y", "x"], lon_rd, {"long_name": "Longitude", "units": "degrees_east"})
    }

    if len(dimensions) == 4:  # Include lev dimension if 4D
        coords["lev"] = (["lev"], ds['lev'].values, {"units": "dimensionless", "long_name": "Hybrid Level"})
    
    ds_var = xr.Dataset(
        {
            var_name: (["time", "y", "x"] if len(dimensions) == 3 else ["time", "lev", "y", "x"], interpolated_var_new, {
                "long_name": var_meta['long_name'], "units": var_meta['units']
            })
        },
        coords=coords
    )

    # Add the Lambert Conformal projection as a coordinate system variable
    ds_var['RD_coordinates'] = xr.DataArray(0, attrs=proj_params_output)

    # Save each interpolated variable to a separate NetCDF file
    output_file_path = os.path.join(output_folder, f"{var_name}_{HA_cycle}_{HA_config}_NETHERLANDS_start_at_{startdate}00.1hr_RD.nc")
    ds_var.to_netcdf(output_file_path)
    print(f"Interpolated dataset for '{var_name}' saved to {output_file_path}")
