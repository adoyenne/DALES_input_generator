"""
Program Description:
This Python script processes hourly meteorological data from GRIB files and converts them into NetCDF format 
needed for preparing input data (on geo lat and lon) for the Dutch Atmospheric Large-Eddy Simulation (DALES) model. 
The program extracts key weather variables from GRIB files, such as temperature, humidity, wind speed, 
and surface fluxes, and save them in a NetCDF format for further use.

Key Components:

1. Reading and Preprocessing GRIB Data:

The script processes multiple GRIB files sequentially, corresponding to different hours of the forecast 
for the specified day. These GRIB files contain meteorological variables like temperature, wind components, 
and specific humidity.

Each GRIB file is opened using the pygrib library, and specific variables are extracted based on 
the parameter mappings defined earlier (e.g., temperature at 2 meters, surface pressure).


2. NetCDF File Creation:

For each GRIB variable, a corresponding NetCDF file is created.
The script handles both 3D variables (e.g., temperature or wind at different pressure levels) and
2D variables (e.g., surf. pressure surface and temperature at 2 meters).
Time information is recorded, and the data is stored in a time-evolving format where each variable can be tracked 
across multiple hours.

Note:
Currently, data are stored only in geo lat and lon coordinates!
Latitude and longitude are saved in 2D arrays and as x and y coordinates which are stored as 1D arrays.

Creator:
Dr. Arseni Doyennel (VU) 2024.
a.doyennel@vu.nl
"""

import pygrib
import numpy as np
from netCDF4 import Dataset, date2num
import datetime
from pyproj import Proj, Transformer
import netCDF4
import pyproj
from pathlib import Path
from scipy.ndimage import gaussian_filter1d
import os

#Load RUNNR from bash script, default = 001

# SETTINGS: Paths and day from GRIB file names
day_str = os.environ.get('STARTDATE', '20000101')  # The start date of HA forecast (in format YYYYMMDD)
HA_config = os.environ.get('HA_CONFIG', 'N20')
HA_cycle= os.environ.get('HA_CYCLE', 'HA43')
pattern = os.environ["HA_grib_file_pattern"]
pattern = pattern.replace("%s", f"{day_str}00", 1)  # replace only the first %s
date_format = '%Y%m%d'  # Date format for parsing
start_date = datetime.datetime.strptime(day_str, date_format)

# Initialize time tracking
forecast_hours = int(os.environ.get('LENGTH_COMP', 3))  # hours in each forecast

# ---- 1. Detect merge hours dynamically ----
first_cycle_len = 28  # 0-27
other_cycle_len = 25  # 3-27

merge_hours = [first_cycle_len - 1]  # first merge at 27
total_hours = forecast_hours
hour_cursor = first_cycle_len

while hour_cursor < total_hours:
    merge_hours.append(hour_cursor)  # merge at start of next cycle +2h offset handled
    hour_cursor += other_cycle_len

print("Automatically detected merge hours:", merge_hours)

path_to_grib_files = Path(os.path.expanduser(os.environ.get('PATH_GRIB_COMP', '~')))
path_to_nc_output = Path(os.path.expanduser(os.environ.get('PATH_NC_COMP', '~')))

# Create folder if it doesn't exist
path_to_nc_output.mkdir(parents=True, exist_ok=True)

grib_path_template = f"{path_to_grib_files}/{pattern.replace('%s', '{hour:03d}')}"

# Define the variable mapping from GRIB to NetCDF
# 3-D variables (time, level, lat, lon) on 'hybrid' levels
three_d_vars = {
    'clw': {'indicatorOfParameter': 76, 'typeOfLevel': 'hybrid', 'nc_name': 'clw'},
    'tke': {'indicatorOfParameter': 200, 'typeOfLevel': 'hybrid', 'nc_name': 'tke'},
    'ta':  {'indicatorOfParameter': 11, 'typeOfLevel': 'hybrid', 'nc_name': 'ta'},
    'ua':  {'indicatorOfParameter': 33, 'typeOfLevel': 'hybrid', 'nc_name': 'ua'},
    'va':  {'indicatorOfParameter': 34, 'typeOfLevel': 'hybrid', 'nc_name': 'va'},
    'hus': {'indicatorOfParameter': 51, 'typeOfLevel': 'hybrid', 'nc_name': 'hus'},
    'wa':  {'indicatorOfParameter': 40, 'typeOfLevel': 'hybrid', 'nc_name': 'wa'},
    'qr':  {'indicatorOfParameter': 181, 'typeOfLevel': 'hybrid', 'nc_name': 'qr'},
}

# 2-D variables (time, lat, lon)
two_d_vars = {
    'hfss': {'shortName': 'shf', 'nc_name': 'hfss'},
    'tas': {'indicatorOfParameter': 11, 'nc_name': 'tas', 'typeOfLevel': 'heightAboveGround','level': 2},
    'huss': {'indicatorOfParameter': 51, 'nc_name': 'huss', 'typeOfLevel': 'heightAboveGround','level': 2},
    'tauu': {'shortName': 'uflx', 'nc_name': 'tauu'},
    'tauv': {'shortName': 'vflx', 'nc_name': 'tauv'},
    'ps':   {'indicatorOfParameter': 1, 'nc_name': 'ps', 'typeOfLevel': 'heightAboveGround', 'level': 0},  # surface pressure
    'cb':   {'indicatorOfParameter': 186, 'nc_name': 'cb', 'level': 0}  # cloud base (masked array)
}

grib_to_nc_var = {**three_d_vars, **two_d_vars}

# Mapping of variable names to descriptions and units 
#(taken from https://hirlam.github.io/HarmonieSystemDocumentation/dev/ForecastModel/Outputlist/#Harmonie-GRIB1-code-table-2-version-253-Indicator-of-parameter)

variable_metadata = {
    'tke': {'long_name': 'Turbulent Kinetic Energy', 'units': 'J kg-1'},
    'hfss': {'long_name': 'Sensible heat flux', 'units': 'J m-2'},
    'lhe':   {'long_name': 'Latent heat flux', 'units': 'J m-2'},
    'huss': {'long_name': 'Specific humidity at 2m', 'units': 'kg kg-1'},
    'ps': {'long_name': 'Surface Air Pressure', 'units': 'Pa'},
    'tas': {'long_name': 'Temperature at 2m', 'units': 'K'},
    'tauu': {'long_name': 'Momentum flux, u-component', 'units': 'N m-2'},  #NOTE: previously we used tauu:long_name = "Accumulated Surface Downward Eastward Stress" ; tauu:units = "kg m-1 s-1" ;
    'tauv': {'long_name': 'Momentum flux, v-component', 'units': 'N m-2'},  #NOTE: previously we used tauv:long_name = "Accumulated Surface Downward Northward Stress" ;tauv:units = "kg m-1 s-1" ;
    'u10':   {'long_name': 'u-component of wind at 10m, relative to model coordinates', 'units': 'm s-1'},
    'v10':   {'long_name': 'v-component of wind at 10m, relative to model coordinates', 'units': 'm s-1'},
    'tcc': {'long_name': 'Total cloud cover', 'units': '0-1'},
    'rain':   {'long_name': 'Accumulated rain', 'units': 'kg m-2'},
    'mld':   {'long_name': 'Mixed layer depth', 'units': 'm'},
    'rh':   {'long_name': 'Relative humidity at 2m', 'units': '0-1'},
}

# Create the NetCDF file once for each variable
nc_files = {}
for var_name in grib_to_nc_var.keys():
    nc_file_path = f'{path_to_nc_output}/{var_name}_{HA_cycle}_{HA_config}_NETHERLANDS_start_at_{day_str}00.1hr.nc'
    nc_files[var_name] = Dataset(nc_file_path, 'w', format='NETCDF4')

time_values = []

# ---- Gaussian smoothing around merge edges ----
def smooth_edges_gaussian(var_array, merge_hours, window=3, sigma=2.0):
    """
    Smooths 2D or 3D variable array along the time axis only around merge_hours.
    var_array: shape (time, y, x) or (time, lev, y, x)
    """
    var_smooth = var_array.copy()
    half_window = window // 2
    for h in merge_hours:
        start = max(h - half_window, 0)
        end = min(h + half_window + 1, var_array.shape[0])
        # Apply Gaussian smoothing along time axis only
        var_smooth[start:end] = gaussian_filter1d(var_array[start:end], sigma=sigma, axis=0, mode='nearest')
    return var_smooth

# Function to process and store data from one GRIB file
def process_grib_file(grib_file, time_idx):
    stored_ps = False
    grbs = pygrib.open(grib_file)

    data_store = {var: {} for var in grib_to_nc_var}
    lats, lons = None, None

    for grb in grbs:
        for var_name, var_info in grib_to_nc_var.items():

            # Handle cloud base first (special case: masked array)
            if var_name == 'cb' and grb.indicatorOfParameter == 186:
                if 'typeOfLevel' in var_info and grb.typeOfLevel != var_info['typeOfLevel']:
                    continue
                if 'level' in var_info and grb.level != var_info['level']:
                    continue
                data_store[var_name] = np.where(grb.values.mask, np.nan, grb.values.data)
                if lats is None:
                    lats, lons = grb.latlons()
                continue

            # Handle surface pressure
            elif var_name == 'ps' and grb.indicatorOfParameter == 1 and grb.typeOfLevel == 'heightAboveGround' and grb.level == 0 and not stored_ps:
                data_store[var_name] = grb.values
                stored_ps = True
                if lats is None:
                    lats, lons = grb.latlons()
                continue

            # Handle 3D variables
            elif var_name in three_d_vars:
                if grb.indicatorOfParameter == var_info['indicatorOfParameter'] and grb.typeOfLevel == var_info['typeOfLevel']:
                    level = grb.level
                    if lats is None:
                        lats, lons = grb.latlons()
                    data_store[var_name][level] = grb.values
                continue

            # Handle other 2D variables
            # For 2D variables (except ps, which is handled above)
            elif var_name in two_d_vars:
                if 'indicatorOfParameter' in var_info:
                    if grb.indicatorOfParameter != var_info['indicatorOfParameter']:
                        continue

                if 'typeOfLevel' in var_info:
                    if grb.typeOfLevel != var_info['typeOfLevel']:
                        continue

                if 'level' in var_info:
                    if grb.level != var_info['level']:
                        continue

                data_store[var_name] = grb.values

                if lats is None:
                    lats, lons = grb.latlons()

    return data_store, lats, lons

# Create NetCDF files for each variable dynamically based on the date
def save_nc_files(ncfile, variable_data, lats, lons, variable_name, levels=None, time_idx=None, current_time=None):
    if time_idx == 0:
        # Create coordinate variables and dimensions only during the first call
        ncfile.createDimension('time', None)  # Unlimited dimension for time
        ncfile.createDimension('y', lats.shape[0])  # Number of latitude points
        ncfile.createDimension('x', lons.shape[1])  # Number of longitude points
        
        # Create 2D latitude variable
        latitudes = ncfile.createVariable('lat', 'f4', ('y', 'x'), zlib=True)
        latitudes.standard_name = "latitude"
        latitudes.long_name = "Latitude"
        latitudes.units = "degrees_north"
        latitudes[:] = lats   # Assuming lats is 2D array with shape (y_points, x_points)

        # Create 2D longitude variable
        longitudes = ncfile.createVariable('lon', 'f4', ('y', 'x'), zlib=True)
        longitudes.standard_name = "longitude"
        longitudes.long_name = "Longitude"
        longitudes.units = "degrees_east"
        longitudes[:] = lons  # Assuming lons is 2D array with shape (y_points, x_points)
        
        
        # Create 1D x variable
        xs = ncfile.createVariable('x', 'f4', ('x',), zlib=True)
        xs.units = "degrees_east"  # Adjust according to your definition
        xs.axis = "X"
        xs.standard_name = "longitude"
        xs.long_name = "Longitude"
        xs[:] = lons[0, :]  # Assign values from the 2D lon array (taking first row for x values)

        # Create 1D y variable
        ys = ncfile.createVariable('y', 'f4', ('y',), zlib=True)
        ys.units = "degrees_north"  # Adjust according to your definition
        ys.axis = "Y"
        ys.standard_name = "latitude"
        ys.long_name = "Latitude"
        ys[:] = lats[:, 0]  # Assign values from the 2D lat array (taking first column for y values)


        # Create a time variable
        time_var = ncfile.createVariable('time', 'f4', ('time',))
        time_var.units = 'hours since {0}'.format(start_date)

    # Add the current time to the time dimension
    time_var = ncfile.variables['time']
    time_var[time_idx] = current_time

    # Check if metadata is available for this variable
    metadata = variable_metadata.get(variable_name, {'long_name': 'Unknown', 'units': 'Unknown'})

    if levels is not None:
        # Handle 3D variables (time, lev, y, x)
        if 'lev' not in ncfile.dimensions:
            ncfile.createDimension('lev', len(levels))  # Create lev dimension here
            hybrid_levels = ncfile.createVariable('lev', 'f4', ('lev',))
            hybrid_levels.units = 'dimensionless'
            hybrid_levels.long_name = 'Hybrid Level'
            hybrid_levels[:] = levels

        # Get or create the main variable
        if variable_name not in ncfile.variables:
            var = ncfile.createVariable(variable_name, 'f4', ('time', 'lev', 'y', 'x'), zlib=True)
            var.long_name = metadata['long_name']
            var.units = metadata['units']
        else:
            var = ncfile.variables[variable_name]

        # Add data for the current time
        var[time_idx, :, :, :] = variable_data
    else:
        # Handle 2D variables (time, y, x)
        if variable_name not in ncfile.variables:
            var = ncfile.createVariable(variable_name, 'f4', ('time', 'y', 'x'), zlib=True)
            var.long_name = metadata['long_name']
            var.units = metadata['units']
        else:
            var = ncfile.variables[variable_name]

        # Add data for the current time
        var[time_idx, :, :] = variable_data

    print(f'NetCDF data updated for {variable_name} at time index {time_idx}')
    
for hour in range(forecast_hours):
    # Now, use the int of string in the file path
    grib_file = grib_path_template.format(hour=hour)  

    print(f'Processing GRIB file: {grib_file}')
    
    # Calculate current time in hours since the start of the day
    current_time = hour  # This would be in hours since start_date
    time_values.append(current_time)
    
    # Get data for this hour
    variable_data, lats, lons = process_grib_file(grib_file, hour)
    
    for var_name, data in variable_data.items():
        if isinstance(data, dict) and data:  # 3D variable with data
            if var_name == 'ps':
                print(f"Warning: ps found as dict, expected 2D array")
            
            levels = sorted(data.keys())
            levels = sorted(data.keys())
            n_levels = len(levels)
            n_lats, n_lons = data[levels[0]].shape

            # Initialize 3D array
            parameter_3d = np.zeros((n_levels, n_lats, n_lons))

            # Populate 3D array
            for i, level in enumerate(levels):
                parameter_3d[i, :, :] = data[level]
            
            # stack time axis first if multiple levels
            parameter_4d = np.expand_dims(parameter_3d, axis=0)  # (time=1, lev, y, x)
            parameter_4d_smooth = smooth_edges_gaussian(parameter_4d, merge_hours)
            parameter_3d = parameter_4d_smooth[0]

            # Append 3D variable data to the respective NetCDF file
            save_nc_files(nc_files[var_name], parameter_3d, lats, lons, var_name, levels, time_idx=hour, current_time=current_time)


        elif isinstance(data, np.ndarray):  # 2D variable
            data_smooth = smooth_edges_gaussian(np.expand_dims(data, axis=0), merge_hours)[0]
            # Append 2D variable data to the respective NetCDF file
            save_nc_files(nc_files[var_name], data_smooth, lats, lons, var_name, time_idx=hour, current_time=current_time)

        else:
            print(f"Skipping {var_name} at hour {hour}: no data")     
              
# Close all the NetCDF files after processing
for nc_file in nc_files.values():
    nc_file.close()

print("Processing and creation of nc-files completed.")
