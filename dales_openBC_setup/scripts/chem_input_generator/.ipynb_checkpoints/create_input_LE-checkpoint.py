# Script to create DALES chemical boundaries from LOTOS-EUROS output

#NOTE: this script should be run after the meteo boundaries preparation is completed, as the program opens 
#initfield and openboundaries nc-files, reads transform projection, and saves chemical boundaries.

#Creator: Arseni Doyennel VUA/IHS EUR., Feb. 2025

import netCDF4 as nc4
import xarray as xr
import numpy as np
import datetime
from pyproj import Proj, CRS, Transformer, Geod
import glob
import os
import json
import sys
import pandas as pd
from scipy.ndimage import gaussian_filter
from scipy.interpolate import interp1d
import interpolate_LE as ip
from GridDales import GridDales
from dateutil import parser

#TO-DO: add setup parameters to the main json script to read from there..
class Setup:

    def __init__(self, input_file='input_coarse.json'): #../
        # Load setup data from JSON
        with open(input_file, 'r') as f:
            input_data = json.load(f)

        self.input = input_data
        
        # File paths
        self.input_type = 'CAMS'
        self.input_dir = '/Users/darseni/Documents/work/from_vu_mac/CAMS_data/'      #'/Users/darseni/Downloads/transfer_3223116_files_3631045d/' 
        self.input_LE =  'cams73_co2_conc_region_1hr_final_05_2022.nc'        #'LE_cams_v8_1_GrETA_conc-3d_final_modified.nc'
        self.input_coarse = input_data['coarse']
        self.output_dir = self.input_coarse['outpath']
        self.iexpnr = f"{self.input_coarse['iexpnr']:03d}"
        # Target grid setup
        self.target_grid = GridDales(self.input_coarse['grid'])
        
        # To calulate sigma for gausian giltering of data at target resolution
        
        # --- Open LOTOS-EUROS dataset (once!) ---
        ds = xr.open_dataset(f"{self.input_dir}{self.input_LE}")

        # --- Estimate native horizontal resolution (in meters) ---
        geod = Geod(ellps="WGS84")
        self.lat_org = ds['latitude'].values
        self.lon_org = ds['longitude'].values

        # Approximate mean Δlat/Δlon resolution at domain center
        lat_mid = float(np.median(self.lat_org))
        lon_mid = float(np.median(self.lon_org))
        dy = geod.inv(lon_mid, lat_mid, lon_mid, lat_mid + (self.lat_org[1] - self.lat_org[0]))[2]
        dx = geod.inv(lon_mid, lat_mid, lon_mid + (self.lon_org[1] - self.lon_org[0]), lat_mid)[2]

        self.original_resolution_x = dx
        self.original_resolution_y = dy

        # --- Define Gaussian filter sigmas dynamically ---
        # You can tune the multiplier if you want stronger smoothing
        self.sigma_x_original = 2.5 * dx / 1000.0  # in km
        self.sigma_y_original = 2.5 * dy / 1000.0  # in km

        print(f"Approx. LOTOS-EUROS horisontal resolution in meters: dx={dx:.1f} m, dy={dy:.1f} m")
        print(f"Using Gaussian sigmas: σx={self.sigma_x_original:.2f} km, σy={self.sigma_y_original:.2f} km")

        # --- Set up tracer names ---
        self.tracer_names = [
            'co2'
        ]

        # --- Build dataset dictionary ---
        # Instead of re-opening files, reuse `ds`
        self.datasets = {name: ds for name in self.tracer_names}
        
        if self.input_type == 'LE':
            self.datasets['altitude'] = ds['altitude']
        elif self.input_type == 'CAMS':
            self.datasets['altitude'] = ds['height_above_reference_ellipsoid']
        
        # DALES constants (modglobal.f90)
        self.cd = dict(p0=1e5, Rd=287.04, Rv=461.5, cp=1004., Lv=2.53e6)
        self.cd['eps'] = self.cd['Rv'] / self.cd['Rd'] - 1.0

        # Projection settings
        self.proj4 = self.get_proj4()
        self.parameters = {'proj4': CRS.from_proj4(self.proj4)}
        
         
        
    def get_proj4(self):
        """
        Try to read the proj4 string from an existing initfields NetCDF file.
        If the file or the information is missing, return the default projection.
        """
        # Define default projection (Rijksdriehoek / Amersfoort)
        default_proj4 = (
            "+proj=sterea +lat_0=52.15616055555555 +lon_0=5.38763888888889 "
            "+k=0.9999079 +x_0=155000 +y_0=463000 +ellps=bessel "
            "+towgs84=565.417,50.3319,465.552,-0.398957,0.343988,-1.8774,4.0725 "
            "+units=m +no_defs"
        )

        file_path = os.path.join(self.output_dir, f"initfields.inp.{self.iexpnr}.nc")

        # If file doesn’t exist — fallback to default
        if not os.path.exists(file_path):
            print(f" File not found: {file_path}. Using default projection.")
            return default_proj4

        try:
            with nc4.Dataset(file_path, "r") as dataset:
                # Check for a variable named "transform"
                transform = dataset.variables.get("transform", None)

                if transform is not None:
                    # Try direct proj4 attribute
                    proj4_str = getattr(transform, "proj4", None)
                    if proj4_str:
                        return proj4_str

                    # Try any other attribute containing "proj4"
                    for attr, value in transform.__dict__.items():
                        if "proj4" in attr.lower() and isinstance(value, str):
                            return value

                # If not found in variable, look at global attributes
                for attr_name in dataset.ncattrs():
                    value = getattr(dataset, attr_name)
                    if "proj4" in attr_name.lower() and isinstance(value, str):
                        return value

                # Nothing found → fallback
                print(f" No proj4 information found in {file_path}. Using default projection.")
                return default_proj4

        except Exception as e:
            print(f" Error reading {file_path}: {e}. Using default projection.")
            return default_proj4


class Transform:
    def __init__(self, parameters):
        self.parameters = parameters
        self.crs_latlon = 'epsg:4326'
        # Construct transformation objects
        self.latlon_to_xy_transform = Transformer.from_crs(self.crs_latlon, self.parameters['proj4'])
        self.xy_to_latlon_transform = Transformer.from_crs(self.parameters['proj4'], self.crs_latlon)

    def latlon_to_xy(self, lat, lon):
        return self.latlon_to_xy_transform.transform(lat, lon)

    def xy_to_latlon(self, x, y):
        return self.xy_to_latlon_transform.transform(x, y)



def calculate_sigma_for_target_resolution(sigma_x_original, sigma_y_original, original_res_x, original_res_y, target_resolution_x, target_resolution_y):
    scaling_factor_x = target_resolution_x / original_res_x
    scaling_factor_y = target_resolution_y / original_res_y
    sigma_x_target = sigma_x_original * scaling_factor_x
    sigma_y_target = sigma_y_original * scaling_factor_y
    return sigma_x_target, sigma_y_target


def get_time_index(dataset, date_run, hour_offset):
    try:
        time_var = dataset['time']
    except KeyError:
        raise KeyError("The dataset does not have a 'time' variable. Available variables: "
                       f"{list(dataset.variables.keys())}")
    
    time_values = pd.to_datetime(time_var.values)
    target_time = date_run + datetime.timedelta(hours=hour_offset)
    if target_time in time_values:
        return list(time_values).index(target_time)
    else:
        raise ValueError(f"Time {target_time} not found in the dataset. Available times: {time_values}")


def write_initfield_nc(dataset, tracer_name, target_grid, hour, minutes, iexpnr,
                       output_dir, transform, sigma_x, sigma_y):
    """
    Create or update initfields NetCDF file for DALES initial conditions.
    """
    name = os.path.join(output_dir, f"initfields.inp.{iexpnr}.nc")

    # Create the file if it doesn't exist
    if not os.path.exists(name):
        print(f"File not found — creating new NetCDF: {name}")
        with nc4.Dataset(name, "w") as f:
            # Create dimensions
            f.createDimension("xt", len(target_grid.xt))
            f.createDimension("yt", len(target_grid.yt))
            f.createDimension("zt", len(target_grid.zt))

            # Create coordinate variables
            f.createVariable("xt", "f4", ("xt",))[:] = target_grid.xt
            f.createVariable("yt", "f4", ("yt",))[:] = target_grid.yt
            f.createVariable("zt", "f4", ("zt",))[:] = target_grid.zt

    # Open in append mode
    with nc4.Dataset(name, "a") as f:
        # Check dimensions
        assert len(f.dimensions["xt"]) == len(target_grid.xt), "Dimension mismatch for 'xt'"
        assert len(f.dimensions["yt"]) == len(target_grid.yt), "Dimension mismatch for 'yt'"
        assert len(f.dimensions["zt"]) == len(target_grid.zt), "Dimension mismatch for 'zt'"

        # Apply Gaussian smoothing (horizontal)
        sv_smooth = gaussian_filter(dataset, sigma=(sigma_x, sigma_y, 1))
        sv_rotated = np.transpose(sv_smooth, axes=(2, 1, 0))  # reorder for DALES (zt, yt, xt)

        # Create or update tracer variable
        if tracer_name not in f.variables:
            var = f.createVariable(tracer_name, "f8", ("zt", "yt", "xt"), fill_value=np.nan)
            var.longname = f"Chemical scalar tracer {tracer_name}"
            var.units = (
                "ppm" if tracer_name == "co2"
                else "kg m-3" if tracer_name in ["pm25", "pm10", "bc"]
                else "ppb"
            )
        else:
            var = f.variables[tracer_name]

        var[:, :, :] = sv_rotated

    print(f"NetCDF file updated: {name}")


def save_boundaries_to_nc(boundary_data, elapsed_times, output_dir, iexpnr, tracer_name,
                          initial_time_DALES, sigma_x, sigma_y, target_grid=None):
    """
    Create or update open boundaries NetCDF file for DALES boundary input.
    """
    file_path = os.path.join(output_dir, f"openboundaries.inp.{iexpnr}.nc")

    # Create the file if it doesn't exist
    if not os.path.exists(file_path):
        print(f"File not found — creating new NetCDF: {file_path}")
        with nc4.Dataset(file_path, "w") as f:
            # Create all standard DALES boundary dimensions
            f.createDimension("time", None)
            if target_grid is not None:
                f.createDimension("xt", len(target_grid.xt))
                f.createDimension("yt", len(target_grid.yt))
                f.createDimension("zt", len(target_grid.zt))

                # Create coordinate variables (optional but useful)
                f.createVariable("xt", "f4", ("xt",))[:] = target_grid.xt
                f.createVariable("yt", "f4", ("yt",))[:] = target_grid.yt
                f.createVariable("zt", "f4", ("zt",))[:] = target_grid.zt
            else:
                print("target_grid not provided — cannot define xt, yt, zt dimensions!")

            # Create the time variable
            time_var = f.createVariable("time", "f8", ("time",), fill_value=np.nan)
            time_var.longname = "Time"
            time_var.units = f"seconds since {initial_time_DALES.isoformat()}"

    # Open in read/write mode
    with nc4.Dataset(file_path, "r+") as f:
        # Ensure the time variable exists
        if "time" not in f.variables:
            time_var = f.createVariable("time", "f8", ("time",), fill_value=np.nan)
            time_var.longname = "Time"
            time_var.units = f"seconds since {initial_time_DALES.isoformat()}"

        # Write time data
        f.variables["time"][:] = elapsed_times

        # Write boundary variables
        for border, data in boundary_data.items():
            var_name = f"{tracer_name}{border}"

            if var_name not in f.variables:
                # Define shape depending on boundary orientation
                if border == "top":
                    shape = ("time", "yt", "xt")
                elif border in ["east", "west"]:
                    shape = ("time", "zt", "yt")
                else:  # north or south
                    shape = ("time", "zt", "xt")

                # Ensure all required dimensions exist
                for dim in shape:
                    if dim not in f.dimensions:
                        raise ValueError(f"Missing dimension '{dim}' in {file_path}")

                # Create variable
                var = f.createVariable(var_name, "f4", shape, fill_value=np.nan)
                var.longname = f"{tracer_name} at {border} boundary"
                var.units = (
                    "ppm" if tracer_name == "co2"
                    else "kg m-3" if tracer_name in ["pm25", "pm10", "bc"]
                    else "ppb"
                )
            else:
                var = f.variables[var_name]

            # Write data
            f[var_name][:] = data

    print(f"NetCDF file updated: {file_path}")


# Main execution
if __name__ == "__main__":
    setup = Setup()

    # Initialize transformation
    transform = Transform(setup.parameters)
    
    x_sw, y_sw = transform.latlon_to_xy(setup.input_coarse['lat_sw'], setup.input_coarse['lon_sw'])
    x_sw = np.round(x_sw, 0)
    y_sw = np.round(y_sw, 0)
    
    target_resolution_x = setup.target_grid.xt[1] - setup.target_grid.xt[0]
    target_resolution_y = setup.target_grid.yt[1] - setup.target_grid.yt[0]
        
    sigma_x, sigma_y = calculate_sigma_for_target_resolution(setup.sigma_x_original, setup.sigma_y_original, setup.original_resolution_x, setup.original_resolution_y,
                                                              target_resolution_x, target_resolution_y)
    
    
    #############################################################
    
    # Create 2D coordinate arrays for all x/y points
    xx, yy = np.meshgrid(
       x_sw + setup.target_grid.xt,
       y_sw + setup.target_grid.yt,
       indexing='ij'
    )

    # Convert all LES grid coordinate to geo lat and lon at once (vectorized)
    les_lats, les_lons  = transform.xy_to_latlon(xx, yy)

    start_date = datetime.datetime.fromisoformat(setup.input['coarse']["start"])
    end_date = datetime.datetime.fromisoformat(setup.input['coarse']["end"])
    time0_date = datetime.datetime.fromisoformat(setup.input['coarse']["time0"])
    

    # Extract t0 as the hour of start_date
    t0 = start_date.hour

    # Calculate the total number of hours in your simulation
    elapsed_hours = int((end_date - start_date).total_seconds() / 3600)

    # Compute t1
    t1 = t0 + elapsed_hours
    
    print(f"Start hour (t0): {t0}")
    print(f"End hour (t1): {t1}")
    
    #t1=3 #for testing

    
    year_of_int = start_date.year
    month_of_int = start_date.month
    day_of_int = start_date.day

    year_of_int_path = str(year_of_int)
    month_of_int_path = f"{month_of_int:02d}"
    day_of_int_path = f"{day_of_int:02d}"

    print("year/month/day of interest:", year_of_int_path, month_of_int_path, day_of_int_path)

    # Step 2: Open the NetCDF file and extract time attributes
    file_path = f'{setup.input_dir}{setup.input_LE}'

    with nc4.Dataset(file_path, 'r') as dataset:
        time_var = dataset.variables.get('time', None)
        if time_var is None:
            raise ValueError("Time variable not found in the NetCDF file.")

        # Extract relevant time attributes
        time_units = time_var.units  # Example: "hours since 1900-01-01 00:00:00.0"
        #time_calendar = getattr(time_var, 'calendar', 'gregorian')
        time_calendar = getattr(time_var, 'calendar', 'standard')

    # Step 3: Parse the units string to extract the reference date
    #Extract the time reference string (remove "UTC")
    time_ref_str = time_units.split("since")[1].strip()  # "2018-01-01 00:00:00 UTC"

    #Remove the 'UTC' part
    time_ref_str = time_ref_str.replace(" UTC", "")  # "2018-01-01 00:00:00"

    #Parse the time reference without microseconds
    #time_ref = datetime.datetime.strptime(time_ref_str, "%Y-%m-%d %H:%M:%S")
    time_ref = parser.isoparse(time_ref_str)


    # Calculate the initial date of the NetCDF time variable
    initial_time = time_ref + datetime.timedelta(hours=0)

    initial_year_cams = initial_time.year
    initial_month_cams = initial_time.month
    initial_day_cams = initial_time.day

    # Step 4: Adjust date of interest based on initial reference day
    date_run = datetime.datetime(year=year_of_int, month=month_of_int, day=day_of_int)
    cams_tref = initial_time

    print("Date run:", date_run)
    print("LE reference date:", cams_tref)    
    
    if setup.input_type == 'LE':
        altitude_xr = setup.datasets["altitude"] #extract altitudes from LE dataset
        altitude = altitude_xr[0, :, :, :].values
    elif setup.input_type == 'CAMS':
        # Height at interfaces
        z_interfaces = setup.datasets["altitude"][0,:,:,:].values
        z_levels = 0.5 * (z_interfaces[:-1] + z_interfaces[1:])
        altitude = z_levels
        

    intp = ip.GridInterpolatorLE(
                        setup.lon_org,  # original LE longitudes x_LS
                        setup.lat_org,   # original LE latitudes y_LS
                        setup.target_grid.zt,  # z_LS (no vertical interpolation here!)
                        les_lons,   # 2-d lon corresponding to DALES x cooridnates (target grid)
                        les_lats,   # 2-d lat corresponding to DALES y cooridnates (target grid)
                        setup.target_grid.zt,  # z (same as z_LS)
                        x_sw, y_sw
                    )
    
    # Determine the number of timesteps
    num_timesteps = t1 - t0+1
    
    
    for tracer_name in setup.tracer_names:
    
    
        # Preallocate arrays for boundary data
        boundary_data = {
            'top': np.empty((num_timesteps, len(setup.target_grid.yt), len(setup.target_grid.xt))),  # (time, yt, xt)
            'east': np.empty((num_timesteps, len(setup.target_grid.zt), len(setup.target_grid.yt))), # (time, zt, yt)
            'west': np.empty((num_timesteps, len(setup.target_grid.zt), len(setup.target_grid.yt))), # (time, zt, yt)
            'north': np.empty((num_timesteps, len(setup.target_grid.zt), len(setup.target_grid.xt))), # (time, zt, xt)
            'south': np.empty((num_timesteps, len(setup.target_grid.zt), len(setup.target_grid.xt)))  # (time, zt, xt)
        }

        #Preallocate array for elapsed times
        elapsed_times = np.empty(num_timesteps)


        # Iterate over timesteps
        for i, t in enumerate(range(t0, t1+1)):  # i tracks the timestep index  t1+1
    
    
            print(f'Processing t={t:>2d}:00 UTC')
            start_time = datetime.datetime.now()
    
            # -- Only execute if timestep exists
            try: 
                # Extract dataset variable
                var_ds = setup.datasets[tracer_name][tracer_name]
                
                # Extract the tracer data at the current time index
                time_index = get_time_index(var_ds, date_run, t)
                raw_data = var_ds[time_index, :, :, :].values

                print(f"Indices for t={t}: index={time_index}")
            
                # Extract data for the given time index and convert from mole moleˆ-1 to ppm/ppb
                if tracer_name=='co2':
                    tracer_data = raw_data * 1e6  # convert to ppm
                elif tracer_name in ['c2h6', 'pm25', 'bc']: #they are in kg m-3 in LE..
                    tracer_data = raw_data  # kg/m3
                else:
                    tracer_data = raw_data * 1e9  # ppb
                
                print(f"Vertical interpolation to DALES grid")
        
                # Perform the vertical interpolation first
                if setup.input_type == 'LE':
                    interpolated_vert = ip.interpolate_vert_mass_conserving(tracer_data, altitude, setup.target_grid.zt)
                elif setup.input_type == 'CAMS':
                    interpolated_vert = ip.interpolate_vert_mass_conserving_CAMS(tracer_data, altitude, setup.target_grid.zt)
                             
                #interpolated_vert is (lat, lon, z)
            
                print(f"Horizontal interpolation to DALES grid")

                interpolated_vert_rotated = np.transpose(interpolated_vert, axes=(1, 0, 2)) #we need lon,lat,z demension order
            
                # Interpolate LE data onto LES grid
                LES_interpolated = intp.interpolate_2d(interpolated_vert_rotated)
                                
                if t == t0:
            
                    print(f'Write initial field to nc for t={t}')
                    #3-D initial fields for chem. are saved to initfields file with meteo input ):
                    write_initfield_nc(LES_interpolated, tracer_name, setup.target_grid, t - t0, 0., setup.iexpnr, setup.output_dir, transform, sigma_x, sigma_y) 
        


                # Calculate elapsed time in seconds, considering hours > 24 as next days
                total_hours = t0 + t  # Accumulate hours beyond a single day
                days = total_hours // 24  # Determine the number of days
                hours = total_hours % 24  # Determine the hours within the current day
                current_time = datetime.datetime(year_of_int, month_of_int, day_of_int + days, hours)

                elapsed_seconds = (current_time -  date_run).total_seconds()
                
                
                elapsed_times[i] = elapsed_seconds

                # Smooth and rotate scalar field
                sv_smooth = gaussian_filter(LES_interpolated, sigma=(sigma_x, sigma_y, 1))
                sv_rotated = np.transpose(sv_smooth, axes=(2, 1, 0))  # ('zt', 'yt', 'xt')

                # Store boundary data in preallocated arrays
                boundary_data['top'][i, :, :] = sv_rotated[-1, :, :]  # (yt, xt)
                boundary_data['east'][i, :, :] = sv_rotated[:, :, -1]  # (zt, yt)
                boundary_data['west'][i, :, :] = sv_rotated[:, :, 0]   # (zt, yt)
                boundary_data['north'][i, :, :] = sv_rotated[:, -1, :] # (zt, xt)
                boundary_data['south'][i, :, :] = sv_rotated[:, 0, :]  # (zt, xt)

                print(f'Processed timestep t={t}')
                
                
                get_time_index(setup.datasets[tracer_name], date_run, t)
        

                # Statistics
                end_time = datetime.datetime.now()
                print('Elapsed = {}'.format(end_time-start_time))
      
            
            except Exception as e:
                print(f"Error processing t={t} e={e}")                                        
        try: 
        
            print(f'Write boundaries to nc')   
         
            save_boundaries_to_nc(boundary_data, elapsed_times, setup.output_dir, setup.iexpnr, tracer_name, date_run, sigma_x, sigma_y, setup.target_grid)
        
        
        except ValueError as e:
            
            print(f"Error: {e}")
