
#To run this code, use:
# python -m DALES_domain_setup

from datetime import datetime, timedelta
from pathlib import Path

import numpy as np

import json

import matplotlib.pyplot as plt
from matplotlib.axes import Axes

from shapely.geometry import box

from pyproj import Transformer, transform

import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.geoaxes import GeoAxes
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER

from GridDales import GridDales

GeoAxes._pcolormesh_patched = Axes.pcolormesh

# This script's main folder
script_folder = Path(__file__).parent


# Define transformation class
class Transform:
    def __init__(self, parameters):
        self.parameters = parameters
        self.crs_latlon = 'epsg:4326'
        self.crs_rd = 'epsg:28992'
        # Construct transformation objects
        self.latlon_to_xy_transform = Transformer.from_crs(self.crs_latlon, self.parameters['proj4'])
        self.xy_to_latlon_transform = Transformer.from_crs(self.parameters['proj4'], self.crs_latlon)
        self.rd_to_lcc_transform = Transformer.from_crs(self.crs_rd, self.parameters['proj4'])

    def latlon_to_xy(self, lat, lon):
        return self.latlon_to_xy_transform.transform(lat, lon)

    def xy_to_latlon(self, x, y):
        return self.xy_to_latlon_transform.transform(x, y)

    def rd_to_lcc(self, x, y):
        return self.rd_to_lcc_transform.transform(x, y)


##################################################
#Main settings:
# Define projection parameters
proj_params = {
    'proj4': '+proj=sterea +lat_0=52.15616055555555 +lon_0=5.38763888888889 +k=0.9999079 +x_0=155000 +y_0=463000 +ellps=bessel '
             '+towgs84=565.417,50.3319,465.552,-0.398957,0.343988,-1.8774,4.0725 +units=m +no_defs'
}
transform = Transform(proj_params)

NESTING=True

author = ""

main_input_generator_folder_path = Path()

source_meteo = "harmonie_rd"

inpath_coarse = Path("/.../RD/")
outpath_coarse = Path("/.../")

inpath_fine= outpath_coarse
outpath_fine = Path("/.../")

ERA5_path = Path("/.../ERA5_soil/")
spatial_data_path = Path("/.../spatial_data/")

lsm_kind="TNO"

lwrite_ags=False

hybrid_lev_file="H43_90lev.txt"

time_start="2025-01-01T00:00"
time0=time_start #change if needed!
time_end="2025-01-02T00:00"

# convert string → datetime
time_start_dt = datetime.fromisoformat(time_start)

# add 2 hours
time_end_filter_dt = time_start_dt + timedelta(hours=2)

# convert back to string if needed
time_end_filter = time_end_filter_dt.isoformat()

print("time_start:", time_start)
print("time0:", time0)
print("time_end:", time_end)
print("time_end_filter:", time_end_filter)

#Check settings (at least the paths) & prep folders
##########################################################
outpath_coarse.mkdir(parents=True, exist_ok=True)
outpath_fine.mkdir(parents=True, exist_ok=True)
paths = [main_input_generator_folder_path, inpath_coarse, inpath_fine, ERA5_path, spatial_data_path]
for path in paths:
    if not path.exists():
        raise OSError(f"Path {path} does not exist. Correct the settings and run again.")


#External domain
#########################################
#West lower corner of the domain in approx. lat and lon:
sw_lat_ext = 51.5
sw_lon_ext = 4.0

x_offset_ext = 0
y_offset_ext = 0

nprocx_coarse = 24
nprocy_coarse = 24

grid_params_external = {
    'xsize': 102400,
    'ysize': 102400,
    'itot': 256,
    'jtot': 256,
    'kmax': 128,
    'dz0': 20,         #if resolution is 400m we can make this also 100m or so??
    'alpha': 0.012,
    'southwest_x': 0,  # Placeholder for coordinates of the southwest corner
    'southwest_y': 0    # Placeholder for coordinates of the southwest corner
}

if NESTING:

    #Internal (nesting) domain
    #Note: resolution is advised to be no smaller than 4 times higher than coarse one
    #########################################
    #West lower corner of the domain in approx. lat and lon:
    sw_lat_nested = 52.0
    sw_lon_nested = 4.5
    
    nprocx_fine = 24
    nprocy_fine = 24

    grid_params_nested = {
        'xsize': 25600,
        'ysize': 25600,
        'itot': 256,
        'jtot': 256,
        'kmax': 128,
        'dz0': 20,
        'alpha': 0.012,
        'southwest_x': 0,  # Placeholder for coordinates of the southwest corner
        'southwest_y': 0   # Placeholder for coordinates of the southwest corner
    }

    ################################################

#External domain creation:
lcc_start_x_ext, lcc_start_y_ext = transform.latlon_to_xy(sw_lat_ext, sw_lon_ext)
# Round to nearest multiple of 5
lcc_start_x_ext = round(lcc_start_x_ext / 5) * 5
lcc_start_y_ext = round(lcc_start_y_ext / 5) * 5

# Define the target grid's southwest corner in flat coordinates
grid_params_external['southwest_x'] = lcc_start_x_ext+x_offset_ext
grid_params_external['southwest_y'] = lcc_start_y_ext+y_offset_ext

les_grid_ext= GridDales(grid_params_external)

dx_ext = les_grid_ext.dx
dy_ext = les_grid_ext.dy

x0_true_ext=lcc_start_x_ext+les_grid_ext.xt[0]
y0_true_ext=lcc_start_y_ext+les_grid_ext.yt[0]

xmax_true_ext=lcc_start_x_ext+les_grid_ext.xt[-1]
ymax_true_ext=lcc_start_y_ext+les_grid_ext.yt[-1]

lat_min_ext, lon_min_ext = transform.xy_to_latlon(x0_true_ext, y0_true_ext)
lat_max_ext, lon_max_ext = transform.xy_to_latlon(xmax_true_ext, ymax_true_ext)

xt_list_ext=les_grid_ext.xt+lcc_start_x_ext
yt_list_ext=les_grid_ext.yt+lcc_start_y_ext


x_mesh_ext, y_mesh_ext = np.meshgrid(xt_list_ext, yt_list_ext)
lat_ext, lon_ext = transform.xy_to_latlon(x_mesh_ext, y_mesh_ext)

center_lat_ext = 0.5 * (lat_min_ext + lat_max_ext)
center_lon_ext = 0.5 * (lon_min_ext + lon_max_ext)


# Compute created domain size

xsize_created = (xmax_true_ext - x0_true_ext) + dx_ext
ysize_created = (ymax_true_ext - y0_true_ext) + dy_ext

print("external domain:", 
      "lat min/lon min:", lat_min_ext, lon_min_ext, 
      "lat max/lon max:", lat_max_ext, lon_max_ext)

print('external domain: center xlat /xlon:',center_lat_ext, center_lon_ext)

print("external domain: South-west corner in xy coordinates (edge!):", 
      "x0/y0:", lcc_start_x_ext,lcc_start_y_ext)

print("external domain in xy coordinates (cell centers):", 
      "x min/y min:", x0_true_ext, y0_true_ext, 
      "x max/y max:", xmax_true_ext, ymax_true_ext)

print("external domain size:",
      "x =", xsize_created, 
      "y =", ysize_created)

print("external domain horizontal resolution:",
      "dx =", dx_ext, 
      "dy =", dy_ext)


print("external domain sanity check (matches grid_params?):",
      "x:", xsize_created == grid_params_external['xsize'],
      "y:", ysize_created == grid_params_external['ysize'])

print('=====================================================')

if NESTING:
    #Internal (nested) domain creation:

    x0_nested, y0_nested = transform.latlon_to_xy(sw_lat_nested,sw_lon_nested)

    x0_nested = round(x0_nested / 5) * 5
    y0_nested = round(y0_nested / 5) * 5

    # Compute the closest index (offset) in external domain to the nested domain's SW corner
    x_offset_nested_index = np.abs(xt_list_ext - x0_nested).argmin()
    y_offset_nested_index = np.abs(yt_list_ext - y0_nested).argmin()

    x_offset_nested=les_grid_ext.xt[x_offset_nested_index]
    y_offset_nested=les_grid_ext.yt[y_offset_nested_index]

    lat_real_nested, lon_real_nested = transform.xy_to_latlon(xt_list_ext[x_offset_nested_index], yt_list_ext[y_offset_nested_index])

    les_grid_nested = GridDales(grid_params_nested)

    xt_list_nested=les_grid_nested.xt+xt_list_ext[x_offset_nested_index]
    yt_list_nested=les_grid_nested.yt+yt_list_ext[y_offset_nested_index]

    lat_max_real_nested, lon_max_real_nested = transform.xy_to_latlon(xt_list_nested[-1], yt_list_nested[-1])

    x0_true_nested=xt_list_nested[0]
    y0_true_nested=yt_list_nested[0]

    xmax_true_nested=xt_list_nested[-1]
    ymax_true_nested=yt_list_nested[-1]

    dx_nested = les_grid_nested.dx
    dy_nested = les_grid_nested.dy

    # Compute created domain size
    xsize_created_nest = (xmax_true_nested - x0_true_nested) + dx_nested
    ysize_created_nest = (ymax_true_nested - y0_true_nested) + dy_nested

    x_mesh_nested, y_mesh_nested = np.meshgrid(xt_list_nested, yt_list_nested)
    lat_nested, lon_nested = transform.xy_to_latlon(x_mesh_nested, y_mesh_nested)

    center_lat_nested = 0.5 * (lat_real_nested + lat_max_real_nested)
    center_lon_nested = 0.5 * (lon_real_nested + lon_max_real_nested)

    xmax_nested_index = np.abs(xt_list_ext - (x0_true_nested + grid_params_nested['xsize'] - dx_nested/2)).argmin()
    ymax_nested_index = np.abs(yt_list_ext - (y0_true_nested + grid_params_nested['ysize'] - dy_nested/2)).argmin()

    print('Nested domain: lat min/lon min:',lat_real_nested, lon_real_nested)
    print('Nested domain: lat max/lon max:',lat_max_real_nested, lon_max_real_nested)
    print('Nested domain: center xlat /xlon:',center_lat_nested, center_lon_nested)

    print("Nested domain: South-west corner in xy coordinates (edge!):", 
          "x0/y0:", xt_list_ext[x_offset_nested_index],yt_list_ext[y_offset_nested_index] )

    print("Nested domain in xy coordinates (cell centers):", 
          "x min/y min:", x0_true_nested, y0_true_nested, 
          "x max/y max:", xmax_true_nested, ymax_true_nested)

    print("Nested domain size:",
          "x =", xsize_created_nest, 
          "y =", ysize_created_nest)

    print("Nested domain horizontal resolution:",
          "dx =", dx_nested, 
          "dy =", dy_nested)


    print("Nested domain sanity check (matches grid_params?):",
          "x:", xsize_created_nest == grid_params_nested['xsize'],
          "y:", ysize_created_nest == grid_params_nested['ysize'])


    print('Nested domain: x and y offsets',x_offset_nested, y_offset_nested)
    print('Nested domain: x and y offsets indexes',x_offset_nested_index, y_offset_nested_index)
    print('crossplane:',y_offset_nested_index+2, ymax_nested_index+2)
    print('crossortho:',x_offset_nested_index+2, xmax_nested_index+2)


# Function to format and plot labels
def annotate_corners(ax, points, color, label=''):
    for i, (lon, lat) in enumerate(points):
        ax.text(
            lon, lat,
            f"{lat:.2f}, {lon:.2f}",
            transform=ccrs.PlateCarree(),
            fontsize=5,
            fontweight='bold',
            color='black',
            bbox=dict(boxstyle='round,pad=0.2', fc='white', ec=color, lw=2),
            ha='center',
            va='center',
            zorder=20
        )

def place_domain_label(ax, lon, lat, text, color='black'):
    lon_min, lon_max = np.nanmin(lon), np.nanmax(lon)
    lat_min, lat_max = np.nanmin(lat), np.nanmax(lat)

    # Margin = small fraction of domain size
    dlon = lon_max - lon_min
    dlat = lat_max - lat_min

    x = lon_max - 0.02 * dlon   # 2% inside from right
    y = lat_min + 0.02 * dlat   # 2% inside from bottom

    ax.text(
        x, y, text,
        transform=ccrs.PlateCarree(),
        ha="right", va="bottom",
        fontsize=8, weight="bold",
        color=color,
        zorder=20
    )


lat_min=np.min(lat_ext)-0.5
lat_max=np.max(lat_ext)+0.5
lon_min=np.min(lon_ext)-1
lon_max=np.max(lon_ext)+1



# Setup figure
fig = plt.figure(figsize=(8, 16), dpi=300)
ax = fig.add_subplot(
    1, 1, 1,
    projection=ccrs.Orthographic(
        central_longitude=((lon_max - lon_min) / 2) + lon_min,
        central_latitude=((lat_max - lat_min) / 2) + lat_min
    )
)


# Top-left corner
corner1 = (lon_ext[0, 0], lat_ext[0, 0])

# Top-right corner
corner2 = (lon_ext[0, -1], lat_ext[0, -1])

# Bottom-right corner
corner3 = (lon_ext[-1, -1], lat_ext[-1, -1])

# Bottom-left corner
corner4 = (lon_ext[-1, 0], lat_ext[-1, 0])

# Combine them into a polygon (either clockwise or counterclockwise)
points_outer = [corner1, corner2, corner3, corner4]


# Line plot using the corners
lon_ext_plt, lat_ext_plt = zip(*points_outer)
lon_ext_plt += (lon_ext_plt[0],)
lat_ext_plt += (lat_ext_plt[0],)
ax.plot(lon_ext_plt, lat_ext_plt, color='blue', linewidth=2, transform=ccrs.PlateCarree(), zorder=10)

if NESTING:
    # Top-left corner
    corner1_nested = (lon_nested[0, 0], lat_nested[0, 0])

    # Top-right corner
    corner2_nested = (lon_nested[0, -1], lat_nested[0, -1])

    # Bottom-right corner
    corner3_nested = (lon_nested[-1, -1], lat_nested[-1, -1])

    # Bottom-left corner
    corner4_nested = (lon_nested[-1, 0], lat_nested[-1, 0])

    # Combine them into a polygon (either clockwise or counterclockwise)
    points_nested = [corner1_nested, corner2_nested, corner3_nested, corner4_nested]


    # Line plot using the corners
    lon_nested_plt, lat_nested_plt = zip(*points_nested)
    lon_nested_plt += (lon_nested_plt[0],)
    lat_nested_plt += (lat_nested_plt[0],)
    ax.plot(lon_nested_plt, lat_nested_plt , color='red', linewidth=2, transform=ccrs.PlateCarree(), zorder=10)

    # Annotate corners
    annotate_corners(ax, points_nested, color='red', label='Nested')

# Annotate corners
annotate_corners(ax, points_outer, color='blue', label='Outer')

# Add geographic features
ax.add_feature(cfeature.COASTLINE, linewidth=1, edgecolor='black')
ax.add_feature(cfeature.LAKES, linewidth=1, edgecolor='black', facecolor='None')
ax.add_feature(cfeature.BORDERS, linewidth=1, linestyle='--', edgecolor='black')

# Gridlines and labels
gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                  linewidth=0.5, color='gray', alpha=0.5, linestyle='--')
gl.top_labels = False
gl.right_labels = False
gl.left_labels = True
gl.bottom_labels = True
gl.xformatter = LONGITUDE_FORMATTER
gl.yformatter = LATITUDE_FORMATTER
gl.xlabel_style = {'color': 'grey', 'weight': 'bold'}
gl.ylabel_style = {'color': 'grey', 'weight': 'bold'}


# Coordinates for Cabauw and Loobos
cabauw_lat, cabauw_lon = 51.9703, 4.9264
loobos_lat, loobos_lon = 52.166447, 5.74355

# Coordinates for additional RITA-2022 stations
slufter_lat, slufter_lon = 51.9575, 4.0031
geulhaven_lat, geulhaven_lon = 51.9050, 4.3860
de_zweth_lat, de_zweth_lon = 51.9600, 4.3900
westmaas_lat, westmaas_lon = 51.7800, 4.4500

# Plot the locations as black dots
ax.plot(cabauw_lon, cabauw_lat, 'ko', markersize=5, transform=ccrs.PlateCarree(), zorder=15)
ax.plot(loobos_lon, loobos_lat, 'ko', markersize=5, transform=ccrs.PlateCarree(), zorder=15)
ax.plot(slufter_lon, slufter_lat, 'ko', markersize=5, transform=ccrs.PlateCarree(), zorder=15)
ax.plot(geulhaven_lon, geulhaven_lat, 'ko', markersize=5, transform=ccrs.PlateCarree(), zorder=15)
ax.plot(de_zweth_lon, de_zweth_lat, 'ko', markersize=5, transform=ccrs.PlateCarree(), zorder=15)
ax.plot(westmaas_lon, westmaas_lat, 'ko', markersize=5, transform=ccrs.PlateCarree(), zorder=15)

# Add coordinate labels
ax.text(cabauw_lon, cabauw_lat + 0.03, "Cabauw\n51.9703, 4.9264",
        transform=ccrs.PlateCarree(), fontsize=2, ha='center', va='bottom',
        bbox=dict(boxstyle='round,pad=0.2', fc='white', ec='black', lw=1),
        zorder=20)

ax.text(loobos_lon, loobos_lat + 0.03, "Loobos\n52.1664, 5.7436",
        transform=ccrs.PlateCarree(), fontsize=2, ha='center', va='bottom',
        bbox=dict(boxstyle='round,pad=0.2', fc='white', ec='black', lw=1),
        zorder=20)

ax.text(slufter_lon, slufter_lat + 0.03, "Slufter\n51.9575, 4.0031",
        transform=ccrs.PlateCarree(), fontsize=2, ha='center', va='bottom',
        bbox=dict(boxstyle='round,pad=0.2', fc='white', ec='black', lw=1),
        zorder=20)

ax.text(geulhaven_lon, geulhaven_lat + 0.03, "Geulhaven\n51.9050, 4.3860",
        transform=ccrs.PlateCarree(), fontsize=2, ha='center', va='bottom',
        bbox=dict(boxstyle='round,pad=0.2', fc='white', ec='black', lw=1),
        zorder=20)

ax.text(de_zweth_lon, de_zweth_lat + 0.03, "De Zweth\n51.9600, 4.3900",
        transform=ccrs.PlateCarree(), fontsize=2, ha='center', va='bottom',
        bbox=dict(boxstyle='round,pad=0.2', fc='white', ec='black', lw=1),
        zorder=20)

ax.text(westmaas_lon, westmaas_lat + 0.03, "Westmaas\n51.7800, 4.4500",
        transform=ccrs.PlateCarree(), fontsize=2, ha='center', va='bottom',
        bbox=dict(boxstyle='round,pad=0.2', fc='white', ec='black', lw=1),
        zorder=20)

place_domain_label(ax, lon_ext, lat_ext, "DALES outer domain", color='blue')

if NESTING:
    place_domain_label(ax, lon_nested, lat_nested, "DALES nested domain", color='red')

# Set the map extent (you can adjust this as needed)
ax.set_extent([lon_min, lon_max, lat_min, lat_max], crs=ccrs.PlateCarree())

plt.tight_layout()

timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
filename = f"your_dales_sim_domains_{timestamp}.png"

plt.savefig(filename, dpi=300, bbox_inches='tight')

plt.show()

#Save domain settings to jSONs (but if the file doesn't exist yet, load a template):
if not (main_input_generator_folder_path / "input_coarse.json").exists():
    infile = script_folder / "input_coarse.json"
else:
    infile = main_input_generator_folder_path / "input_coarse.json"
with open(infile, "r") as f:
    data = json.load(f)

data["author"] = author
data["coarse"]["source"] = source_meteo
data["coarse"]["inpath"] = str(inpath_coarse)
data["coarse"]["outpath"] = str(outpath_coarse)
  
data["coarse"]["start"] = time_start
data["coarse"]["time0"] = time0
data["coarse"]["end"] = time_end

data["coarse"]["lat_sw"] = sw_lat_ext
data["coarse"]["lon_sw"] = sw_lon_ext

data["coarse"]["nprocx"] = int(nprocx_coarse)
data["coarse"]["nprocy"] = int(nprocy_coarse)

data["coarse"]["xlat"] = round(center_lat_ext, 3)
data["coarse"]["xlon"] = round(center_lon_ext, 3)

data["coarse"]["hybrid_file"] = hybrid_lev_file

data["coarse"]["grid"]["xsize"] = grid_params_external["xsize"]
data["coarse"]["grid"]["ysize"] = grid_params_external["ysize"]

data["coarse"]["grid"]["itot"] = grid_params_external["itot"]
data["coarse"]["grid"]["jtot"] = grid_params_external["jtot"]
data["coarse"]["grid"]["kmax"] = grid_params_external["kmax"]
data["coarse"]["grid"]["dz0"] = grid_params_external["dz0"]
data["coarse"]["grid"]["alpha"] = grid_params_external["alpha"]

data["coarse"]["_tskin"]["ERA5_path"] = str(ERA5_path)
data["coarse"]["data_path"] = str(spatial_data_path)
data["coarse"]["LSM"]["ERA5_path"] = str(ERA5_path)
data["coarse"]["LSM"]["spatial_data_path"] = str(spatial_data_path)
data["coarse"]["LSM"]["lsm_kind"] = lsm_kind
data["coarse"]["LSM"]["lwrite_ags"] = lwrite_ags


if NESTING:
    data["coarse"]["grid"]["x_offset"] = x_offset_nested
    data["coarse"]["grid"]["y_offset"] = y_offset_nested
    data["coarse"]["grid"]["xsize_fine"] = grid_params_nested["xsize"]
    data["coarse"]["grid"]["ysize_fine"] = grid_params_nested["ysize"]
else:
    data["coarse"]["grid"]["x_offset"] = None
    data["coarse"]["grid"]["y_offset"] = None
    data["coarse"]["grid"]["xsize_fine"] = None
    data["coarse"]["grid"]["ysize_fine"] = None

data["coarse"]["synturb"]["dx"] = grid_params_external["xsize"]
data["coarse"]["synturb"]["dy"] = grid_params_external["ysize"]

data["coarse"]["filter"]["time_start"] = time_start
data["coarse"]["filter"]["time_end"] = time_end_filter

#Write BACK to the SAME file
with open(main_input_generator_folder_path / "input_coarse.json", "w") as f:
    json.dump(data, f, indent=4)

if NESTING:
    if not (main_input_generator_folder_path / "input_fine.json").exists():
        infile = script_folder / "input_fine.json"
    else:
        infile = main_input_generator_folder_path / "input_fine.json"
    with open(infile, "r") as f:
        data_nest = json.load(f)
        
    data_nest["author"] = author
    data_nest["fine"]["inpath"] = str(inpath_fine)
    data_nest["fine"]["outpath"] = str(outpath_fine)
  
    data_nest["fine"]["start"] = time_start
    data_nest["fine"]["time0"] = time0
    data_nest["fine"]["end"] = time_end


    data_nest["fine"]["xlat"] = round(center_lat_nested, 3)
    data_nest["fine"]["xlon"] = round(center_lon_nested, 3)

    data_nest["fine"]["x_fine"] = xt_list_ext[x_offset_nested_index]
    data_nest["fine"]["y_fine"] = yt_list_ext[y_offset_nested_index]

    data_nest["fine"]["x_offset"] = x_offset_nested
    data_nest["fine"]["y_offset"] = y_offset_nested

    data_nest["fine"]["dx_coarse"] = dx_ext
    data_nest["fine"]["dy_coarse"] = dy_ext
    
    data_nest["fine"]["nprocx"] = int(nprocx_fine)
    data_nest["fine"]["nprocy"] = int(nprocy_fine)

    data_nest["fine"]["grid"]["xsize"] = grid_params_nested["xsize"]
    data_nest["fine"]["grid"]["ysize"] = grid_params_nested["ysize"]
    data_nest["fine"]["grid"]["itot"] = grid_params_nested["itot"]
    data_nest["fine"]["grid"]["jtot"] = grid_params_nested["jtot"]
    data_nest["fine"]["grid"]["kmax"] = grid_params_nested["kmax"]
    data_nest["fine"]["grid"]["dz0"] = grid_params_nested["dz0"]
    data_nest["fine"]["grid"]["alpha"] = grid_params_nested["alpha"]
    
    data_nest["fine"]["_tskin"]["ERA5_path"] = str(ERA5_path)
    data_nest["fine"]["data_path"] = str(spatial_data_path)
    data_nest["fine"]["LSM"]["ERA5_path"] = str(ERA5_path)
    data_nest["fine"]["LSM"]["spatial_data_path"] = str(spatial_data_path)
    data_nest["fine"]["LSM"]["lsm_kind"] = lsm_kind
    data_nest["fine"]["LSM"]["lwrite_ags"] = lwrite_ags

    data_nest["fine"]["synturb"]["dx"] = grid_params_nested["xsize"]
    data_nest["fine"]["synturb"]["dy"] = grid_params_nested["ysize"]
    
    data_nest["fine"]["filter"]["time_start"] = time_start
    data_nest["fine"]["filter"]["time_end"] = time_end_filter
    
    data_nest["fine"]["filter"]["time_start"] = time_start
    data_nest["fine"]["filter"]["time_end"] = time_end_filter

    # 3. Write BACK to the SAME file
    with open(main_input_generator_folder_path / "input_fine.json", "w") as f:
        json.dump(data_nest, f, indent=4)


print("domain settings in jSONs have been updated")