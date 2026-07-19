# DALES meteorological input and open boundary generation

After the HARMONIE-AROME preprocessing step is completed, the generated RD-projected NetCDF files contain the meteorological variables required for DALES.

The next step is to create:

- DALES initial fields
- meteorological open boundary conditions
- surface fields

This procedure is located in:

```
dales_openBC_setup/scripts/
```

The main configuration file is:

```
DALES_domain_setup.py
```

---

# DALES domain setup

`DALES_domain_setup.py` controls:

- input/output paths,
- coordinate projection,
- simulation period,
- coarse domain settings,
- optional nested fine domain settings,
- surface data,
- vertical grid configuration.

---

## Coordinate system

The DALES domains are defined in the Dutch RD coordinate system.

The projection is specified using:

```python
proj_params = {
    'proj4':
    '+proj=sterea +lat_0=52.15616055555555 '
    '+lon_0=5.38763888888889 +k=0.9999079 '
    '+x_0=155000 +y_0=463000 '
    '+ellps=bessel '
    '+towgs84=565.417,50.3319,465.552,'
    '-0.398957,0.343988,-1.8774,4.0725 '
    '+units=m +no_defs'
}
```

---

# General settings

Important options:

```python
NESTING=True
```

Enables a nested fine domain.

Input meteorology:

```python
source_meteo = "harmonie_rd"
```

Paths:

```python
inpath_coarse = Path("/.../RD/")
outpath_coarse = Path("/.../")
```

For nested simulations:

```python
inpath_fine = outpath_coarse
outpath_fine = Path("/...")
```

Additional required data:

```python
ERA5_path = Path("/.../ERA5_soil/")
spatial_data_path = Path("/.../spatial_data/")
```

---

# Simulation period

The DALES simulation period is defined by:

```python
time_start="2025-01-01T00:00"
time_end="2025-01-02T00:00"
```

The meteorological input must cover the complete DALES simulation period.

---

# Vertical grid

The HARMONIE vertical grid is specified using:

```python
hybrid_lev_file="H43_90lev.txt"
```

The vertical level definition must correspond to the HARMONIE input data.

---

# Domain configuration

The setup supports two domains:

- external/coarse domain,
- internal/fine nested domain.

---

## Coarse domain

The southwest corner is specified approximately in latitude/longitude:

```python
sw_lat_ext = 51.5
sw_lon_ext = 4.0
```

Parallel decomposition:

```python
nprocx_coarse = 24
nprocy_coarse = 24
```

Grid parameters:

```python
grid_params_external = {
    'xsize': 102400,
    'ysize': 102400,
    'itot': 256,
    'jtot': 256,
    'kmax': 128,
    'dz0': 20,
    'alpha': 0.012
}
```

---

## Fine nested domain

If:

```python
NESTING=True
```

the internal domain is created.

Example:

```python
sw_lat_nested = 52.0
sw_lon_nested = 4.5
```

Parallel decomposition:

```python
nprocx_fine = 24
nprocy_fine = 24
```

Grid:

```python
grid_params_nested = {
    'xsize': 25600,
    'ysize': 25600,
    'itot': 256,
    'jtot': 256,
    'kmax': 128,
    'dz0': 20,
    'alpha': 0.012
}
```

The nested domain resolution should normally be at least a factor of four higher than the coarse domain.

---

# Optional DALES namelist settings

After defining the domain, general DALES settings can be modified in:

```
setting_files/
```

Available files:

```
namoptions_coarse.001
namoptions_fine.001
```

These files can be adjusted before running DALES, although many settings can also be modified after input generation.

---

# Running DALES meteorological input generation

The execution script is:

```
dales_openBC_setup/scripts/run_script_input.csh
```

Before running, specify the domain:

```bash
export DOMAIN_NAME='coarse'
```

or

```bash
export DOMAIN_NAME='fine'
```

depending on which input is required.

The output experiment directory can also be changed:

```bash
rm -r /.../test_exp_${DOMAIN_NAME}/
mkdir -p /.../test_exp_${DOMAIN_NAME}/
```

Run:

```bash
./run_script_input.csh
```

The generated files contain:

- DALES initial meteorological fields,
- meteorological open boundary conditions,
- surface input.

---

# Chemistry open boundary generation

Chemistry boundary generation is performed separately:

```
chem_input_generator/
```

# Complete workflow summary

The complete DALES input generation workflow is therefore:

```
HARMONIE-AROME archive
        |
        v
harmonie_preprocess
        |
        v
RD-projected meteorological NetCDF
        |
        v
dales_openBC_setup
        |
        v
DALES meteorological initial fields
and open boundaries
        |
        +----------------+
                         |
                         v
             chem_input_generator
                         |
                         v
             DALES chemistry boundaries
```

After these steps, all required DALES input files are prepared for running a meteorological-chemical LES simulation.