# Chemistry open boundary generation
```

The main configuration file is:

```
config.py
```

Currently supported chemistry sources:

- LOTOS-EUROS
- CAMS (IFS)

---

# Chemistry configuration

Example:

```python
input_type = 'LE'
```

Available options:

```python
input_type = 'LE'
input_type = 'CAMS' 
```

Input directory:

```python
input_dir = '/.../'
```

Input file:

```python
input_file = '....nc'
```

Important:

The LOTOS-EUROS/CAMS/IFS input file must cover the **entire DALES simulation period**.

A single NetCDF file containing the complete simulation period is required.

---

# Chemical tracer variables

Tracer names must be provided in lowercase:

```python
tracer_names = [
    'hno3',
    'n2o5',
    'no',
    'no2',
    'o3',
    'n2',
    'iso',
    'ch2o',
    'co'
]
```

The input NetCDF file must contain:

- tracer variables,
- altitude variable above ground,
- complete time period.

---

# Chemistry units

Required units:

| Variable type | Unit |
|---|---|
| CO2 | ppm |
| Other chemical gases | ppb |
| Aerosols | kg/m3 |

Incorrect units will result in incorrect DALES chemistry forcing.

---

# Running chemistry boundary generation

The main script is:

```
Create_chem_input_DALES.csh
```

Internally it calls:

```
create_chem_input.py
```

which generates DALES-compatible chemistry boundaries.

The chemistry fields are added to the same DALES files created during meteorological preprocessing:

- initial field NetCDF files,
- open boundary NetCDF files.

---