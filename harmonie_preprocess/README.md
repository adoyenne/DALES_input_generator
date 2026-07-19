# DALES Input Generator

This repository contains preprocessing tools for generating meteorological and chemical input files for **DALES (Dutch Atmospheric Large-Eddy Simulation)**.
for open boundaries setup.

Currently, the workflow supports:

- **Meteorology:** HARMONIE-AROME operational forecasts
- **Chemistry:** LOTOS-EUROS / CAMS products

The preprocessing converts the original forecast products into DALES-ready NetCDF files on the **Dutch Rijksdriehoekscoördinaten (RD)** projection.

---

# Workflow

The complete HA preprocessing consists of several stages:

```
HARMONIE archive (.tar/.zip)
            │
            ▼
Extract forecast GRIB files
            │
            ▼
Compose continuous forecast
            │
            ▼
Convert GRIB → NetCDF (lat/lon)
            │
            ▼
Reproject NetCDF → RD coordinates
            │
            ▼
DALES meteorological input
```

The chemistry preprocessing follows a similar philosophy but uses LOTOS-EUROS and CAMS datasets.

---

# HARMONIE-AROME preprocessing

The HARMONIE preprocessing is located in

```
harmonie_preprocess/
```

The main entry point is

```
main_program_HRM_grib_nc_convert.csh
```

This script controls the complete preprocessing workflow.

---

# Configuration

The main configuration options are

```bash
ROOT_DIR=''
```

Root directory where all temporary and generated files are stored.

Simulation period

```bash
export STARTDATE='20250101'
export ENDDATE='20250102'
```

Model configuration

```bash
export HA_CYCLE='HA43'
export HA_CONFIG='N20'
```

The archive and GRIB filename patterns are configurable:

```bash
export HA_archive_pattern='${HA_CYCLE}_${HA_CONFIG}_%s.tar'
export HA_grib_file_pattern='fc%s+%sCONTROL_GB_UWCW01_${HA_CONFIG}'
```

---

# Forecast composition

DALES requires one continuous meteorological time series.

Operational HARMONIE forecasts overlap in time because a new forecast cycle is produced every 24 hours.

The preprocessing therefore constructs a continuous forecast composite using the following strategy:

- **First forecast cycle**
  - forecast hours **0–27**

- **Every subsequent cycle**
  - forecast hours **3–27**

This approach

- preserves the complete spin-up period for the first forecast,
- removes overlapping spin-up hours from subsequent forecasts,
- produces one continuous meteorological forcing for DALES.

The composed GRIB files are written to

```bash
${PATH_GRIB_COMP}
```

---

# GRIB → NetCDF conversion

The composite GRIB files are converted to NetCDF on the native latitude/longitude grid using

```bash
python3 grib_to_nc.py
```

The resulting files contain only the variables required by DALES.

---

# Reprojection to RD coordinates

The latitude/longitude NetCDF files are subsequently reprojected to the Dutch RD coordinate system using

```bash
python3 harmonie_latlon_to_RD.py
```

The resulting NetCDF files constitute the meteorological forcing used by the DALES input generator.

---

# Output

The generated NetCDF files are written to

```bash
${PATH_NC_COMP}
```

The number of output hours is controlled by

```bash
export LENGTH_COMP=30
```

where

```
LENGTH_COMP = DALES simulation length + 1 hour
```

For example,

- 24-hour DALES simulation → `LENGTH_COMP=25`
- 48-hour DALES simulation → `LENGTH_COMP=49`

---

# Running the preprocessing

Execute

```bash
cd harmonie_preprocess

csh main_program_HRM_grib_nc_convert.csh
```

The workflow automatically performs:

1. extraction of archived HARMONIE forecasts,
2. construction of the forecast composite,
3. GRIB to NetCDF conversion,
4. reprojection from latitude/longitude to RD coordinates.

After successful completion, the RD-projected NetCDF files are ready for generating DALES meteorological input.