# DALES Input Generator

This repository contains a complete preprocessing workflow to generate input files for **DALES (Dutch Atmospheric Large-Eddy Simulation)**.

The workflow prepares:

- **meteorological forcing**
- **surface fields**
- **meteorological open boundary conditions**
- **chemical tracer boundary conditions**

for DALES simulations using operational weather forecast data and atmospheric chemistry products.

Currently supported input sources:

- **Meteorology:** HARMONIE-AROME forecasts
- **Chemistry:** LOTOS-EUROS and CAMS datasets

The generated files are compatible with DALES simulations including nested domains and coupled meteorology-chemistry experiments.

---

# Workflow overview

The complete workflow consists of two main stages:

```
HARMONIE-AROME forecast archives
              |
              v
   harmonie_preprocess
              |
              v
RD-projected meteorological NetCDF files
              |
              v
      dales_openBC_setup
              |
              v
DALES meteorological + surface input
and open boundary conditions
              |
              v
DALES chemistry boundaries
(LOTOS-EUROS / CAMS)
```

---

# Repository structure

```
DALES_input_generator/

│
├── harmonie_preprocess/
│
│   HARMONIE-AROME preprocessing:
│   - extracts forecast archives
│   - creates continuous forecast composites
│   - converts GRIB files to NetCDF
│   - interpolates meteorological fields
│   - transforms data to the RD coordinate system
│
│
└── dales_openBC_setup/
    scripts/
        Main DALES input generation workflow
    cases/
        Example configurations
    simulations/
        Example DALES experiments
```

---

# 1. HARMONIE-AROME preprocessing

The folder

```
harmonie_preprocess/
```

contains the tools required to prepare HARMONIE-AROME meteorological input.

The main tasks are:

1. Reading archived HARMONIE-AROME forecast files
2. Creating a continuous forecast composite suitable for DALES
3. Selecting and converting required meteorological variables
4. Creating NetCDF files
5. Reprojecting the data to the Dutch RD coordinate system

The output of this stage is a set of RD-projected meteorological NetCDF files that are used as input for DALES initialization and boundary generation.

For details, see:

```
harmonie_preprocess/README.md
```

---

# 2. DALES input generation

The folder

```
dales_openBC_setup/
```

contains the main DALES input generation workflow.

It creates:

- DALES meteorological initial conditions
- surface fields
- lateral open boundary conditions
- nested-domain input
- chemistry boundary conditions

Supported chemistry sources:

- LOTOS-EUROS
- CAMS

The generated NetCDF files are directly used for DALES simulations.

For details, see:

```
dales_openBC_setup/scripts/README.md
```

---

# Required datasets

Depending on the simulation setup, the following datasets are required:

## Meteorology

- HARMONIE-AROME forecast archives

## Surface information

- ERA5 soil data
- land-surface and spatial datasets

## Chemistry (optional)

One of:

- LOTOS-EUROS output
- CAMS atmospheric composition data

Chemistry input files must contain the required tracer variables and cover the complete DALES simulation period.

---

# Typical usage workflow

A typical simulation preparation consists of:

### Step 1 — Prepare meteorological input

Run the HARMONIE preprocessing:

```
harmonie_preprocess/
```

This produces RD-projected meteorological NetCDF files.

---

### Step 2 — Generate DALES input

Run the DALES input generator:

```
dales_openBC_setup/
```

This creates:

- initial fields,
- surface fields,
- meteorological open boundaries.

---

### Step 3 — Add chemistry (optional)

If a chemistry simulation is required:

```
dales_openBC_setup/chem_input_generator/
```

creates chemistry initial and boundary conditions from LOTOS-EUROS or CAMS/IFS data.

---

# Notes

- The workflow is currently designed for HARMONIE-AROME meteorological forcing.
- Other meteorological models are not supported without additional preprocessing.
- Chemistry input is currently limited to LOTOS-EUROS and CAMS/IFS formats.
- All input datasets must cover the complete DALES simulation period.

---
# Contact

For questions, improvements, or bug reports, please contact the repository maintainers.