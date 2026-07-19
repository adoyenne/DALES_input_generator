#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=128
#SBATCH --ntasks-per-node=128
#SBATCH --partition=rome
#SBATCH --time=09:00:00

# Other useful SBATCH options
# #SBATCH --ntasks-per-node=16

#SYST=gnu-fast

# one-time installation of Python modules
#pip install cdsapi
#pip install numpy xarray netcdf4 matplotlib # some are already present
#pip install numba "dask[complete]" progress
#pip install cartopy
############################################################

python3 create_chem_input.py ../input_coarse.json | tee -a out_create_chem_input.txt

echo 'DALES chemical boundaries have been prepared and merged to initfield and openboundaries'