#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --partition=genoa
#SBATCH --time=72:00:00

# Other useful SBATCH options
# #SBATCH --ntasks-per-node=16

#SYST=gnu-fast

# one-time installation of Python modules
#pip install cdsapi
#pip install numpy xarray netcdf4 matplotlib # some are already present
#pip install numba "dask[complete]" progress
#pip install cartopy
############################################################

export DOMAIN_NAME='coarse' #'fine'

#Update JSON setting scripts before running it:
#python -m DALES_domain_setup

rm -r /.../test_exp_${DOMAIN_NAME}/
mkdir -p /.../test_exp_${DOMAIN_NAME}/
    
rm -rf log_create_input_${DOMAIN_NAME}.txt
    
python3 create_input.py input_${DOMAIN_NAME}.json | tee -a log_create_input_${DOMAIN_NAME}.txt
