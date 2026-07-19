#!/bin/bash

### Setup on Snellius 2022

#SBATCH -n 1
#SBATCH --partition=rome
#SBATCH -t 06:00:00

# Other useful SBATCH options
# #SBATCH --ntasks-per-node=16

#SYST=gnu-fast
############################################################

# Input configuration
ROOT_DIR='' #Root dir for the whole program and where all gribs/nc files will be stored

export STARTDATE='20250101'
export ENDDATE="20250102"

export HA_CYCLE='HA43'
export HA_CONFIG='N20'

#Current usage of HARMONIE-AROME forecast in DALES:
#Unpack HARMONIE-AROME grib files, rolling short-range forecast composite from the 60h HA forecast:
# First cycle: 0–27h (DALES spin-up OK)
# Other cycles: take 3–27 (25 hours)

# Folder with unarchived grib HA files:
export SRC_DIR_gribs="${ROOT_DIR}/harmonie-grib"

#Composed composite folder:
export PATH_GRIB_COMP="${SRC_DIR_gribs}/composed_${STARTDATE}0000"

#Pattern of HA archive name:
export HA_archive_pattern='${HA_CYCLE}_N20_%s.tar'

#Pattern of HA grib file name:
export HA_grib_file_pattern='fc%s+%sCONTROL_GB_UWCW01_N20'

source compose_forecast_N20.sh

##########################################

#Export folder for nc-files derived from HA grib files: 
export PATH_NC_COMP="${ROOT_DIR}/harm_ncdfs/Jan_01_2025/"

if [ ! -d "$PATH_NC_COMP" ]; then
    mkdir -p "$PATH_NC_COMP"
fi

export LENGTH_COMP=30  #the length of the nc- composite (specify here the length of your DALES simulation in hours +1 hour)

#Convert grib to nc on lat lon grid:
python3 grib_to_nc.py  | tee -a log_grib_to_nc.txt
#Convert nc on lat lon grid to nc on RD xy grid:
python3 harmonie_latlon_to_RD.py  | tee -a log_latlon_to_RD.txt

echo 'HARMONIE ncdf composites have been prepared!'