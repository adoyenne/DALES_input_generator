#Setting file to create emissions
#Note: only LOTOS-EUROS and CAMS/IFS input is accepted for now!
#If you run DALES simulation with several days period, LE/CAMS/IFS input file should cover the entire period!

input_type = 'LE' #'LE', 'CAMS' or 'IFS' only!
input_dir = '/.../'
input_LE = '....nc'

#units of input must be ppm for co2, ppb for other chemical tracer gases, and kg/m3 for aerosols!
tracer_names = ['hno3','n2o5','no','no2','o3','n2','iso','ch2o','co']

