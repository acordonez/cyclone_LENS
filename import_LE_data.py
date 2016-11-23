"""
import_LE_data.py

Create timeseries for selected variables and time periods 
in the CESM LENS

For the Arctic

Eras:
1980-1999
2000-2019
2020-2039
2040-2059
2060-2079
2080-2099

"""

import numpy as np
from netCDF4 import Dataset

def read_ice_data(varname,num):
    """read_ice_data(varname,num)
    varname: string
        name of the variable to read from file
    num: string
        ID number of the ensemble member to read
    """
    t1 = ((1980-1920)*365) - 1

    if len(num) < 3:
        num = num.zfill(3)

    fdir = '/glade/p/cesm0005/CESM-CAM5-BGC-LE/ice/proc/tseries/daily/' 
           + varname + '_d/'
    if num == '001':
        fname1 = fdir + 'b.e11.B20TRC5CNBDRD.f09_g16.'
                 + num + '.cice.h1.' + varname 
                 + '_d_nh.18500101-20051231.nc'
    else:
        fname1 = fdir + 'b.e11.B20TRC5CNBDRD.f09_g16.' 
                 + num + '.cice.h1.' + varname 
                 + '_d_nh.19200101-20051231.nc' 
    fname1 = fdir + 'b.e11.BRCP85C5CNBDRD.f09_g16.' 
             + num + '.cice.h1.' + varname 
             + '_d_nh.20060101-20801231.nc' 
    fname3 = fdir + 'b.e11.BRCP85C5CNBDRD.f09_g16.'
             + num + '.cice.h1.' + varname 
             + '_d_nh.20810101-21001231.nc'
    
    data = Dataset(fname1)
    var1 = data.variables[varname][t1:,:,:]

    data = Dataset(fname2)
    var2 = data.variables[varname][:]

    data = Dataset(fname3)
    var3 = data.variables[varname][:(365*20),:,:]

    var = np.concatenate((var1,var2), axis = 0)
    var = np.concatenate((var,var3), axis = 0)
    return var

def read_atm_data(varname,num):
    """read_atm_data(varname,num)
    varname: string
        name of the variable to read from file
    num: string
        ID number of the ensemble member to read
    """
    t1 = ((1980-1920)*365) - 1

    if len(num) < 3:
        num = num.zfill(3)

    fdir = '/glade/p/cesm0005/CESM-CAM5-BGC-LE/atm/proc/tseries/daily/' 
           + varname + '/'
    if num == '001':
        fname1 = fdir + 'b.e11.B20TRC5CNBDRD.f09_g16.' 
                 + num + '.cam.h1.' + varname
                 + '.18500101-20051231.nc'
    else:
        fname1 = fdir + 'b.e11.B20TRC5CNBDRD.f09_g16.' 
                 + num + '.cam.h1.' + varname 
                 + '.19200101-20051231.nc'
    fname2 = fdir + 'b.e11.BRCP85C5CNBDRD.f09_g16.' 
             + num + '.cam.h1.' + varname 
             + '.20060101-20801231.nc'
    fname3 = fdir + 'b.e11.BRCP85C5CNBDRD.f09_g16.' 
             + num + '.cam.h1.' + varname 
             + '.20810101-21001231.nc'

    data = Dataset(fname1)
    var1 = data.variables[varname][t1:,:,:]

    data = Dataset(fname2)
    var2 = data.variables[varname][:]

    data = Dataset(fname3)
    var3 = data.variables[varname][:(365*19),:,:]

    var = np.concatenate((var1,var2), axis = 0)
    var = np.concatenate((var,var3),axis = 0)
    return var

def get_lat_lon():
    """THIS IS A COPY; FIX IT
    """
    name_dict = {"hires":"tx0.1v2", "lowres":"gx1v6"}
    ncf = '/glade/scratch/aordonez/SOM1850_proc/area_' + name_dict[res] + '.nc'
    area_data = Dataset(ncf)
    tarea = area_data.variables['tarea'][0:eq,:]
    return tarea

def do_by_era(var):
    nyrs = [0,20,40,60,80]
    for yr in nyrs:
        start = yr * 365
        end = start + 20 * 365
        vartmp = var[start:end,:,:]

def get_storm_variables():
    PSL = read_atm_data('PSL','001')
    TAUX = read_atm_data('TAUX','001')
    TAUY = read_atm_data('TAUY','001')
    TS = read_atm_data('TS','001')



