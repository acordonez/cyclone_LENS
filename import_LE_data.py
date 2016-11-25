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
    var: numpy array
        3-D array of variable from 1980-2100 in NH
    """

    if len(num) < 3:
        num = num.zfill(3)

    fdir = ('/glade/p/cesm0005/CESM-CAM5-BGC-LE/ice/proc/tseries/daily/' + 
           varname + '_d/')
    if num == '001':
        t1 = ((1980-1850)*365) - 1
        fname1 = (fdir + 'b.e11.B20TRC5CNBDRD.f09_g16.' +
                 num + '.cice.h1.' + varname  +
                 '_d_nh.18500101-20051231.nc') 
    else:
        t1 = ((1980-1920)*365) - 1
        fname1 = (fdir + 'b.e11.B20TRC5CNBDRD.f09_g16.' 
                 + num + '.cice.h1.' + varname 
                 + '_d_nh.19200101-20051231.nc' )

    fname2 = (fdir + 'b.e11.BRCP85C5CNBDRD.f09_g16.' 
             + num + '.cice.h1.' + varname 
             + '_d_nh.20060101-20801231.nc') 

    fname3 = (fdir + 'b.e11.BRCP85C5CNBDRD.f09_g16.'
             + num + '.cice.h1.' + varname 
             + '_d_nh.20810101-21001231.nc')
   
    data = Dataset(fname1)
    var1 = data.variables[varname+'_d'][t1:,:,:]

    data = Dataset(fname2)
    var2 = data.variables[varname+'_d'][:]

    data = Dataset(fname3)
    var3 = data.variables[varname+'_d'][:(365*20),:,:]

    var = np.concatenate((var1,var2), axis = 0)
    var = np.concatenate((var,var3), axis = 0)
    #recode fill values as nan
    #var[var > 10000] = np.nan
    return var

def read_atm_data(varname,num):
    """read_atm_data(varname,num)
    varname: string
        name of the variable to read from file
    num: string
        ID number of the ensemble member to read
    var: numpy array
        3-D array of variable from 1980-2100 in NH
    """

    if len(num) < 3:
        num = num.zfill(3)

    fdir = ('/glade/p/cesm0005/CESM-CAM5-BGC-LE/atm/proc/tseries/daily/' 
           + varname + '/')
    if num == '001':
        t1 = ((1980-1850)*365) - 1
        fname1 = (fdir + 'b.e11.B20TRC5CNBDRD.f09_g16.' 
                 + num + '.cam.h1.' + varname
                 + '.18500101-20051231.nc')
    else:
        t1 = ((1980-1920)*365) - 1
        fname1 = (fdir + 'b.e11.B20TRC5CNBDRD.f09_g16.' 
                 + num + '.cam.h1.' + varname 
                 + '.19200101-20051231.nc')

    fname2 = (fdir + 'b.e11.BRCP85C5CNBDRD.f09_g16.' 
             + num + '.cam.h1.' + varname 
             + '.20060101-20801231.nc')
    fname3 = (fdir + 'b.e11.BRCP85C5CNBDRD.f09_g16.' 
             + num + '.cam.h1.' + varname 
             + '.20810101-21001231.nc')

    data = Dataset(fname1)
    var1 = data.variables[varname][t1:,:,:]

    data = Dataset(fname2)
    var2 = data.variables[varname][:]

    data = Dataset(fname3)
    var3 = data.variables[varname][:(365*19),:,:]

    var = np.concatenate((var1,var2), axis = 0)
    var = np.concatenate((var,var3),axis = 0)
    return var

def save_ice_vars_as_stereo():
    """save_ice_vars_as_stereo
    Regrids ice data from gx1v6 native grid
    to NH stereo polar grid and saves as .npy
    """
    from grid1togrid2 import *
    varlist = ['aice','daidtd','daidtt','frazil','congel','snoice','melts','meltb','meltt']
    ncfile =  '/glade/p/work/aordonez/cesm_mapping/map_gx1v6NH_TO_stereo25km_blin.161123.nc'
    project_and_save_stereo(varlist,ncfile)

def save_atm_vars_as_stereo():
    """save_atm_vars_as_stereo
    Regrids atmosphere data from fv0.9x1.25 native grid
    to NH stereo polar grid and saves as .npy
    """
    from grid1togrid2 import *
    varlist = ['TS','PSL','TAUX','TAUY']
    ncfile = '/glade/p/work/aordonez/cesm_mapping/map_fv0.9x1.25_TO_stereo25km_blin.161123.nc'
    project_and_save_stereo(varlist,ncfile)

def project_and_save_stereo(varlist,ncfile):
    for varname in varlist:
        var = read_atm_data(varname,'001')
        varproj = grid1togrid2(var,ncfile)
        fname = '/glade/scratch/aordonez/LE_tmp/' + varname + '_stereo_001.npy'
        np.save(fname,varproj)

if __name__ = "__main__":
    save_atm_vars_as_stereo()
