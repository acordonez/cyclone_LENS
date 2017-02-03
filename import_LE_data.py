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
from grid1togrid2 import *

"""
Northern Hemisphere Daily
"""

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
    var3 = data.variables[varname][:(365*20),:,:]

    var = np.concatenate((var1,var2), axis = 0)
    var = np.concatenate((var,var3),axis = 0)
    return var

def read_ocn_data(varname,num):
    """read_ocn_data(varname,num)
    varname: string
        name of the variable to read from file
    num: string
        ID number of the ensemble member to read
    var: numpy array
        3-D array of variable from 1980-2100 in NH
    """

    if len(num) < 3:
        num = num.zfill(3)

    fdir = ('/glade/p/cesm0005/CESM-CAM5-BGC-LE/ocn/proc/tseries/daily/' 
           + varname + '/')
    if num == '001':
        t1 = ((1980-1850)*365) - 1
        fname1 = (fdir + 'b.e11.B20TRC5CNBDRD.f09_g16.' 
                 + num + '.pop.h.nday1.' + varname
                 + '.18500102-20051231.nc')
    else:
        t1 = ((1980-1920)*365) - 1
        fname1 = (fdir + 'b.e11.B20TRC5CNBDRD.f09_g16.' 
                 + num + '.pop.h.nday1.' + varname 
                 + '.19200102-20051231.nc')

    fname2 = (fdir + 'b.e11.BRCP85C5CNBDRD.f09_g16.' 
             + num + '.pop.h.nday1.' + varname 
             + '.20060102-20801231.nc')
    fname3 = (fdir + 'b.e11.BRCP85C5CNBDRD.f09_g16.' 
             + num + '.pop.h.nday1.' + varname 
             + '.20810101-21001231.nc')

    data = Dataset(fname1)
    var1 = data.variables[varname][t1:,280:384,:]
    data = Dataset(fname2)
    var2 = data.variables[varname][:,280:384,:]
    data = Dataset(fname3)
    var3 = data.variables[varname][:(365*20),280:384,:]

    var = np.concatenate((var1,var2), axis = 0)
    var = np.concatenate((var,var3),axis = 0)
    return var

def read_native_lat_lon_atm():
    fname = '/glade/p/cesm0005/CESM-CAM5-BGC-LE/atm/proc/tseries/daily/TS/b.e11.B20TRC5CNBDRD.f09_g16.001.cam.h1.TS.18500101-20051231.nc'
    data = Dataset(fname)
    lat = data.variables['lat'][:]
    lon = data.variables['lon'][:]
    lon,lat = np.meshgrid(lon,lat)
    return lat, lon

def read_native_lat_lon_ice():
    fname = '/glade/p/cesm0005/CESM-CAM5-BGC-LE/ice/proc/tseries/daily/aice_d/b.e11.B20TRC5CNBDRD.f09_g16.002.cice.h1.aice_d_nh.19200101-20051231.nc'
    data = Dataset(fname)
    lat = data.variables['TLAT'][:]
    lon = data.variables['TLON'][:]
    lat = lat.data
    lon = lon.data
    lat[lat > 1e20] = np.nan
    lon[lon > 1e20] = np.nan
    return lat, lon
   

def read_stereo_lat_lon():
    """read_stereo_lat_lon
    Returns the latitude and longitude grids for the 
    EASE Northern Hemisphere stereo grid
    """
    ncfile = '/glade/p/work/aordonez/cesm_mapping/stereo_gridinfo.nc'
    data = Dataset(ncfile)
    lat = data.variables['grid_center_lat'][:]
    lon = data.variables['grid_center_lon'][:]
    lat = lat * (180./np.pi)
    lon = lon * (180./np.pi)
    lat = np.reshape(lat,(304,448))
    lon = np.reshape(lon,(304,448))
    return lat,lon

def read_area_ice():
    ncfile = '/glade/p/cesm0005/CESM-CAM5-BGC-LE/ice/proc/tseries/daily/aice_d/b.e11.B20TRC5CNBDRD.f09_g16.002.cice.h1.aice_d_nh.19200101-20051231.nc'
    data = Dataset(ncfile)
    area = data.variables['tarea'][:]
    area = area.data
    area[area > 1e20] = np.nan
    return area

def save_ice_vars_as_stereo(num):
    """save_ice_vars_as_stereo
    Regrids ice data from gx1v6 native grid
    to NH stereo polar grid and saves as .npy
    """
    aice = np.zeros((1,1))
    daidtd = np.zeros((1,1))
    daidtt = np.zeros((1,1))
    frazil = np.zeros((1,1))
    congel = np.zeros((1,1)) 
    snoice = np.zeros((1,1))
    melts = np.zeros((1,1))
    meltb = np.zeros((1,1))
    meltt = np.zeros((1,1))
    sst = np.zeros((1,1))
    #varlist = {'aice':aice,'daidtd':daidtd,'daidtt':daidtt}
    #           'frazil':frazil,'congel':congel,'snoice':snoice,
    #           'melts':melts,'meltb':meltb,'meltt':meltt}
    varlist = {'daidtd':daidtd,'daidtt':daidtt,'aice':aice,'meltb':meltb}
    ncfile =  '/glade/p/work/aordonez/cesm_mapping/map_gx1v6NH_TO_stereo25km_blin.161123.nc'
    for varname in varlist:
        print "reading data from file"
        var = read_ice_data(varname,num)
        var = np.concatenate((var,var[:,:,0:1]),axis = 2)
        print "regridding to stereo"
        tmp= grid1togrid2(var,ncfile)
        print "reshaping"
        tmp = np.transpose(tmp,(2,0,1))
        varlist[varname] = tmp
        np.save('/glade/scratch/aordonez/'+varname+'_'+num+'_proj.npy',tmp)
    return varlist

def save_ocn_vars_as_stereo(num):
    """save_ocn_vars_as_stereo
    Regrids ocn data from gx1v6 native grid
    to NH stereo polar grid and saves as .npy
    """
    SST = np.zeros((1,1))
    varlist = {'SST':SST}
    ncfile =  '/glade/p/work/aordonez/cesm_mapping/map_gx1v6NH_TO_stereo25km_blin.161123.nc'
    for varname in varlist:
        print "reading data from file"
        var = read_ocn_data(varname,num)
        var = np.concatenate((var,var[:,:,0:1]),axis = 2)
        print "regridding to stereo"
        tmp= grid1togrid2(var,ncfile)
        print "reshaping"
        tmp = np.transpose(tmp,(2,0,1))
        varlist[varname] = tmp
        np.save('/glade/scratch/aordonez/'+varname+'_'+num+'_proj.npy',tmp)
    return varlist


def save_atm_vars_as_stereo(num):
    """save_atm_vars_as_stereo
    Regrids atmosphere data from fv0.9x1.25 native grid
    to NH stereo polar grid and saves as .npy
    """
    TS = np.zeros((1,1))
    PSL = np.zeros((1,1))
    TAUX = np.zeros((1,1,))
    TAUY = np.zeros((1,1))
    ICE = np.zeros((1,1))
    varlist = {'TS':TS,'PSL':PSL,'TAUX':TAUX,'TAUY':TAUY,'ICEFRAC':ICE}
    #varlist = {'TAUX':TAUX,'TAUY':TAUY}
    ncfile = '/glade/p/work/aordonez/cesm_mapping/map_fv0.9x1.25_TO_stereo25km_blin.161123.nc'
    for varname in varlist:
        print "reading data from file"
        var = read_atm_data(varname,num)
        var = np.concatenate((var,var[:,:,0:1]),axis = 2)
        print "regridding to stereo"
        tmp = grid1togrid2(var,ncfile)
        print "reshaping"
        tmp = np.transpose(tmp,(2,0,1))
        varlist[varname] = tmp
        np.save('/glade/scratch/aordonez/' + varname + '_' + num +'_proj.npy',tmp)
    return varlist

"""
SOUTHERN HEMISPHERE FUNCTIONS
"""

def read_ice_data_SH(varname,num):
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
                 '_d_sh.18500101-20051231.nc') 
    else:
        t1 = ((1980-1920)*365) - 1
        fname1 = (fdir + 'b.e11.B20TRC5CNBDRD.f09_g16.' 
                 + num + '.cice.h1.' + varname 
                 + '_d_sh.19200101-20051231.nc' )

    fname2 = (fdir + 'b.e11.BRCP85C5CNBDRD.f09_g16.' 
             + num + '.cice.h1.' + varname 
             + '_d_sh.20060101-20801231.nc') 

    fname3 = (fdir + 'b.e11.BRCP85C5CNBDRD.f09_g16.'
             + num + '.cice.h1.' + varname 
             + '_d_sh.20810101-21001231.nc')
   
    data = Dataset(fname1)
    var1 = data.variables[varname+'_d'][t1:,:,:]
    data = Dataset(fname2)
    var2 = data.variables[varname+'_d'][:]
    data = Dataset(fname3)
    var3 = data.variables[varname+'_d'][:(365*20),:,:]

    var = np.concatenate((var1,var2), axis = 0)
    var = np.concatenate((var,var3), axis = 0)
    return var
    
def read_ocn_data_SH(varname,num):
    """read_ocn_data(varname,num)
    varname: string
        name of the variable to read from file
    num: string
        ID number of the ensemble member to read
    var: numpy array
        3-D array of variable from 1980-2100 in NH
    """

    if len(num) < 3:
        num = num.zfill(3)

    fdir = ('/glade/p/cesm0005/CESM-CAM5-BGC-LE/ocn/proc/tseries/daily/' 
           + varname + '/')
    if num == '001':
        t1 = ((1980-1850)*365) - 1
        fname1 = (fdir + 'b.e11.B20TRC5CNBDRD.f09_g16.' 
                 + num + '.pop.h.nday1.' + varname
                 + '.18500102-20051231.nc')
    else:
        t1 = ((1980-1920)*365) - 1
        fname1 = (fdir + 'b.e11.B20TRC5CNBDRD.f09_g16.' 
                 + num + '.pop.h.nday1.' + varname 
                 + '.19200102-20051231.nc')

    fname2 = (fdir + 'b.e11.BRCP85C5CNBDRD.f09_g16.' 
             + num + '.pop.h.nday1.' + varname 
             + '.20060102-20801231.nc')
    fname3 = (fdir + 'b.e11.BRCP85C5CNBDRD.f09_g16.' 
             + num + '.pop.h.nday1.' + varname 
             + '.20810101-21001231.nc')

    data = Dataset(fname1)
    var1 = data.variables[varname][t1:,0:76,:]
    data = Dataset(fname2)
    var2 = data.variables[varname][:,0:76,:]
    data = Dataset(fname3)
    var3 = data.variables[varname][:(365*20),0:76,:]

    var = np.concatenate((var1,var2), axis = 0)
    var = np.concatenate((var,var3),axis = 0)
    return var


def read_native_lat_lon_ice_SH():
    fname = '/glade/p/cesm0005/CESM-CAM5-BGC-LE/ice/proc/tseries/daily/aice_d/b.e11.B20TRC5CNBDRD.f09_g16.002.cice.h1.aice_d_sh.19200101-20051231.nc'
    data = Dataset(fname)
    lat = data.variables['TLAT'][:]
    lon = data.variables['TLON'][:]
    lat = lat.data
    lon = lon.data
    lat[lat > 1e20] = np.nan
    lon[lon > 1e20] = np.nan
    return lat, lon

def save_ice_vars_as_SHstereo(num):
    """save_ice_vars_as_stereo
    Regrids ice data from gx1v6 native grid
    to SH stereo polar grid and saves as .npy
    """
    aice = np.zeros((1,1))
    daidtd = np.zeros((1,1))
    daidtt = np.zeros((1,1))
    frazil = np.zeros((1,1))
    congel = np.zeros((1,1)) 
    snoice = np.zeros((1,1))
    melts = np.zeros((1,1))
    meltb = np.zeros((1,1))
    meltt = np.zeros((1,1))
    SST = np.zeros((1,1))
    #varlist = {'aice':aice,'daidtd':daidtd,'daidtt':daidtt,
    #           'frazil':frazil,'congel':congel,'snoice':snoice,
    #           'melts':melts,'meltb':meltb,'meltt':meltt,'sst':SST}
    varlist = {'aice':aice,'daidtd':daidtd,'daidtt':daidtt,'meltb':meltb}
    varlist = {'SST':SST}
    ncfile =  '/glade/p/work/aordonez/cesm_mapping/map_gx1v6SH_TO_SHstereo25km_blin.161213.nc'
    for varname in varlist:
        print "reading data from file"
        var = read_ice_data_SH(varname,num)
        var = np.concatenate((var,var[:,:,0]),axis = 2)
        print "regridding to stereo"
        tmp= grid1togrid2(var,ncfile)
        print "reshaping"
        tmp = np.transpose(tmp,(2,0,1))
        varlist[varname] = tmp
        np.save('/glade/scratch/aordonez/'+varname+'_' + num + '_SHproj.npy',tmp)
    return varlist

def save_ocn_vars_as_SHstereo():
    """save_ocn_vars_as_stereo
    Regrids ocn data from gx1v6 native grid
    to SH stereo polar grid and saves as .npy
    """
    SST = np.zeros((1,1))
    varlist = {'SST':SST}
    ncfile =  '/glade/p/work/aordonez/cesm_mapping/map_gx1v6SH_TO_SHstereo25km_blin.161213.nc'
    for varname in varlist:
        print "reading data from file"
        var = read_ocn_data(varname,'001')
        print "regridding to stereo"
        tmp= grid1togrid2(var,ncfile)
        print "reshaping"
        tmp = np.transpose(tmp,(2,0,1))
        varlist[varname] = tmp
        np.save('/glade/scratch/aordonez/'+varname+'SHproj.npy',tmp)
    return varlist

def save_atm_vars_as_SHstereo():
    """save_atm_vars_as_stereo
    Regrids atmosphere data from fv0.9x1.25 native grid
    to SH stereo polar grid and saves as .npy
    """
    TS = np.zeros((1,1))
    PSL = np.zeros((1,1))
    UBOT = np.zeros((1,1,))
    VBOT = np.zeros((1,1))
    ICE = np.zeros((1,1))
    #varlist = {'TS':TS,'PSL':PSL,'UBOT':UBOT,'VBOT':VBOT}
    varlist = {'ICEFRAC':ICE}
    ncfile = '/glade/p/work/aordonez/cesm_mapping/map_fv0.9x1.25_TO_SHstereo25km_blin.161213.nc'
    for varname in varlist:
        print "reading data from file"
        var = read_atm_data(varname,'001')
        print "regridding to stereo"
        tmp = grid1togrid2(var,ncfile)
        print "reshaping"
        tmp = np.transpose(tmp,(2,0,1))
        varlist[varname] = tmp
        np.save('/glade/scratch/aordonez/'+varname+'SHproj.npy',tmp)
    return varlist

"""
NORTHERN HEMISPHERE MONTHLY
"""
def read_ice_data_monthly(varname,num):
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

    fdir = ('/glade/p/cesm0005/CESM-CAM5-BGC-LE/ice/proc/tseries/monthly/' + 
           varname + '/')
    if num == '001':
        t1 = ((1980-1850)*12) - 1
        fname1 = (fdir + 'b.e11.B20TRC5CNBDRD.f09_g16.' +
                 num + '.cice.h.' + varname  +
                 '_nh.185001-200512.nc') 
    else:
        t1 = ((1980-1920)*12) - 1
        fname1 = (fdir + 'b.e11.B20TRC5CNBDRD.f09_g16.' 
                 + num + '.cice.h.' + varname 
                 + '_nh.192001-200512.nc' )

    fname2 = (fdir + 'b.e11.BRCP85C5CNBDRD.f09_g16.' 
             + num + '.cice.h.' + varname 
             + '_nh.200601-208012.nc') 

    fname3 = (fdir + 'b.e11.BRCP85C5CNBDRD.f09_g16.'
             + num + '.cice.h.' + varname 
             + '_nh.208101-210012.nc')
   
    data = Dataset(fname1)
    var1 = data.variables[varname][t1:,:,:]
    data = Dataset(fname2)
    var2 = data.variables[varname][:]
    data = Dataset(fname3)
    var3 = data.variables[varname][:(12*20),:,:]

    var = np.concatenate((var1,var2), axis = 0)
    var = np.concatenate((var,var3), axis = 0)
    return var

def read_atm_data_monthly(varname,num):
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

    fdir = ('/glade/p/cesm0005/CESM-CAM5-BGC-LE/atm/proc/tseries/monthly/' 
           + varname + '/')
    if num == '001':
        t1 = ((1980-1850)*12) - 1
        fname1 = (fdir + 'b.e11.B20TRC5CNBDRD.f09_g16.' 
                 + num + '.cam.h0.' + varname
                 + '.185001-200512.nc')
    else:
        t1 = ((1980-1920)*12) - 1
        fname1 = (fdir + 'b.e11.B20TRC5CNBDRD.f09_g16.' 
                 + num + '.cam.h0.' + varname 
                 + '.192001-200512.nc')

    fname2 = (fdir + 'b.e11.BRCP85C5CNBDRD.f09_g16.' 
             + num + '.cam.h0.' + varname 
             + '.200601-208012.nc')
    fname3 = (fdir + 'b.e11.BRCP85C5CNBDRD.f09_g16.' 
             + num + '.cam.h0.' + varname 
             + '.208101-210012.nc')

    data = Dataset(fname1)
    var1 = data.variables[varname][t1:,:,:]
    data = Dataset(fname2)
    var2 = data.variables[varname][:]
    data = Dataset(fname3)
    var3 = data.variables[varname][:(12*20),:,:]

    var = np.concatenate((var1,var2), axis = 0)
    var = np.concatenate((var,var3),axis = 0)
    return var

