"""
import_LE_data.py

Create timeseries for selected variables and time periods 
in the CESM LENS. Functions for both the Arctic and Antarctic,
and monthly or daily time series. 

Ana Ordonez 10/2018
"""

import numpy as np
from netCDF4 import Dataset
from grid1togrid2 import grid1togrid2
import os.path

"""
Northern Hemisphere Daily
"""

def read_atm_data_stereo(varname,num):
    """read_atm_data_stereo()
    load previously generated numpy array from scratch space
    data is NH on stereo grid
    """
    if len(num) < 3:
        num = num.zfill(3)
    fname = '/glade/scratch/aordonez/' + varname + '_' + num + '_proj_test_global_blin.npy'
    data =  np.load(fname)
    return data

def read_ice_data_stereo(varname,num):
    if len(num) < 3:
        num = num.zfill(3)
    fname = '/glade/scratch/aordonez/' + varname + '_' + num + '_proj_blin.npy'
    data =  np.load(fname)
    return data

def read_ocn_data_stereo(varname,num):
    if len(num) < 3:
        num = num.zfill(3)
    fname = '/glade/scratch/aordonez/' + varname + '_' + num + '_proj_blin.npy'
    data =  np.load(fname)
    return data

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
    fname4 = (fdir + 'b.e11.BRCP85C5CNBDRD.f09_g16.'
             + num + '.cice.h1.' + varname 
             + '_d_nh.20060101-21001231.nc')
    
    data = Dataset(fname1)
    var1 = data.variables[varname+'_d'][t1:,:,:]
    if os.path.isfile(fname2):
        data = Dataset(fname2)
        var2 = data.variables[varname+'_d'][:]
        data = Dataset(fname3)
        var3 = data.variables[varname+'_d'][:(365*40),:,:]

        var = np.concatenate((var1,var2), axis = 0)
        var = np.concatenate((var,var3), axis = 0)
    elif os.path.isfile(fname4):
        data = Dataset(fname4)
        var2 = data.variables[varname + '_d'][:(94*365),:,:]
        var = np.concatenate((var1,var2),axis = 0)

    return var

def read_ice_data_hr6(varname,num):
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

    fdir = ('/glade/p/cesm0005/CESM-CAM5-BGC-LE/ice/proc/tseries/hourly6/' + 
           varname + '_h/')

    t1 = ((1980-1850)*365) - 1
    fname1 = (fdir + 'b.e11.B20TRC5CNBDRD.f09_g16.' +
             num + '.cice.h2_06h.' + varname  +
             '_h_nh.1990010100Z-2005123118Z.nc') 
    fname2 = (fdir + 'b.e11.BRCP85C5CNBDRD.f09_g16.' 
             + num + '.cice.h2_06h.' + varname 
             + '_h_nh.2026010100Z-2035123118Z.nc') 
    fname3 = (fdir + 'b.e11.BRCP85C5CNBDRD.f09_g16.'
             + num + '.cice.h2_06h.' + varname 
             + '_h_nh.2071010100Z-2080123118Z.nc')
    
    data = Dataset(fname1)
    var1 = data.variables[varname+'_h'][:]
    data = Dataset(fname2)
    var2 = data.variables[varname+'_h'][:]
    data = Dataset(fname3)
    var3 = data.variables[varname+'_h'][:]

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
    fname4 = (fdir + 'b.e11.BRCP85C5CNBDRD.f09_g16.' 
             + num + '.cam.h1.' + varname 
             + '.20060101-21001231.nc')
    data = Dataset(fname1)
    var1 = data.variables[varname][t1:,:,:]
    if os.path.isfile(fname2):
        data = Dataset(fname2)
        var2 = data.variables[varname][:]
        data = Dataset(fname3)
        var3 = data.variables[varname][:(365*40),:,:]
        var = np.concatenate((var1,var2), axis = 0)
        var = np.concatenate((var,var3),axis = 0)
    elif os.path.isfile(fname4):
        data = Dataset(fname4)
        var2 = data.variables[varname][:(94*365),:,:]
        var = np.concatenate((var1,var2), axis = 0)
    return var

def read_ocn_data(varname,num):
    """read_ocn_data(varname,num)
    varname: string
        name of the variable to read from file
    num: strings
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
    if (num == '034') | (num == '035'):
        t1 = ((1980-1920)*365) - 1
        fname1 = (fdir + 'b.e11.B20TRC5CNBDRD.f09_g16.' 
                 + num + '.pop.h.nday1.' + varname 
                 + '.19200101-20051231.nc')
        fdir1 = '/glade/scratch/aordonez/'
        fname2 = (fdir1 + 'b.e11.BRCP85C5CNBDRD.f09_g16.' 
                 + num + '.pop.h.nday1.' + varname 
                 + '.20060101-21001231.nc')
        data = Dataset(fname2)
        var2 = data.variables[varname][:,280:384,:]
    else:
        fname2 = (fdir + 'b.e11.BRCP85C5CNBDRD.f09_g16.' 
                 + num + '.pop.h.nday1.' + varname 
                 + '.20060102-20801231.nc')
        fname3 = (fdir + 'b.e11.BRCP85C5CNBDRD.f09_g16.' 
                 + num + '.pop.h.nday1.' + varname 
                 + '.20810101-21001231.nc')
        data = Dataset(fname2)
        var2 = data.variables[varname][:,280:384,:]
        data = Dataset(fname3)
        var3 = data.variables[varname][:(365*20),280:384,:] 
        var2 = np.concatenate((var2,var3),axis = 0)

    data = Dataset(fname1)
    var1 = data.variables[varname][t1:,280:384,:]

    var = np.concatenate((var1,var2), axis = 0)
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

def read_ocean_depth():
    ncfile = '/glade/p/cesm0005/CESM-CAM5-BGC-LE/ocn/proc/tseries/monthly/TEMP/b.e11.BRCP85C5CNBDRD.f09_g16.101.pop.h.TEMP.200601-210012.nc'
    data = Dataset(ncfile)
    z_t = data.variables['z_t'][:]
    return z_t

def read_area_ice():
    ncfile = '/glade/p/cesm0005/CESM-CAM5-BGC-LE/ice/proc/tseries/daily/aice_d/b.e11.B20TRC5CNBDRD.f09_g16.002.cice.h1.aice_d_nh.19200101-20051231.nc'
    data = Dataset(ncfile)
    area = data.variables['tarea'][:]
    area = area.data
    area[area > 1e20] = np.nan
    return area

def read_region_mask_arctic():
    """Returns the mask of regions and list of region names
    for the Arctic from the 2016 Sea ice outlook region mask
    on the EASE grid.
    """
    ncfile = '/glade/u/home/aordonez/sio_2016_mask.nc'
    data = Dataset(ncfile)
    mask = data.variables['mask'][:]
    names = data.variables['region_names'][:]
    return mask, names

def read_area_stereo():
    fname = '/glade/p/work/aordonez/psn25area_v3.dat'
    area = np.fromfile(fname,dtype = np.int32) / 1000
    area = np.reshape(area,(448,304))
    area = np.transpose(area)
    return area

def get_atm_var_as_stereo(var):
    ncfile = '/glade/p/work/aordonez/cesm_mapping/map_fv0.9x1.25_TO_stereo25km_blin.170209.nc'
    stereovar = grid1togrid2(var,ncfile)
    #if len(stereovar.shape) > 2:
        #stereovar = np.transpose(stereovar,(2,0,1))
    return stereovar

def save_ice_vars_as_stereo(num):
    """save_ice_vars_as_stereo
    Regrids ice data from gx1v6 native grid
    to NH stereo polar grid and saves as .npy
    """
    aice = np.zeros((1,1))
    daidtd = np.zeros((1,1))
    daidtt = np.zeros((1,1))
    dvidtd = np.zeros((1,1))
    dvidtt = np.zeros((1,1))
    frazil = np.zeros((1,1))
    congel = np.zeros((1,1)) 
    snoice = np.zeros((1,1))
    melts = np.zeros((1,1))
    meltb = np.zeros((1,1))
    meltt = np.zeros((1,1))
    #a1 = np.zeros((1,1))
    #a2 = np.zeros((1,1))
    #a3 = np.zeros((1,1))
    #a4 = np.zeros((1,1))
    #a5 = np.zeros((1,1))
    hi = np.zeros((1,1))
    #varlist = {'frazil':frazil,'snoice':snoice,'hi':hi}
    #varlist = {'aice':aice,'daidtd':daidtd,'daidtt':daidtt,'dvidtd':dvidtd,'dvidtt':dvidtt,'hi':hi}
    #           'frazil':frazil,'congel':congel,'snoice':snoice,
    #           'melts':melts,'meltb':meltb,'meltt':meltt}
    varlist = {'meltb':meltb,'meltt':meltt,'melts':melts,'frazil':frazil,'congel':congel,'snoice':snoice}
    #varlist = {'aicen001':a1,'aicen002':a2,'aicen003':a3,'aicen004':a4,'aicen005':a5}
    #varlist = {'daidtd':daidtd,'daidtt':daidtt,'meltb':meltb,'dvidtd':dvidtd,'dvidtt':dvidtt}
    #varlist = {'aice':aice}
    ncfile =  '/glade/p/work/aordonez/cesm_mapping/map_gx1v6_TO_stereo25km_blin.170210.nc'
    for varname in varlist:
        print("reading data from file")
        print(varname)
        var = read_ice_data(varname,num)
        print("regridding to stereo")
        tmp= grid1togrid2(var,ncfile)
        print("reshaping")
        tmp = np.transpose(tmp,(2,0,1))
        varlist[varname] = tmp
        np.save('/glade/scratch/aordonez/'+varname+'_'+num+'_proj_blin.npy',tmp)
    return varlist

def save_ocn_vars_as_stereo(num):
    """save_ocn_vars_as_stereo
    Regrids ocn data from gx1v6 native grid
    to NH stereo polar grid and saves as .npy
    """
    SST = np.zeros((1,1))
    HBLT_2 = np.zeros((1,1))
    #varlist = {'HBLT_2':HBLT_2}
    varlist = {'SST':SST}
    ncfile =  '/glade/p/work/aordonez/cesm_mapping/map_gx1v6_TO_stereo25km_blin.170210.nc'
    for varname in varlist:
        print("reading data from file")
        var = read_ocn_data(varname,num)
        print("regridding to stereo")
        tmp= grid1togrid2(var,ncfile)
        print("reshaping")
        #tmp = np.transpose(tmp,(2,0,1))
        varlist[varname] = tmp
        np.save('/glade/scratch/aordonez/'+varname+'_'+num+'_proj_blin.npy',tmp)
    return varlist


def save_atm_vars_as_stereo(num):
    """save_atm_vars_as_stereo
    Regrids atmosphere data from fv0.9x1.25 native grid
    to NH stereo polar grid and saves as .npy
    """
    TS = np.zeros((1,1))
    TREFHT = np.zeros((1,1))
    PSL = np.zeros((1,1))
    TAUX = np.zeros((1,1,))
    TAUY = np.zeros((1,1))
    ICEFRAC = np.zeros((1,1))
    SHFLX = np.zeros((1,1))
    LHFLX = np.zeros((1,1))
    FSNS = np.zeros((1,1))
    FLNS = np.zeros((1,1))
    TMQ = np.zeros((1,1))
    #varlist = {'PSL':PSL,'TS':TS,'PSL':PSL,'TAUX':TAUX,'TAUY':TAUY,'ICEFRAC':ICEFRAC}
    #varlist = {'PSL':PSL,'ICEFRAC':ICEFRAC}
    #varlist = {'TAUX':TAUX,'TAUY':TAUY}
    #varlist = {'SHFLX':SHFLX,'LHFLX':LHFLX,'FSNS':FSNS,'FLNS':FLNS}
    varlist = {'TMQ':TMQ}
    #ncfile = '/glade/p/work/aordonez/cesm_mapping/map_fv0.9x1.25_TO_stereo25km_blin.170209.nc'
    #attempt to regrid with land:
    ncfile = '/glade/p/work/aordonez/cesm_mapping/map_fv0.9x1.25_TO_stereo25km_blin.180405.nc'
    for varname in varlist:
        print("reading data from file")
        var = read_atm_data(varname,num)  
        print("regridding to stereo")
        tmp = grid1togrid2(var,ncfile)
        print("reshaping")
        varlist[varname] = tmp
        np.save('/glade/scratch/aordonez/' + varname + '_' + num +'_proj_test_global_blin_lnd.npy',tmp)
    return varlist

def save_psl_laplacian_as_stereo(num):
    lapl_stereo = np.zeros((1,1))
    ncfile = '/glade/p/work/aordonez/cesm_mapping/map_fv0.9x1.25_TO_stereo25km_blin.180405.nc'
    print("reading data from file")
    fname = '/glade/scratch/aordonez/PSL_laplacian_'+str(num)+'_1980-2079.nc'
    dataset = Dataset(fname)
    lapl_fv = dataset.variables['laplacian'][:]
    print("regridding to stereo")
    lapl_stereo = grid1togrid2(lapl_fv,ncfile)
    np.save('/glade/scratch/aordonez/PSL_laplacian_' + num +'_proj_blin_lnd.npy',lapl_stereo)
    return lapl_stereo

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
    fname4 = (fdir + 'b.e11.BRCP85C5CNBDRD.f09_g16.'
             + num + '.cice.h1.' + varname 
             + '_d_sh.20060101-21001231.nc')
    
    data = Dataset(fname1)
    var1 = data.variables[varname+'_d'][t1:,:,:]
    if os.path.isfile(fname2):
        data = Dataset(fname2)
        var2 = data.variables[varname+'_d'][:]
        data = Dataset(fname3)
        var3 = data.variables[varname+'_d'][:(365*20),:,:]

        var = np.concatenate((var1,var2), axis = 0)
        var = np.concatenate((var,var3), axis = 0)
    elif os.path.isfile(fname4):
        data = Dataset(fname4)
        var2 = data.variables[varname + '_d'][:(94*365),:,:]
        var = np.concatenate((var1,var2),axis = 0)

    return var

def read_ice_data_hr6_SH(varname,num):
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

    fdir = ('/glade/p/cesm0005/CESM-CAM5-BGC-LE/ice/proc/tseries/hourly6/' + 
           varname + '_h/')

    fname1 = (fdir + 'b.e11.B20TRC5CNBDRD.f09_g16.' +
             num + '.cice.h2_06h.' + varname  +
             '_h_sh.1990010100Z-2005123118Z.nc') 
    fname2 = (fdir + 'b.e11.BRCP85C5CNBDRD.f09_g16.' 
             + num + '.cice.h2_06h.' + varname 
             + '_h_sh.2026010100Z-2035123118Z.nc') 
    fname3 = (fdir + 'b.e11.BRCP85C5CNBDRD.f09_g16.'
             + num + '.cice.h2_06h.' + varname 
             + '_h_sh.2071010100Z-2080123118Z.nc')
    
    data = Dataset(fname1)
    var1 = data.variables[varname+'_h'][0:23354,:,:]
    data = Dataset(fname2)
    var2 = data.variables[varname+'_h'][0:14595,:,:]
    data = Dataset(fname3)
    var3 = data.variables[varname+'_h'][0:14595,:,:] # leaving out last half day

    var = np.concatenate((var1,var2), axis = 0)
    var = np.concatenate((var,var3), axis = 0)

    # convert to daily mean
    time,x,y = var.shape
    var = np.reshape(var,(4,time/4,x,y))
    var = np.nanmean(var,axis = 0)

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
    if (num == '034') | (num == '035'):
        t1 = ((1980-1920)*365) - 1
        fname1 = (fdir + 'b.e11.B20TRC5CNBDRD.f09_g16.' 
                 + num + '.pop.h.nday1.' + varname 
                 + '.19200101-20051231.nc')
        fdir1 = '/glade/scratch/aordonez/'
        fname2 = (fdir1 + 'b.e11.BRCP85C5CNBDRD.f09_g16.' 
                 + num + '.pop.h.nday1.' + varname 
                 + '.20060101-21001231.nc')
        data = Dataset(fname2)
        var2 = data.variables[varname][:,0:76,:]
    else:
        fname2 = (fdir + 'b.e11.BRCP85C5CNBDRD.f09_g16.' 
                 + num + '.pop.h.nday1.' + varname 
                 + '.20060102-20801231.nc')
        fname3 = (fdir + 'b.e11.BRCP85C5CNBDRD.f09_g16.' 
                 + num + '.pop.h.nday1.' + varname 
                 + '.20810101-21001231.nc')
        data = Dataset(fname2)
        var2 = data.variables[varname][:,0:76,:]
        data = Dataset(fname3)
        var3 = data.variables[varname][:(365*20),0:76,:] 
        var2 = np.concatenate((var2,var3),axis = 0)

    data = Dataset(fname1)
    var1 = data.variables[varname][t1:,0:76,:]

    var = np.concatenate((var1,var2), axis = 0)
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

def read_native_area_ice_SH():
    fname = '/glade/p/cesm0005/CESM-CAM5-BGC-LE/ice/proc/tseries/daily/aice_d/b.e11.B20TRC5CNBDRD.f09_g16.002.cice.h1.aice_d_sh.19200101-20051231.nc'
    data = Dataset(fname)
    area = data.variables['tarea'][:]
    area = area.data
    area[area > 1e20] = np.nan
    return area

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
    varlist = {           'frazil':frazil,'congel':congel,'snoice':snoice,
               'melts':melts,'meltb':meltb,'meltt':meltt,'sst':SST}
    #varlist = {'aice':aice,'daidtd':daidtd,'daidtt':daidtt,'meltb':meltb}
    varlist = {'SST':SST}
    ncfile =  '/glade/p/work/aordonez/cesm_mapping/map_gx1v6SH_TO_SHstereo25km_blin.161213.nc' 
    for varname in varlist:
        print("reading data from file")
        var = read_ice_data_SH(varname,num)
        x,y,z = var.shape
        var = np.concatenate((var[:,:,z-1:z],var),axis = 2)
        var = np.concatenate((var,var[:,:,1:2]),axis = 2)
        print("regridding to stereo")
        tmp= grid1togrid2(var,ncfile)
        print("reshaping")
        #tmp = np.transpose(tmp,(2,0,1))
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
        print("reading data from file")
        var = read_ocn_data(varname,'001')
        print("regridding to stereo")
        tmp= grid1togrid2(var,ncfile)
        print("reshaping")
        #tmp = np.transpose(tmp,(2,0,1))
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
        print("reading data from file")
        var = read_atm_data(varname,'001')
        print("regridding to stereo")
        tmp = grid1togrid2(var,ncfile)
        print("reshaping")
        #tmp = np.transpose(tmp,(2,0,1))
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
    fname4 = (fdir + 'b.e11.BRCP85C5CNBDRD.f09_g16.'
             + num + '.cice.h.' + varname 
             + '_nh.200601-210012.nc')
    
    data = Dataset(fname1)
    var1 = data.variables[varname][t1:,:,:]
    if os.path.isfile(fname2):
        data = Dataset(fname2)
        var2 = data.variables[varname][:]
        data = Dataset(fname3)
        var3 = data.variables[varname][:(12*20),:,:]

        var = np.concatenate((var1,var2), axis = 0)
        var = np.concatenate((var,var3), axis = 0)
    elif os.path.isfile(fname4):
        data = Dataset(fname4)
        var2 = data.variables[varname][:(74 / 12 * 365),:,:]
        var = np.concatenate((var1,var2),axis = 0)

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

def read_ocn_data_monthly(varname,num):
    """read_ocn_data_monthly(varname,num)
    varname: string
        name of the variable to read from file
    num: string
        ID number of the ensemble member to read
    var: numpy array
        3-D array of variable from 1980-2100 in NH

    if variable has depth dimensions, only the first 10 levels 
    are returned
    """

    if len(num) < 3:
        num = num.zfill(3)

    fdir = ('/glade/p/cesm0005/CESM-CAM5-BGC-LE/ocn/proc/tseries/monthly/' 
           + varname + '/')
    if num == '001':
        t1 = ((1980-1850)*12) - 1
        fname1 = (fdir + 'b.e11.B20TRC5CNBDRD.f09_g16.' 
                 + num + '.pop.h.' + varname
                 + '.185001-200512.nc')
    else:
        t1 = ((1980-1920)*12) - 1
        fname1 = (fdir + 'b.e11.B20TRC5CNBDRD.f09_g16.' 
                 + num + '.pop.h.' + varname 
                 + '.192001-200512.nc')

    data = Dataset(fname1)
    dims = data.variables[varname].shape
    if len(dims) == 4:
      var1 = data.variables[varname][t1:,0:10,280:384,:]
    else:
      var1 = data.variables[varname][t1:,280:384,:]

    if (num == '034') | (num == '035'):
        t1 = ((1980-1920)*365) - 1
        fname1 = (fdir + 'b.e11.B20TRC5CNBDRD.f09_g16.' 
                 + num + '.pop.h.' + varname 
                 + '.192001-200512.nc')
        fname2 = (fdir + 'b.e11.BRCP85C5CNBDRD.f09_g16.' 
                 + num + '.pop.h.' + varname 
                 + '.200601-210012.nc')
        data = Dataset(fname2)
        if len(dims) == 4:
            var2 = data.variables[varname][:,0:10,280:384,:]
        else:
            var2 = data.variables[varname][:,280:384,:]
    else:
        fname2 = (fdir + 'b.e11.BRCP85C5CNBDRD.f09_g16.' 
                 + num + '.pop.h.' + varname 
                 + '.200601-208012.nc')
        fname3 = (fdir + 'b.e11.BRCP85C5CNBDRD.f09_g16.' 
                 + num + '.pop.h.' + varname 
                 + '.208101-210012.nc')
        data2 = Dataset(fname2)
        data3 = Dataset(fname3)
        if len(dims) == 4:
            var2 = data2.variables[varname][:,0:10,280:384,:]
            var3 = data3.variables[varname][:(12*20),0:10,280:384,:]
            var2 = np.concatenate((var2,var3),axis = 0)
        else:
            var2 = data2.variables[varname][:,280:384,:]
            var3 = data3.variables[varname][:(12*20),280:384,:]
            var2 = np.concatenate((var2,var3),axis = 0)

    var = np.concatenate((var1,var2), axis = 0)
    return var


def read_atm_hyb_coord_monthly(num):
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
           + 'TS' + '/')
    if num == '001':
        t1 = ((1980-1850)*12) - 1
        fname1 = (fdir + 'b.e11.B20TRC5CNBDRD.f09_g16.' 
                 + num + '.cam.h0.' + 'TS'
                 + '.185001-200512.nc')
    else:
        t1 = ((1980-1920)*12) - 1
        fname1 = (fdir + 'b.e11.B20TRC5CNBDRD.f09_g16.' 
                 + num + '.cam.h0.' + 'TS' 
                 + '.192001-200512.nc')

    data = Dataset(fname1)
    P0 = data.variables['P0'][:]
    hyam = data.variables['hyam'][:]
    hybm = data.variables['hybm'][:]

    return P0, hyam, hybm

"""
SOUTHERN HEMISPHERE MONTHLY
"""
def read_ice_data_monthly_SH(varname,num):
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
                 '_sh.185001-200512.nc') 
    else:
        t1 = ((1980-1920)*12) - 1
        fname1 = (fdir + 'b.e11.B20TRC5CNBDRD.f09_g16.' 
                 + num + '.cice.h.' + varname 
                 + '_sh.192001-200512.nc' )

    fname2 = (fdir + 'b.e11.BRCP85C5CNBDRD.f09_g16.' 
             + num + '.cice.h.' + varname 
             + '_sh.200601-208012.nc') 

    fname3 = (fdir + 'b.e11.BRCP85C5CNBDRD.f09_g16.'
             + num + '.cice.h.' + varname 
             + '_sh.208101-210012.nc')
   
    data = Dataset(fname1)
    var1 = data.variables[varname][t1:,:,:]
    data = Dataset(fname2)
    var2 = data.variables[varname][:]
    data = Dataset(fname3)
    var3 = data.variables[varname][:(12*20),:,:]

    var = np.concatenate((var1,var2), axis = 0)
    var = np.concatenate((var,var3), axis = 0)
    return var
