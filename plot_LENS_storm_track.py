"""
plot_LENS_storm_track.py

Makes plots of ensemble mean values for different indicators 
of the storm track. Uses bandpass filtered sea level pressure and 500 mb
heights along with the eady growth rate. Plots compare values 
in different decades.
"""
import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from import_LE_data import *
from cyclone_composite_LENS import *
from scipy.signal import butter, lfilter


def calculate_g(height):
    G = 6.67e-11 #N m**2 / kg**2
    m = 5.98e24 # kg
    r = 6.38e6 + height # m
    g = (G * m) / r**2
    return g

def convert_geopotential_to_meters(Z):
    # height relative to sea level
    # Got this from integrating Z = 1/g0 * integral(gdz) from 0 to z2
    # by hand, so there could be mistakes
    G = 6.67e-11
    m = 5.98e24
    re = 6.38e6
   
    z = -1 * ( (Z*9.81 /(G*m) - 1./re)**(-1) + re)
    return z

def calculate_eady(lat,Z500,T500,TS):
    # eady_growth_rate = 0.3098*g*abs(f)*abs(du/dz)/brunt_vaisala  
    # def from NCL page
    # estimated at 500 mb
    time,rown,coln = TS.shape
    g = 9.81
    w = 2. * np.pi / (60.*60.*24.)
    f = 2. * w * np.sin(lat[:,1:,:] * 2. * np.pi / 180.)
    z500 = convert_geopotential_to_meters(Z500)
    N = np.sqrt(g / TS[:,1:,:] * (TS[:,1:,:] - T500[:,1:,:])/z500[:,1:,:])
    
    dT = TS[:,1:rown,:] - TS[:,0:rown-1,:]
    dy = (lat[:,1:rown,:] - lat[:,0:rown-1,:]) * 110000
    dT_dy = dT / dy
    du_dz = -1 * g/(f*TS[:,1:,:]) * dT_dy #(approximate)
    eady = 0.3098 * g * abs(f) * abs(du_dz) / N
    eady[eady > 1e100] = np.nan
    return eady * (60 * 60 * 24) #per day

def calculate_eady_surface(lat,num):
    # eady_growth_rate = 0.3098*g*abs(f)*abs(du/dz)/brunt_vaisala  
    # def from NCL page
    # estimated at 500 mb
    T,dT_dz,Z = read_T_and_Z(num)
    time,_,rown,coln = T.shape
    g = 9.81
    w = 2. * np.pi / (60.*60.*24.)
    f = 2. * w * np.sin(lat[:,1:,:] * 2. * np.pi / 180.)
    #Z = convert_geopotential_to_meters(dZ)
    #N = np.sqrt(g / T1[:,1:,:] * (T1[:,1:,:] - T2[:,1:,:])/Z[:,1:,:])
    N = np.sqrt(g / T[:,1,:,:] * dT_dZ
    dT = T1[:,1,1:rown,:] - T1[:,1,0:rown-1,:]
    dy = (lat[:,1:rown,:] - lat[:,0:rown-1,:]) * 110000
    dT_dy = dT / dy
    du_dz = -1 * g/(f*T[:,1,1:,:]) * dT_dy #(approximate)
    eady = 0.3098 * g * abs(f) * abs(du_dz) / N
    eady[eady > 1e100] = np.nan
    return eady * (60 * 60 * 24) #per day
  
# filter functions from https://scipy.github.io/old-wiki/pages/Cookbook/ButterworthBandpass
def butter_bandpass(lowcut, highcut, fs, order=5):
    nyq = 0.5 * fs
    low = lowcut / nyq #Wn scaled where nyquist == 1
    high = highcut / nyq
    b, a = butter(order, [low, high], btype='band')
    return b, a
 
def butter_bandpass_filter(data, lowcut, highcut, fs, order=5):
    b, a = butter_bandpass(lowcut, highcut, fs, order=order)
    u,v,w = data.shape
    data = np.reshape(data,(u,v*w))
    y = np.zeros((u,v*w))
    for point in range(0,v*w):
        y[:,point] = lfilter(b, a, data[:,point])
    return np.reshape(y,(u,v,w))

def get_jja_monthly(data):
    nmnths,x,y = data.shape
    nyrs = nmnths / 12.
    data = np.reshape(data,(12,nyrs,x,y))
    data = data[6:9,:,:,:]
    return np.reshape(data,(3*nyrs,x,y))

def read_T_and_Z(num):
    ncfile = '/glade/scratch/aordonez/LENS_T_850_' + str(num) + '_1920-2005.nc'
    data = Dataset(ncfile)
    T = data.variables['T'][:]
    dT = abs(T[2]-T[0])

    Z = data.variables['Z'][:]
    h1 = convert_geopotential_to_meters(Z[0])
    h2 = convert_geopotential_to_meters(Z[2])
    dH = abs(h2-h1)
    dT_dz = dT / dH
    return T,dT_dz,dH


nlens = 30
pstd = np.zeros((nlens,192,288))
pstd2 = np.zeros((nlens,192,288))
icemean = np.zeros((nlens,192,288))
icemean2 = np.zeros((nlens,192,288))
zstd = np.zeros(np.shape(pstd))
zstd2 = np.zeros(np.shape(pstd2))
eady1 = np.zeros((nlens,191,288))
eady2 = np.zeros((nlens,191,288))
lat,lon = read_native_lat_lon_atm()

for num in range(1,nlens+1):
    test_num = str(num)
    psl = read_atm_data('PSL',test_num) / 100.
    icefrac = read_atm_data('ICEFRAC',test_num)
    Z500 = read_atm_data('Z500',test_num)
    T = read_atm_data_monthly('T',test_num)
    T1 = T[:,28,:,:]
    T2 = T[:,29,:,:]
    Z = read_atm_data_monthly('Z3',test_num)
    dZ = Z[:,29,:,:]-Z[:,28,:,:]
    lattile = np.tile(lat,(dZ.shape[0],1,1))

    start = 0
    end = 365*20
    end1 = 12*20

    ptmp = psl[start:end,:,:]
    ptmp_filt = butter_bandpass_filter(ptmp, 1./6., 1./2.5, 1, order=5)
    ptmp_filt = get_jja(ptmp_filt)
    pstd[num-1,:,:,] = np.var(ptmp_filt,axis = 0)
    ztmp = Z500[start:end,:,:]
    ztmp_filt = butter_bandpass_filter(ztmp, 1./6., 1./2.5, 1, order=5)
    ztmp_filt = get_jja(ztmp_filt)
    zstd[num-1,:,:,] = np.var(ztmp_filt,axis = 0)
    icemean[num-1,:,:] = np.nanmean(icefrac[start:end,:,:],axis = 0)
    etmp = calculate_eady_surface(get_jja_monthly(lattile[start:end1,:,:]),
                                  get_jja_monthly(dZ[start:end1,:,:]),
                                  get_jja_monthly(T1[start:end1,:,:]),
                                  get_jja_monthly(T2[start:end1,:,:]))
    etmp = np.nanmean(etmp,axis = 0)
    eady1[num-1,:,:] = etmp

    start = 365*79
    end = 365*100
    start2 = 12*79
    end2 = 12*100

    ptmp2 = psl[start:end,:,:]
    ptmp2_filt = butter_bandpass_filter(ptmp2, 1./6., 1./2.5, 1, order=5)
    ptmp2_filt = get_jja(ptmp2_filt)
    pstd2[num-1,:,:] = np.var(ptmp2_filt, axis = 0)
    ztmp2 = Z500[start:end,:,:]
    ztmp2_filt = butter_bandpass_filter(ztmp2, 1./6., 1./2.5, 1, order=5)
    ztmp2_filt = get_jja(ztmp2_filt)
    zstd2[num-1,:,:,] = np.var(ztmp2_filt,axis = 0)
    icemean2[num-1,:,:] = np.nanmean(icefrac[start:end,:,:],axis = 0)
    etmp = calculate_eady_surface(get_jja_monthly(lattile[start2:end2,:,:]), 
                          get_jja_monthly(dZ[start2:end2,:,:]),
                          get_jja_monthly(T1[start2:end2,:,:]),
                          get_jja_monthly(T2[start2:end2,:,:]))
    etmp = np.nanmean(etmp,axis = 0)
    eady2[num-1,:,:] = etmp

end

pstd = np.nanmean(pstd,axis = 0)
pstd2 = np.nanmean(pstd2,axis = 0)
icemean = np.nanmean(icemean,axis = 0)
icemean2 = np.nanmean(icemean2,axis = 0)
eady1 = np.nanmean(eady1,axis = 0)
eady2 = np.nanmean(eady2,axis = 0)

zstd = np.nanmean(zstd,axis = 0)
zstd2 = np.nanmean(zstd2,axis = 0)
proj = 'npstere'

#----------------------
# PSL
#----------------------
f,axs=plt.subplots(1,3)
map = Basemap(projection = proj, lat_0 = 90, lon_0 = 180, boundinglat = 40, round = True, ax = axs[0])
x,y = map(lon,lat)
m = map.pcolormesh(x,y,pstd,vmin = 0,vmax = 20,cmap = 'YlOrRd')
m3 = map.contour(x,y,icemean,1,colors = 'k',linewidths = 1)
map.drawcoastlines(color = 'white')
axs[0].set_title('1980-1999',fontsize = 8)
cbar = map.colorbar(m,location = 'bottom',pad = '5%')
cbar.ax.tick_params(labelsize=6) 


map = Basemap(projection = proj, lat_0 = 90, lon_0 = 180, boundinglat = 40, round = True, ax = axs[1])
x,y = map(lon,lat)
m=map.pcolormesh(x,y,pstd2,vmin = 0,vmax = 20,cmap = 'YlOrRd')
m3 = map.contour(x,y,icemean2,1,colors = 'r',linewidths = 1)
map.drawcoastlines(color = 'white')
axs[1].set_title('2060-2079',fontsize = 8)
cbar = map.colorbar(m,location = 'bottom',pad = '5%')
cbar.ax.tick_params(labelsize=6) 

map = Basemap(projection = proj, lat_0 = 90, lon_0 = 180, boundinglat = 40, round = True, ax = axs[2])
x,y = map(lon,lat)
m=map.pcolormesh(x,y,(pstd2-pstd)/pstd * 100,vmin = -15,vmax = 15,cmap = 'RdBu')
axs[2].set_title('%difference',fontsize = 8)
m3 = map.contour(x,y,icemean,1,colors = 'k',linewidths = 1)
m3 = map.contour(x,y,icemean2,1,colors = 'r',linewidths = 1)
map.drawcoastlines(color = 'white')
cbar = map.colorbar(m,location = 'bottom',pad = '5%')
cbar.ax.tick_params(labelsize=6) 
f.savefig('LENS_stormtrack_jja_T850.png')

#----------------------
# 500 mb
#----------------------

f,axs=plt.subplots(1,3)
map = Basemap(projection = proj, lat_0 = 90, lon_0 = 180, boundinglat = 40, round = True, ax = axs[0])
x,y = map(lon,lat)
m = map.pcolormesh(x,y,zstd,vmin = 0,vmax = 800,cmap = 'YlOrRd') #,vmin = 30,vmax = 70)
m3 = map.contour(x,y,icemean,1,colors = 'k',linewidths = 1)
map.drawcoastlines(color = 'white')
axs[0].set_title('1980-1999',fontsize = 8)
cbar = map.colorbar(m,location = 'bottom',pad = '5%')
cbar.ax.tick_params(labelsize=6) 


map = Basemap(projection = proj, lat_0 = 90, lon_0 = 180, boundinglat = 40, round = True, ax = axs[1])
x,y = map(lon,lat)
m=map.pcolormesh(x,y,zstd2,vmin = 0,vmax = 800,cmap = 'YlOrRd') #,vmin = 30,vmax = 70)
m3 = map.contour(x,y,icemean2,1,colors = 'r',linewidths = 1)
map.drawcoastlines(color = 'white')
axs[1].set_title('2060-2079',fontsize = 8)
cbar = map.colorbar(m,location = 'bottom',pad = '5%')
cbar.ax.tick_params(labelsize=6) 

map = Basemap(projection = proj, lat_0 = 90, lon_0 = 180, boundinglat = 40, round = True, ax = axs[2])
x,y = map(lon,lat)
m=map.pcolormesh(x,y,(zstd2-zstd)/zstd * 100,vmin = -15,vmax = 15,cmap = 'RdBu_r')
axs[2].set_title('%difference',fontsize = 8)
m3 = map.contour(x,y,icemean,1,colors = 'k',linewidths = 1)
m3 = map.contour(x,y,icemean2,1,colors = 'r',linewidths = 1)
map.drawcoastlines(color = 'white')
cbar = map.colorbar(m,location = 'bottom',pad = '5%')
cbar.ax.tick_params(labelsize=6) 
f.savefig('LENS_stormtrack_Z500_jja_T850.png')


#----------------------
# Eady
#----------------------

f,axs=plt.subplots(1,3)
map = Basemap(projection = proj, lat_0 = 90, lon_0 = 180, boundinglat = 40, round = True, ax = axs[0])
x,y = map(lon,lat)
x1 = x[1:,:]
y1 = y[1:,:]
m = map.pcolormesh(x1,y1,eady1,vmin = 0,vmax = 10,cmap = 'YlOrRd') #,vmin = 30,vmax = 70)
m3 = map.contour(x,y,icemean,1,colors = 'k',linewidths = 1)
map.drawcoastlines(color = 'white')
axs[0].set_title('1980-1999',fontsize = 8)
cbar = map.colorbar(m,location = 'bottom',pad = '5%')
cbar.ax.tick_params(labelsize=6) 


map = Basemap(projection = proj, lat_0 = 90, lon_0 = 180, boundinglat = 40, round = True, ax = axs[1])
x,y = map(lon,lat)
m=map.pcolormesh(x1,y1,eady2,vmin = 0,vmax = 10,cmap = 'YlOrRd') #,vmin = 30,vmax = 70)
m3 = map.contour(x,y,icemean2,1,colors = 'r',linewidths = 1)
map.drawcoastlines(color = 'white')
axs[1].set_title('2060-2079',fontsize = 8)
cbar = map.colorbar(m,location = 'bottom',pad = '5%')
cbar.ax.tick_params(labelsize=6) 

map = Basemap(projection = proj, lat_0 = 90, lon_0 = 180, boundinglat = 40, round = True, ax = axs[2])
x,y = map(lon,lat)
m=map.pcolormesh(x1,y1,(eady2-eady1)/eady1 * 100,vmin = -100,vmax = 100,cmap = 'RdBu_r')
axs[2].set_title('%difference',fontsize = 8)
m3 = map.contour(x,y,icemean,1,colors = 'k',linewidths = 1)
m3 = map.contour(x,y,icemean2,1,colors = 'r',linewidths = 1)
map.drawcoastlines(color = 'white')
cbar = map.colorbar(m,location = 'bottom',pad = '5%')
cbar.ax.tick_params(labelsize=6) 
f.savefig('LENS_stormtrack_eady_jja_monthly.png')
f.show()
