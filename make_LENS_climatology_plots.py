from import_LE_data import *
from cyclone_composite_LENS import *
import matplotlib.pyplot as plt


"""
Northern Hemisphere
"""
print 'loading data'
pproj = np.load('/glade/scratch/aordonez/pproj.npy')
iproj = np.load('/glade/scratch/aordonez/iproj.npy')
tproj = np.load('/glade/scratch/aordonez/tproj.npy') 
lat,lon = read_stereo_lat_lon()
years = ['1980-1999','2000-2019','2020-2039',
         '2040-2059','2060-2079']


for n in range(0,5):
    start = 20*365*n 
    end = start + (20*365)
    yrs = years[n]
    # min and max ice extent contours over mean extent pcolor
    # 1980-1999
    icemean = np.nanmean(iproj[start:end,:,:],axis = 0)
    icemax = np.nanmean(iproj[range(start+73,end,365),:,:],axis = 0)
    icemin = np.nanmean(iproj[range(start+257,end,365),:,:],axis = 0)
    icemax = np.select([icemax >= 0.15, icemax < 0.15, np.isnan(icemax)],[1, 0, np.nan])
    icemin = np.select([icemin >= 0.15, icemin < 0.15, np.isnan(icemin)],[1, 0, np.nan])

    #    map = Basemap(projection = proj, lat_0 = -90, lon_0 = 180,boundinglat = -60, round = True, ax = ax2)
#    m22 = map.contourf(x,y, pol_freq_mam, ax = ax2)
    blat = np.arange(40,91,1)
    blon = np.arange(-180,181,1)
    blon,blat = np.meshgrid(blon,blat)
    bground = np.zeros(blon.shape)

    proj = 'npstere'
    f,axs=plt.subplots(1,1)
    map = Basemap(projection = proj, lat_0 = 90, lon_0 = 180, boundinglat = 40, round = True, ax = axs)
    x,y = map(lon,lat)
    bx,by = map(blon,blat)
    m1 = map.pcolormesh(bx,by,bground,vmin = 0,vmax = 1,cmap = 'PuBu_r',edgecolors = 'none')
    m2 = map.pcolormesh(x,y,icemean,vmin = 0,vmax = 1,cmap = 'PuBu_r',edgecolors = 'none')
    m3 = map.contour(x,y,icemax,1,colors = 'b',linewidths = 2)
    m4 = map.contour(x,y,icemin,1,colors = 'r',linewidths = 2)
    map.drawcoastlines(color = 'none')
    map.fillcontinents(color='white')
    axs.set_title(yrs + ' ice extent')
    f.savefig('ice_extent_' + yrs + '.png')

    lows = find_cyclone_center(pproj[start:end,:,:],iproj[start:end,:,:],104000,90000)
    lowsmap = np.nansum(lows,axis = 0)
    f,axs = plt.subplots(1,1)
    map = Basemap(projection = proj, lat_0 = 90, lon_0 = 180, boundinglat = 40, round = True, ax = axs)
    map.pcolormesh(x,y,lowsmap)
    map.drawcoastlines(color = 'none')
    map.fillcontinents(color='white')
    axs.set_title(yrs + ' storm count')
    f.savefig('storm_count_' + yrs + '.png')
    
plt.show()
# mean storm location map


# seasonal storm count bar chart

"""
Southern Hemisphere
"""
print 'loading data'
pproj = np.load('/glade/scratch/aordonez/PSLSHproj.npy')
iproj = np.load('/glade/scratch/aordonez/ICEFRACSHproj.npy')
tproj = np.load('/glade/scratch/aordonez/TSSHproj.npy') 
lat,lon = read_SHstereo_lat_lon()
years = ['1980-1999','2000-2019','2020-2039',
         '2040-2059','2060-2079']
    
for n in range(0,5):
    start = 20*365*n 
    end = start + (20*365)
    yrs = years[n]
    # min and max ice extent contours over mean extent pcolor
    # 1980-1999
    icemean = np.nanmean(iproj[start:end,:,:],axis = 0)
    icemax = np.nanmean(iproj[range(start+73,end,365),:,:],axis = 0)
    icemin = np.nanmean(iproj[range(start+257,end,365),:,:],axis = 0)
    icemax = np.select([icemax >= 0.15, icemax < 0.15, np.isnan(icemax)],[1, 0, np.nan])
    icemin = np.select([icemin >= 0.15, icemin < 0.15, np.isnan(icemin)],[1, 0, np.nan])

    #    map = Basemap(projection = proj, lat_0 = -90, lon_0 = 180,boundinglat = -60, round = True, ax = ax2)
#    m22 = map.contourf(x,y, pol_freq_mam, ax = ax2)
    blat = np.arange(40,91,1)
    blon = np.arange(-180,181,1)
    blon,blat = np.meshgrid(blon,blat)
    bground = np.zeros(blon.shape)

    proj = 'spstere'
    f,axs=plt.subplots(1,1)
    map = Basemap(projection = proj, lat_0 = -90, lon_0 = 180, boundinglat = -40, round = True, ax = axs)
    x,y = map(lon,lat)
    bx,by = map(blon,blat)
    m1 = map.pcolormesh(bx,by,bground,vmin = 0,vmax = 1,cmap = 'PuBu_r',edgecolors = 'none')
    m2 = map.pcolormesh(x,y,icemean,vmin = 0,vmax = 1,cmap = 'PuBu_r',edgecolors = 'none')
    m3 = map.contour(x,y,icemax,1,colors = 'b',linewidths = 2)
    m4 = map.contour(x,y,icemin,1,colors = 'r',linewidths = 2)
    map.drawcoastlines(color = 'none')
    map.fillcontinents(color='white')
    axs.set_title(yrs + ' ice extent')
    f.savefig('ice_extent_SH_' + yrs + '.png')

plt.show()
