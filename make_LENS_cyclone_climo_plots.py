from import_LE_data import *
from cyclone_composite_LENS import *
import matplotlib.pyplot as plt
import matplotlib.colors as colors

# more or bigger storms?
# make bullet list of points - have storms, number of storms change?
# lance bosart - papers on big storm fueled by sea ice loss NAAS workshop video

# True: Northern Hemisphere, False: Southern Hemisphere
nh = False

print 'loading data'
#psl = read_atm_data('PSL','001') 
#icefrac = read_atm_data('ICEFRAC','001')
pproj = np.load('/glade/scratch/aordonez/PSL_001_proj_test_global_blin.npy')
iproj = np.load('/glade/scratch/aordonez/ICEFRAC_001_test_proj_global_blin.npy')
if nh:
    ice_func = read_ice_data
    ocn_func = read_ocn_data 
    latlon_func = read_native_lat_lon_ice
    lows_func = find_cyclone_center
else:
    ice_func = read_ice_data_SH
    ocn_func = read_ocn_data_SH
    latlon_func = read_native_lat_lon_ice_SH
    lows_func = find_cyclone_center_SH
#lat,lon = read_native_lat_lon_atm()
#tlat,tlon = latlon_func()
lat,lon = read_stereo_lat_lon()

print "finding lows"
start = 0
end = 20*365
lows1 = find_cyclone_center(pproj[start:end,:,:],iproj[start:end,:,:],lat,107000,90000)
klows = np.where(lows1 == 1.0)
#p1 = psl[klows] / 100.

box,_,_ = get_boxes(lows1,pproj[start:end,:,:],50,lat,lon,90000) 
# replace no data and absurd values with nan:
box[box == 0.0] = np.nan
box[box > 108000] = np.nan   
box[box < 90000] = np.nan
boxslice = box[:,41:60,41:60]
pmean = np.nanmean(np.nanmean(boxslice,axis = 2),axis = 1)
k = np.where(box[:,50,50] < pmean)
box = box[k[0],:,:]
kx = klows[0][k]
ky = klows[1][k]
kz = klows[2][k]
lows1_filt = np.zeros(lows1.shape)
lows1_filt[kx,ky,kz] = 1
#lows = find_cyclone_center_SH(psl[start:end,:,:],icefrac[start:end,:,:],lat,107000,90000)
#k = np.where(lows == 1.0)
#psh1 = psl[k] / 100.
"""
start = 365*20*4
end = 365*20*5
lows2 = find_cyclone_center(pproj[start:end,:,:],iproj[start:end,:,:],lat,107000,90000)
klows = np.where(lows2 == 1.0)
#ptmp = psl[start:end,:,:]
#p2 = ptmp[k] / 100.
box,_,_ = get_boxes(lows2,pproj[start:end,:,:],50,lat,lon,90000) 
# replace no data and absurd values with nan:
box[box == 0.0] = np.nan
box[box > 108000] = np.nan   
box[box < 90000] = np.nan
boxslice = box[:,41:60,41:60]
pmean = np.nanmean(np.nanmean(boxslice,axis = 2),axis = 1)
k = np.where(box[:,50,50] < pmean)
kx = klows[0][k]
ky = klows[1][k]
kz = klows[2][k]
lows2_filt = np.zeros(lows2.shape)
lows2_filt[kx,ky,kz]=1
"""
#lows = find_cyclone_center_SH(psl[start:end,:,:],icefrac[start:end,:,:],lat,107000,90000)
#k = np.where(lows == 1.0)
#psh2 = ptmp[k] / 100.

print "making plots"

#histogram of central pressure
"""f3,axs = plt.subplots(2,1)
ax1 = axs[0].hist(lows1_filt, bins = 20, histtype = 'stepfilled', color = 'r', label = '1980-1999', normed = True)
ax2 = axs[0].hist(lows2_filt, bins = 20, histtype = 'stepfilled', alpha = 0.5, color = 'b', label = '2060-2079', normed = True)
axs[0].set_ylabel('Probability')
axs[0].set_xlabel('Central pressure (mb)')
axs[0].set_title('Northern Hemisphere')
axs[0].legend(loc = 'upper left')"""

lowsum = np.nansum(lows1, axis = 0)
k = np.where(lowsum == 1)
newx = np.arange(-30,-90,-0.5)
newy = np.arange(0,360,0.5)
X,Y = np.meshgrid(newx,newy)
s = interpolate.griddata((lon[k],lat[k]),lowsum[k],(X,Y),method = 'linear')

time = kx
lowrow = ky
lowcol = kz
lows = lows1

#k=np.where(lows_edit == 1.)
latplot = np.tile(lat,(lows.shape[0],1,1))
lonplot = np.tile(lon,(lows.shape[0],1,1))
latplot = latplot[time,lowrow,lowcol]
lonplot = lonplot[time,lowrow,lowcol]
latplot = np.reshape(latplot,(latplot.shape[0]))
lonplot = np.reshape(lonplot,(lonplot.shape[0]))
#lonplot = np.select([lows1 == 1.][lonplot])
proj = 'npstere'
f,axs = plt.subplots(1,1)
map = Basemap(projection = proj, lat_0 = -90, lon_0 = 180, boundinglat = -40, round = True, ax = axs)
map.drawcoastlines()
m  = map.plot(lonplot,latplot,'bo',markersize = 1, latlon = True)
plt.show()
f.savefig('sh_cyclone_count.png')



"""f,axs=plt.subplots(1,1)
xedges = np.arange(0,365.,1)
yedges = np.arange(30.,90.,1)
map = Basemap(projection = proj, lat_0 = 90, lon_0 = 180, boundinglat = 40, round = True, ax = axs)
H,xedges,yedges = np.histogram2d(latplot,lonplot,bins = [yedges,xedges])
Y,X = np.meshgrid(yedges[:-1],xedges[:-1])
m = map.pcolormesh(X,Y,H,latlon = True)
plt.colorbar(m,ax = axs)
f.show()"""

