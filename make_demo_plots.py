from import_LE_data import *
from cyclone_composite_LENS import *
import matplotlib.pyplot as plt

psl = read_atm_data('PSL','001')
psl = psl[36:38,:,:]
pproj = np.load('/glade/scratch/aordonez/pproj_small.npy')
iproj = np.load('/glade/scratch/aordonez/iproj_small.npy')
lat,lon = read_stereo_lat_lon()

#np.save('/glade/scratch/aordonez/pproj_small.npy',pproj[0:3650,:,:])
#np.save('/glade/scratch/aordonez/iproj_small.npy',iproj[0:3650,:,:])

pproj = pproj[36:38,:,:]
iproj = iproj[36:38,:,:]
lows = find_cyclone_center(pproj,iproj,104000,90000)
lon[lon < 0.] = lon[lon < 0.] + 360.
data = pproj
size = 50
long_size = ((size *2) + 1)
mylow = np.where(lows == 1)
nlows = mylow[0].shape[0]
data_box = np.zeros((nlows,long_size,long_size))
(tmax, ymax, xmax) = data.shape
# get lon where north is up
lon0 = lon[0,(xmax/2)-1]
count = 0

for ind in range(2,-1,-1):
    time = mylow[0][ind]
    lowrow = mylow[1][ind]
    lowcol = mylow[2][ind]
    # -----------------
    # rotation to north
    # -----------------
    mylon = lon[lowrow,lowcol]
    low_mask = np.zeros((ymax,xmax))
    low_mask[lowrow,lowcol] = 1
    if lon0 < mylon:
        deg = mylon - lon0
    elif lon0 >= mylon:
        deg = (360 + mylon) - lon0
    low_rotated = interpolation.rotate(low_mask,deg)
    # because of interpolation, lows != 1
    ynew,xnew = np.where(low_rotated == low_rotated.max())   
    data_rotated = interpolation.rotate(data[time,:,:],deg)
    # take out noisy grid cells near coast
    coast = buffer_coast(data_rotated, buf = (8,8), edgedif = 90000.)
    data_rotated = data_rotated * coast 
    # -----------------
    # extracting box
    # -----------------
    y1 = ynew - size
    y2 = ynew + size + 1
    x1 = xnew - size
    x2 = xnew + size + 1
    if (y1 < 0) | (x1 < 0) | (y2 > ymax) | (x2 > xmax):
        # too close to edge of map
        continue
    else:
        data_box[count,:,:] = data_rotated[y1:y2,x1:x2]
        count += 1

f1,ax1 = plt.subplots(1,1)
a1 = ax1.pcolormesh(psl[0,:,:],vmin = 94000,vmax = 104000)
f1.colorbar(a1,ax = ax1)
ax1.set_title('Model output sea level pressure')
f1.savefig('demo_1_model.png')

f2,ax2 = plt.subplots(1,1)
a2 = ax2.pcolormesh(pproj[0,:,:],vmin = 94000,vmax = 104000)
f2.colorbar(a2,ax = ax2)
ax2.set_title('NH Stereo sea level pressure')
f2.savefig('demo_2_stereo.png')

f3,ax3 = plt.subplots(1,1)
a3 = ax3.pcolormesh(pproj[0,:,:],vmin = 94000,vmax = 104000)
k = np.where(lows[0,:,:] == 1)
ax3.scatter(k[1], k[0],color = 'k') 
f3.colorbar(a3,ax = ax3)
ax3.set_title('NH Stereo sea level pressure with lows')
f3.savefig('demo_3_lows.png')

f4,ax4 = plt.subplots(1,1)
a4 = ax4.pcolormesh(data_rotated,vmin = 94000,vmax = 104000)
f4.colorbar(a4, ax = ax4)
ax4.set_title('Rotated NH')
f4.savefig('demo_4_rotated.png')

f5,ax5 = plt.subplots(1,1)
a5 = ax5.pcolormesh(data_box[2,:,:],vmin = 94000,vmax = 104000)
f5.colorbar(a5,ax = ax5)
ax5.set_title('Clipped rotated box')
f5.savefig('demo_5_box.png')
