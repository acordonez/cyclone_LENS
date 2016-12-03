from import_LE_data import *
from cyclone_composite_LENS import *
import matplotlib.pyplot as plt

pproj = np.load('/glade/scratch/aordonez/pproj.npy')
iproj = np.load('/glade/scratch/aordonez/iproj.npy')
tproj = np.load('/glade/scratch/aordonez/tproj.npy')
tproj = tproj[0:3650 * 3,:,:]
iproj = iproj[0:3650 * 3,:,:]
pproj = pproj[0:3650 * 3,:,:]

lat,lon = read_stereo_lat_lon()

lows = find_cyclone_center(pproj,iproj,104000,90000)
if np.max(lows) == 1.0:
    box,_,_ = get_boxes(lows,pproj,50,lat,lon,90000) 
    # replace no data and absurd values with nan:
    box[box == 0.0] = np.nan
    box[box > 108000] = np.nan   
    box[box < 90000] = np.nan

    tbox,_,_=get_boxes(lows,tproj,50,lat,lon,90000)
    tbox[tbox == 0.0] = np.nan
    tbox[tbox < 180.] = np.nan
    tbox[tbox > 360.] = np.nan

    # remove samples where SLP at 'low' is greater than box average:
    pmean = np.nanmean(np.nanmean(box,axis = 2),axis = 1)
    k = np.where(box[:,50,50] < pmean)
    box = box[k[0],:,:]
    tbox = tbox[k[0],:,:]

    plt.pcolormesh(np.nanmean(tbox,axis = 0))
    plt.contour(np.nanmean(box,axis = 0),range(96000,103000,100))
    
else:
    print "Max of lows is ",np.max(lows)
    print "Check that pmax and pmin are in correct order"


