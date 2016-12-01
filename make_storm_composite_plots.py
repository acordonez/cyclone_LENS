from import_LE_data import *
from cyclone_composite_LENS import *
import matplotlib.pyplot as plt

pproj = np.load('/glade/scratch/aordonez/pproj_small.npy')
iproj = np.load('/glade/scratch/aordonez/iproj_small.npy')
#tproj = np.load('/glade/scratch/aordonez/tproj.npy')

lat,lon = read_stereo_lat_lon()

#tproj = tproj[0:(365*10),:,:]

lows = find_cyclone_center(pproj,iproj,98000,90000)
if np.max(lows) == 1.0:
    box = get_boxes(lows,pproj,50,lon,90000) 
    box[box == 0.0] = np.nan
    #tbox = get_boxes(lows,tproj,50,lon,200)
else:
    print "Max of lows is ",np.max(lows)
    print "Check that pmax and pmin are in correct order"

plot_mean(box)

20yrs = 20 * 365 
n = 5

for era in range(0,n):
    start = era * 20yrs
    end = (era + 1) * 20yrs
    boxes = box[start:end,:,:]
