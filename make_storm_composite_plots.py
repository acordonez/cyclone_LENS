from import_LE_data import *
from cyclone_composite_LENS import *
import matplotlib.pyplot as plt

pproj = np.load('/glade/scratch/aordonez/pproj.npy')
iproj = np.load('/glade/scratch/aordonez/iproj.npy')

lat,lon = read_stereo_lat_lon()

lows = find_cyclone_center(pproj,iproj,0,0)
box = get_boxes(lows,pproj,15,lon)

20yrs = 20 * 365 
n = 5

for era in range(0,n):
    start = era * 20yrs
    end = (era + 1) * 20yrs
    boxes = box[start:end,:,:]
