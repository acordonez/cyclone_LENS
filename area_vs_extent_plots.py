"""area_vs_extent_plots.py

Creates scatter plot of the ice area:extent ratio in 
all 30 members of the CESM LENS for years 1980-2100
"""

from import_LE_data import read_area_ice,read_ice_data_monthly
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.gridspec as gridspec

area = read_area_ice()
num_list = [str(x) for x in range(1,31)]
ratio = np.zeros((30,1453))
for ind,model_num in enumerate(num_list):
    aice = read_ice_data_monthly('aice',model_num)
    # remove fill values and mask
    aice[aice > 110] = 0
    aice = aice.data
    aice = aice / 100.
    aice[aice < 0.15] = 0
    ice_area = np.nansum(np.nansum(aice * area, axis = 2),axis = 1)
    ice_extent = np.select([aice>=0.15],[1])
    ice_extent = np.nansum(np.nansum(ice_extent * area, axis = 2),axis = 1)
    ratio[ind,:] = ice_area / ice_extent

# plot annual ratios with +- 2 std (based on 1980-1999)
ratio_ann = np.reshape(ratio[:,0:1452],(30,121,12))
ratio_ann = np.mean(ratio_ann,axis = 2)
ratio_std = np.std(ratio_ann[:,0:30].flatten()) * 2.
ratio_mean = np.nanmean(ratio_ann[:,0:30].flatten())
ratio_ann = np.transpose(ratio_ann)

years = np.transpose(np.tile(np.arange(1980,2101),(30,1)))

#1980-2100 subplot
gs = gridspec.GridSpec(2,2)
ax1 = plt.subplot(gs[0,:])
ax1.set_autoscale_on(False)
ax1.fill_between(years[:,0],ratio_mean + ratio_std,ratio_mean - ratio_std, color = [0.7,0.7,0.7])
colors = iter(cm.rainbow(np.linspace(0,1,30)))
for y in range(0,30):
    ax1.scatter(years[:,y],ratio_ann[:,y],color = next(colors))
ax1.set_title('1980-2100')
ax1.axis([1980,2100,0.5,1.0])

#1980-2030 subplot
ax2 = plt.subplot(gs[1,:])
ax2.set_autoscale_on(False)
ax2.fill_between(years[:,0],ratio_mean + ratio_std,ratio_mean - ratio_std, color = [0.7,0.7,0.7])
colors = iter(cm.rainbow(np.linspace(0,1,30)))
for y in range(0,30):
    ax2.scatter(years[0:50,y],ratio_ann[0:50,y],color = next(colors))
ax2.set_title('1980-2030')
ax2.axis([1980,2030,0.75,0.95])

#2060-2100 subplot
#ax3 = plt.subplot(gs[1,1])
#ax3.set_autoscale_on(False)
#ax3.fill_between(years[:,0],ratio_mean + ratio_std,ratio_mean - ratio_std, color = [0.7,0.7,0.7])
#colors = iter(cm.rainbow(np.linspace(0,1,30)))
#for y in range(0,30):
#    ax3.scatter(years[80:120,y],ratio_ann[80:120,y],color = next(colors))
#ax3.set_title('2060-2100')
#ax3.axis([2060,2100,0.6,3.1])

plt.suptitle('Ratio of ice area:extent in 30 LENS members')
#plt.show()
plt.savefig('LENS_area_extent_ratio_21stcentury.png')
