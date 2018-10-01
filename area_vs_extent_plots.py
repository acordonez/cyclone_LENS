"""area_vs_extent_plots.py

Creates scatter plot of the ice area:extent ratio in 
all 30 members of the CESM LENS for years 1980-2100
"""

from import_LE_data import read_area_ice,read_ice_data_monthly
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.gridspec as gridspec
import matplotlib.lines as mlines

area = read_area_ice()
num_list = [str(x) for x in range(1,36)]
ratio = np.zeros((35,121))
ratio2 = np.zeros((35,121))
areamean = np.zeros((35,121))
extentmean = np.zeros((35,121))
volmean = np.zeros((35,121))
for ind,model_num in enumerate(num_list):
    print model_num
    aice = read_ice_data_monthly('aice',model_num)
    hi = read_ice_data_monthly('hi',model_num)
    time = aice.shape[0]
    aice = aice[range(8,time,12),:,:]
    hi = hi[range(8,time,12),:,:]
    # remove fill values and mask
    aice[aice > 110] = 0
    hi[hi > 100] = 0
    aice = aice.data
    hi = hi.data
    aice = aice / 100.
    aice[aice < 0.15] = 0
    hi[aice < 0.15] = 0
    ice_area = np.nansum(np.nansum(aice * area, axis = 2),axis = 1)
    ice_extent = np.select([aice>=0.15],[1])
    ice_extent = np.nansum(np.nansum(ice_extent * area, axis = 2),axis = 1)
    ratio[ind,:] = ice_area / ice_extent
    ice_vol = np.nansum(np.nansum(hi*area, axis  = 2), axis = 1)
    ratio2[ind] = ice_area / ice_vol
    areamean[ind,:] = ice_area
    extentmean[ind,:] = ice_extent
    volmean[ind,:] = ice_vol

# plot annual ratios with +- 2 std (based on 1980-1999)
"""ratio_ann = np.reshape(ratio[:,0:1452],(30,121,12))
ratio_ann = np.mean(ratio_ann,axis = 2)
ratio_std = np.std(ratio_ann[:,0:30].flatten()) * 2.
ratio_mean = np.nanmean(ratio_ann[:,0:30].flatten())
ratio_ann = np.transpose(ratio_ann)
"""
# get standard deviation and mean from 1980-2010
ratio_ann = ratio[:,0:120]
ratio_std = np.std(ratio_ann[:,0:30].flatten()) * 2
ratio_mean = np.nanmean(ratio_ann[:,0:30].flatten())
ratio_ann = np.transpose(ratio_ann)

ratio_ann2 = ratio2[:,0:120]
ratio_ann2 = np.transpose(ratio_ann2)

areamean = np.transpose(np.nanmean(areamean[:,0:120],axis = 0))
extentmean = np.transpose(np.nanmean(extentmean[:,0:120],axis = 0))
volmean = np.transpose(np.nanmean(volmean[:,0:120],axis = 0))

obsfile = open("area_extent_ratio_obs_111978-122015_nasateam.txt",'r')
obs = []
for line in obsfile:
    data = line.split()
    if data == []:
        continue
    obs.append(float(data[0]))
obs = np.array(obs)
# get september
obs = obs[10:12:445]

years = np.transpose(np.tile(np.arange(1980,2100),(35,1)))

#1980-2100 subplot
gs = gridspec.GridSpec(2,2)
ax1 = plt.subplot(gs[0,:])
#ax1.set_autoscale_on(False)
ax1.fill_between(years[:,0],ratio_mean + ratio_std,ratio_mean - ratio_std, color = [0.7,0.7,0.7])
colors = iter(cm.rainbow(np.linspace(0,1,35)))
for y in range(0,35):
    ax1.scatter(years[:,y],ratio_ann[:,y],color = next(colors))
ax1.set_title('area:extent')
ax1.axis([1980,2100,0.0,1.0])

"""ax2= plt.subplot(gs[1,:])
ax3 = ax2.twinx()
a = ax2.plot(years[:,0],areamean,color = 'b',linewidth = 2)
b = ax2.plot(years[:,0],extentmean,color = 'r',linewidth = 2)
c = ax3.plot(years[:,0],volmean,color = 'g',linewidth = 2)
ax2.set_title('Ensemble mean')
ax2.set_xlim([1980,2100])
ax2.set_ylabel('area m**2')
ax3.set_ylabel('volume m**3')
aline = mlines.Line2D([],[],color = 'b',label = 'area')
bline = mlines.Line2D([],[],color = 'r',label = 'extent')
cline = mlines.Line2D([],[],color = 'g',label = 'volume')
mylines = [aline,bline,cline]
ax2.legend(handles = mylines,labels = [h.get_label() for h in mylines])"""

"""#1980-2100 area/volume subplot
gs = gridspec.GridSpec(2,2)
ax2 = plt.subplot(gs[1,:])
#ax1.set_autoscale_on(False)
#ax1.fill_between(years[:,0],ratio_mean + ratio_std,ratio_mean - ratio_std, color = [0.7,0.7,0.7])
colors = iter(cm.rainbow(np.linspace(0,1,35)))
for y in range(0,35):
    ax2.scatter(years[:,y],ratio_ann2[:,y],color = next(colors))
ax2.set_title('area:volume')
ax2.set_xlim([1980,2100])"""

#1980-2030 subplot
ax2 = plt.subplot(gs[1,:])
ax2.set_autoscale_on(False)
ax2.fill_between(years[:,0],ratio_mean + ratio_std,ratio_mean - ratio_std, color = [0.7,0.7,0.7])
colors = iter(cm.rainbow(np.linspace(0,1,35)))
for y in range(0,35):
    ax2.scatter(years[0:50,y],ratio_ann[0:50,y],color = next(colors))
print years.shape, obs.shape
ax2.scatter(years[0:38,y],obs,color = 'k')
ax2.set_title('1980-2030')
ax2.axis([1980,2030,0.25,0.95])

#2060-2100 subplot
#ax3 = plt.subplot(gs[1,1])
#ax3.set_autoscale_on(False)
#ax3.fill_between(years[:,0],ratio_mean + ratio_std,ratio_mean - ratio_std, color = [0.7,0.7,0.7])
#colors = iter(cm.rainbow(np.linspace(0,1,30)))
#for y in range(0,30):
#    ax3.scatter(years[80:120,y],ratio_ann[80:120,y],color = next(colors))
#ax3.set_title('2060-2100')
#ax3.axis([2060,2100,0.6,3.1])

plt.suptitle('September area:extent ratio in 35 LENS members',fontsize = 14)
plt.savefig('LENS35_area_extent_ratio_21stcentury_sept_withobs.png')
