import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from mpl_toolkits.basemap import Basemap
from scipy import ndimage
import scipy.ndimage.filters as filters
import scipy.ndimage.morphology as morphology
import scipy.ndimage.interpolation as interpolation
from netCDF4 import Dataset

#------------------------------------------------
# functions to read data and get low positions
#------------------------------------------------

def detect_local_minima(arr):
    # http://stackoverflow.com/questions/3684484/peak-detection-in-a-2d-array/3689710#3689710
    """
    Takes an array and detects the troughs using the local minimum filter.
    Returns a boolean mask of the troughs (i.e. 1 when
    the pixel's value is the neighborhood minimum, 0 otherwise)
    """
    # define an connected neighborhood
    neighborhood = morphology.generate_binary_structure(2,2)
    neighborhood = morphology.binary_dilation(neighborhood,iterations = 5)
    # apply the local minimum filter; all locations of minimum value 
    # in their neighborhood are set to 1
    # filter multiple times to get just one point per cyclone(A.O.)
    tmp = filters.minimum_filter(arr, footprint=neighborhood)
    local_min = (filters.minimum_filter(tmp, footprint=neighborhood)==arr)
    background = (arr==0)
    # 
    # a little technicality: we must erode the background in order to 
    # successfully subtract it from local_min, otherwise a line will 
    # appear along the background border (artifact of the local minimum filter)
    eroded_background = morphology.binary_erosion(
        background, structure=neighborhood, border_value=1)
    #
    # we obtain the final mask, containing only peaks, 
    # by removing the background from the local_min mask
    detected_minima = local_min - eroded_background 
    return detected_minima 

def find_cyclone_center(psl,icefrac,pmax,pmin):
    """
    find_cyclone_center
    
    Returns a matrix (time x lon x lat). Cells with 
    a "1" indicate a low pressure center; cells equal "0" 
    otherwise.

    psl: numpy array of sea level pressure. land areas masked with np.nan
    icefrac: numpy array of sea ice concentration on atmosphere grid, max = 1.
    pmax: numeric, maximum allowed value of central pressure
    pmin: numeric, minimum allowed value of central pressure
    """
    time,rows,cols = psl.shape
    lows = np.zeros((psl.shape)) 
    # find the lows
    for n in range(0,psl.shape[0]):
        low_n = detect_local_minima(psl[n,:,:])
        lows[n,:,:] = np.select([(low_n == True) & (icefrac[n,:,:] > 0) &
                                 (psl[n,:,:] != np.nan)],[low_n])
    return lows


def get_boxes(lows,data,size,lat,lon):
    """
    box = get_boxes(lows, data, size)
    lows: binary matrix where 1 = low pressure center
    data: the data to subset
    size: half the length of the 2D subset box
    box:  flattened array of data around low pressure centers
    """
    long_size = ((size *2) + 1)
    mylow = np.where(lows == 1)
    nlows = mylow[0].shape[0]
    data_box = np.zeros((nlows,long_size,long_size))
    (tmax, ymax, xmax) = data.shape
    count = 0

    for ind in range(0,100):
        # Rotate grid so that north is always up
        # and grab box of data around each low
        time = mylow[0][ind]
        lowcol = mylow[1][ind]
        lowrow = mylow[2][ind]
        mylon = lon[lowrow,lowcol]
        low_mask = np.zeros((ymax,xmax))
        low_mask[lowrow,lowcol] = 1
        deg = mylon - lon[0,(xmax/2)-1]
        low_rotated = interpolation.rotate(low_mask,deg)
        # because of interpolation, lows != 1
        ynew,xnew = np.where(low_rotated == low_rotated.max())
        data_rotated = interpolation.rotate(data[time,:,:],deg)
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
    return data_box[0:count,:,:]

def get_mam(data):
    """Pulls out a timeseries only containing days in 
    March, April, and May. Input timeseries may not contain
    partial years of data.
    """
    nyrs = data.shape[0] / 365 
    time = len(range(59,151))
    for yr in range(0,nyrs):
        if yr == 0:
            data_mam = data[59:151,:,:]
        else:
            data_mam = np.concatenate((data_mam,data[59+(yr*365):151+(yr*365),:]),axis = 0)
    return data_mam

def get_jja(data):
    """Pulls out a timeseries only containing days in 
    June, July, and August
    """
    nyrs = data.shape[0] / 365 
    time = len(range(151,243))
    for yr in range(0,nyrs):
        if yr == 0:
            data_son = data[151:243,:,:]
        else:
            data_son = np.concatenate((data_son,data[151+(yr*365):243+(yr*365),:,:]),axis = 0)
    return data_son

def get_son(data):
    """Pulls out a timeseries only containing days in 
    September, October, and November
    """
    nyrs = data.shape[0] / 365 
    time = len(range(243,334))
    for yr in range(0,nyrs):
        if yr == 0:
            data_son = data[243:334,:,:]
        else:
            data_son = np.concatenate((data_son,data[243+(yr*365):334+(yr*365),:,:]),axis = 0)
    return data_son

def get_djf(data):
    """Pulls out a timeseries only containing days in 
    December, January, and February
    """
    nyrs = data.shape[0] / 365 
    time = len(range(0,60)+range(334,365))
    for yr in range(0,nyrs):
        if yr == 0:
            data_djf = data[0:60,:,:]
            data_djf = np.concatenate((data_djf,data[334:365,:,:]),axis = 0)
        else:
            data_djf = np.concatenate((data_djf,data[0+(yr*365):60+(yr*365),:]),
                       axis = 0)
            data_djf = np.concatenate((data_djf,data[334+(yr*365):365+(yr*365),:]),
                       axis = 0)
    return data_djf


def get_difference_from_mean(icedata,vardata,tarea):
    """Calculates the daily anomalies relative to values
    over the ice pack where concentration >= 15%
    """
    if np.max(icedata) > 1:
        icedata[icedata < 15] = 0.
        icedata[icedata > 0] = 1.
    else:
        icedata[icedata < 0.15] = 0.
        icedata[icedata > 0] = 1.
    tarea = np.tile(tarea,(icedata.shape[0],1,1))
    # compute area weights over sea ice
    wgt_total = np.nansum(np.nansum((tarea * icedata),axis =2), axis = 1)
    wgt_total = np.tile(wgt_total,(vardata.shape[1],vardata.shape[2],1))
    wgt_total = np.transpose(wgt_total,(2,0,1))
    wgt_total[icedata != 1] = 0.
    wgt = tarea * icedata / wgt_total
    # get time series of mean values
    varmean = np.nansum(np.nansum(wgt * vardata, axis = 2),axis = 1)
    varmean_overall = np.nanmean(varmean,axis = 0)
    #print varmean[0]
    varmean = np.tile(varmean, (vardata.shape[1],vardata.shape[2],1)) 
    varmean = np.transpose(varmean,(2,0,1))
    vardata = vardata - varmean
    return vardata, varmean_overall

def plot_composites(res,storm_set,season,X1,Y1,X5,Y5,varlist,
                    titles,sname,Pmean,Tmean,Umean,Vmean):
    """Makes a composite plot of the input data field with 
    temperature and pressure contours and wind streamlines. 
    Saves the plot.
    """
    # figure out rows and columns
    if len(varlist) < 3:
        row_n = 1
        col_n = len(varlist)
    else:
        row_n = 2
        col_n = int(np.ceil(len(varlist)/2.))
    # create plot
    f,axs = plt.subplots(nrows = row_n, ncols = col_n, figsize=(row_n * 3., col_n * 2.5))
    axs = np.reshape(axs,(row_n*col_n,1))
    Tlevels = np.arange(200,310,5)

    for ind in range(0,len(varlist)):
        cmax = np.max([abs(np.nanmax(varlist[ind])),abs(np.nanmin(varlist[ind]))])
        ax1 = axs[ind,0].pcolormesh(X1,Y1,varlist[ind],cmap='PuOr',norm=colors.SymLogNorm(linthresh=0.1,linscale=0.1,vmin=-1*cmax,vmax=cmax))
        c1 = axs[ind,0].contour(X5,Y5,Pmean,colors='k',linewidths = 1)
        c2 = axs[ind,0].contour(X5,Y5,Tmean,Tlevels,colors='r',linewidths = 1)
        #axs[ind,0].streamplot(X5,Y5,Umean,Vmean,linewidth = 1)
        plot_spaced_quivers(axs[ind,0],X5,Y5,Umean,Vmean,spacing = 5);
        axs[ind,0].clabel(c1, inline=1, fontsize=8)
        axs[ind,0].clabel(c2, inline=1, fontsize=8)
        axs[ind,0].set_title(titles[ind], fontsize = 10)    
        cb1 = f.colorbar(ax1, ax = axs[ind,0], ticks = [-1.0*cmax,-0.5*cmax,-0.2*cmax,0.,0.2*cmax,0.5*cmax,1.0*cmax])    
        cb1.ax.tick_params(labelsize = 8)   
        #plt.setp(axs[ind,0].get_xticklabels(), fontsize = 8)
        #plt.setp(axs[ind,0].get_yticklabels(), fontsize = 8)
        axs[ind,0].get_xaxis().set_ticks([])
        axs[ind,0].get_yaxis().set_ticks([])
    savename = "cyclone_composite_" + sname + "_" + res + "_c_" + storm_set + "_" + season + ".png"
    f.savefig(savename)

def get_p_limits(storm_set):
    """Returns the minimum and maximum pressure 
    for the input storm category
    """
    storm_dict = {'strong':[94000,90000],'medium':[98000,94000],'bulk':[97000,90000]}
    pmax = storm_dict[storm_set][0]
    pmin = storm_dict[storm_set][1]
    return pmax, pmin


