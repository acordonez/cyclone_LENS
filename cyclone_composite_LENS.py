import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from mpl_toolkits.basemap import Basemap
from scipy import ndimage
import scipy.ndimage.filters as filters
import scipy.ndimage.morphology as morphology
import scipy.ndimage.interpolation as interpolation
import scipy.interpolate as interpolate
import scipy.misc as misc
from netCDF4 import Dataset

#------------------------------------------------
# functions to read data and get low positions
#------------------------------------------------

def detect_local_minima(arr):
    # http://stackoverflow.com/questions/3684484/peak-detection-in-a-2d-array/3689710#3689710
    """
    Takes an array and detects the troughs using the 
    local minimum filter.Returns a boolean mask of the troughs 
    (i.e. 1 when the pixel's value is the neighborhood 
    minimum, 0 otherwise)

    Parameters:
    --------------------
    arr: numpy array where the land is masked with 0
    detected_minima: numpy array 
    """
    # define an connected neighborhood
    neighborhood = morphology.generate_binary_structure(2,2)
    neighborhood = morphology.binary_dilation(neighborhood,iterations = 20)
    # apply the local minimum filter; all locations of minimum value
    # in their neighborhood are set to 1
    # filter multiple times to get just one point per cyclone(A.O.)
    tmp = filters.minimum_filter(arr, footprint=neighborhood)
    tmp = filters.minimum_filter(tmp, footprint=neighborhood)
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

def buffer_coast(pdata,buf = (5,5), edgedif = 90000.):
    """buffer_coast()
  
    Returns a mask that can be used to select only those
    grid cells which are not adjacent to land within
    a distance given by 'buf'
    
    Parameters:
    --------------------
    pdata: numpy array of pressure data
    buf: 2d buffer shape. (5,5) minimum for continuous coast outline
    edgedif: expected difference between water and land values
    mask: numpy array. 0 at coast, 1 away from coast
    """
    edge = morphology.morphological_gradient(pdata,\
    size = buf)
    edge[edge < edgedif] = 1.
    edge[edge >= edgedif] = 0.
    bi = morphology.generate_binary_structure(2,2)
    mask = morphology.binary_dilation(edge,\
    structure = bi, iterations = 1)
    return mask.astype(int)


def find_cyclone_center(psl,icefrac,lat,pmax,pmin):        
    """
    find_cyclone_center
    
    Returns a matrix (time x lon x lat). Cells with 
    a "1" indicate a low pressure center; cells equal "0" 
    otherwise.

    For a pixel to be counted as a low, it must meet these criteria:
    There must be a local minima in the sea level pressure (SLP),  
    a local maxima in the laplacian of SLP, greater than 15% ice 
    cover, and the SLP must be between the bounds 'pmin' and 'pmax'.
    Grid cells near the coast are not included due to noise from 
    the stereo regridding and rotation done in get_boxes()

    Parameters:
    --------------------
    psl: numpy array of sea level pressure. land areas masked with 0
    icefrac: numpy array of sea ice concentration on atmosphere grid, 
             max = 1
    pmax: numeric, maximum allowed value of central pressure
    pmin: numeric, minimum allowed value of central pressure
    lows: numpy array
    """
    time,rows,cols = psl.shape
    lows = np.zeros((psl.shape)) 
    # find the lows
    for n in range(0,psl.shape[0]):
        lap = filters.laplace(psl[n,:,:])
        lapmax = detect_local_minima(lap*-1.)
        coast = buffer_coast(psl[n,:,:],buf = (5,5))
        ptmp = psl[n,:,:] * coast
        low_n = detect_local_minima(ptmp)
        lows[n,:,:] = np.select([(low_n == True) & 
                                 #(icefrac[n,:,:] > 0.15) &
                                 (lat > 65) &
                                 (psl[n,:, :] <= pmax) & 
                                 (psl[n,:,:] >= pmin) &
                                 (lapmax ==1) & 
                                 (coast == 1)],[low_n])
    return lows


def get_boxes(lows,data,size,lat,lon,edgedif):
    """
    box = get_boxes(lows, data, size)
   
    Clips a square of length(2 x size) + 1 around each low
    pressure center in lows and returns an array with all the
    boxes.

    Parameters:
    --------------------
    lows: binary matrix where 1 = low pressure center
    data: numpy array, land masked with 0
    size: numeric, half the length of the 2D subset box
    edgedif: numeric, roughly the difference 
             in value between data and land grid cells
    box:  numpy array of data around low pressure centers
    """
    lon[lon < 0.] = lon[lon < 0.] + 360.

    long_size = ((size *2) + 1)
    mylow = np.where(lows == 1)
    nlows = mylow[0].shape[0]
    data_box = np.zeros((nlows,long_size,long_size))
    lat_box = np.zeros(data_box.shape)
    lon_box = np.zeros(data_box.shape)
    (tmax, ymax, xmax) = data.shape
    # get lon where north is up
    lon0 = lon[0,(xmax/2)-1]
    count = 0

    for ind in range(0,nlows):
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
        coast = buffer_coast(data_rotated, buf = (8,8), edgedif = edgedif)
        data_rotated = data_rotated * coast 
        #ynew,xnew = lowrow,lowcol
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
            #data_box[count,:,:] = data[ind,y1:y2,x1:x2]
            lat_box[count,:,:] = lat[y1:y2,x1:x2]
            lon_box[count,:,:] = lon[y1:y2,x1:x2]
            count += 1
    return data_box[0:count,:,:], lat_box[0:count,:,:], lon_box[0:count,:,:]

def regrid_to_conic(lat,lon,lat_ref,lon_ref,lat_stnd1,lat_stnd2):
    # regrid to conformal conic
    # equations from https://en.wikipedia.org/wiki/Lambert_conformal_conic_projection
    row,col = lat.shape
    lat = lat * (np.pi / 180.)
    lon = lon * (np.pi / 180.)
    lon_ref = lon_ref * (np.pi / 180.)
    lat_ref = lat_ref * (np.pi / 180.)
    lat_stnd1 = np.complex(lat_stnd1 * (np.pi / 180.))
    lat_stnd2 = np.complex(lat_stnd2 * (np.pi / 180.))
    #lon_ref = lon[0,col/2]
    #lat_ref = lat[row/2,col/2]
    #lat_stnd2 = np.complex(lat[row/2,col/2] - 0.01)
    #lat_stnd1 = np.complex(lat[row/2,col/2] + 0.01)
    
    n_top = np.log(np.cos(lat_stnd1) * 1./np.cos(lat_stnd2))
    n_bottom = np.log(np.tan(0.25 * np.pi + 0.5 * lat_stnd2) *
                   1./np.tan(0.25 * np.pi + 0.5 * lat_stnd1))
    n = n_top / n_bottom
    F = (np.cos(lat_stnd1) * np.power(np.tan(0.25 * np.pi + 0.5 * lat_stnd1),n)) / n
    rho = F * 1./np.tan(0.25 * np.pi + 0.5 * lat)**n
    rho_0 = F * 1./np.tan(0.25 * np.pi + 0.5 * lat_ref)**n
    x = rho * np.sin(n * (lon - lon_ref))
    y = rho_0 - rho * np.cos(n * (lon - lon_ref))

    return np.real(x),np.real(y)

def resample_and_composite(data,lat,lon):
    """IN PROGRESS
    """
    x = np.zeros(data.shape)
    y = np.zerointers(data.shape)
    for ind in range(0,data.shape[0]):
        xtmp,ytmp = regrid_to_conic(lat[ind,:,:],lon[ind,:,:])
        x[ind,:,:] = xtmp
        y[ind,:,:] = ytmp
    X,Y = np.meshgrid(np.arange(-0.5,0.7,0.001),np.arange(-0.5,0.7,0.001))
    new_data = np.zeros((data.shape[0],X.shape[0],X.shape[1]))
    for ind in range(0,data.shape[0]):
        xnew,ynew = x[ind,:,:].flatten(), y[ind,:,:].flatten()
        data_tmp = data[ind,:,:].flatten()
        new_data[ind,:,:] = interpolate.griddata((xnew,ynew),data_tmp,(X,Y),method = 'linear')
    return xnew,ynew,new_data

def test_regrid_methods():
    psl = read_atm_data('PSL','002')
    icefrac = read_atm_data('ICEFRAC','002')
    lat,lon = read_native_lat_lon_atm()
    end = 20*365
    lows = find_cyclone_center(psl[0:end,:,:],icefrac[0:end,:,:],lat,104000,90000)

    t = 200
    n = 0
    l = lows[t,:,:]
    c = np.where(l == 1.0)
    lattest = lat[c[0][n],c[1][n]]
    lontest = lon[c[0][n],c[1][n]]
    x,y=regrid_to_conic(lat,lon,lattest,lontest,lattest+5,lattest - 5)
    k = np.where(np.isnan(x) == False)
    p = psl[t,:,:]
    xnew,ynew = np.meshgrid(np.arange(-0.05,.05,0.005),np.arange(-0.05,0.05,0.005))
    s = interpolate.griddata((x[k],y[k]),p[k],(xnew,ynew),method = 'linear')
    f1 = plt.figure()
    plt.pcolormesh(x,y,p,vmin = 96000,vmax = 102000)
    f1.show()
    f2 = plt.figure()
    plt.pcolormesh(xnew,ynew,s,vmin = 96000, vmax = 102000)
    f2.show()

def get_conic_boxes(lows,data,lat,lon):
    """like get_boxes, but clips from a 
    small, conic-projected area
    """
    mylow = np.where(lows == 1)
    nlows = mylow[0].shape[0]
    dataregrid = np.zeros((nlows,21,21))
    count = 0
    for ind in range(0,nlows):
        time = mylow[0][ind]
        lowrow = mylow[1][ind]
        lowcol = mylow[2][ind]
        lattest = lat[lowrow,lowcol]
        lontest = lon[lowrow,lowcol]
        x,y=regrid_to_conic(lat,lon,lattest,lontest,lattest+5,lattest - 5)
        # area of interest is small region at center of regrid:
        xnew,ynew = np.meshgrid(np.arange(-0.05,.055,0.005),np.arange(-0.05,0.055,0.005))
        d = data[time,:,:]
        k = np.where(np.isnan(x) == False)
        s = interpolate.griddata((x[k],y[k]),d[k],(xnew,ynew),method = 'linear')
        # extra low filter for quality cont
        if s[10,10] < np.nanmean(s):
            dataregrid[count,:,:] = s
            count += count
    return dataregrid[0:count,:,:]

def plot_lows_on_map(lows,psl,time = 230):
    """plot_lows_on_map
    Function for quickly assessing the find_cyclone_center
    results. Makes a plot of sea level pressure with 
    identified low pressure centers marked

    Parameters:
    --------------------
    lows: 3D numpy array
    psl: 3D numpy array
    time (optional): numeric
    """
    lowsmap = lows[time,:,:]
    pslmap = psl[time,:,:]
    k = np.where(lowsmap == 1)
    f,ax = plt.subplots(1,1)
    ax.pcolormesh(pslmap, vmin = 90000, vmax = 104000)
    ax.scatter(k[1], k[0],color = 'k') 
    f.show()

def plot_box(box,time = 0):
    boxplot = box[time,:,:]
    f,axs = plt.subplots(1,1)
    h = axs.pcolormesh(boxplot,vmin = 90000, vmax = 104000)
    f.colorbar(h,ax = axs)
    f.show()

def plot_mean(data,cmin = 90000,cmax = 101000):
    f,axs = plt.subplots(1,1)
    h = axs.pcolormesh(np.nanmean(data,axis = 0),vmin = cmin, vmax = cmax)
    f.colorbar(h,ax = axs)
    f.show()


def get_mam(data):
    """Pulls out a timeseries only containing days in 
    March, April, and May. Input timeseries may not contain
    partial years of data.
    """
    nyrs = data.shape[0] / 365 
    mam = len(range(59,151))
    data_mam = np.zeros((nyrs*mam,data.shape[1],data.shape[2]))
    for yr in range(0,nyrs):
        if yr == 0:
            data_mam[0:mam,:,:] = data[59:151,:,:]
            last_ind = mam
        else:
            data_mam[last_ind:last_ind + mam,:,:] = data[59+(yr*365):151+(yr*365),:]
            last_ind = last_ind + mam
    return data_mam

def get_jja(data):
    """Pulls out a timeseries only containing days in 
    June, July, and August
    """
    nyrs = data.shape[0] / 365 
    jja = len(range(151,243))
    data_jja = np.zeros((nyrs*jja,data.shape[1],data.shape[2]))
    for yr in range(0,nyrs):
        if yr == 0:
            data_jja[0:jja,:,:] = data[151:243,:,:]
            last_ind = jja
        else:
            data_jja[last_ind:last_ind + jja,:,:] = data[151+(yr*365):243+(yr*365),:,:]
            last_ind = last_ind + jja
    return data_jja

def get_son(data):
    """Pulls out a timeseries only containing days in 
    September, October, and November
    """
    nyrs = data.shape[0] / 365 
    son = len(range(243,334))
    data_son = np.zeros((nyrs*son,data.shape[1],data.shape[2]))
    for yr in range(0,nyrs):
        if yr == 0:
            data_son[0:son,:,:] = data[243:334,:,:]
            last_ind = son
        else:
            data_son[last_ind:last_ind + son,:,:] = data[243+(yr*365):334+(yr*365),:,:]
            last_ind = last_ind + son
    return data_son

def get_djf(data):
    """Pulls out a timeseries only containing days in 
    December, January, and February
    """
    nyrs = data.shape[0] / 365 
    jf = len(range(0,60))
    d = len(range(334,365))
    data_djf = np.zeros((nyrs*(d+jf),data.shape[1],data.shape[2]))
    for yr in range(0,nyrs):
        if yr == 0:
            data_djf[0:60,:,:] = data[0:60,:,:]
            data_djf[60:60+d,:,:] =data[334:365,:,:]
            last_ind = 60+d
        else:
            data_djf[last_ind:last_ind + jf,:,:] = data[0+(yr*365):60+(yr*365),:]
            data_djf[last_ind + jf:last_ind + d + jf,:,:] = data[334+(yr*365):365+(yr*365),:]
            last_ind = last_ind + d + jf
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


