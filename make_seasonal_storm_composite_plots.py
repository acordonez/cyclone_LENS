from import_LE_data import *
from cyclone_composite_LENS import *
import matplotlib.pyplot as plt

print 'loading data'
pproj = np.load('/glade/scratch/aordonez/pproj.npy')
iproj = np.load('/glade/scratch/aordonez/iproj.npy')
tproj = np.load('/glade/scratch/aordonez/tproj.npy') 
uproj = np.load('/glade/scratch/aordonez/UBOTproj.npy')
vproj = np.load('/glade/scratch/aordonez/VBOTproj.npy')   
dtdproj = np.load('/glade/scratch/aordonez/daidtdproj.npy') 
dttproj = np.load('/glade/scratch/aordonez/daidttproj.npy')
lat,lon = read_stereo_lat_lon()

# get ice as anomaly from 5-day mean
wgts = [1./5,1./5,1./5,1./5,1./5,]
def get_anomaly_from_ma(data,wgts):
    n = len(wgts)
    half = n/2
    datama = np.zeros(data.shape)
    for ind in range(half,data.shape[0]-2):
        n0 = ind - half
        nf = ind + half + 1
        datama[ind,:,:] = np.average(data[n0:nf,:,:],axis = 0,weights = wgts)
    datama = datama - data
    data[0:2,:,:] = 0
    data[(data.shape[0]-2):,:,:] = 0
    return datama
dtdma = get_anomaly_from_ma(dtdproj,wgts)
dttma = get_anomaly_from_ma(dttproj,wgts)
iaprojma = get_anomaly_from_ma(iproj,wgts)

seasons = ['djf','mam','jja','son']
season_functions = {'djf':get_djf,'mam':get_mam,'jja':get_jja,'son':get_son}

print 'making composites'
for n in range(4,5):
    print 'set = ', str(n)
    start = 20*365*n 
    end = start + (20*365)

    tprojn = tproj[start:end,:,:]
    iprojn = iproj[start:end,:,:]
    iaprojn = iaprojma[start:end,:,:]
    pprojn = pproj[start:end,:,:]
    uprojn = uproj[start:end,:,:]
    vprojn = vproj[start:end,:,:]
    dtdn = dtdma[start:end,:,:]
    dttn = dttma[start:end,:,:]

    for season in seasons:
        get_season = season_functions[season]
        pprojseas = get_season(pprojn)
        iprojseas = get_season(iprojn)
        iaprojseas = get_season(iaprojn)
        tprojseas = get_season(tprojn)
        uprojseas = get_season(uprojn)
        vprojseas = get_season(vprojn)
        dttseas = get_season(dttn)
        dtdseas = get_season(dtdn)
        lows = find_cyclone_center(pprojseas,iprojseas,lat,104000,90000)
        print lows.shape
        if (len(lows) < 1):
            continue
        elif (np.max(lows) == 1.0):
            box,_,_ = get_boxes(lows,pprojseas,50,lat,lon,90000) 
            # replace no data and absurd values with nan:
            box[box == 0.0] = np.nan
            box[box > 108000] = np.nan   
            box[box < 90000] = np.nan

            tbox,_,_=get_boxes(lows,tprojseas,50,lat,lon,90000)
            tbox[tbox == 0.0] = np.nan
            tbox[tbox < 180.] = np.nan
            tbox[tbox > 360.] = np.nan

            ibox,_,_  = get_boxes(lows,iprojseas,50,lat,lon,10)
            ibox[ibox <= 0.0] = np.nan
            ibox[ibox > 110] = np.nan

            iabox,_,_ = get_boxes(lows,iaprojseas,50,lat,lon,10)
            iabox[iabox > 100] = np.nan
            iabox[iabox == 0.0] = np.nan
 
            dttbox,_,_ = get_boxes(lows,dttseas,50,lat,lon,5)
            dttbox[dttbox == 0.0] = np.nan
            dttbox[abs(dttbox) > 100] = np.nan

            dtdbox,_,_ = get_boxes(lows,dtdseas,50,lat,lon,5)
            dtdbox[dtdbox == 0.0] = np.nan
            dtdbox[abs(dtdbox) > 100] = np.nan

            ubox,_,_ = get_boxes(lows,uprojseas,50,lat,lon,500)
            vbox,_,_ = get_boxes(lows,vprojseas,50,lat,lon,500)
            ubox[ubox == 0.0] = np.nan
            vbox[vbox == 0.0] = np.nan
            ubox[abs(ubox) > 500] = np.nan
            vbox[abs(vbox) > 500] = np.nan

            # remove samples where SLP at 'low' is greater than box average:
            pmean = np.nanmean(np.nanmean(box,axis = 2),axis = 1)
            k = np.where(box[:,50,50] < pmean)
            box = box[k[0],:,:]
            tbox = tbox[k[0],:,:]
            ibox = ibox[k[0],:,:]
            iabox = iabox[k[0],:,:]
            ubox = ubox[k[0],:,:]
            vbox = vbox[k[0],:,:]
            dtdbox = dtdbox[k[0],:,:]
            dttbox = dttbox[k[0],:,:]

            X,Y = np.meshgrid(range(0,101),range(0,101))

            ubox = np.nanmean(ubox, axis = 0)
            vbox = np.nanmean(vbox, axis = 0)

            f,axs = plt.subplots(1,1)
            h = axs.pcolormesh(X,Y,np.nanmean(ibox,axis = 0),cmap='PuBu_r',
                                              vmin = 0,vmax = 1)
            axs.streamplot(X,Y,ubox,vbox,linewidth = 1)
            axs.contour(X,Y,np.nanmean(tbox,axis = 0),range(240,300,5),colors = 'r')
            axs.contour(X,Y,np.nanmean(box,axis = 0),range(96000,103000,100),colors = 'k')
            f.colorbar(h,ax = axs)
            f.savefig('test2_iceanom_' + str(n) + '_' + season + 'png')

            f,axs = plt.subplots(1,1)
            h = axs.pcolormesh(X,Y,np.nanmean(dtdbox,axis = 0),cmap='PuOr',
                                              norm=colors.SymLogNorm(linthresh=0.01,
                                              linscale=0.01,vmin=-0.5,vmax=0.5))
            axs.streamplot(X,Y,ubox,vbox,linewidth = 1)
            axs.contour(X,Y,np.nanmean(tbox,axis = 0),range(240,300,5),colors = 'r')
            axs.contour(X,Y,np.nanmean(box,axis = 0),range(96000,103000,100),colors = 'k')
            f.colorbar(h,ax = axs)
            f.savefig('test2_dtd_' + str(n) + '_' + season + 'png')

            f,axs = plt.subplots(1,1)
            h = axs.pcolormesh(X,Y,np.nanmean(dttbox,axis = 0),cmap='PuOr',
                                              norm=colors.SymLogNorm(linthresh=0.01,
                                              linscale=0.01,vmin=-0.5,vmax=0.5))
            axs.streamplot(X,Y,ubox,vbox,linewidth = 1)
            axs.contour(X,Y,np.nanmean(tbox,axis = 0),range(240,300,5),colors = 'r')
            axs.contour(X,Y,np.nanmean(box,axis = 0),range(96000,103000,100),colors = 'k')
            f.colorbar(h,ax = axs)
            f.savefig('test2_dtt_' + str(n) + '_' + season + 'png')





