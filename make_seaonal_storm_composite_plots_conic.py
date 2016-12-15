from import_LE_data import *
from cyclone_composite_LENS import *
import matplotlib.pyplot as plt

print 'loading data'
psl = read_atm_data('PSL','001')
icefrac = read_atm_data('ICEFRAC','001')
ts = read_atm_data('TS','001')
u = read_atm_data('UBOT','001')
v = read_atm_data('VBOT','001')
#daidtd = read_ice_data('daidtd','001')
#daidtt = read_ice_data('daidtt','001')
lat,lon = read_native_lat_lon_atm()
tlat,tlon = read_native_lat_lon_ice()

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
#dtdma = get_anomaly_from_ma(daidtd,wgts)
#dttma = get_anomaly_from_ma(daidtt,wgts)


seasons = ['djf','mam','jja','son']
season_functions = {'djf':get_djf,'mam':get_mam,'jja':get_jja,'son':get_son}


print 'making composites'
for n in range(4,5):
    print 'set = ', str(n)
    start = 20*365*n 
    end = start + (20*365)

    tsn = ts[start:end,:,:]
    icen = icefrac[start:end,:,:]
    psln = psl[start:end,:,:]
    un = u[start:end,:,:]
    vn = v[start:end,:,:]
    #dtdn = dtdma[start:end,:,:]
    #dttn = dttma[start:end,:,:]

    for season in seasons:
        get_season = season_functions[season]
        pseas = get_season(psln)
        iseas = get_season(icen)
        tseas = get_season(tsn)
        useas = get_season(un)
        vseas = get_season(vn)
        #dttseas = get_season(dttn)
        #dtdseas = get_season(dtdn)
        lows = find_cyclone_center(pprojseas,iprojseas,lat,104000,90000)

        print lows.shape
        if (len(lows) < 1):
            continue
        elif (np.max(lows) == 1.0):
            data = {'psl':pseas,'icefrac':iseas,'ts':tseas,
                    'u':useas,'v':vseas,'dtd':dtdseas,'dtt':dttseas}
            latlist = {'psl':lat,'icefrac':lat,'ts':lat,
                       'u':lat,'v':lat,'dtd':tlat,'dtt':tlat}
                       'ice':tlat}
            lonlist = {'psl':lon,'icefrac':lon,'ts':lon,
                       'u':lon,'v':lon,'dtd':tlon,'dtt':tlon,
                       'ice':tlon}
            types = {'psl':'atm','icefrac':'atm','ts':'atm',
                     'u':'atm','v':'atm','dtd':'ice','dtt':'ice'}
            boxes = get_conic_boxes(lows,data,types,latliat,lonlist) 
            pbox = np.nanmean(boxes['psl'],axis = 0)
            ibox = np.nanmean(boxes['icefrac'],axis = 0)
            tbox = np.nanmean(boxes['ts'],axis = 0)
            ubox = np.nanmean(boxes['u'],axis = 0)
            vbox = np.nanmean(boxes['v'],axis = 0)
            
            
            f,axs = plt.subplots(1,1)
            h = axs.pcolormesh(X,Y,ibox,cmap='PuBu_r',
                                              vmin = 0,vmax = 1)
            axs.streamplot(X,Y,ubox,vbox,linewidth = 1)
            axs.contour(X,Y,tbox,range(240,300,5),colors = 'r')
            axs.contour(X,Y,pbox,range(96000,103000,100),colors = 'k')
            f.colorbar(h,ax = axs)
            f.savefig('test3_iceanom_' + str(n) + '_' + season + 'png')

            """f,axs = plt.subplots(1,1)
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
            f.savefig('test2_dtt_' + str(n) + '_' + season + 'png')"""




