import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from import_LE_data import *
from cyclone_composite_LENS import *

# more or bigger storms?
# make bullet list of points - have storms, number of storms change?
# lance bosart - papers on big storm fueled by sea ice loss NAAS workshop video

# True: Northern Hemisphere, False: Southern Hemisphere
nh = True
print 'NH = ',nh
print 'loading data'

psl = read_atm_data('PSL','002')
icefrac = read_atm_data('ICEFRAC','002')
ts = read_atm_data('TS','002')
u = read_atm_data('TAUX','002') * -1
v = read_atm_data('TAUY','002') * -1
if nh:
    ice_func = read_ice_data
    ocn_func = read_ocn_data # FIX OCN COORD STUFF
    latlon_func = read_native_lat_lon_ice
    lows_func = find_cyclone_center
else:
    ice_func = read_ice_data_SH
    ocn_func = read_ocn_data_SH
    latlon_func = read_native_lat_lon_ice_SH
    lows_func = find_cyclone_center_SH
daidtd = ice_func('daidtd','002')
daidtt = ice_func('daidtt','002')
daidtd[daidtd  > 1e20] = np.nan
daidtt[daidtt > 1e20] = np.nan
sst = ocn_func('SST','002')
sst[sst > 1e20] = np.nan
lat,lon = read_native_lat_lon_atm()
tlat,tlon = latlon_func()

wgts = [1./5,1./5,1./5,1./5,1./5,]
def get_anomaly_from_ma(data,wgts):
    n = len(wgts)
    half = n/2
    datama = np.zeros(data.shape)
    for ind in range(half,data.shape[0]-half):
        n0 = ind - half
        nf = ind + half + 1
        datama[ind,:,:] = np.average(data[n0:nf,:,:],axis = 0,weights = wgts)
    datama = datama - data
    data[0:half,:,:] = 0
    data[(data.shape[0]-half):,:,:] = 0
    return datama
dtdma = get_anomaly_from_ma(daidtd,wgts)
dttma = get_anomaly_from_ma(daidtt,wgts)
#dttma = daidtt
#dtdma = daidtd
iceanom = get_anomaly_from_ma(icefrac,wgts)
sstanom = get_anomaly_from_ma(sst,wgts)


seasons = ['jja']
season_functions = {'djf':get_djf,'mam':get_mam,'jja':get_jja,'son':get_son}
plows = np.zeros((5,4,300))
numlist = {'djf':0,'mam':1,'jja':2,'son':3}

print 'making composites'
for n in [4]:
    print 'set = ', str(n)
    start = 0 + 365 * 20 * n 
    end = start + 20 *365

    tsn = ts[start:end,:,:]
    ian = iceanom[start:end,:,:]
    #ishiftn = iceanom[start + 1:end + 1,:,:]
    icen = icefrac[start:end,:,:]
    psln = psl[start:end,:,:]
    un = u[start:end,:,:]
    vn = v[start:end,:,:]
    dtdn = dtdma[start:end,:,:]
    dttn = dttma[start:end,:,:]
    sstn = sstanom[start:end,:,:]
    print 'A'
    for season in seasons:
        get_season = season_functions[season]
        pseas = get_season(psln)
        iseas = get_season(icen)
        iaseas = get_season(ian)
        #isseas = get_season(ishiftn)
        tseas = get_season(tsn)
        useas = get_season(un)
        vseas = get_season(vn)
        dttseas = get_season(dttn)
        dtdseas = get_season(dtdn)
        sstseas = get_season(sstn)
        
        # get small random sample of low pressure centers to speed analysis
        lows = find_cyclone_center(pseas,iseas,lat,104000,90000)
        k = np.where(lows == 1)
        kind = np.random.choice(k[0].shape[0],size = 2000)
        lows_new = np.zeros(lows.shape)
        lows_new[k[0][kind],k[1][kind],k[2][kind]] = 1.0
        print 'B'
        
        data = {'psl':pseas,'icefrac':iseas,'ts':tseas,
                'u':useas,'v':vseas,'dtd':dtdseas,'dtt':dttseas,
                'iceanom':iaseas,'ice':dtdseas,'sst':sstseas}
                #'ishift':isseas}
        latlist = {'psl':lat,'icefrac':lat,'ts':lat,
                   'u':lat,'v':lat,'dtd':tlat,'dtt':tlat,
                   'iceanom':lat,'ice':tlat,'sst':tlat}
                   #'ishift':tlat}
        lonlist = {'psl':lon,'icefrac':lon,'ts':lon,
                   'u':lon,'v':lon,'dtd':tlon,'dtt':tlon,
                   'iceanom':lon,'ice':tlon,'sst':tlon}
                   #'ishift':tlon}
        types = {'psl':'atm','icefrac':'atm','ts':'atm',
                 'u':'atm','v':'atm','dtd':'ice','dtt':'ice',
                 'iceanom':'atm','ice':'ice','sst':'ice'}
                 #'ishift':'ice'}
        boxes,_ = get_conic_boxes(lows_new,data,types,latlist,lonlist) 
        print 'C'
        pbox = np.nanmean(boxes['psl'],axis = 0)     
        ibox = np.nanmean(boxes['icefrac'],axis = 0)
        iabox = np.nanmean(boxes['iceanom'],axis = 0)
        #isbox = np.nanmean(boxes['ishift'],axis = 0)
        tbox = np.nanmean(boxes['ts'],axis = 0)
        ubox = np.nanmean(boxes['u'],axis = 0)
        vbox = np.nanmean(boxes['v'],axis = 0)
        sbox = np.nanmean(boxes['sst'],axis = 0)
        dtdbox = np.nanmean(boxes['dtd'],axis = 0)
        dttbox = np.nanmean(boxes['dtt'],axis = 0)

        #p = boxes['psl']
        #t,x,y = p.shape
        #num = numlist[season]
        #plows[n,num,0:t] = boxes['psl'][:,x/2,y/2] / 100.
            
        size = pbox.shape[0]
        X,Y = np.meshgrid(range(0,size),range(0,size))
        if nh:
            pre = 'marg_conic_002_'
        else:
            pre = 'marg_conic_002_SH_'
           
        # ice concentration
        f,axs = plt.subplots(1,1)
        h = axs.pcolormesh(X,Y,ibox,cmap='PuBu_r',
                                          vmin = 0,vmax = 1)
        axs.streamplot(X,Y,ubox,vbox,linewidth = 1)
        cs = axs.contour(X,Y,tbox,np.arange(240,300,2.5),colors = 'r')
        plt.clabel(cs,inline = 1,fontsize = 8)
        cs = axs.contour(X,Y,pbox/100.,range(96000,103000,100),colors = 'k')
        plt.clabel(cs,inline = 1,fontsize = 12)
        plt.title(season + ' ice concentration')
        f.colorbar(h,ax = axs)
        f.savefig(pre + 'icefrac_' + str(n) + '_' + season + 'png')
        plt.close() 

        # ice concentration anomalies
        f,axs = plt.subplots(1,1)
        h = axs.pcolormesh(X,Y,iabox,cmap='PuOr_r',vmin = -0.005,vmax = 0.005)
        axs.streamplot(X,Y,ubox*-1,vbox*-1,linewidth = 1)
        cs= axs.contour(X,Y,tbox,np.arange(240,300,2.5),colors = 'r',linewidth = 1.5)
        plt.clabel(cs,inline = 1,fontsize = 12)
        cs = axs.contour(X,Y,pbox,range(96000,103000,100),colors = 'k',linewidth = 1.5)
        plt.clabel(cs,inline = 1,fontsize = 12)
        f.colorbar(h,ax = axs)
        plt.title(season + ' ice concentration anomaly')
        f.savefig(pre + 'iceanom_' + str(n) + '_' + season + 'png')
        plt.close()

        # ice concentration anomalies day + 1
        """f,axs = plt.subplots(1,1)
        h = axs.pcolormesh(X,Y,isbox,cmap='PuOr_r',vmin = -0.005,vmax = 0.005)
        axs.streamplot(X,Y,ubox*-1,vbox*-1,linewidth = 1)
        cs= axs.contour(X,Y,tbox,np.arange(240,300,2.5),colors = 'r',linewidth = 1.5)
        plt.clabel(cs,inline = 1,fontsize = 12)
        cs = axs.contour(X,Y,pbox,range(96000,103000,100),colors = 'k',linewidth = 1.5)
        plt.clabel(cs,inline = 1,fontsize = 12)
        f.colorbar(h,ax = axs)
        plt.title(season + ' ice concentration anomaly')
        f.savefig(pre + 'iceanom_n+1_' + str(n) + '_' + season + 'png')
        """
        # dynamic tendency anomalies
        f,axs = plt.subplots(1,1)
        h = axs.pcolormesh(X,Y,dtdbox,cmap='PuOr',
                           norm=colors.SymLogNorm(linthresh=0.01,
                           linscale=0.01,vmin=-0.2,vmax=0.2))
        axs.streamplot(X,Y,ubox*-1,vbox*-1,linewidth = 1)
        cs = axs.contour(X,Y,tbox,np.arange(240,300,2.5),colors = 'r',linewidth = 1.5)
        plt.clabel(cs,inline = 1,fontsize = 12)
        cs = axs.contour(X,Y,pbox,range(96000,103000,100),colors = 'k',linewidth = 1.5)
        plt.clabel(cs,inline = 1,fontsize = 12)
        f.colorbar(h,ax = axs)
        plt.title(season + ' dynamic tendency anomaly')
        f.savefig(pre + 'dtd_' + str(n) + '_' + season + 'png')
        plt.close()

        # thermodynamic tendency anomalies
        f,axs = plt.subplots(1,1)
        h = axs.pcolormesh(X,Y,dttbox,cmap='PuOr',
                           norm=colors.SymLogNorm(linthresh=0.01,
                          linscale=0.01,vmin=-0.2,vmax=0.2))
        axs.streamplot(X,Y,ubox,vbox,linewidth = 1)
        cs = axs.contour(X,Y,tbox,np.arange(240,300,2.5),colors = 'r',linewidth = 1.5)
        plt.clabel(cs,inline = 1,fontsize = 12)
        cs = axs.contour(X,Y,pbox,range(96000,103000,100),colors = 'k',linewidth = 1.5)
        plt.clabel(cs,inline = 1,fontsize = 12)
        f.colorbar(h,ax = axs)
        plt.title(season + ' thermodynamic tendency anomaly')
        f.savefig(pre + 'dtt_' + str(n) + '_' + season + 'png')
        plt.close()

        # total tendency anomaly
        f,axs = plt.subplots(1,1)
        h = axs.pcolormesh(X,Y,dttbox+dtdbox,cmap='PuOr',
                           norm=colors.SymLogNorm(linthresh=0.01,
                           linscale=0.01,vmin=-0.2,vmax=0.2))
        axs.streamplot(X,Y,ubox,vbox,linewidth = 1)
        cs = axs.contour(X,Y,tbox,np.arange(240,300,2.5),colors = 'r',linewidth = 1.5)
        plt.clabel(cs,inline = 1,fontsize = 12)
        cs = axs.contour(X,Y,pbox,range(96000,103000,100),colors = 'k',linewidth = 1.5)
        plt.clabel(cs,inline = 1,fontsize = 12)
        f.colorbar(h,ax = axs)
        plt.title(season + ' total tendency anomaly')
        f.savefig(pre + 'tendency_' + str(n) + '_' + season + 'png')
        plt.close()

        # SST
        f,axs = plt.subplots(1,1)
        h = axs.pcolormesh(X,Y,sbox,cmap='PuOr',vmin = -0.03,vmax = 0.03)
        axs.streamplot(X,Y,ubox,vbox,linewidth = 1)
        cs = axs.contour(X,Y,tbox,np.arange(240,300,2.5),colors = 'r',linewidth = 1.5)
        plt.clabel(cs,inline = 1,fontsize = 12)
        cs = axs.contour(X,Y,pbox,range(96000,103000,100),colors = 'k',linewidth = 1.5)
        plt.clabel(cs,inline = 1,fontsize = 12)
        f.colorbar(h,ax = axs)
        plt.title(season + ' sea surface temperature anomaly')
        f.savefig(pre + 'sst_' + str(n) + '_' + season + 'png')
        plt.close()

