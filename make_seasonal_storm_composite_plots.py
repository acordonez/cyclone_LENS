"""make_seasonal_storm_composite_plots.py
Makes storm composite plots for DJF,MAM,JJA, and SON in the Arctic. 
Should double check low ID criteria in cyclone_composite_LENS to make sure 
it is compositing over the right areas before running.

Try to filter storms that seem to be the same one as in a previous time step to get
parked storms - want day before storm arrived, day after it left
"""
import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from import_LE_data import *
from cyclone_composite_LENS import *
from SH_plots import get_difference, remove_daily_mean
import scipy.ndimage.filters as filters


tmpmask, names = read_region_mask_arctic()
rmask = np.zeros(tmpmask.shape)
rmask[(tmpmask > 10) & (tmpmask < 14)] = 1 #Beaufort, Chukchi, E Siberian
rmask[(tmpmask > 5) & (tmpmask < 9)] = 2 #Baffin, Greenland, Barents
rmask[(tmpmask == 15)] = 3 # central Arctic
rmask  = np.reshape(rmask,(1,rmask.shape[0],rmask.shape[1]))
reg = ['Chukchi','NA','Central']
rnum = 1


print 'loading data'
pproj = np.load('/glade/scratch/aordonez/PSL_034_proj_test_global_blin.npy')
iproj = np.load('/glade/scratch/aordonez/ICEFRAC_034_proj_test_global_blin.npy')
tproj = np.load('/glade/scratch/aordonez/TS_034_proj_test_global_blin.npy') 
uproj = np.load('/glade/scratch/aordonez/TAUX_034_proj_test_global_blin.npy') * -1
vproj = np.load('/glade/scratch/aordonez/TAUY_034_proj_test_global_blin.npy')  * -1 
shflx = np.load('/glade/scratch/aordonez/SHFLX_034_proj_test_global_blin.npy')
lhflx = np.load('/glade/scratch/aordonez/LHFLX_034_proj_test_global_blin.npy')
fsns = np.load('/glade/scratch/aordonez/FSNS_034_proj_test_global_blin.npy')
dtdproj = np.load('/glade/scratch/aordonez/daidtd_034_proj_blin.npy') 
dttproj = np.load('/glade/scratch/aordonez/daidtt_034_proj_blin.npy')
dvtproj = np.load('/glade/scratch/aordonez/dvidtd_034_proj_blin.npy') 
dvdproj = np.load('/glade/scratch/aordonez/dvidtt_034_proj_blin.npy')
mproj = np.load('/glade/scratch/aordonez/meltt_034_proj_blin.npy')
fproj = np.load('/glade/scratch/aordonez/frazil_034_proj_blin.npy')
cproj = np.load('/glade/scratch/aordonez/congel_034_proj_blin.npy')
hblt = np.load('/glade/scratch/aordonez/HBLT_2_034_proj_blin.npy')
sproj = np.load('/glade/scratch/aordonez/SST_005_proj_blin.npy')
#dvtproj = 0 * hblt
#dvdproj = 0 * hblt
#dtdproj = 0*hblt
#dttproj = 0*hblt
#mproj = 0 * hblt
lat,lon = read_stereo_lat_lon()


print 'Calculating anomalies'
# get ice as anomaly from 5-day mean
wgts = [1./5,1./5,1./5,1./5,1./5,]
def get_anomaly_from_ma(data,wgts):
    n = len(wgts)
    half = n
    datama = np.zeros(data.shape)
    for ind in range(half,data.shape[0]-half):
        n0 = ind - half
        nf = ind
        datama[ind,:,:] = np.average(data[n0:nf,:,:],axis = 0,weights = wgts)
    datama = data- datama
    data[0:half,:,:] = 0
    data[(data.shape[0]-half):,:,:] = 0
    return datama
#dtdma = get_anomaly_from_ma(dtdproj,wgts)
#dttma = get_anomaly_from_ma(dttproj,wgts)
#dvdma = get_anomaly_from_ma(dvdproj,wgts)
#dvtma = get_anomaly_from_ma(dvtproj,wgts)
#mltma = get_anomaly_from_ma(mproj,wgts)
#iaprojma = get_anomaly_from_ma(iproj,wgts)
dtdma = (dtdproj)
dttma = (dttproj)
dvdma = (dvdproj)
dvtma = (dvtproj)
mltma = (mproj)
iaprojma = get_difference(iproj)
#tma = get_diff_withseasonalerence(tproj)
#hbltma = get_diff_withseasonalerence(hblt)
tma = tproj
hbltma = get_difference(hblt)

seasons = ['djf','mam','jja','son']
season_functions = {'djf':get_djf,'mam':get_mam,'jja':get_jja,'son':get_son}

print 'making composites'
for n in [0,4]:
    print 'set = ', str(n)
    start = 20*365*n 
    end = start + (20*365)

    tprojn = tproj[start:end,:,:]
    iprojn = iproj[start:end,:,:]
    iaprojn = iaprojma[start:end,:,:] - get_nday_trend(iprojn,5)
    print 'got iceanom'
    pprojn = pproj[start:end,:,:]
    #panom = pproj[start:end,:,:] - trend_predict(pproj[start:end,:,:],5)
    uprojn = uproj[start:end,:,:]
    vprojn = vproj[start:end,:,:]
    shn = remove_daily_mean(shflx[start:end,:,:])
    lhn = remove_daily_mean(lhflx[start:end,:,:])
    fsn = remove_daily_mean(fsns[start:end,:,:])
    print 'fsn'
    dtdn = dtdma[start:end,:,:] - trend_predict(dtdma[start:end,:,:],5)
    dttn = dttma[start:end,:,:] - trend_predict(dttma[start:end,:,:],5)
    dvdn = dvdma[start:end,:,:] - trend_predict(dvdma[start:end,:,:],5)
    dvtn = dvtma[start:end,:,:] - trend_predict(dvdma[start:end,:,:],5)
    fn = fproj[start:end,:,:] - trend_predict(fproj[start:end,:,:],5)
    cn = cproj[start:end,:,:] - trend_predict(cproj[start:end,:,:],5)
    print 'cn'
    mltn = mltma[start:end,:,:] - trend_predict(mltma[start:end,:,:],5)
    print 1
    tn = tma[start:end,:,:] - trend_predict(tma[start:end,:,:],5)
    print 2
    #hbltn = hbltma[start:end,:,:] - trend_predict(hbltma[start:end,:,:],5)
    print 3
    sn = sproj[start:end,:,:] - trend_predict(sproj[start:end,:,:],5)
    print 4
    #iaprojn = iaprojma[start:end,:,:]
    #dtdn = dtdma[start:end,:,:]
    #dttn = dttma[start:end,:,:]
    #dvdn = dvdma[start:end,:,:]
    #dvtn = dvtma[start:end,:,:]
    #mltn = mltma[start:end,:,:]
    #tn = tma[start:end,:,:]
    #hbltn = hbltma[start:end,:,:]

    for season in ['jja']:
        print season
        get_season = season_functions[season]
        pprojseas = get_season(pprojn)
        #panomseas = get_season(panom)
        iprojseas = get_season(iprojn)
        iaprojseas = get_season(iaprojn)
        tprojseas = get_season(tprojn)
        uprojseas = get_season(uprojn)
        vprojseas = get_season(vprojn)
        shseas = get_season(shn)
        lhseas = get_season(lhn)
        fsseas = get_season(fsn)
        dttseas = get_season(dttn)
        dtdseas = get_season(dtdn)
        dvtseas = get_season(dvtn)
        dvdseas = get_season(dvdn) 
        fseas = get_season(fn)
        cseas = get_season(cn)
        mseas = get_season(mltn)
        tdifseas = get_season(tn)
        #hbltseas = get_season(hbltn)
        sstseas = get_season(sn)

        lowss = find_cyclone_center(pprojseas,iprojseas,107000,90000)
        lowss = remove_parked_lows(lowss)
        #lowss[panomseas > -5] = 0
        print lowss.shape

        if (len(lowss) < 1):
            continue
        elif (np.max(lowss) == 1.0):
            #k = np.where(rmask == rnum)
            #mask = np.zeros(rmask.shape)
            #mask[k] = 1
            #lows = lowss * mask
            lows = lowss

            pprojseas[np.isnan(pprojseas)]= -1000
            tprojseas[np.isnan(tprojseas)] = -1000
            iprojseas[np.isnan(iprojseas)]= -1000
            iaprojseas[np.isnan(iaprojseas)] = -1000
            uprojseas[np.isnan(uprojseas)] = -1000
            vprojseas[np.isnan(vprojseas)] = -1000
            shseas[np.isnan(shseas)]= -1000
            lhseas[np.isnan(lhseas)] = -1000
            fsseas[np.isnan(fsseas)] = -1000
            dtdseas[np.isnan(dtdseas)]= -1000
            dttseas[np.isnan(dttseas)] = -1000
            dvdseas[np.isnan(dvdseas)] = -1000
            dvtseas[np.isnan(dvtseas)] = -1000
            fseas[np.isnan(fseas)] = -1000
            cseas[np.isnan(cseas)] = -1000
            mseas[np.isnan(mseas)] = -1000
            #hbltseas[np.isnan(hbltseas)] = -1000
            tdifseas[np.isnan(tdifseas)] = -1000
            sstseas[np.isnan(sstseas)] = -1000

            box,_ = get_boxes(lows,pprojseas,50,lat,lon,90000) 
            # replace no data and absurd values with nan:
            #box[box == 0.0] = np.nan
            box[box > 108000] = np.nan   
            box[box < 90000] = np.nan

            tbox,_=get_boxes(lows,tprojseas,50,lat,lon,900)
            #tbox[tbox == 0.0] = np.nan
            tbox[tbox < 180.] = np.nan
            tbox[tbox > 360.] = np.nan

            ibox,_  = get_boxes(lows,iprojseas,50,lat,lon,900)
            ibox[ibox <= 0.0] = np.nan
            ibox[ibox > 110] = np.nan
            ibox[ibox > 100] = 100.

            iabox,_ = get_boxes(lows,iaprojseas,50,lat,lon,900)
            iabox[iabox >= 100] = np.nan
            iabox[iabox <= -100] = np.nan

            shbox,_  = get_boxes(lows,shseas,50,lat,lon,500)
            shbox[shbox > 500] = np.nan
            shbox[shbox < -500] = np.nan
            lhbox,_  = get_boxes(lows,lhseas,50,lat,lon,500) 
            lhbox[lhbox > 500] = np.nan
            lhbox[lhbox < -500] = np.nan
            fsbox,_  = get_boxes(lows,fsseas,50,lat,lon,500)
            fsbox[fsbox > 500] = np.nan
            fsbox[fsbox < -500] = np.nan
            dttbox,_ = get_boxes(lows,dttseas,50,lat,lon,500)
            dttbox[abs(dttbox) > 500] = np.nan

            dtdbox,_ = get_boxes(lows,dtdseas,50,lat,lon,500)
            dtdbox[abs(dtdbox) > 500] = np.nan

            dvdbox,_ = get_boxes(lows,dvdseas,50,lat,lon,500)
            dvdbox[abs(dvdbox) > 500] = np.nan

            dvtbox,_ = get_boxes(lows,dvtseas,50,lat,lon,500)
            dvtbox[abs(dvtbox) > 500] = np.nan

            fbox,_ = get_boxes(lows,fseas,50,lat,lon,500)
            fbox[abs(fbox) > 500] = np.nan
    
            cbox,_= get_boxes(lows,cseas,50,lat,lon,500)
            cbox[abs(cbox) > 500] = np.nan

            mbox,_ = get_boxes(lows,mseas,50,lat,lon,500)
            mbox[abs(mbox) > 500] = np.nan

            #hbox,_ = get_boxes(lows,hbltseas,50,lat,lon,500)
            #hbox[abs(hbox) > 900] = np.nan

            sstbox,_ = get_boxes(lows,sstseas,50,lat,lon,500)
            sstbox[abs(sstbox) > 100] = np.nan

            tdifbox,_=get_boxes(lows,tdifseas,50,lat,lon,5000)
            tdifbox[abs(tdifbox) > 360.] = np.nan

            ubox,_ = get_boxes(lows,uprojseas,50,lat,lon,500)
            vbox,_ = get_boxes(lows,vprojseas,50,lat,lon,500)
            #ubox[ubox == 0.0] = np.nan
            #vbox[vbox == 0.0] = np.nan
            ubox[abs(ubox) > 500] = np.nan
            vbox[abs(vbox) > 500] = np.nan

            # remove samples where SLP at 'low' is greater than box average:
            boxslice = box[:,41:60,41:60]
            pmean = np.nanmean(np.nanmean(boxslice,axis = 2),axis = 1)
            k = np.where(box[:,50,50] < pmean)
            box = box[k[0],:,:]
            tbox = tbox[k[0],:,:]
            ibox = ibox[k[0],:,:]
            iabox = iabox[k[0],:,:]
            ubox = ubox[k[0],:,:]
            vbox = vbox[k[0],:,:]
            shbox = shbox[k[0],:,:]
            lhbox = lhbox[k[0],:,:]
            fsbox = fsbox[k[0],:,:]
            dtdbox = dtdbox[k[0],:,:]
            dttbox = dttbox[k[0],:,:]
            dvtbox = dvtbox[k[0],:,:]
            dvdbox = dvdbox[k[0],:,:]
            fbox = fbox[k[0],:,:]
            cbox = cbox[k[0],:,:]
            mbox = mbox[k[0],:,:]
            sstbox = sstbox[k[0],:,:]
            #hbox[k[0],:,:]
            tboxdif = tdifbox[k[0],:,:]

            ubox = np.nanmean(ubox, axis = 0)
            vbox = np.nanmean(vbox, axis = 0)

            tbox_filt = filters.gaussian_filter(np.nanmean(tbox,axis = 0),3)
            pbox_filt = filters.gaussian_filter(np.nanmean(box, axis = 0),3)

            X,Y = np.meshgrid(range(0,101),range(0,101))

            #plt.close('all')
            # ice concentration
            f,axs = plt.subplots(1,1)
            h = axs.pcolormesh(X,Y,np.nanmean(ibox*100,axis = 0),cmap='PuBu_r',
                                              vmin = 0,vmax = 100)
            axs.streamplot(X,Y,ubox,vbox,linewidth = 1)
            cs = axs.contour(X,Y,tbox_filt,range(240,300,1),colors = 'r')
            plt.clabel(cs,inline = 1,fontsize = 12)
            cs = axs.contour(X,Y,pbox_filt,range(90000,107000,100),colors = 'k')
            plt.clabel(cs,inline = 1,fontsize = 12)
            f.colorbar(h,ax = axs)
            f.savefig('NH_034_ice_' + str(n) + '_' + season + '_dailyanom_nopark_detrend.png')

            # ice concentration anomaly
            cmax = np.nanmax(abs(np.nanmean(iabox*100, axis = 0)))
            f,axs = plt.subplots(1,1)
            h = axs.pcolormesh(X,Y,np.nanmean(iabox*100,axis = 0),cmap='PuOr',
                                              vmin = -1*cmax,vmax = cmax)
            axs.streamplot(X,Y,ubox,vbox,linewidth = 1)
            cs = axs.contour(X,Y,tbox_filt,range(240,300,1),colors = 'r')
            plt.clabel(cs,inline = 1,fontsize = 12)
            cs = axs.contour(X,Y,pbox_filt,range(90000,107000,100),colors = 'k')
            plt.clabel(cs,inline = 1,fontsize = 12)
            f.colorbar(h,ax = axs)
            f.savefig('NH_034_iceanom_' + str(n) + '_' + season + '_diff_nopark_detrend.png')

            # mean zonal wind stress
            cmax = np.nanmax(abs(ubox))
            f,axs = plt.subplots(1,1)
            h = axs.pcolormesh(X,Y,ubox,cmap='PuOr',vmin = -1*cmax,vmax = cmax)
            axs.streamplot(X,Y,ubox,vbox,linewidth = 1)
            cs = axs.contour(X,Y,tbox_filt,range(240,300,1),colors = 'r')
            plt.clabel(cs,inline = 1,fontsize = 12)
            cs = axs.contour(X,Y,pbox_filt,range(90000,107000,100),colors = 'k')
            plt.clabel(cs,inline = 1,fontsize = 12)
            f.colorbar(h,ax = axs)
            f.savefig('NH_034_taux_' + str(n) + '_' + season + '_mean_nopark_detrend.png')


            # mean meridional wind stress
            cmax = np.nanmax(abs(vbox))
            f,axs = plt.subplots(1,1)
            h = axs.pcolormesh(X,Y,vbox,cmap='PuOr',vmin = -1*cmax,vmax = cmax)
            axs.streamplot(X,Y,ubox,vbox,linewidth = 1)
            cs = axs.contour(X,Y,tbox_filt,range(240,300,1),colors = 'r')
            plt.clabel(cs,inline = 1,fontsize = 12)
            cs = axs.contour(X,Y,pbox_filt,range(90000,107000,100),colors = 'k')
            plt.clabel(cs,inline = 1,fontsize = 12)
            f.colorbar(h,ax = axs)
            f.savefig('NH_034_tauy_' + str(n) + '_' + season + '_mean_nopark_detrend.png')


            # sensible heat flux anomaly
            cmax = np.nanmax(abs(np.nanmean(shbox, axis = 0)))
            f,axs = plt.subplots(1,1)
            h = axs.pcolormesh(X,Y,np.nanmean(shbox,axis = 0),cmap='PuOr',
                                              vmin = -1*cmax,vmax = cmax)
            axs.streamplot(X,Y,ubox,vbox,linewidth = 1)
            cs = axs.contour(X,Y,tbox_filt,range(240,300,1),colors = 'r')
            plt.clabel(cs,inline = 1,fontsize = 12)
            cs = axs.contour(X,Y,pbox_filt,range(90000,107000,100),colors = 'k')
            plt.clabel(cs,inline = 1,fontsize = 12)
            f.colorbar(h,ax = axs)
            f.savefig('NH_034_shflx_' + str(n) + '_' + season + '_diff_nopark_detrend.png')

            # latent heat flux anomaly
            cmax = np.nanmax(abs(np.nanmean(lhbox, axis = 0)))
            f,axs = plt.subplots(1,1)
            h = axs.pcolormesh(X,Y,np.nanmean(lhbox,axis = 0),cmap='PuOr',
                                              vmin = -1*cmax,vmax = cmax)
            axs.streamplot(X,Y,ubox,vbox,linewidth = 1)
            cs = axs.contour(X,Y,tbox_filt,range(240,300,1),colors = 'r')
            plt.clabel(cs,inline = 1,fontsize = 12)
            cs = axs.contour(X,Y,pbox_filt,range(90000,107000,100),colors = 'k')
            plt.clabel(cs,inline = 1,fontsize = 12)
            f.colorbar(h,ax = axs)
            f.savefig('NH_034_lhflx_' + str(n) + '_' + season + '_diff_nopark_detrend.png')

            # total heat flux anomaly
            cmax = np.nanmax(abs(np.nanmean(shbox+lhbox, axis = 0)))
            f,axs = plt.subplots(1,1)
            h = axs.pcolormesh(X,Y,np.nanmean(shbox+lhbox,axis = 0),cmap='PuOr',
                                              vmin = -1*cmax,vmax = cmax)
            axs.streamplot(X,Y,ubox,vbox,linewidth = 1)
            cs = axs.contour(X,Y,tbox_filt,range(240,300,1),colors = 'r')
            plt.clabel(cs,inline = 1,fontsize = 12)
            cs = axs.contour(X,Y,pbox_filt,range(90000,107000,100),colors = 'k')
            plt.clabel(cs,inline = 1,fontsize = 12)
            f.colorbar(h,ax = axs)
            f.savefig('NH_034_sfcflx_' + str(n) + '_' + season + '_diff_nopark_detrend.png')

            # net surface shortwave flux anomaly
            cmax = np.nanmax(abs(np.nanmean(fsbox, axis = 0)))
            f,axs = plt.subplots(1,1)
            h = axs.pcolormesh(X,Y,np.nanmean(fsbox,axis = 0),cmap='PuOr',
                                              vmin = -1*cmax,vmax = cmax)
            axs.streamplot(X,Y,ubox,vbox,linewidth = 1)
            cs = axs.contour(X,Y,tbox_filt,range(240,300,1),colors = 'r')
            plt.clabel(cs,inline = 1,fontsize = 12)
            cs = axs.contour(X,Y,pbox_filt,range(90000,107000,100),colors = 'k')
            plt.clabel(cs,inline = 1,fontsize = 12)
            f.colorbar(h,ax = axs)
            f.savefig('NH_034_fsns_' + str(n) + '_' + season + '_diff_nopark_detrend.png')

            """# mixed layer depth anomaly
            cmax = np.nanmax(abs(np.nanmean(hbox / 100, axis = 0)))
            f,axs = plt.subplots(1,1)
            h = axs.pcolormesh(X,Y,np.nanmean(hbox/100,axis = 0),cmap='PuOr',
                                              vmin = -1*cmax,vmax = cmax)
            axs.streamplot(X,Y,ubox,vbox,linewidth = 1)
            cs = axs.contour(X,Y,tbox_filt,range(240,300,1),colors = 'r')
            plt.clabel(cs,inline = 1,fontsize = 12)
            cs = axs.contour(X,Y,pbox_filt,range(90000,107000,100),colors = 'k')
            plt.clabel(cs,inline = 1,fontsize = 12)
            f.colorbar(h,ax = axs)
            f.savefig('NH_034_hblt_' + str(n) + '_' + season + '_diff_nopark_detrend.png')"""

            # dynamic area tendency anomaly
            cmax = np.nanmax(abs(np.nanmean(dtdbox, axis = 0)))
            f,axs = plt.subplots(1,1)
            h = axs.pcolormesh(X,Y,np.nanmean(dtdbox,axis = 0),cmap='PuOr',
                                              norm=colors.SymLogNorm(linthresh=0.01,
                                              linscale=0.1,vmin=-1*cmax,vmax=cmax))
            axs.streamplot(X,Y,ubox,vbox,linewidth = 1)
            cs = axs.contour(X,Y,tbox_filt,range(240,300,1),colors = 'r')
            plt.clabel(cs,inline = 1,fontsize = 12)
            cs = axs.contour(X,Y,pbox_filt,range(90000,107000,100),colors = 'k')
            plt.clabel(cs,inline = 1,fontsize = 12)
            f.colorbar(h,ax = axs)
            f.savefig('NH_034_dtd_' + str(n) + '_' + season + '_dailyanom_nopark_detrend.png')

            # thermodynamic area tendency anomaly
            cmax = np.nanmax(abs(np.nanmean(dttbox, axis = 0)))
            f,axs = plt.subplots(1,1)
            h = axs.pcolormesh(X,Y,np.nanmean(dttbox,axis = 0),cmap='PuOr',
                                              norm=colors.SymLogNorm(linthresh=0.01,
                                              linscale=0.1,vmin=-1*cmax,vmax=cmax))
            axs.streamplot(X,Y,ubox,vbox,linewidth = 1)
            cs = axs.contour(X,Y,tbox_filt,range(240,300,1),colors = 'r')

            plt.clabel(cs,inline = 1,fontsize = 12)
            cs = axs.contour(X,Y,pbox_filt,range(90000,107000,100),colors = 'k')
            plt.clabel(cs,inline = 1,fontsize = 12)
            f.colorbar(h,ax = axs)
            f.savefig('NH_034_dtt_' + str(n) + '_' + season + '_dailyanom_nopark_detrend.png')
 
            # dynamic volume tendency anomaly
            cmax = np.nanmax(abs(np.nanmean(dvdbox, axis = 0)))
            f,axs = plt.subplots(1,1)
            h = axs.pcolormesh(X,Y,np.nanmean(dvdbox,axis = 0),cmap='PuOr',
                                              norm=colors.SymLogNorm(linthresh=0.01,
                                              linscale=0.1,vmin=-1*cmax,vmax=cmax))
            axs.streamplot(X,Y,ubox,vbox,linewidth = 1)
            cs = axs.contour(X,Y,tbox_filt,range(240,300,1),colors = 'r')
            plt.clabel(cs,inline = 1,fontsize = 12)
            cs = axs.contour(X,Y,pbox_filt,range(90000,107000,100),colors = 'k')
            plt.clabel(cs,inline = 1,fontsize = 12)
            f.colorbar(h,ax = axs)
            f.savefig('NH_034_dvidtd_' + str(n) + '_' + season + '_dailyanom_nopark_detrend.png')

            # thermodynamic volume tendency anomaly
            cmax = np.nanmax(abs(np.nanmean(dvtbox, axis = 0)))
            f,axs = plt.subplots(1,1)
            h = axs.pcolormesh(X,Y,np.nanmean(dvtbox,axis = 0),cmap='PuOr',
                                              norm=colors.SymLogNorm(linthresh=0.01,
                                              linscale=0.1,vmin=-1*cmax,vmax=cmax))
            axs.streamplot(X,Y,ubox,vbox,linewidth = 1)
            cs = axs.contour(X,Y,tbox_filt,range(240,300,1),colors = 'r')
            plt.clabel(cs,inline = 1,fontsize = 12)
            cs = axs.contour(X,Y,pbox_filt,range(90000,107000,100),colors = 'k')
            plt.clabel(cs,inline = 1,fontsize = 12)
            f.colorbar(h,ax = axs)
            f.savefig('NH_034_dvidtt_' + str(n) + '_' + season + '_dailyanom_nopark_detrend.png')

            # bottom melt anomaly
            cmax = np.nanmax(abs(np.nanmean(mbox, axis = 0)))
            f,axs = plt.subplots(1,1)
            h = axs.pcolormesh(X,Y,np.nanmean(mbox,axis = 0),cmap='PuOr',
                                              norm=colors.SymLogNorm(linthresh=0.01,
                                              linscale=0.1,vmin=-1*cmax,vmax=1*cmax))
            axs.streamplot(X,Y,ubox,vbox,linewidth = 1)
            cs = axs.contour(X,Y,tbox_filt,range(240,300,5),colors = 'r')
            plt.clabel(cs,inline = 1,fontsize = 12)
            cs = axs.contour(X,Y,pbox_filt,range(90000,107000,100),colors = 'k')
            plt.clabel(cs,inline = 1,fontsize = 12)
            f.colorbar(h,ax = axs)
            f.savefig('NH_034_meltt_' + str(n) + '_' + season + '_dailyanom_nopark_detrend.png')

            # frazil anomaly
            cmax = np.nanmax(abs(np.nanmean(fbox, axis = 0)))
            f,axs = plt.subplots(1,1)
            h = axs.pcolormesh(X,Y,np.nanmean(fbox,axis = 0),cmap='PuOr',
                                              norm=colors.SymLogNorm(linthresh=0.01,
                                              linscale=0.1,vmin=-1*cmax,vmax=1*cmax))
            axs.streamplot(X,Y,ubox,vbox,linewidth = 1)
            cs = axs.contour(X,Y,tbox_filt,range(240,300,5),colors = 'r')
            plt.clabel(cs,inline = 1,fontsize = 12)
            cs = axs.contour(X,Y,pbox_filt,range(90000,107000,100),colors = 'k')
            plt.clabel(cs,inline = 1,fontsize = 12)
            f.colorbar(h,ax = axs)
            f.savefig('NH_034_frazil_' + str(n) + '_' + season + '_dailyanom_nopark_detrend.png')

            # congelation anomaly
            cmax = np.nanmax(abs(np.nanmean(cbox, axis = 0)))
            f,axs = plt.subplots(1,1)
            h = axs.pcolormesh(X,Y,np.nanmean(cbox,axis = 0),cmap='PuOr',
                                              norm=colors.SymLogNorm(linthresh=0.01,
                                              linscale=0.1,vmin=-1*cmax,vmax=1*cmax))
            axs.streamplot(X,Y,ubox,vbox,linewidth = 1)
            cs = axs.contour(X,Y,tbox_filt,range(240,300,5),colors = 'r')
            plt.clabel(cs,inline = 1,fontsize = 12)
            cs = axs.contour(X,Y,pbox_filt,range(90000,107000,100),colors = 'k')
            plt.clabel(cs,inline = 1,fontsize = 12)
            f.colorbar(h,ax = axs)
            f.savefig('NH_034_congel_' + str(n) + '_' + season + '_dailyanom_nopark_detrend.png')


            #surface temperature anomaly
            f,axs = plt.subplots(1,1)
            h = axs.pcolormesh(X,Y,np.nanmean(tdifbox,axis = 0),cmap='PuBu_r',
                                              vmin = -5,vmax = 5)
            axs.streamplot(X,Y,ubox,vbox,linewidth = 1)
            cs = axs.contour(X,Y,tbox_filt,range(240,300,1),colors = 'r')
            plt.clabel(cs,inline = 1,fontsize = 12)
            cs = axs.contour(X,Y,pbox_filt,range(90000,107000,100),colors = 'k')
            plt.clabel(cs,inline = 1,fontsize = 12)
            f.colorbar(h,ax = axs)
            f.savefig('NH_034_ts_' + str(n) + '_' + season + '_dailyanom_nopark_detrend.png')


            #surface temperature anomaly
            f,axs = plt.subplots(1,1)
            cmax = np.nanmax(abs(np.nanmean(sstbox, axis = 0)))
            h = axs.pcolormesh(X,Y,np.nanmean(sstbox,axis = 0),cmap='PuBu_r',
                                              vmin = -1*cmax,vmax = cmax)
            axs.streamplot(X,Y,ubox,vbox,linewidth = 1)
            cs = axs.contour(X,Y,tbox_filt,range(240,300,1),colors = 'r')
            plt.clabel(cs,inline = 1,fontsize = 12)
            cs = axs.contour(X,Y,pbox_filt,range(90000,107000,100),colors = 'k')
            plt.clabel(cs,inline = 1,fontsize = 12)
            f.colorbar(h,ax = axs)
            f.savefig('NH_005_sst_' + str(n) + '_' + season + '_dailyanom_nopark_detrend.png')

            plt.close('all')




