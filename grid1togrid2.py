def grid1togrid2(ongrid1,ncfile):
    """
    ncfile = '/glade/p/work/aordonez/cesm_mapping/map_gx1v6NH_TO_stereo25km_blin.161123.nc'

 b = np.concatenate((np.zeros((280,320)),a),axis = 0)

    """    

    griddata = Dataset(ncfile)
    N=griddata.variables['dst_grid_dims'][:]
    Ntot=N[0]*N[1]

    # nasa stereo gridcell location
    dst=griddata.variables['row'][:] 
    # gx1v6 pop ocn gridcell location
    src=griddata.variables['col'][:] 
    # map weight
    S=griddata.variables['S'][:] 
    S = np.reshape(S,(len(S),1))

    # convert to python 0-indexing
    dst = dst - 1
    src = src - 1

    totwts=np.zeros((Ntot,1))
    totwts[dst]=totwts[dst]+S

    NN=ongrid1.shape

    if len(NN)<3: # just 2 D input 
        # weights are upside-down relative to how python reads netcdf
        tmp = np.flipud(ongrid1).flatten()
        tmp = np.reshape(tmp,(len(tmp),1))
        ongrid2=np.zeros((Ntot,1))
        # multiply each element of the matrix (not matrix multiplication)
        for ind in np.unique(dst):
            k = np.where(dst == ind)
            totwt = np.sum(S[k])
            ongrid2[ind]=np.sum(S[k]*tmp[src[k]]) / totwt
        ongrid2=np.reshape(ongrid2,(N[1],N[0]))
        plt.pcolormesh(ongrid2)
        plt.colorbar()
        plt.show()
    else:
        ongrid2=np.zeros((NN[0],Ntot))
        for n=1:NN[0]:
            ongrid2[n,dst]=S*ongrid1[n,src]
        end
        ongrid2=ongrid2/(np.ones((NN([0],1))*totwts)
        ongrid2=np.reshape(ongrid2,(NN[0],N[0],N[1])) 
        ongrid2=np.transpose(ongrid2,(0 2 1)) 


