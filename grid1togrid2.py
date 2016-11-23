def grid1togrid2(ongrid1,ncfile):

    griddata = Dataset(ncfile)
    N=griddata['dst_grid_dims'][:]
    Ntot=N[0]*N[1]

    dst=griddata['row'][:] # nasa stereo gridcell location
    src=griddata['col'][:] # gx1v6 pop ocn gridcell location
    S=griddata['S'][:] # map weight

    totwts=np.zeros((Ntot,1))
    totwts[dst]=totwts[dst]+S

    NN=ongrid1.shape

    if len(NN)<3: # just 2 D input 
        tmp = ongrid1.flatten()
        ongrid2=np.zeros((Ntot,1))
        ongrid2[dst]=S.*tmp[src] # multiply each element of the matrix (not matrix multiplication)
        ongrid2=ongrid2/totwts;  # divide each element of the matrix
        ongrid2=np.reshape(ongrid2,(N[0],N[1])) 
    else:
        ongrid2=np.zeros((NN[0],Ntot))
        for n=1:NN[0]:
            ongrid2[n,dst]=S*ongrid1[n,src]
        end
        ongrid2=ongrid2/(np.ones((NN([0],1))*totwts)
        ongrid2=np.reshape(ongrid2,(NN[0],N[0],N[1])) 
        ongrid2=np.transpose(ongrid2,(0 2 1)) 


