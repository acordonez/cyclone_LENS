import numpy as np
from netCDF4 import Dataset

def grid1togrid2(ongrid1,ncfile):
    """
    ncfile = '/glade/p/work/aordonez/cesm_mapping/map_gx1v6NH_TO_stereo25km_blin.161123.nc'

    ongrid1: numpy array on NH ice/ocean grid

    """    

    griddata = Dataset(ncfile)
    N=griddata.variables['dst_grid_dims'][:]

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

    Ntot=N[0]*N[1]
    NN=ongrid1.shape
    
    if len(NN) < 3:
        # weights are upside-down relative to how python reads netcdf
        tmp = ongrid1.flatten()
        tmp = np.reshape(tmp,(len(tmp),1))
        ongrid2=np.zeros((Ntot,1))
        # multiply each element of the matrix (not matrix multiplication)
        for ind in np.unique(dst):
            k = np.where(dst == ind)
            totwt = np.sum(S[k])
            ongrid2[ind]=np.sum(S[k]*tmp[src[k]]) / totwt
        ongrid2=np.reshape(ongrid2,(N[1],N[0]))
    else:
       time,y,x = ongrid1.shape
       tmp = np.reshape(ongrid1,(time,y*x))
       tmp = np.transpose(tmp,(1,0))
       ongrid2 = np.zeros((Ntot,time))
       for ind in np.unique(dst):
           k = np.where(dst == ind)
           totwt = np.sum(S[k])
           ongrid2[ind,:] = np.sum(S[k] * tmp[src[k],time]) / totwt
       ongrid2 = np.reshape(ongrid2,(N[1],N[0],time))

    return ongrid2






