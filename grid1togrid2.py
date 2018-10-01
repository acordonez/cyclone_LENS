from netCDF4 import Dataset
import numpy as np

def test_grid1togrid2_atm():
    from import_LE_data import read_atm_data
    psl_on_latlon = read_atm_data('PSL','001')
    # SCRIP regridding file with weights
    ncfile = '/glade/p/work/aordonez/cesm_mapping/map_fv0.9x1.25_TO_stereo25km_blin.170209.nc'
    psl_on_stereo = grid1togrid2(psl_on_latlon, ncfile)

def test_grid1togrid2_ice():
    from import_LE_data import read_ice_data
    aice_on_latlon = read_ice_data('aice','001')[0:4,:,:]
    # SCRIP regridding file with weights
    ncfile =  '/glade/p/work/aordonez/cesm_mapping/map_gx1v6_TO_stereo25km_blin.170210.nc'
    aice_on_stereo = grid1togrid2(aice_on_latlon, ncfile)

def grid1togrid2(ongrid1,ncfile):
    # based on matlab script by Cecilia Bitz

    if len(ongrid1.shape) < 3:
        ongrid1 = np.tile(ongrid1,(1,1,1))
    time,dim1,dim2 = ongrid1.shape

    griddata = Dataset(ncfile)
    N=griddata.variables['dst_grid_dims'][:]
    Ntot=N[0]*N[1]

    # nasa stereo gridcell location
    dst=griddata.variables['row'][:] 
    # gx1v6 pop ocn gridcell location
    src=griddata.variables['col'][:] 
    # map weight
    S=griddata.variables['S'][:] 

    # convert to 0 based index
    dst -= 1
    src -= 1

    totwts=np.zeros((Ntot))
    totwts[dst]=totwts[dst]+S

    NN=ongrid1.shape
 
    tmp = np.reshape(ongrid1,(time,dim1*dim2))
    ongrid2=np.zeros((time,Ntot))

    # get unique destination points
    dunique,dindices = np.unique(dst,return_index = True)
    # get the last valid index
    indN = np.where(dst == dunique[-1])[0][-1]
    # multiply each element of the matrix (not matrix multiplication)
    for n in range(0,len(dindices)):
        ind1 = dindices[n]
        if n == len(dindices)-1:
            ind2 = indN + 1
        else:
            ind2 = dindices[n+1]
        ongrid2[:,dunique[n]] = np.nansum(S[ind1:ind2]*tmp[:,src[ind1:ind2]],axis = 1)

    ongrid2 = np.reshape(ongrid2,(time,N[1],N[0]))

    return ongrid2


