"""save_stereo_vars.py
calls functions in import_LE_data to project arctic LE data
to stereo grid and save the resulting ouput. Used to 
make input data for make_seasonal_storm_composite_plots.py
"""
from import_LE_data import *


num = ['034','035']
#num = ['005']
#num = ['001','002','003','004','005','006','007','008','009','010']
#num = ['011','012','013','014','015','016','017','018','019','020']
#num = ['021','022','023','024','025','026','027','028','029','030']
for n in num:
    save_psl_laplacian_as_stereo(n)
    #save_atm_vars_as_stereo(n)
    #save_ice_vars_as_stereo(n)
    #save_ocn_vars_as_stereo(n)



