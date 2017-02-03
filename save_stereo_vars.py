"""save_stereo_vars.py
calls functions in import_LE_data to project arctic LE data
to stereo grid and save the resulting ouput. Used to 
make input data for make_seasonal_storm_composite_plots.py
"""
from import_LE_data import *

num = '011'
save_atm_vars_as_stereo(num)
save_ice_vars_as_stereo(num)
save_ocn_vars_as_stereo(num)



