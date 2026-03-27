## copy bkg files to all observations

import numpy as np
import glob
import os
from astropy.io import fits 
import pandas as pd
from astropy.time import Time

obsid_list = []
with open('./obsid_list.txt', 'r') as f:
    for line in f:
        obsid_list.append(line.strip())

#fxt_observations = glob.glob('./swift/*/xrt')
gal_list = ['NGC1313', 'NGC5128', 'NGC247', 'NGC7793', 'NGC55', 'NGC0224', 'NGC1316', 'NGC4038', 'NGC0925', 'IC0342', 'M82']

for gal_name in gal_list:
    gal_obsid_temp = []
    gal_obstime_temp = []
    gal_obsexp_temp = []
    bkg_exist_list = []
    bkg_not_exist_list = []
    for i in range(len(obsid_list)):
        obsid = obsid_list[i]
        fits_file = fits.open('./data/' + obsid + '/a/a.img')
        object = fits_file[0].header['OBJECT']
        obstime = fits_file[0].header['DATE-OBS']
        exposure = fits_file[0].header['EXPOSURE']
        
        if object=='2XMM J034615.6+681112':
            object = 'IC0342'
        
        if object==gal_name:
            print(gal_name, obsid)
            if os.path.exists('./data/' + obsid + '/bkg.reg'):
                bkg_exist_list.append(obsid)
            else:
                bkg_not_exist_list.append(obsid)

    for j in range(len(bkg_not_exist_list)):
        obsid = bkg_not_exist_list[j]
        os.system('cp ' + './data/' + bkg_exist_list[0] + '/bkg.reg' + ' ./data/' + obsid)
        print(gal_name + ', copy bkg file from ' + bkg_exist_list[0] + ' to ' + obsid)
            


