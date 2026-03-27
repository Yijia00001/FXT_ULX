## generate spec files for FXT observations

import os, glob
import matplotlib.pyplot as plt
import numpy as np
from astropy import units as u
from astropy.io import fits

## 
def get_folders_only(path_pattern):
    all_items = glob.glob(path_pattern)
    folders = [item for item in all_items if os.path.isdir(item)]
    return folders

obsid_list = []
with open('./obsid_list.txt', 'r') as f:
    for line in f:
        obsid_list.append(line.strip())

#fxt_observations = glob.glob('./swift/*/xrt')
gal_list = ['NGC1313', 'NGC5128', 'NGC247', 'NGC7793', 'NGC55', 'NGC0224', 'NGC1316', 'NGC4038', 'NGC0925', 'IC0342', 'M82']

with open('run_gen_bkg_rate_flux', 'w') as f:
    f.write('## The script to generate bkg rate flux files for all observations\n\n')

for i in range(len(obsid_list)):
    obsid = obsid_list[i]
    #xrt_observation = xrt_observations[i]
    #gal_name = gal_list[i]
    x_source_path = get_folders_only('data/' + obsid + '/bkg_spec')[0]
    ## search for only directory
    with open('run_gen_bkg_rate_flux', 'a') as f:
        ## copy fit_spec.py to each xrt observation folder
        f.write('cp bkg_rate_flux.py ' + x_source_path + '\n')
        f.write('cd ' + x_source_path + '\n')
        f.write('python3 bkg_rate_flux.py\n')
        f.write('cd /home/yijia/work/EP_ULX\n\n')
    print('Spec fitting files generated for observation ' + x_source_path)
        
    #print(len(filter_file))
