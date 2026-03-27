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

with open('run_fit_fxt_spec', 'w') as f:
    f.write('## The script to fit FXT spec files for all observations\n\n')

for i in range(len(obsid_list)):
    obsid = obsid_list[i]
    #xrt_observation = xrt_observations[i]
    #gal_name = gal_list[i]
    x_source_path_list = get_folders_only('data/' + obsid + '/x*')
    ## search for only directory

    ## The download file has been combined 
    if len(x_source_path_list)==0:
        print('No FXT source folder found for observation ' + obsid)
        continue
    for j in range(len(x_source_path_list)):
        x_source_path = x_source_path_list[j]

        with open('run_fit_fxt_spec', 'a') as f:
            ## copy fit_spec.py to each xrt observation folder
            f.write('cp fit_spec.py ' + x_source_path + '\n')
            f.write('cd ' + x_source_path + '\n')
            f.write('python3 fit_spec.py\n')
            f.write('cd ../../.. \n\n')

        print('Spec fitting files generated for observation ' + x_source_path)
    #print(len(filter_file))
