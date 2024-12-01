## check ULXs in galaxies, using Walton 2022

import numpy as np
from astropy.table import Table
from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy import units as u
import csv, re
import pandas as pd

## Walton 2022
hdulist = fits.open('/home/yijia/work/radio-X-ray/files/J_MNRAS_509_1587_master.dat.gz.fits')
evt_ULX = hdulist[1].data  # EVENT extension
hdr_ULX = hdulist[1].header # HEADER
name_ULX = evt_ULX['Name']
ra_ULX = evt_ULX['RABdeg']
dec_ULX = evt_ULX['DEBdeg']
dist_ULX = evt_ULX['Dist']
Lpeak_2SXPS = evt_ULX['Lpeak2SXPS']
Lpeak_4XMM = evt_ULX['Lpeak4XMM']
Lpeak_CSC2 = evt_ULX['LpeakCSC2']
ULX_index = np.array(range(len(name_ULX)))

gal_list = ['NGC1313', 'NGC5128', 'NGC0247', 'NGC7793', 'NGC0055', 'NGC0224', 'NGC1316']

with open('ULX_list.py', 'w') as f:
    f.write('## The coordinate of ULXs \n')
    f.write('from astropy.coordinates import SkyCoord\n')
    f.write('from astropy import units as u\n')
    f.write('source_dic = {\\\n')
    
    for gal_name in gal_list:
        source_index = 1
        for i in range(len(evt_ULX)):
            if gal_name in name_ULX[i]:
                f.write('\'' + gal_name + '_X' + str(source_index) + '\': ' + 'SkyCoord(ra=' + str(ra_ULX[i]) + '*u.degree, dec=' + str(dec_ULX[i]) + '*u.degree, frame=\'icrs\'),\\\n')
                source_index += 1
                print(name_ULX[i], Lpeak_2SXPS[i], Lpeak_4XMM[i])
    
    f.write('\'NGC55_X1\': SkyCoord(ra=3.870375*u.degree, dec=-39.2218889*u.degree, frame=\'icrs\'),\\\n')
    f.write('\'NGC1313_X0\': SkyCoord(ra=49.5916667*u.degree, dec=-66.6011944*u.degree, frame=\'icrs\'),\\\n')
    f.write('}\n')


