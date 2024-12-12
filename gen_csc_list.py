## generate source list from CSC2.0 for transients searching

import numpy as np
from astropy.table import Table
from astropy.coordinates import SkyCoord
from astropy import units as u
import matplotlib.pyplot as plt
import csv

csc2_table = Table.read('/home/yijia/work/radio-X-ray/files/csc2_mrt.txt', format='ascii.cds')
csc2_index = np.array(range(len(csc2_table['2CXO'])))
CXO_IDs = np.array(csc2_table['2CXO'])
ra_csc = np.array(csc2_table['RAdeg'])
dec_csc = np.array(csc2_table['DEdeg'])
flux_b = np.array(csc2_table['Fluxb'])
coord_csc2 = SkyCoord(ra=ra_csc*u.degree, dec=dec_csc*u.degree, frame='icrs')

gal_coord_dic = {'NGC1313':SkyCoord(ra=49.56686*u.degree, dec=-66.49826*u.degree, frame='icrs'),\
                'NGC5128':SkyCoord(ra=201.36506*u.degree, dec=-43.01911*u.degree, frame='icrs'),\
                'NGC247':SkyCoord(ra=11.78564*u.degree, dec=-20.7604*u.degree, frame='icrs'),\
                'NGC7793':SkyCoord(ra=358.81392*u.degree, dec=-32.59103*u.degree, frame='icrs'),\
                'NGC55':SkyCoord(ra=3.723333*u.degree, dec=-39.47996*u.degree, frame='icrs'),\
                'NGC0224':SkyCoord(ra=10.68471*u.degree, dec=41.26875*u.degree, frame='icrs'),\
                'NGC1316':SkyCoord(ra=50.67412*u.degree, dec=-37.2082*u.degree, frame='icrs')}

gal_name_list = list(gal_coord_dic.keys())

r = 2600 ## arcsec
content = []
content.append(['galaxy_name', 'ra', 'dec'])

for i in range(len(gal_name_list)):
    gal_name = gal_name_list[i]
    source_gal = gal_coord_dic[gal_name]
    sep = coord_csc2.separation(source_gal).arcsec
    q = (sep < r) & (flux_b > 1e-14)
    print(gal_name, len(coord_csc2[q]))
    for j in range(len(coord_csc2[q])):
        ra, dec = ra_csc[q][j], dec_csc[q][j]
        content.append([gal_name, ra, dec])

with open('./files/csc_source.csv', 'w') as f:
    writer = csv.writer(f)
    writer.writerows(content)
