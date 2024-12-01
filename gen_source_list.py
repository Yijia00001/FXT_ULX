## generate source list for a given source

import os, glob, csv
import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.time import Time

def get_coord_wcs(reg_file):
    with open(reg_file, 'r') as f:
        for line in f:
            reg = line
    temp = reg.split('(')[1].split(')')[0].split(',')
    ra, dec = temp[0], temp[1]
    return SkyCoord(ra+' '+dec, frame='icrs', unit=(u.hourangle, u.degree)), temp[2], temp[3], temp[4]

source_dic = {\
'NGC1313_X1': SkyCoord(ra=49.4477949806808*u.degree, dec=-66.5026721661653*u.degree, frame='icrs'),\
'NGC1313_X2': SkyCoord(ra=49.5766173515857*u.degree, dec=-66.5007637891666*u.degree, frame='icrs'),\
'NGC1313_X3': SkyCoord(ra=49.5830052739027*u.degree, dec=-66.4863726381794*u.degree, frame='icrs'),\
'NGC1313_X4': SkyCoord(ra=49.5871657179149*u.degree, dec=-66.5096654560853*u.degree, frame='icrs'),\
'NGC1313_X5': SkyCoord(ra=50.0206123443106*u.degree, dec=-66.7037929296161*u.degree, frame='icrs'),\
'NGC5128_X1': SkyCoord(ra=201.25684*u.degree, dec=-43.19227*u.degree, frame='icrs'),\
'NGC5128_X2': SkyCoord(ra=201.32695*u.degree, dec=-43.05168*u.degree, frame='icrs'),\
'NGC5128_X3': SkyCoord(ra=201.34156*u.degree, dec=-43.02377*u.degree, frame='icrs'),\
'NGC5128_X4': SkyCoord(ra=201.35207*u.degree, dec=-43.0267*u.degree, frame='icrs'),\
'NGC5128_X5': SkyCoord(ra=201.36183*u.degree, dec=-43.00543*u.degree, frame='icrs'),\
'NGC5128_X6': SkyCoord(ra=201.362665634207*u.degree, dec=-43.0196541961934*u.degree, frame='icrs'),\
'NGC5128_X7': SkyCoord(ra=201.37518*u.degree, dec=-43.01027*u.degree, frame='icrs'),\
'NGC5128_X8': SkyCoord(ra=201.41*u.degree, dec=-43.0348*u.degree, frame='icrs'),\
'NGC5128_X9': SkyCoord(ra=201.426147049614*u.degree, dec=-42.9954472965537*u.degree, frame='icrs'),\
'NGC0247_X1': SkyCoord(ra=11.7661604923803*u.degree, dec=-20.7955512941367*u.degree, frame='icrs'),\
'NGC7793_X1': SkyCoord(ra=359.44958837373*u.degree, dec=-32.582575627892*u.degree, frame='icrs'),\
'NGC7793_X2': SkyCoord(ra=359.462573618803*u.degree, dec=-32.6241577824994*u.degree, frame='icrs'),\
'NGC7793_X3': SkyCoord(ra=359.536758427298*u.degree, dec=-32.5675827208591*u.degree, frame='icrs'),\
'NGC0224_X1': SkyCoord(ra=10.58937*u.degree, dec=41.26675*u.degree, frame='icrs'),\
'NGC0224_X2': SkyCoord(ra=10.69618*u.degree, dec=41.27075*u.degree, frame='icrs'),\
'NGC0224_X3': SkyCoord(ra=10.71669*u.degree, dec=41.51873*u.degree, frame='icrs'),\
'NGC0224_X4': SkyCoord(ra=10.7213490711122*u.degree, dec=41.2398919286776*u.degree, frame='icrs'),\
'NGC0224_X5': SkyCoord(ra=11.43997*u.degree, dec=41.66142*u.degree, frame='icrs'),\
'NGC1316_X1': SkyCoord(ra=50.657205410073*u.degree, dec=-37.223112408439*u.degree, frame='icrs'),\
'NGC1316_X2': SkyCoord(ra=50.66906*u.degree, dec=-37.19327*u.degree, frame='icrs'),\
'NGC1316_X3': SkyCoord(ra=50.67174*u.degree, dec=-37.21405*u.degree, frame='icrs'),\
'NGC1316_X4': SkyCoord(ra=50.67992*u.degree, dec=-37.18529*u.degree, frame='icrs'),\
'NGC1316_X5': SkyCoord(ra=50.7134966484313*u.degree, dec=-37.1635953695803*u.degree, frame='icrs'),\
'NGC1316_X6': SkyCoord(ra=50.773866178381*u.degree, dec=-37.1864054985313*u.degree, frame='icrs'),\
'NGC55_X1': SkyCoord(ra=3.870375*u.degree, dec=-39.2218889*u.degree, frame='icrs'),\
'NGC1313_X0': SkyCoord(ra=49.5916667*u.degree, dec=-66.6011944*u.degree, frame='icrs'),\
}


# source_dic = {'NGC1313_X1': SkyCoord(ra=49.583333*u.degree, dec=-66.486361*u.degree, frame='icrs'),\
#               'NGC1313_X2': SkyCoord(ra=49.5916667*u.degree, dec=-66.6011944*u.degree, frame='icrs'),\
#               'NGC5128_X1': SkyCoord(ra=201.3325*u.degree, dec=-43.053333*u.degree, frame='icrs'),\
#               'NGC5128_X2': SkyCoord(ra=201.3258333*u.degree, dec=-43.051111*u.degree, frame='icrs'),
#               'NGC247_X1': SkyCoord(ra=11.77*u.degree, dec=-20.7952778*u.degree, frame='icrs'),\
#               'NGC7793_X1': SkyCoord(ra=359.4620833*u.degree, dec=-32.6238889*u.degree, frame='icrs'),\
#               'NGC55_X1': SkyCoord(ra=3.870375*u.degree, dec=-39.2218889*u.degree, frame='icrs'),\
#               'NGC224_X1': SkyCoord(ra=10.72125*u.degree, dec=41.239444*u.degree, frame='icrs'),\
#               'NGC224_X2': SkyCoord(ra=10.6816667*u.degree, dec=41.421944*u.degree, frame='icrs'),\
#               'NGC1316_X1': SkyCoord(ra=50.67125*u.degree, dec=-37.209722*u.degree, frame='icrs'),\
#               'NGC1316_X2': SkyCoord(ra=50.67375*u.degree, dec=-37.2080556*u.degree, frame='icrs'),\
#               'NGC1316_X3': SkyCoord(ra=50.67416667*u.degree, dec=-37.206944*u.degree, frame='icrs'),\
#               'NGC1316_X4': SkyCoord(ra=50.713333*u.degree, dec=-37.1636111*u.degree, frame='icrs'),\
#               }


obsid_list = []
with open('obsid_list.txt', 'r') as f:
    for line in f:
        obsid_list.append(line.strip())

source_names = list(source_dic.keys())

for i in range(len(source_dic.keys())):
    source0 = source_dic[source_names[i]]
    content = source_names[i]
    content_head = 'source_name'
    for j in range(len(obsid_list)):
        source_obs_list = glob.glob('data/' + obsid_list[j] + '/a/x*.reg')
        content_head = content_head + ',' + obsid_list[j]
        #print(content_head)
        if_in_obs = False
        for k in range(len(source_obs_list)):
            source1, ra_axis, dec_axis, angle = get_coord_wcs(source_obs_list[k])
            sep = source1.separation(source0).arcsec
            if sep < 20:
                content = content + ',' + source_obs_list[k].split('/')[-1].split('.')[0]
                if_in_obs = True
                break
        if not if_in_obs:
            content = content + ',-' 
    if i==0:
        with open('source_list.txt', 'w') as f:
            f.write(content_head + '\n')
    
    with open('source_list.txt', 'a') as f:
        f.write(content + '\n')

