##check the flux of ULXs from FXT observation
## generate region file for different sources

import os, glob, csv
import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from astropy import units as u

class fxt_source_info:
    def __init__(self, source_names, ra, ra_err, dec, dec_err):#, flux, flux_err, SNR):
        self.source_names = source_names
        self.ra = ra
        self.dec = dec
        self.ra_err = ra_err
        self.dec = dec
        self.dec_err = dec_err
        # self.flux = flux
        # self.flux_err = flux_err
        # self.SNR = SNR

def read_fxt_source(fxt_table):
    source_names = np.array(fxt_table['Name'])
    ra = np.array(fxt_table['RA'])
    ra_err = np.array(fxt_table['RA_Err (arcsec)'])
    dec = np.array(fxt_table['Dec'])
    dec_err = np.array(fxt_table['Dec_Err (arcsec)'])
    flux = np.array(fxt_table['Flux (erg/cm^2/s)'])
    flux_err = np.array(fxt_table['Flux Err'])
    SNR = np.array(fxt_table['SNR'])
    return fxt_source_info(source_names, ra, ra_err, dec, dec_err, flux, flux_err, SNR)

## read reg file line by line
def read_reg_file(reg_file):
    output = []
    with open(reg_file, 'r') as f:
        for line in f:
            output.append(line)
    return output

def get_coord(reg):
    temp = reg.split('(')[1].split(')')[0].split(',')
    x, y, x_axis, y_axis = float(temp[0]), float(temp[1]), float(temp[2]), float(temp[3])
    angle = float(temp[4])
    return x, y, x_axis, y_axis, angle
    #return ra, dec, ra_err, dec_err
    # return fxt_source_info(source_names='x0', ra=ra, dec=dec, ra_err=ra_err, dec_err=dec_err)

def write_reg_file(path, ra, dec, ra_axis, dec_axis, angle):
    with open(path, 'w') as f:
        f.write('# Region file format: DS9 version 4.1\n')
        f.write('global color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n')
        f.write('ellipse(' + ra + ',' + dec + ',' + ra_axis + ',' + dec_axis + ',' + angle + ')\n')

obs_id_path = glob.glob('./data/*')
pixel_size = 9.66866 ## arcsec per pixel

for i in range(len(obs_id_path)):
    src_a = fits.open(obs_id_path[i] + '/a/wav_src.fits')[1].data
    src_b = fits.open(obs_id_path[i] + '/b/wav_src.fits')[1].data
    ra_list_a, dec_list_a, x_list_a, y_list_a = src_a['RA'], src_a['DEC'], src_a['X'], src_a['Y']
    ra_list_b, dec_list_b, x_list_b, y_list_b = src_b['RA'], src_b['DEC'], src_b['X'], src_b['Y']

    a_reg = read_reg_file(obs_id_path[i] + '/a/wavdetect.reg')
    b_reg = read_reg_file(obs_id_path[i] + '/b/wavdetect.reg')
    source_num = 1
    for k in range(len(a_reg)):
        x0_a, y0_a, x0_axis_a, y0_axis_a, angle_a = get_coord(a_reg[k])
        r_a = ((x_list_a-x0_a)**2 + (y_list_a-y0_a)**2)**0.5
        ra_a, dec_a = ra_list_a[r_a<1][0], dec_list_a[r_a<1][0]
        #print(source_a.ra)
        source_a = SkyCoord(ra=ra_a*u.degree, dec=dec_a*u.degree, frame='icrs')
        #source_a = w_a.pixel_to_world(x0_a, y0_a)
        for j in range(len(b_reg)):
            x0_b, y0_b, x0_axis_b, y0_axis_b, angle_b = get_coord(b_reg[j])
            #source_b = w_b.pixel_to_world(x0_b, y0_b)
            r_b = ((x_list_b-x0_b)**2 + (y_list_b-y0_b)**2)**0.5
            ra_b, dec_b = ra_list_b[r_b<1][0], dec_list_b[r_b<1][0]
            source_b = SkyCoord(ra=ra_b*u.degree, dec=dec_b*u.degree, frame='icrs')
            sep = source_a.separation(source_b)
            if sep.arcsec<10:
                coord_temp_a = source_a.to_string('hmsdms').split()
                coord_temp_b = source_b.to_string('hmsdms').split()
                ra_hms_a, dec_dms_a = coord_temp_a[0].replace('h', ':').replace('m', ':')[:-1], coord_temp_a[1].replace('d', ':').replace('m', ':')[:-1]
                ra_hms_b, dec_dms_b = coord_temp_b[0].replace('h', ':').replace('m', ':')[:-1], coord_temp_b[1].replace('d', ':').replace('m', ':')[:-1]
                write_reg_file(path=obs_id_path[i] + '/a/x' + str(source_num) + '.reg', ra=ra_hms_a, dec=dec_dms_a,\
                                ra_axis=str(x0_axis_a*pixel_size)+'"', dec_axis=str(y0_axis_a*pixel_size)+'"', angle=str(angle_a))
                write_reg_file(path=obs_id_path[i] + '/b/x' + str(source_num) + '.reg', ra=ra_hms_b, dec=dec_dms_b,\
                                ra_axis=str(x0_axis_b*pixel_size)+'"', dec_axis=str(y0_axis_b*pixel_size)+'"', angle=str(angle_b))
                print(x0_a, y0_a, source_a.ra, source_a.dec, sep.arcsec, obs_id_path[i], source_num)
                source_num = source_num + 1
        #break

    #break
    #print(source_list_A.ra)


