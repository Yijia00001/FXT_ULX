## plot the images of the searching galaxies

import os, glob, csv
import numpy as np
import pandas as pd
from astropy.io import fits
#from pdf2image import convert_from_path
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.time import Time
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import matplotlib.colors as colors
import matplotlib.gridspec as gridspec
from matplotlib.patches import Ellipse

def get_coord_wcs(reg_file):
    with open(reg_file, 'r') as f:
        for line in f:
            reg = line
    temp = reg.split('(')[1].split(')')[0].split(',')
    ra, dec = temp[0], temp[1]
    return SkyCoord(ra+' '+dec, frame='icrs', unit=(u.hourangle, u.degree)), temp[2][:-1], temp[3][:-1], temp[4]

source_info = pd.read_csv('source_list.txt')
source_info_array = np.array(source_info)
source_name_list = np.array(source_info['source_name'])
obsid_list = np.array(source_info.keys())[1:]

sxps_source = pd.read_csv('./files/sxps_source.csv')
sxps_source_gal = np.array(sxps_source['galaxy_name'])
sxps_source_ra = np.array(sxps_source['ra'])
sxps_source_dec = np.array(sxps_source['dec'])

gal_list = ['NGC1313', 'NGC5128', 'NGC247', 'NGC7793', 'NGC55', 'NGC0224', 'NGC1316']
gal_obsid_dic = {}

plot_num = 0

## count the observation number of each galaxy
for gal_name in gal_list:
    gal_obsid_temp = []
    for i in range(len(obsid_list)):
        obsid = obsid_list[i]
        fits_file = fits.open('./data/' + obsid + '/a/a.img')
        object = fits_file[0].header['OBJECT']
        if object==gal_name:
            #print(obsid, gal_name)
            #print(plot_num)
            plot_num += 1
            gal_obsid_temp.append(obsid)
    gal_obsid_dic[gal_name] = gal_obsid_temp

print(gal_obsid_dic)

pixel_size = 9.66866
v_range = [0.1, 10]

for gal_name in gal_list:
    fig = plt.figure(figsize=(6*2*(len(gal_obsid_dic[gal_name])//4+1.2), 2+4*6))
    gs = gridspec.GridSpec(4, 2*(len(gal_obsid_dic[gal_name])//4+1), wspace=0.1, hspace=0.2)

    gal_obsid_list = gal_obsid_dic[gal_name]
    for i in range(len(gal_obsid_list)):
        img_file = fits.open('data/' + gal_obsid_list[i] + '/a/a.img')
        hdr, evt = img_file[0].header, img_file[0].data
        exposure = hdr['EXPOSURE']
        w_a = WCS(hdr)
        #x0, y0 = w_a.world_to_pixel(source0)
        ax_a = fig.add_subplot(gs[i%4, i//4*2], projection=w_a)#, projection=w)
        im = ax_a.imshow(evt, cmap='gist_yarg', norm=colors.LogNorm(vmin=v_range[0], vmax=v_range[1]))
        ax_a.set_title(hdr['DATE-OBS'].split('T')[0] + ' A, exp=' + format(exposure, '.0f') + 's', fontsize=20)
        ax_a.coords[0].set_axislabel('R.A.')
        ax_a.coords[1].set_axislabel('Dec.')
        reg_file_list = glob.glob('./data/' + gal_obsid_list[i] + '/a/x*.reg')
        for reg_file in reg_file_list:
            source0, ra_axis0, dec_axis0, angle0  = get_coord_wcs(reg_file)
            x0, y0 = w_a.world_to_pixel(source0)
            e = Ellipse(xy=(x0, y0), width=2*float(ra_axis0)/pixel_size, height=2*float(dec_axis0)/pixel_size, angle=float(angle0),\
                facecolor='None', edgecolor='green', linewidth=1.0)
            ax_a.add_artist(e)

        img_file = fits.open('data/' + gal_obsid_list[i] + '/b/b.img')
        hdr, evt = img_file[0].header, img_file[0].data
        exposure = hdr['EXPOSURE']
        w_b = WCS(hdr)
        #x0, y0 = w_a.world_to_pixel(source0)
        ax_b = fig.add_subplot(gs[i%4, i//4*2+1], projection=w_b)#, projection=w)
        im = ax_b.imshow(evt, cmap='gist_yarg', norm=colors.LogNorm(vmin=v_range[0], vmax=v_range[1]))
        ax_b.set_title(hdr['DATE-OBS'].split('T')[0] + ' B, exp=' + format(exposure, '.0f') + 's', fontsize=20)
        ax_b.coords[0].set_axislabel('R.A.')
        ax_b.coords[1].set_axislabel('Dec.')
        reg_file_list = glob.glob('./data/' + gal_obsid_list[i] + '/b/x*.reg')
        for reg_file in reg_file_list:
            source0, ra_axis0, dec_axis0, angle0  = get_coord_wcs(reg_file)
            x0, y0 = w_b.world_to_pixel(source0)
            e = Ellipse(xy=(x0, y0), width=2*float(ra_axis0)/pixel_size, height=2*float(dec_axis0)/pixel_size, angle=float(angle0),\
                facecolor='None', edgecolor='green', linewidth=1.0)
            ax_b.add_artist(e)
        
        q_sxps = (sxps_source_gal==gal_name)
        for j in range(len(sxps_source_gal[q_sxps])):
            ra_sxps, dec_sxps = sxps_source_ra[q_sxps][j], sxps_source_dec[q_sxps][j]
            source_sxps = SkyCoord(ra=ra_sxps, dec=dec_sxps, unit=(u.degree, u.degree), frame='icrs')
            x_sxps_a, y_sxps_a = w_a.world_to_pixel(source_sxps)
            ax_a.scatter(x_sxps_a, y_sxps_a, color='red', s=3, marker='.')
            x_sxps_b, y_sxps_b = w_b.world_to_pixel(source_sxps)
            ax_b.scatter(x_sxps_b, y_sxps_b, color='red', s=3, marker='.')
            
    plt.rcParams['savefig.bbox'] = 'tight'
    plt.savefig('./figure/whole_img/' + gal_name + '.pdf')
    #break

#plt.show()
