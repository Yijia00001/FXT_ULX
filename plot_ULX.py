## plot the basic info of ULXs

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
source_index = 0

for source_name in source_name_list:
    ## aimed target
    #source_name = 'IC342_X1'
    #source0 = SkyCoord(ra=56.481717*u.degree, dec=68.082025*u.degree, frame='icrs') ## X1
    #source0 = SkyCoord(ra=56.565*u.degree, dec=68.186667*u.degree, frame='icrs') ## X2
    source_info_temp = source_info_array[source_index][1:]
    if len(source_info_temp[source_info_temp!='-'])==0:
        source_index += 1
        continue
    else:
        obsid_list_target = obsid_list[source_info_temp!='-']
        print(obsid_list_target)
    fig = plt.figure(figsize=(1+8*len(obsid_list_target), 4 + 4*6))
    gs = gridspec.GridSpec(6, 2*len(obsid_list_target), wspace=0.1, hspace=0.2)

    width = 800
    pixel_size = 9.66866
    ax_flux = fig.add_subplot(gs[3, :])
    ax_index = fig.add_subplot(gs[4, :])
    ax_NH = fig.add_subplot(gs[5, :])

    for i in range(len(obsid_list_target)):
        print(np.array(source_info[obsid_list_target[i]])[source_name_list==source_name][0])
        source_obs_id = np.array(source_info[obsid_list_target[i]])[source_name_list==source_name][0]

        source0, ra_axis0, dec_axis0, angle0  = get_coord_wcs('data/' + obsid_list_target[i] + '/a/' + source_obs_id + '.reg')
        print(source0, ra_axis0, dec_axis0, angle0)
        ## for A 
        img_file = fits.open('data/' + obsid_list_target[i] + '/a/a.img')
        hdr, evt = img_file[0].header, img_file[0].data
        w_a = WCS(hdr)
        x0, y0 = w_a.world_to_pixel(source0)
        ax_a = fig.add_subplot(gs[0, 2*i])#, projection=w)
        im = ax_a.imshow(evt, cmap='gist_yarg', norm=colors.LogNorm(vmin=0.1, vmax=30))
        e = Ellipse(xy=(x0, y0), width=2*float(ra_axis0)/pixel_size, height=2*float(dec_axis0)/pixel_size, angle=float(angle0),\
                facecolor='None', edgecolor='green', linewidth=2)
        ax_a.add_artist(e)
        x1, x2 = x0 - width/2/pixel_size, x0 + width/2/pixel_size
        y1, y2 = y0 - width/2/pixel_size, y0 + width/2/pixel_size
        ax_a.set_xlim(x1, x2)
        ax_a.set_ylim(y1, y2)
        ax_a.set_title(hdr['DATE-OBS'].split('T')[0] + ' A', fontsize=20)
        ax_a.set_xticks([])
        ax_a.set_yticks([])

        ## for B 
        img_file = fits.open('data/' + obsid_list_target[i] + '/b/b.img')
        source0, ra_axis0, dec_axis0, angle0 = get_coord_wcs('data/' + obsid_list_target[i] + '/b/' + source_obs_id + '.reg')
        hdr, evt = img_file[0].header, img_file[0].data
        w_b = WCS(hdr)
        ax_b = fig.add_subplot(gs[0, 2*i+1])#, projection=w)
        x0, y0 = w_b.world_to_pixel(source0)
        im = ax_b.imshow(evt, cmap='gist_yarg', norm=colors.LogNorm(vmin=0.1, vmax=30))
        e = Ellipse(xy=(x0, y0), width=2*float(ra_axis0[0:5])/pixel_size, height=2*float(dec_axis0[0:5])/pixel_size, angle=float(angle0),\
                facecolor='None', edgecolor='green', linewidth=2)
        ax_b.add_artist(e)
        x1, x2 = x0 - width/2/pixel_size, x0 + width/2/pixel_size
        y1, y2 = y0 - width/2/pixel_size, y0 + width/2/pixel_size
        ax_b.set_xlim(x1, x2)
        ax_b.set_ylim(y1, y2)
        ax_b.set_title(hdr['DATE-OBS'].split('T')[0] + ' B', fontsize=20)
        ax_b.set_xticks([])
        ax_b.set_yticks([])

        ## plot the fitting results

        fit_result_file = 'data/' + obsid_list_target[i] + '/' + source_obs_id + '_fit_results.csv'
        fit_result = pd.read_csv(fit_result_file)
        para = np.array(fit_result['para'])
        value = np.array(fit_result['value'])
        value_down = np.array(fit_result['value_down'])
        value_up = np.array(fit_result['value_up'])

        fit_img_path = 'data/' + obsid_list_target[i] + '/fit_image/' + source_obs_id + '.png'
        if not os.path.exists(fit_img_path):
            continue

        t0 = Time(hdr['DATE-OBS'], format='isot', scale='utc')
        fitting_result = mpimg.imread(fit_img_path)
        ax_fitting = fig.add_subplot(gs[1:3, 2*i:(2*i+2)])
        ax_fitting.imshow(fitting_result, aspect='auto')
        ax_fitting.set_xticks([])
        ax_fitting.set_yticks([])
        ax_fitting.set_title(r'$\rm \chi^2/dof$ ' + format(value[para=='dof'][0], '.1f') + '/' + str(int(value_up[para=='dof'][0])), fontsize=20)

        ax_flux.plot([t0.mjd, t0.mjd], [value_down[para=='Flux'][0], value_up[para=='Flux'][0]], color='blue')
        ax_flux.scatter(t0.mjd, value[para=='Flux'][0], s=20, color='k')
        ax_flux.tick_params(which='both', labelsize=20)
        ax_flux.tick_params(which='both', direction='in')
        ax_flux.tick_params(which='major', length=6)
        ax_flux.tick_params(which='both', top='on', right='on')
        ax_flux.set_xlabel('MJD', fontsize=20)
        ax_flux.set_ylabel(r'log($\rm Flux / erg~s^{-1}~cm^{-2}$)', fontsize=20)

        ax_index.plot([t0.mjd, t0.mjd], [value_down[para=='PhotonIndex'][0], value_up[para=='PhotonIndex'][0]], color='blue')
        ax_index.scatter(t0.mjd, value[para=='PhotonIndex'][0], s=20, color='k')
        ax_index.tick_params(which='both', labelsize=20)
        ax_index.tick_params(which='both', direction='in')
        ax_index.tick_params(which='major', length=6)
        ax_index.tick_params(which='both', top='on', right='on')
        ax_index.set_xlabel('MJD', fontsize=20)
        ax_index.set_ylabel(r'$\Gamma$', fontsize=20)

        ax_NH.plot([t0.mjd, t0.mjd], [value_down[para=='NH'][0], value_up[para=='NH'][0]], color='blue')
        ax_NH.scatter(t0.mjd, value[para=='NH'][0], s=20, color='k')
        ax_NH.tick_params(which='both', labelsize=20)
        ax_NH.tick_params(which='both', direction='in')
        ax_NH.tick_params(which='major', length=6)
        ax_NH.tick_params(which='both', top='on', right='on')
        ax_NH.set_xlabel('MJD', fontsize=20)
        ax_NH.set_ylabel(r'${\rm N_H\times 10^{22}~cm^{-2}}$', fontsize=20)
        #break

    source_index += 1
    plt.rcParams['savefig.bbox'] = 'tight'
    plt.savefig('figure/' + source_name + '.pdf')
    #plt.show()