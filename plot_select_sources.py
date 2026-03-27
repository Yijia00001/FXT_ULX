## plot those selected sources

import os, glob, csv
import numpy as np
import pandas as pd
from astropy.io import fits
#from pdf2image import convert_from_path
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from astropy import units as u
from photutils.aperture import CircularAperture, aperture_photometry
from astropy.time import Time
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import matplotlib.colors as colors
import matplotlib.gridspec as gridspec
from matplotlib.patches import Ellipse
from scipy.stats import poisson
from scipy.optimize import root_scalar

def find_lambda_for_poisson_cdf(n, prob=0.9): 
    low = 1e-6
    high = max(1, n)
    def f(lam):
        return poisson.cdf(n, lam) - prob
    while f(high) > 0:
        high *= 2
    sol = root_scalar(f, bracket=[low, high], method='bisect')
    return sol.root

def find_n_for_poisson_cdf(lam, prob=0.9):
    high = 0
    def f(n):
        return poisson.cdf(n, lam) - prob
    while f(high) < 0:
        high += 1
    return high

def get_coord_wcs(reg_file):
    with open(reg_file, 'r') as f:
        for line in f:
            reg = line
    temp = reg.split('(')[1].split(')')[0].split(',')
    ra, dec = temp[0], temp[1]
    if len(temp)==5:
        return SkyCoord(ra+' '+dec, frame='icrs', unit=(u.hourangle, u.degree)), temp[2][:-1], temp[3][:-1], temp[4]
    elif len(temp)==3:
        return SkyCoord(ra+' '+dec, frame='icrs', unit=(u.hourangle, u.degree)), temp[2][:-1]

# ## generate effective number for a err
def output_err_str(err):
    err_str = format(err, '.2g')
    first_eff_num = '0'
    for i in err_str:
        if i!='0' and i!='.':
            first_eff_num = i
            break
    if first_eff_num!='1':
        if first_eff_num!='0':
            err_str = format(err, '.1g')
        else:
            err_str = '0'
    return err_str

# The value and err are float 
def output_sci_results(value, err):
    if type(err)!=list:
        err_str = output_err_str(err)
        ## with . in the error 
        if '.' in err_str:
            decimal_places = len(err_str.split('.')[1])
        ## int error
        else:
            decimal_places = 1
        value_str = format(value, '.' + str(decimal_places) + 'f')
        return '$' + value_str + '\pm' + err_str + '$'
    elif len(err)==2:
        err_down_str, err_up_str = output_err_str(err[0]), output_err_str(err[1])
        if ('.' in err_down_str) and ('.' in err_up_str):
            decimal_places_down, decimal_places_up = len(err_down_str.split('.')[1]), len(err_up_str.split('.')[1])
            decimal_places = np.max([decimal_places_down, decimal_places_up])
        else:
            decimal_places = 1
        value_str = format(value, '.' + str(decimal_places) + 'f')
        err_down_str, err_up_str = format(err[0], '.' + str(decimal_places) + 'f'), format(err[1], '.' + str(decimal_places) + 'f')
        return '$' + value_str + '_{-' + err_down_str + '}^{+' + err_up_str + '}$'

galaxy_NH_file = pd.read_csv('./combined_data/galaxt_NH.csv')
galaxy_NH_galaxy_list = np.array(galaxy_NH_file['galaxy'])
galaxy_NH_value_list = np.array(galaxy_NH_file['NH'])

obsid_list = []
with open('obsid_list.txt', 'r') as f:
    for line in f:
        obsid_list.append(line.strip())
## selected sources
source_dic = {\
'NGC55_X1': SkyCoord(ra=3.8706903*u.degree, dec=-39.22260014*u.degree, frame='icrs'),\
'NGC55_X2': SkyCoord(ra=3.69416667*u.degree, dec=-39.1905556*u.degree, frame='icrs'),\
#'NGC55_ULX3': SkyCoord(ra=3.7778257*u.degree, dec=-39.2430267*u.degree, frame='icrs'),\
'NGC247_X1': SkyCoord(ra=11.7661604923803*u.degree, dec=-20.7955512941367*u.degree, frame='icrs'),\
'NGC0925_ULX1': SkyCoord(ra=36.86483*u.degree, dec=33.57850*u.degree, frame='icrs'),\
'NGC0925_ULX2': SkyCoord(ra=36.8391667*u.degree, dec=33.5827778*u.degree, frame='icrs'),\
'NGC0925_ULX3': SkyCoord(ra=36.8345833*u.degree, dec=33.5688889*u.degree, frame='icrs'),\
#'NGC0925_ULX4': SkyCoord(ra=36.7625000*u.degree, dec=33.5972222*u.degree, frame='icrs'),\
'NGC1313_X1':SkyCoord(ra=49.5833333*u.degree, dec=-66.4863610*u.degree, frame='icrs'),\
'NGC1313_X2':SkyCoord(ra=49.5916667*u.degree, dec=-66.6011943*u.degree, frame='icrs'),\
'NGC1316_X3':SkyCoord(ra=50.5037500*u.degree, dec=-37.3213888*u.degree, frame='icrs'),\
'NGC1316_X6':SkyCoord(ra=50.6070833*u.degree, dec=-37.1533332*u.degree, frame='icrs'),\
#'NGC1316_X5':SkyCoord(ra=50.6712500*u.degree, dec=-37.2094443*u.degree, frame='icrs'),\
'NGC1316_X4':SkyCoord(ra=50.7133333*u.degree, dec=-37.2097193*u.degree, frame='icrs'),\
#'NGC1316_X7':SkyCoord(ra=50.670000*u.degree, dec=-37.2066666*u.degree, frame='icrs'),\
'IC0342_X1': SkyCoord(ra=56.4833333*u.degree, dec=+68.0819444*u.degree, frame='icrs'),\
'IC0342_X2': SkyCoord(ra=56.5656674082348*u.degree, dec=68.186911317252*u.degree, frame='icrs'),\
#'IC0342_X3': SkyCoord(ra=56.7033333*u.degree, dec=+68.0961111*u.degree, frame='icrs'),\
'M82_S9': SkyCoord(ra=148.964583*u.degree, dec=+69.6769444*u.degree, frame='icrs'),\
'NGC5128_X1': SkyCoord(ra=201.332500*u.degree, dec=-43.0533332*u.degree, frame='icrs'),\
'NGC5128_X2': SkyCoord(ra=201.325833*u.degree, dec=-43.0511110*u.degree, frame='icrs'),\
#'NGC5128_AGN': SkyCoord(ra=201.364583*u.degree, dec=-43.0188888*u.degree, frame='icrs'),\
'NGC7793_P9': SkyCoord(ra=359.536250*u.degree, dec=-32.5672221*u.degree, frame='icrs'),\
'NGC7793_P13': SkyCoord(ra=359.462083*u.degree, dec=-32.6236110*u.degree, frame='icrs'),\
}

#source_dic = {'NGC55_S30': SkyCoord(ra=3.7778257*u.degree, dec=-39.2430267*u.degree, frame='icrs')}

pixel_size = 9.66866
## The source index in Table 2
## NGC 55 ULX2: NGC 55 S24
## IC 342 X2: IC342 S16 

#gal_list = ['NGC55', 'NGC0224', 'NGC247', 'NGC0925', 'NGC1313', 'NGC1316', 'IC0342', 'NGC4038', 'NGC5128', 'NGC7793', 'M82']

target_info = pd.read_csv('./combined_data/target_list.csv')
gal_list = np.array(target_info['name'])
distance_list = np.array(target_info['distance'])

# ULX transients
# 'NGC0925_ULX4': 'NGC 925 ULX-4'
#, 'NGC55_ULX3': 'NGC 55 ULX-3'

source_name_dic = {'NGC55_X1': 'NGC 55 ULX-1', 'NGC55_X2': 'NGC 55 ULX-2',\
                   'NGC247_X1': 'NGC 247 X-1',\
                   'NGC0925_ULX1': 'NGC 925 ULX-1', 'NGC0925_ULX2': 'NGC 925 ULX-2', 'NGC0925_ULX3': 'NGC 925 ULX-3',\
                   'NGC1313_X1': 'NGC 1313 X-1', 'NGC1313_X2': 'NGC 1313 X-2',\
                   'NGC1316_X3': 'NGC 1316 X-3', 'NGC1316_X4': 'NGC 1316 X-4', 'NGC1316_X5': 'NGC 1316 X-5', 'NGC1316_X6': 'NGC 1316 X-6', 'NGC1316_X7': 'NGC 1316 X-7',\
                   'IC0342_X1': 'IC 342 X-1', 'IC0342_X3': 'IC 342 AGN', 'IC0342_X2': 'IC 342 X-2',\
                   'M82_S9': 'M 82 S9',\
                   'NGC5128_AGN': 'NGC 5128 AGN', 'NGC5128_X1': 'NGC 5128 ULX-1', 'NGC5128_X2': 'NGC 5128 ULX-2',\
                    'NGC7793_P13': 'NGC 7793 P13', 'NGC7793_P9': 'NGC 7793 P9'}

gal_obsid_dic = {}
gal_obstime_dic = {}

for gal_name in gal_list:
    gal_obsid_temp = []
    gal_obstime_temp = []
    for i in range(len(obsid_list)):
        obsid = obsid_list[i]
        fits_file = fits.open('data/' + obsid + '/a/a.img')
        object = fits_file[0].header['OBJECT']
        obstime = fits_file[0].header['DATE-OBS']
        if object=='2XMM J034615.6+681112':
            object = 'IC0342'
        if object==gal_name:
            #print(obsid, gal_name)
            #print(plot_num)
            gal_obsid_temp.append(obsid)
            gal_obstime_temp.append(Time(obstime, scale='tt').to_value('mjd') - 60000)
    gal_obsid_dic[gal_name] = gal_obsid_temp
    gal_obstime_dic[gal_name] = gal_obstime_temp

source_names = list(source_dic.keys())

with open('./results/selected_source_table.txt', 'w') as f:
    f.write('\\tablehead{ \n')
    f.write('\colhead{ULX name} & \colhead{Index} & \colhead{${\\rm log}(f_{\\rm 0.5-10~keV})$} & \colhead{${\\rm log}(L_{\\rm 0.5-10~keV})$} & \colhead{$\Gamma$} & \colhead{$\\rm N_H$} & \colhead{$\chi^2$/dof} \\\ \n')
    f.write('\colhead{ } & \colhead{ } & \colhead{($\\rm log( erg~s^{-1}~cm^{-2})$)} & \colhead{($\\rm log( erg~s^{-1})$)} & \colhead{ } & \colhead{($10^{22}~{\\rm cm^{-2}}$)} & \colhead{} \n')
    f.write('}\n')
    f.write('\startdata \n')

content = []
content.append(['ULX_name', 'ra', 'dec', 'out_name'])

for i in range(len(source_names)):
    source0 = source_dic[source_names[i]]
    gal_name = source_names[i].split('_')[0]
    distance = distance_list[gal_list==gal_name][0]
    galaxy_NH_value = galaxy_NH_value_list[galaxy_NH_galaxy_list==gal_name][0]
    obsid_list_temp = gal_obsid_dic[gal_name]
    obstime_list_temp = gal_obstime_dic[gal_name]
    flux_list_temp = []

    content.append([source_names[i], source0.ra.degree, source0.dec.degree, source_name_dic[source_names[i]]])

    ## HERE we need to load the fitting results of the combined image
    combined_source_list = glob.glob('combined_data/' + gal_name + '/a/x*.reg')
    for j in range(len(combined_source_list)):
        source_combined, ra_axis_combined, dec_axis_combined, angle_combined = get_coord_wcs(combined_source_list[j])
        sep_combined = source0.separation(source_combined).arcsec
        x_index = -1
        if sep_combined < 20:
            x_index = combined_source_list[j].split('/')[-1].split('.')[0]
            break
    if x_index==-1:
        print('No matched source in combined image for ' + source_names[i])
        continue
    
    fit_results = pd.read_csv('combined_data/' + gal_name + '/' + x_index + '/fxt_fit_results.csv')
    # fit_results = pd.read_csv('data/' + obsid_list_temp[j] + '/' + source_obs_list[k].split('/')[-1].split('.')[0] + '/fxt_fit_results.csv')
    log10flux, log10flux_down, log10flux_up = np.array(fit_results['log10flux']), np.array(fit_results['log10flux_down']), np.array(fit_results['log10flux_up'])
    photon_index, photon_index_down, photon_index_up = np.array(fit_results['photon_index']), np.array(fit_results['photon_index_down']), np.array(fit_results['photon_index_up'])
    nh, nh_down, nh_up = np.array(fit_results['nh']), np.array(fit_results['nh_down']), np.array(fit_results['nh_up'])
    chi2, dof, stat_method = np.array(fit_results['chi2']), np.array(fit_results['dof']), np.array(fit_results['stat_method'])
    ## using aprates to estimate the count rate, and utilize the fitted results to estimate the flux
    aprates_results_combine = pd.read_csv('combined_data/' + gal_name + '/' + x_index + '/aprates_rate.csv')
    aprates_inst = np.array(aprates_results_combine['inst'])
    count_rate_net_combine, count_rate_net_up_combine, count_rate_net_down_combine = np.array(aprates_results_combine['src_rate'])[aprates_inst=='a'],\
                                                                        np.array(aprates_results_combine['src_rate_err_up'])[aprates_inst=='a'], np.array(aprates_results_combine['src_rate_err_low'])[aprates_inst=='a']
    #flux_unabsorbed, count_rate_model = np.array(fit_results['flux_unabsorbed']), np.array(fit_results['count_rate_model'])

    fig = plt.figure(figsize=(18, 10))
    gs = gridspec.GridSpec(6, 2, wspace=0.1, hspace=0.0)

    ax_flux = fig.add_subplot(gs[0:6, :])
    # ax_index = fig.add_subplot(gs[2:4, :])
    # ax_NH = fig.add_subplot(gs[4:6, :])
    fs = 36 # font size
    for j in range(len(obsid_list_temp)):
        ## THe x2 sourc ein this observation shows very strange spectrum, skip it
        source_obs_list = glob.glob('data/' + obsid_list_temp[j] + '/a/x*.reg')
        #content_head = content_head + ',' + obsid_list[j]
        #print(content_head)
        if_in_obs = False
        for k in range(len(source_obs_list)):
            source1, ra_axis, dec_axis, angle = get_coord_wcs(source_obs_list[k])
            sep = source1.separation(source0).arcsec
            if sep < 20:
                if_in_obs = True
                obstime = obstime_list_temp[j]
                if glob.glob('data/' + obsid_list_temp[j] + '/' + source_obs_list[k].split('/')[-1].split('.')[0] + '/fxt_fit_results.csv')==[]:
                    print('No fit results for ' + source_names[i] + ' in ' + obsid_list_temp[j])
                    break

                aprates_results = pd.read_csv('data/' + obsid_list_temp[j] + '/' + source_obs_list[k].split('/')[-1].split('.')[0] + '/aprates_rate.csv')
                aprates_inst = np.array(aprates_results['inst'])
                ## inst a is utilize for estimate the flux 
                count_rate_net, count_rate_net_up, count_rate_net_down = np.array(aprates_results['src_rate'])[aprates_inst=='a'],\
                                                                        np.array(aprates_results['src_rate_err_up'])[aprates_inst=='a'], np.array(aprates_results['src_rate_err_low'])[aprates_inst=='a']

                count_rate_model = count_rate_net_combine
                flux_from_count_rate = count_rate_net * 10**(log10flux) / count_rate_model ## unabsorbed flux estimation from count rate, 0.5 - 10 keV
                flux_from_count_rate_up = count_rate_net_up * 10**(log10flux) / count_rate_model
                flux_from_count_rate_down = count_rate_net_down * 10**(log10flux) / count_rate_model

                print(source_names[i], source_obs_list[k], obstime, log10flux[0], photon_index[0], nh[0], chi2[0], dof[0])

                # if np.isnan(log10flux_down[0]) or log10flux_down[0]==-1 or log10flux_up[0]==0:
                    # if log10flux[0] < -16:## hard coding! attention!
                    #     print('This observation have something bad, need checking!')
                    #     continue
                #     ax_flux.scatter(obstime, flux_from_count_rate[0], s=80, color='black', marker='s')
                # else:
                ax_flux.scatter(obstime, flux_from_count_rate[0], s=80, color='blue', marker='s')
                ax_flux.plot([obstime, obstime], [flux_from_count_rate_down[0], flux_from_count_rate_up[0]], color='blue')
                ax_flux.scatter(obstime, log10flux[0], s=80, color='blue')
                flux_list_temp.append(flux_from_count_rate[0])
                # ax_flux.tick_params(which='both', labelsize=fs)
                # ax_flux.tick_params(axis='x', labelsize=0.1)
                # ax_flux.tick_params(which='both', direction='in')
                # ax_flux.tick_params(which='major', length=6)
                # ax_flux.minorticks_on()
                # ax_flux.tick_params(which='minor', length=3)
                # ax_flux.tick_params(which='both', top='on', right='on')
                # ax_flux.set_xlabel('MJD', fontsize=0.1)
                # ax_flux.set_ylabel(r'log($\rm Flux / erg~s^{-1}~cm^{-2}$)', fontsize=fs)

                # if np.isnan(photon_index_down[0]) or photon_index_down[0]==-1 or photon_index_up[0]==0:
                #     ax_index.scatter(obstime, photon_index[0], s=40, color='black', marker='s')
                # else:
                #     ax_index.scatter(obstime, photon_index[0], s=40, color='blue')
                #     ax_index.plot([obstime, obstime], [photon_index_down[0], photon_index_up[0]], color='blue')
                # ax_index.tick_params(which='both', labelsize=fs)
                # ax_index.tick_params(axis='x', labelsize=0.1)
                # ax_index.tick_params(which='both', direction='in')
                # ax_index.tick_params(which='major', length=6)
                # ax_index.minorticks_on()
                # ax_index.tick_params(which='minor', length=3)
                # ax_index.tick_params(which='both', top='on', right='on')
                # ax_index.set_xlabel('MJD', fontsize=0.1)
                # ax_index.set_ylabel(r'$\Gamma$', fontsize=fs)

                # if np.isnan(nh_down[0]) or nh_down[0]==-1 or nh_up[0]==0:
                #     ax_NH.scatter(obstime, nh[0]+galaxy_NH_value, s=40, color='black', marker='s')
                # else:
                #     ax_NH.scatter(obstime, nh[0]+galaxy_NH_value, s=40, color='blue')
                #     ax_NH.plot([obstime, obstime], [nh_down[0]+galaxy_NH_value, nh_up[0]+galaxy_NH_value], color='blue')
                # ax_NH.tick_params(which='both', labelsize=fs)
                # ax_NH.tick_params(which='both', direction='in')
                # ax_NH.tick_params(which='major', length=6)
                # ax_NH.minorticks_on()
                # ax_NH.tick_params(which='minor', length=3)
                # ax_NH.tick_params(which='both', top='on', right='on')
                # ax_NH.set_xlabel('MJD', fontsize=fs)
                # ax_NH.set_ylabel(r'${\rm N_H\times 10^{22}~cm^{-2}}$', fontsize=fs)
                break
        if not if_in_obs:
            # evt_file_path = glob.glob('data/' + obsid_list_temp[j] + '/fxt/products/fxt_a*')[0]
            # evt_file = fits.open(evt_file_path)
            img_file = fits.open('data/' + obsid_list_temp[j] + '/a/a.img')
            hdr, evt = img_file[0].header, img_file[0].data
            source_bkg, radius_bkg = get_coord_wcs('data/' + obsid_list_temp[j] + '/bkg.reg')
            w_a = WCS(hdr)
            x0, y0 = w_a.world_to_pixel(source0)
            x_bkg, y_bkg = w_a.world_to_pixel(source_bkg)
            ## estimate the number of photons
            aperture = CircularAperture((x0, y0), r=30 / pixel_size)
            bkg_radius = float(radius_bkg) / pixel_size
            aperture_bkg = CircularAperture((x_bkg, y_bkg), r=bkg_radius)
            phot_table = aperture_photometry(evt, aperture)
            phot_table_bkg = aperture_photometry(evt, aperture_bkg)
            count_number = phot_table['aperture_sum'][0]
            count_number_bkg = phot_table_bkg['aperture_sum'][0]
            net_count = count_number - count_number_bkg * (30/float(radius_bkg))**2
            
            upper_limit = True

            count_rate_bkg = count_number_bkg * (30/float(radius_bkg))**2 / hdr['EXPOSURE'] ## in a unit of counts/s
            bkg_flux_rate_file = pd.read_csv('data/' + obsid_list_temp[j] + '/bkg_spec/bkg_rate_flux.csv')
            bkg_flux_ref, bkg_count_rate_ref = bkg_flux_rate_file['flux'][0], bkg_flux_rate_file['count_rate'][0]

            ## the total count num of bkg in the reg
            count_num_bkg  = count_rate_bkg * hdr['EXPOSURE']
            bkg_count_rate_upper_limit  = (count_num_bkg**0.5) / hdr['EXPOSURE'] ## 1sigma upper limit
            ## extimate with the count rate of the source region
            net_count_rate_upper_limit = net_count / hdr['EXPOSURE']

            ## max of the two methods
            #flux_upperlimit = find_lambda_for_poisson_cdf(n=max(count_num_bkg, net_count), prob=0.9) / hdr['EXPOSURE'] * bkg_flux_ref / bkg_count_rate_ref ## 90% upper limit
            flux_upperlimit = find_n_for_poisson_cdf(lam=max(count_num_bkg, net_count), prob=0.9) / hdr['EXPOSURE'] * bkg_flux_ref / bkg_count_rate_ref ## 90% upper limit

            ax_flux.scatter(obstime_list_temp[j], flux_upperlimit, s=80, color='red', marker='v')
            flux_list_temp.append(flux_upperlimit)
            print(obsid_list_temp[j], ' No target source ' + source_names[i], obsid_list_temp[j])
            print('Number of photons is ', count_number, count_number_bkg, net_count, count_rate_bkg)
    interval = max(obstime_list_temp) - min(obstime_list_temp)
    ax_flux.set_xlim(min(obstime_list_temp)-interval*0.1, max(obstime_list_temp)+interval*0.1)
    ax_flux.tick_params(which='both', labelsize=fs)
    ax_flux.tick_params(axis='x', labelsize=fs)
    ax_flux.tick_params(which='both', direction='in')
    ax_flux.tick_params(which='major', length=12)
    ax_flux.minorticks_on()
    ax_flux.tick_params(which='minor', length=6)
    ax_flux.tick_params(which='both', top='on', right='on')
    ax_flux.set_xlabel('MJD $-$ 60000', fontsize=fs)
    ax_flux.set_ylabel(r'$\rm f_X / erg~s^{-1}~cm^{-2}$', fontsize=fs)
    ylim = [np.min(flux_list_temp)*10**(-0.3), np.max(flux_list_temp)*10**(0.3)]
    ax_flux.set_ylim(ylim[0], ylim[1])
    ax_flux.set_yscale('log')

    ## set the left y axis
    ax = plt.gca().twinx()
    ax.set_ylabel(r'$\rm L_X / erg~s^{-1}$', fontsize=fs)
    ax.set_ylim(ylim[0] * (float(distance)*3.09*10**18*10**6)**2*4*np.pi, ylim[1] * (float(distance)*3.09*10**18*10**6)**2*4*np.pi)
    ax.tick_params(which='both', labelsize=fs)
    ax.tick_params(which='both', direction='in')
    ax.tick_params(which='major', length=12)
    ax.minorticks_on()
    ax.tick_params(which='minor', length=6)
    ax.set_yscale('log')

        # ax_index.set_xlim(min(obstime_list_temp)-interval*0.1, max(obstime_list_temp)+interval*0.1)
        # ax_NH.set_xlim(min(obstime_list_temp)-interval*0.1, max(obstime_list_temp)+interval*0.1)
    ax_flux.text(0.65, 0.9, source_name_dic[source_names[i]], fontsize=fs, color='black', style='normal', weight='bold', clip_on=True, transform=ax_flux.transAxes)
    plt.rcParams['savefig.bbox'] = 'tight'
    plt.savefig('figure/selected_ULX/' + source_names[i] + '.pdf')
    plt.clf()


####\colhead{ULX name} & \colhead{Index} & \colhead{${\\rm log}(f_{\\rm 0.5-10~keV})$} & \colhead{${\\rm log}(L_{\\rm 0.5-10~keV})$} & \colhead{$\Gamma$} & \colhead{$\\rm N_H$} & \colhead{STAT} 
    ## Here we generate a table for the selected results
    log10flux_out = output_sci_results(log10flux[0], [log10flux[0]-log10flux_down[0], log10flux_up[0]-log10flux[0]])
    log10lumin_out = output_sci_results(log10flux[0] + np.log10((float(distance)*3.09*10**18*10**6)**2*4*np.pi), \
                                        [(log10flux[0]-log10flux_down[0]),\
                                            (log10flux_up[0]-log10flux[0])])
    photon_index_out = output_sci_results(photon_index[0], [photon_index[0]-photon_index_down[0], photon_index_up[0]-photon_index[0]])
    if ~np.isnan(nh[0]):
        if ~np.isnan(nh_down[0]):
            if nh_down[0]!=-1:
                nh_out = format(galaxy_NH_value, '.2f') + ' (fixed)' + ' + ' + output_sci_results(nh[0], [nh[0]-nh_down[0], nh_up[0]-nh[0]])
                ## wrong with the NH estimation in the fitting 
                if (nh_down[0]==0) and (nh_up[0]==0):
                    nh_out = format(galaxy_NH_value, '.2f') + ' (fixed)'
                elif (nh[0]-nh_down[0] < 0.001) and (nh[0]<0.001) and (nh_up[0]>nh[0]):
                    if '.' in output_err_str(nh_up[0]-nh[0]):
                        decimal_places = len(output_err_str(nh_up[0]-nh[0]).split('.')[1])
                    else:
                        decimal_places = 0
                    nh_out = format(galaxy_NH_value, '.2f') + ' (fixed)' + ' + $' + format(nh[0], '.' + str(decimal_places) + 'f') + '^{+' + output_err_str(nh_up[0]-nh[0]) + '}$'
                if ('--' in nh_out):
                    nh_out = format(galaxy_NH_value, '.2f') + ' (fixed)'
            else:
                nh_out = format(galaxy_NH_value, '.2f') + ' (fixed)'
        else:
            nh_out = format(galaxy_NH_value, '.2f') + ' (fixed)'
        if ~np.isnan(chi2[0]):
            chi2_dof_out = format(chi2[0], '.1f') + '/' + str(dof[0])
        else:
            chi2_dof_out = '/' + str(dof[0])
    with open('./results/selected_source_table.txt', 'a') as f:
        f.write(source_name_dic[source_names[i]] + ' & ' + x_index.replace('x', 'S') + ' & ' + log10flux_out + ' & ' + log10lumin_out + ' & ' + photon_index_out + ' & ' +\
                nh_out + ' & ' + chi2_dof_out +  '\\\ \n')

with open('./files/selected_source.csv', 'w') as f:
    writer = csv.writer(f)
    writer.writerows(content)
