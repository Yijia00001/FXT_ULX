## fit the spectrum using pyXSPEC

import glob, os, sys, csv
import numpy as np
import matplotlib.pyplot as plt
import xspec
from astropy.io import fits
import pandas as pd

galaxy_NH_file = pd.read_csv('/home/yijia/work/EP_ULX/combined_data/galaxt_NH.csv')
img_file = glob.glob('../a/a.img')[0]
fits_file = fits.open(img_file)
gal_name = fits_file[0].header['OBJECT']
if gal_name=='2XMM J034615.6+681112':
    gal_name = 'IC0342'

galaxy_NH_galaxy_list = np.array(galaxy_NH_file['galaxy'])
galaxy_NH_value_list = np.array(galaxy_NH_file['NH'])
#script_path = sys.argv[0]
galaxy_NH_value = galaxy_NH_value_list[galaxy_NH_galaxy_list==gal_name][0]

#source_file_path_list = glob.glob('./a/x*.reg')
x_source_index = glob.glob('./x*_a_grp.pi')[0].split('_')[0]
#for i in range(len(source_file_path_list)):
content = []
content.append(['log10flux', 'log10flux_down', 'log10flux_up', 'photon_index', 'photon_index_down', 'photon_index_up',\
                 'nh', 'nh_down', 'nh_up', 'chi2', 'dof', 'stat_method',\
                  'flux_unabsorbed', 'count_rate_flux', 'count_rate_model', 'count_rate_net', 'count_rate_net_err',\
                  'count_rate_flux_low', 'count_rate_net_low', 'count_rate_net_low_err',\
                  'count_rate_flux_high', 'count_rate_net_high', 'count_rate_net_high_err'])

data = xspec.AllData('1:1 ' + x_source_index + '_a_grp.pi' + ' 2:2 ' + x_source_index + '_b_grp.pi')

xspec.AllData.ignore("1:**-0.5 1.-**")
xspec.AllData.ignore("2:**-0.5 1.-**")
xspec.AllData.ignore("bad")
spectrum = xspec.AllData(1)
rate_reuslt = spectrum.rate
#count_rate_model_low = rate_reuslt[3] ## model predict rate
count_rate_net_low = rate_reuslt[0]
count_rate_net_low_err = rate_reuslt[1]
xspec.AllData.clear()

data = xspec.AllData('1:1 ' + x_source_index + '_a_grp.pi' + ' 2:2 ' + x_source_index + '_b_grp.pi')
xspec.AllData.ignore("1:**-1. 10.-**")
xspec.AllData.ignore("2:**-1. 10.-**")
xspec.AllData.ignore("bad")
spectrum = xspec.AllData(1)
rate_reuslt = spectrum.rate
#count_rate_model_high = rate_reuslt[3] ## model predict rate
count_rate_net_high = rate_reuslt[0]
count_rate_net_high_err = rate_reuslt[1]
xspec.AllData.clear()

data = xspec.AllData('1:1 ' + x_source_index + '_a_grp.pi' + ' 2:2 ' + x_source_index + '_b_grp.pi')
xspec.AllData.ignore("1:**-0.5 10.-**")
xspec.AllData.ignore("2:**-0.5 10.-**")
xspec.AllData.ignore("bad")

m1 = xspec.Model("const*TBabs*Tbabs*cflux*powerlaw")
m1 = xspec.AllModels(1)
#print(dir(m1.TBabs_3))
par_norm = m1.powerlaw.norm
par_norm.frozen = True
nh_2 = m1.TBabs_3.nH
nh_2.values = 0.0
nh = m1.TBabs.nH
nh.values = 0.0
photon_index= m1.powerlaw.PhoIndex
photon_index.values = 2
xspec.AllModels.calcFlux("0.5 10.0")
spectrum = xspec.AllData(1)
flux_result = spectrum.flux
## unabsorbed flux
flux_unabsorbed = flux_result[0]
nh_2.values = galaxy_NH_value
nh_2.frozen = True

m2 = xspec.AllModels(2)
const2 = m2.constant.factor
const2.untie()
const2.frozen = True
xspec.Xset.abund = "wilm"
chi2 = xspec.Fit.statistic
dof = xspec.Fit.dof
xspec.Fit.query = "yes"

mincts_a = np.loadtxt('./mincts_a')
mincts_b = np.loadtxt('./mincts_b')
stat_method = 'chi'
if mincts_a<5 or mincts_b<5:
    xspec.Fit.statMethod = 'cstat'
    stat_method = 'cstat'

xspec.AllModels.calcFlux("0.5 10.0")
spectrum = xspec.AllData(1)
flux_result = spectrum.flux
count_rate_flux = flux_result[3]
rate_reuslt = spectrum.rate
count_rate_model = rate_reuslt[3] ## model predict rate
count_rate_net = rate_reuslt[0]
count_rate_net_err = rate_reuslt[1]
## flux_unaborbed, count_rate_flux, count_rate
## erg cm-2 s-1, counts cm-2 s-1, counts s-1
## If one would like to estimate the flux using count rate
## flux = flux_unabsorbed / count_rate_model * count_rate_net

## estimatte the 0.5 - 2.0 keV count flux
xspec.AllModels.calcFlux("0.5 2.0")
spectrum = xspec.AllData(1)
flux_result = spectrum.flux
count_rate_flux_low = flux_result[3]

## estimatte the 2.0 - 10.0 keV count rate flux
xspec.AllModels.calcFlux("2.0 10.0")
spectrum = xspec.AllData(1)
flux_result = spectrum.flux
count_rate_flux_high = flux_result[3]

if dof <= 1 or stat_method=='cstat':
    print("Not enough data points to fit the spectrum.")
    ## fit the first NH component to be 0
    nh = m1.TBabs.nH
    nh.values = 0.0
    nh.frozen = True
    #############################
    # photon_index= m1.powerlaw.PhoIndex
    # photon_index.values = 2
    # photon_index.frozen = True
    const1 = m1.constant.factor
    const1.frozen = True
    chi2 = xspec.Fit.statistic
    dof = xspec.Fit.dof

    if dof >= 0:
        xspec.Fit.perform()
        chi2 = xspec.Fit.statistic
        dof = xspec.Fit.dof
        log10flux = m1.cflux.lg10Flux.values
        log10flux_err = ['nan', 'nan']
        photon_index_err = ['nan', 'nan']
        if (chi2 >= 2*dof) or (dof==0):
            log10flux_err = [-1, -1]
        else:
            xspec.Fit.error('6, 7')
            log10flux_err = m1.cflux.lg10Flux.error
            photon_index_err = m1.powerlaw.PhoIndex.error
        nh = m1.TBabs.nH.values
        photon_index = m1.powerlaw.PhoIndex.values
        nh_err = ['nan', 'nan']
    else:
        photon_index= m1.powerlaw.PhoIndex
        photon_index.values = 2
        photon_index.frozen = True
        xspec.Fit.perform()
        chi2 = xspec.Fit.statistic
        dof = xspec.Fit.dof

        log10flux = m1.cflux.lg10Flux.values
        log10flux_err = ['nan', 'nan']
        xspec.Fit.error('6')
        log10flux_err = m1.cflux.lg10Flux.error
        photon_index_err = ['nan', 'nan']
        nh_err = ['nan', 'nan']
        if dof < 0:
            log10flux == ['nan']
            photon_index = ['nan']
            nh = ['nan']
            chi2 = 'nan'
            dof = xspec.Fit.do
else:
    xspec.Fit.query = "yes"
    #xspec.Fit.nIterations = 5000
    xspec.Fit.perform()
    chi2 = xspec.Fit.statistic
    dof = xspec.Fit.dof
    fitting_time = 0
    while(1):
        if fitting_time > 5:
            nh = m1.TBabs.nH
            nh.values = 0.0
            nh.frozen = True
            xspec.Fit.perform()
            chi2 = xspec.Fit.statistic
            dof = xspec.Fit.dof
            #############################
            if chi2 >= 2*dof:
                log10flux = m1.cflux.lg10Flux.values
                log10flux_err = [-1, -1]
                nh = m1.TBabs.nH.values
                nh_err = [-1, -1]
                photon_index = m1.powerlaw.PhoIndex.values
                photon_index_err = [-1, -1]
                break
            else:
                xspec.Fit.error('6, 7')
                log10flux = m1.cflux.lg10Flux.values
                log10flux_err = m1.cflux.lg10Flux.error
                nh = m1.TBabs.nH.values
                nh_err = [-1, -1]
                photon_index = m1.powerlaw.PhoIndex.values
                photon_index_err = m1.powerlaw.PhoIndex.error
                break
        if chi2 >= 2*dof:
            xspec.Fit.perform()
            chi2 = xspec.Fit.statistic
            dof = xspec.Fit.dof
            fitting_time += 1
        else:
            xspec.Fit.error('2')
            nh = m1.TBabs.nH.values
            nh_err = m1.TBabs.nH.error
            ## set the nh to be frozen if estimation of nh error is wrong
            ## in some fitting, the NH could not be restricted
            if nh_err[0] > nh_err[1]:
                m1.TBabs.nH.frozen = True
                nh_err = m1.TBabs.nH.error
            xspec.Fit.perform()
            xspec.Fit.error('2, 6, 7')
            log10flux = m1.cflux.lg10Flux.values
            log10flux_err = m1.cflux.lg10Flux.error
            #nh = m1.TBabs.nH.values
            photon_index = m1.powerlaw.PhoIndex.values
            photon_index_err = m1.powerlaw.PhoIndex.error
            
            print("Log10 Flux: ", log10flux[0], log10flux_err[0], log10flux_err[1])
            print('Photon Index: ', photon_index[0], photon_index_err[0], photon_index_err[1])
            print('NH: ', nh[0], nh_err[0], nh_err[1])
            break
    chi2 = xspec.Fit.statistic
    dof = xspec.Fit.dof
content.append([log10flux[0], log10flux_err[0], log10flux_err[1], \
                photon_index[0], photon_index_err[0], photon_index_err[1],\
                nh[0], nh_err[0], nh_err[1], \
                chi2, dof, stat_method, flux_unabsorbed, count_rate_flux, count_rate_model, count_rate_net, count_rate_net_err,\
                count_rate_flux_low, count_rate_net_low, count_rate_net_low_err,\
                count_rate_flux_high, count_rate_net_high, count_rate_net_high_err])

## save all
xcm_file_name = x_source_index + ".xcm"
if os.path.exists(xcm_file_name):
    os.remove(xcm_file_name)
xspec.Xset.save(x_source_index + ".xcm", info='a')


## plot
xspec.Plot.device = '/null'
xspec.Plot('eeuf')

xVals = xspec.Plot.x()
yVals = xspec.Plot.y()
xErrs = xspec.Plot.xErr()
yErrs = xspec.Plot.yErr()
unfolded = xspec.Plot.model()
header = 'xVals yVals xErr yErr model'
np.savetxt(f'data_tbabs_pl_A.txt', np.array([xVals, yVals, xErrs, yErrs, unfolded]).T, header=header)

xVals2 = xspec.Plot.x(2)
yVals2 = xspec.Plot.y(2) 
xErrs2 = xspec.Plot.xErr(2)
yErrs2 = xspec.Plot.yErr(2)
unfolded = xspec.Plot.model(2)
np.savetxt(f'data_tbabs_pl_B.txt', np.array([xVals2, yVals2, xErrs2, yErrs2, unfolded]).T, header=header)

with open('fxt_fit_results.csv', 'w') as f:
    writer = csv.writer(f)
    writer.writerows(content)
print(flux_result)
print(rate_reuslt)
# xrt_observations = glob.glob('/home/yijia/work/EP_AGN/mrk1044/swift/*/xrt')

# for i in range(len(xrt_observations)):
#     xrt_observation = xrt_observations[i]
#     spec_file = xrt_observation + '/x1_grp.pha'
#     print(xrt_observation)
#     data = xspec.AllData(f" {spec_file}")

#     break
