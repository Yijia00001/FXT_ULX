# create a text grppha file to rebin the spectrum 
# with at least mincts counts/ch from the energy enlow

from astropy.io import fits
import sys, os, glob
import numpy as np

#exec(open("func.py").read())

#source_name = 'x1'

obsid_list = []
with open('./obsid_list.txt', 'r') as f:
    for line in f:
        obsid_list.append(line.strip())

finish_obsid_list = np.loadtxt('finish.txt', dtype=int)

#gal_list = ['NGC1313', 'NGC5128', 'NGC247', 'NGC7793', 'NGC55', 'NGC0224', 'NGC1316', 'NGC4038', 'NGC0925', 'IC0342', 'M82']

for obsid in obsid_list:
    if int(obsid) in finish_obsid_list:
        print('skip finished obsid: ', obsid)
        continue

    source_name_path = glob.glob('./data/' + obsid + '/a/' + 'x*.reg')
    print(source_name_path)
    for i in range(len(source_name_path)):
        source_name = source_name_path[i].split('/')[-1].split('.')[0]
        for inst in ['a', 'b', 'bkg']:
            mincts = 15
            # fsrc = obsid + '/ener_res/' + 'x1.pi'
            # frmf = obsid + '/ener_res/' + 'x1.rmf'
            if inst=='bkg':
                fsrc = 'data/' + obsid + '/bkg_spec/bkg_a.pi'
                frmf = 'data/' + obsid + '/bkg_spec/bkg_a.rmf'
                fgrp_name = 'data/' + obsid + '/bkg_spec/bkg_a.grp'
            else:
                fsrc = 'data/' + obsid + '/' + source_name + '/' + source_name + '_' + inst + '.pi'
                frmf = 'data/' + obsid + '/' + source_name + '/' + source_name + '_' + inst  + '.rmf'
                fgrp_name = 'data/' + obsid + '/' + source_name + '/' + source_name + '_' + inst + '.grp'
            enlow = 0.51
            
            hdulist = fits.open(fsrc)
            hdu_spec = hdulist[1]
            src_ch = hdu_spec.data['CHANNEL']
            src_cts = hdu_spec.data['COUNTS']
            hdulist.close()
            
            hdulist = fits.open(frmf)
            hdu_rmf = hdulist[2]
            channel = hdu_rmf.data['CHANNEL']
            e_min = hdu_rmf.data['E_MIN']
            e_max = hdu_rmf.data['E_MAX']
            hdulist.close()
            
            chen = (e_min + e_max) / 2.0 # energy at each AD channel
            nch = len(channel) # number of AD channels
            
            q = np.where(e_min >= enlow)
            ch = q[0][0]
            q_max = np.where(e_max <= 10.0)
            ch_max = q_max[-1][-1]
            
            fgrp = open(fgrp_name, "w")
            fgrp.write("%5d %5d %5d\n" % (channel[0], ch - 1, ch - channel[0]))
            
            # first and last channel no.
            ch0 = channel[0] # ch0 may be 0 or 1
            maxch = channel[-1]
            total_counts = src_cts[ch-ch0:ch_max-ch0+1].sum()
            if total_counts<30:
                mincts = 3
            with open('data/' + obsid + '/' + source_name + '/mincts_' + inst, 'w') as f_temp:
                f_temp.write(str(mincts) + '\n')
            
            while ch <= maxch:
                dch = 1
                ch1 = ch + dch - 1
                if ch1 > maxch:
                    ch1 = maxch
                    dch = maxch - ch + 1
                
                while (ch1 < maxch) and (src_cts[ch-ch0:ch1-ch0+1].sum() < mincts):
                    dch = dch + 1
                    ch1 = ch + dch - 1
                    
                fgrp.write("%5d %5d %5d\n" % (ch, ch1, dch))
                # print("%d, %d, %d" % (ch, ch1, src_cts[ch-ch0:ch1-ch0+1].sum()))
                ch = ch + dch
            
            fgrp.close()
            print("Done! please run grppha to apply the binning.")
        
