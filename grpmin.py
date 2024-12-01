# create a text grppha file to rebin the spectrum 
# with at least mincts counts/ch from the energy enlow

from astropy.io import fits
import sys, os, glob
import numpy as np

#exec(open("func.py").read())

#source_name = 'x1'
mincts = 15

# if len(sys.argv)>1:
#     obsid = sys.argv[1]
# if len(sys.argv)>2:
#     source_name = sys.argv[2]
# if len(sys.argv)>3:
#     mincts = int(sys.argv[3])

# obsid_list = get_lines_from_file('source_ra_dec.txt')
# obsid_list = ['893_2_0_3761.reg']
obsid_list = []
with open('./obsid_list.txt', 'r') as f:
    for line in f:
        obsid_list.append(line.strip())

for obsid in obsid_list:
    source_name_path = glob.glob('./data/' + obsid + '/a/' + 'x*.reg')
    print(source_name_path)
    for i in range(len(source_name_path)):
        source_name = source_name_path[i].split('/')[-1].split('.')[0]
        for inst in ['a', 'b']:
            # fsrc = obsid + '/ener_res/' + 'x1.pi'
            # frmf = obsid + '/ener_res/' + 'x1.rmf'
            fsrc = 'data/' + obsid + '/' + inst + '/' + source_name + '/' + source_name + '.pi'
            frmf = 'data/' + obsid + '/' + inst + '/' + source_name + '/' + source_name + '.rmf'
            #mincts = 15
            enlow = 0.51
            
            #~ args = sys.argv
            
            #~ if len(args) == 5:
                #~ fsrc = args[1]
                #~ frmf = args[2]
                #~ mincts = np.float(args[3])
                #~ enlow = np.float(args[4])
                #~ if not os.path.exists(fsrc):
                    #~ print("*** File %s not exist ***" & fsrc)
                #~ if not os.path.exists(frmf):
                    #~ print("*** File %s not exist ***" & frmf)    
            #~ else:
                #~ print("Usage:")
                #~ print("   grpmin f.pi f.rmf 25 0.3")
                #~ sys.exit()
            
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
            
            fgrp_name = fsrc.split('.')[0] + '.grp'
            
            fgrp = open(fgrp_name, "w")
            fgrp.write("%5d %5d %5d\n" % (channel[0], ch - 1, ch - channel[0]))
            
            # first and last channel no.
            ch0 = channel[0] # ch0 may be 0 or 1
            maxch = channel[-1]
            
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
        
