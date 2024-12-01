readcol, 'flux.txt', format='A,A,D,D,D', date0, date1, flux, flux1, flux2
print, 'Observation Date:  ', date0, format='(A32,A)'
print, 'Flux (erg/cm2/s):', 10d^flux, format='(A30, E10.2)'
print, '90% Error Range (erg/cm2/s):', 10d^flux1, ' - ', 10d^flux2, format='(A30, E10.2, A, E8.2)'
