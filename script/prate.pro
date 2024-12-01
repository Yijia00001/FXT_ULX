pro prate, saveps=saveps
; find net count rate and scale it to intrinsic luminosity

include_82125 = 0

if include_82125 then files = file_search('./*/nuc.pi') else files = file_search('./*33469*/nuc.pi')

n_obs = n_elements(files)
date = dblarr(n_obs)
rate = dblarr(n_obs, 3)

; first load the RMF file and get channel range for 0.3-10 keV
; 0.3-10 keV: 30-999 (start from 0) or 31-1000 (natural row number)
; checked that all RMF files are the same

for i = 0L, n_obs - 1L do begin
	fsrc = files[i]
	dir = file_dirname(fsrc)
	fbkg = dir + '/nuc_bkg.pi'
	
	dsrc = mrdfits(fsrc, 1, shdr, /silent)
	dbkg = mrdfits(fbkg, 1, bhdr, /silent)
	exptime = sxpar(shdr, 'EXPOSURE')
	
	scales = sxpar(shdr, 'BACKSCAL')
	scaleb = sxpar(bhdr, 'BACKSCAL')
	f = scales / scaleb
	
	csrc = total((dsrc.counts)[30:999])
	cbkg = total((dbkg.counts)[30:999])
	
	aprates, csrc, cbkg, scales, scaleb, s0, err0, err1, CL=0.68
	; calculate the 90% upper limit for non-detection
	if err0 eq 0 then aprates, csrc, cbkg, scales, scaleb, s0, err0, err1, CL=0.90
	
	rate[i, *] = [s0, err0, err1] / exptime
	; print, 'SRC, BKG:', csrc, cbkg, s0, err0, err1, format='(A10, I4, I4, F10.2, F10.2, F10.2)'
	print, 'Rate:', csrc/exptime, cbkg/exptime, format='(A10, F15.4, F15.4)'
	
	; get MJD at the mid-point of the observation
	date0 = sxpar(shdr, 'DATE-OBS')
	date1 = sxpar(shdr, 'DATE-END')
	date[i] = (date_conv(date0, 'MODIFIED') + date_conv(date1, 'MODIFIED')) / 2d
endfor

if keyword_set(saveps) then begin
	current_device = !d.name ; save current device
	set_plot, 'ps'
	epsname = 'ngc247nuc_lc1.eps'
	device, filename=epsname, /encapsulated
endif

; convert from counts/s to unabsorbed flux erg/cm2/s
rate2flux = 8.2705d-13 / 1.07158d-2
; convert from counts/s to intrinsic erg/s
rate2lum = rate2flux * f2l(3.4d)

; flux = rate * rate2flux
flux = rate * rate2lum

; xr = [0, date[-1] - date[0] + 7] + round((date[0] - 56000)/10.0)*10

if include_82125 then xr = [650, 1200] $ ; include the first point in 2014-02
else xr = [800, 1200] ; [800, 1080]

yr = [0, 3.1]
cs = 2
fnt = 1
co = [cgcolor('red'), cgcolor('orange'), cgcolor('dark green')]

plot, xr, /nodata, charsize=cs, font=fnt, $
	xstyle=1, xrange=xr, ystyle=1, yrange=yr, $
	xmargin=[6, 2], ymargin=[3, 1], $
	xtitle='MJD - 56000', ytitle='Luminosity (10!U39!N erg s!U-1!N)'

xx = date - 56000 ; MJD
yy = flux[*, 0] * 1d-39
y0 = flux[*, 1] * 1d-39
y1 = flux[*, 2] * 1d-39

qdet = where(y0 gt 0)
qund = where(y0 eq 0)

psize = 1.5
xsym = [0, 0.5, 0, -0.5,  0,  0,  0.5,   0, -0.5] * psize
ysym = [0,   0, 0,    0,  0, -3, -2.4,  -3, -2.4] * psize
usersym, xsym, ysym
oplot, xx[qund], y1[qund], psym=8, color=co[1]

plotsym, 0, /fill
oplot, xx[qdet], yy[qdet], psym=8, color=co[1]
errplot, xx[qdet], y0[qdet], y1[qdet], color=co[1]

; stop

; add XMM and Chandra points
; flux corrected for absorption, err 1.0 6
xmm_flux = 10d ^ [-11.8433, -11.8523, -11.8343] * f2l(3.4d)
; DATE-OBS= '2014-07-01T05:01:03' / Start Time (UTC) of exposure
; DATE-END= '2014-07-01T14:24:10' / End Time (UTC) of exposure
xmm_date = (date_conv('2014-07-01T05:01:03', 'MODIFIED') + $
				date_conv('2014-07-01T14:24:10', 'MODIFIED')) / 2d

xx = xmm_date - 56000
yy = xmm_flux[0] * 1d-39
y0 = xmm_flux[1] * 1d-39
y1 = xmm_flux[2] * 1d-39
plots, xx, yy, psym=8, color=co[0]
errplot, xx, y0, y1, color=co[0]

; nuc_po.xcm err 1.0 5
; -11.7468 -11.7899 -11.6944
; nuc_bp.xcm, err 1.0 5
; -11.8443 -11.8739 -11.816
; nuc_bp2.xcm, err 1.0 5
; -11.8839 -11.9055 -11.8633

cha_flux = 10d ^ [-11.8839, -11.9055, -11.8633] * f2l(3.4d)
; DATE-OBS= '2014-11-12T05:34:48' / Observation start date
; DATE-END= '2014-11-12T07:30:13' / Observation end date
cha_date = (date_conv('2014-11-12T05:34:48', 'MODIFIED') + $
				date_conv('2014-11-12T07:30:13', 'MODIFIED')) / 2d

xx = cha_date - 56000
yy = cha_flux[0] * 1d-39
y0 = cha_flux[1] * 1d-39
y1 = cha_flux[2] * 1d-39
plots, xx, yy, psym=8, color=co[2]
errplot, xx, y0, y1, color=co[2]

al_legend, ['XMM-Newton', 'Swift', 'Chandra'], psym=[8,8,8], color=co, $
	charsize=cs-0.5, font=fnt, box=0

if keyword_set(saveps) then begin
	device, /close
	chepslinewidth, epsname, 40
	set_plot, current_device
endif

stop

end
