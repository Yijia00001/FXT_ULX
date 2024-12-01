pro lcforall, saveps=saveps
; plot the lightcurve using all of the data

files = file_search('./*33469*/ulx.pi')

n_obs = n_elements(files)
date = dblarr(n_obs)
rate = dblarr(n_obs)
rate_err = dblarr(n_obs)
flconv = f2l(3.4d) / 1d39

; first load the RMF file and get channel range for 0.3-10 keV
; E_MIN: channel
; 0.2 keV: 20
; 0.3 keV: 30
; 0.75 keV: 75
; 1.0 keV: 100
; 1.5 keV: 150
; 2 keV: 200
; 10 keV: 1000
; checked that all RMF files are the same

ebound = [30, 200] ; energy bounds

for i = 0L, n_obs - 1L do begin
	fsrc = files[i]
	dir = file_dirname(fsrc)
	fbkg = dir + '/ulx_bkg.pi'
	
	dsrc = mrdfits(fsrc, 1, shdr, /silent)
	dbkg = mrdfits(fbkg, 1, bhdr, /silent)
	exptime = sxpar(shdr, 'EXPOSURE')
	
	scales = sxpar(shdr, 'BACKSCAL')
	scaleb = sxpar(bhdr, 'BACKSCAL')
	f = scales / scaleb
	
	ss = dsrc.counts
	bb = dbkg.counts
		
	csrc = total(ss[ebound[0]:ebound[1] - 1L])
	cbkg = total(bb[ebound[0]:ebound[1] - 1L])
	
	rate[i] = (csrc - cbkg * f) / exptime
	var = csrc + cbkg * f^2d
	rate_err[i] = sqrt(var) / exptime
	print, dir + ' Total/back counts: ', csrc, cbkg * f, $
		format='(A, I5, F10.4)'

	; get MJD at the mid-point of the observation
	date0 = sxpar(shdr, 'DATE-OBS')
	date1 = sxpar(shdr, 'DATE-END')
	date[i] = (date_conv(date0, 'MODIFIED') + date_conv(date1, 'MODIFIED')) / 2d
endfor

; get combined count rate
cmfiles = file_search('./combine/*/ulx.pi')

cmn_obs = n_elements(cmfiles)
cmrate = dblarr(cmn_obs)

for i = 0L, cmn_obs - 1L do begin
	fsrc = cmfiles[i]
	dir = file_dirname(fsrc)
	fbkg = dir + '/ulx_bkg.pi'
	
	dsrc = mrdfits(fsrc, 1, shdr, /silent)
	dbkg = mrdfits(fbkg, 1, bhdr, /silent)
	exptime = sxpar(shdr, 'EXPOSURE')
	
	scales = sxpar(shdr, 'BACKSCAL')
	scaleb = sxpar(bhdr, 'BACKSCAL')
	f = scales / scaleb
	
	ss = dsrc.counts
	bb = dbkg.counts
		
	csrc = total(ss[ebound[0]:ebound[1] - 1L])
	cbkg = total(bb[ebound[0]:ebound[1] - 1L])
	
	cmrate[i] = (csrc - cbkg * f) / exptime
	print, dir + ' Total/back counts: ', csrc, cbkg * f, $
		format='(A, I5, F10.4)'
endfor


; convert form observed rate to 0.3-2 keV luminosity
; group, intrinsic luminosity
; 1, 01/02/06, 4.2888e-12
; 2, 03/04/10/26, 4.9029e-12
; 3, 05/15/16/20/22, 2.3082e-12
; 4, 07/08/21/24/25/28, 5.6661e-12
; 5, 09/11/12/13/14/17/19/23/27, 2.1084e-12

xrt_lum = rate * 0.0

grp = [1, 2, 6] - 1
xrt_lum[grp] = 4.2888d-12 * flconv * rate[grp] / total(cmrate[0])

grp = [3, 4, 10, 26] - 1
xrt_lum[grp] = 4.9029d-12 * flconv * rate[grp] / total(cmrate[1])

grp = [05, 15, 16, 20, 22] - 1
xrt_lum[grp] = 2.3082d-12 * flconv * rate[grp] / total(cmrate[2])

grp = [07, 08, 18, 21, 24, 25, 28] - 1
xrt_lum[grp] = 5.6661d-12 * flconv * rate[grp] / total(cmrate[3])

grp = [09, 11, 12, 13, 14, 17, 19, 23, 27] - 1
xrt_lum[grp] = 2.1084d-12 * flconv * rate[grp] / total(cmrate[4])

xrt_lumerr = rate_err / rate * xrt_lum

if keyword_set(saveps) then begin
	current_device = !d.name ; save current device
	set_plot, 'ps'
	epsname = 'ulx_lcall.eps'
	device, filename=epsname, /encapsulated
endif

xr1 = [100, 700d]
xr2 = [1800, 2400d]

dxr1 = xr1[1] - xr1[0]
dxr2 = xr2[1] - xr2[0]

px0 = 0.13
py0 = 0.17
px1 = 0.95
py1 = 0.95

pgap = (px1 - px0) * 0.05
pdx1 = (px1 - px0) * 0.95 * dxr1 / (dxr1 + dxr2)
pdx2 = (px1 - px0) * 0.95 * dxr2 / (dxr1 + dxr2)

pos1 = [px0, py0, px0 + pdx1, py1]
pos2 = [px0 + pdx1 + pgap, py0, px1, py1]

yr = [1e-1, 1e2] & ylog = 1
yr = [0, 16] & ylog = 0
cs = 2
fnt = 1
co = [cgcolor('blue'), cgcolor('red'), cgcolor('dark green')]


plot, xr2, ylog=ylog, /nodata, charsize=cs, font=fnt, position=pos2, $
	xstyle=1, xrange=xr2, ystyle=1, yrange=yr, ytickname=replicate(' ', 30), $
	xminor=2, xtickinterval=200


MJD0 = 55000L
xx = date - MJD0 ; MJD
yy = xrt_lum
yerr = xrt_lumerr
yy0 = yy - yerr
yy1 = yy + yerr

for i = 0, n_elements(xx) - 1 do begin
	oplot, xx[i]*[1,1], [yy0[i], yy1[i]], color=co[0]
	plotsym, 0, 1.1, /fill
	plots, xx[i], yy[i], psym=8, color=co[0]
	plotsym, 0, 0.8, /fill
	plots, xx[i], yy[i], psym=8, color=cgcolor('background')

endfor

;#############################################################
; xmm data
; 07 high: -11.0963 -11.2278     -10.9513
; 07 low: -11.6505 -11.6892     -11.6132
; 06: -11.4282 -11.6324     -11.1896
; DATE-OBS= '2009-12-27T19:41:35' / Start Time (UTC) of exposure
; DATE-END= '2009-12-28T04:50:30' / End Time (UTC) of exposure
; DATE-OBS= '2014-07-01T05:01:03' / Start Time (UTC) of exposure
; DATE-END= '2014-07-01T14:24:10' / End Time (UTC) of exposure

xmm_date = dblarr(3)
date0 = '2014-07-01T05:01:03'
date1 = '2014-07-01T14:24:10'
xmm_date[0] = (date_conv(date0, 'MODIFIED') + date_conv(date1, 'MODIFIED')) / 2d
xmm_date[1] = xmm_date[0]
date0 = '2009-12-27T19:41:35'
date1 = '2009-12-28T04:50:30'
xmm_date[2] = (date_conv(date0, 'MODIFIED') + date_conv(date1, 'MODIFIED')) / 2d

xmm_lum = 10d ^ [-11.0963, -11.6505, -11.4282] * flconv
xmm_lum0 = 10d ^ [-11.2278, -11.6892, -11.6324] * flconv
xmm_lum1 = 10d ^ [-10.9513, -11.6132, -11.1896] * flconv

xx = xmm_date - MJD0
yy = xmm_lum
yy0 = xmm_lum0
yy1 = xmm_lum1

for i = 0, n_elements(xx) - 1 do begin
	oplot, xx[i]*[1,1], [yy0[i], yy1[i]], color=co[1]
endfor
plotsym, 3, 1.1, /fill
oplot, xx, yy, psym=8, color=co[1]
plotsym, 3, 0.8, /fill
oplot, xx, yy, psym=8, color=cgcolor('background')

;#############################################################
; Chandra data
; 12437, -11.2747 -11.3451     -11.2094
; 17547, -11.4559 -11.5425     -11.3751 
; DATE-OBS= '2011-02-01T14:32:15' / Observation start date
; DATE-END= '2011-02-01T16:34:31' / Observation end date
; DATE-OBS= '2014-11-12T05:34:48' / Observation start date
; DATE-END= '2014-11-12T07:30:13' / Observation end date

cha_date = dblarr(2)
date0 = '2011-02-01T14:32:15'
date1 = '2011-02-01T16:34:31'
cha_date[0] = (date_conv(date0, 'MODIFIED') + date_conv(date1, 'MODIFIED')) / 2d
date0 = '2014-11-12T05:34:48'
date1 = '2014-11-12T07:30:13'
cha_date[1] = (date_conv(date0, 'MODIFIED') + date_conv(date1, 'MODIFIED')) / 2d

cha_lum = 10d ^ [-11.2747, -11.4559] * flconv
cha_lum0 = 10d ^ [-11.3451, -11.5425] * flconv
cha_lum1 = 10d ^ [-11.2094, -11.3751] * flconv

xx = cha_date - MJD0
yy = cha_lum
yy0 = cha_lum0
yy1 = cha_lum1

for i = 0, n_elements(xx) - 1 do begin
	oplot, xx[i]*[1,1], [yy0[i], yy1[i]], color=co[2]
endfor
plotsym, 3, 1.1, /fill
oplot, xx, yy, psym=8, color=co[2]
plotsym, 3, 0.8, /fill
oplot, xx, yy, psym=8, color=cgcolor('background')

plot, xr1, ylog=ylog, /nodata, charsize=cs, font=fnt, position=pos1, /noerase, $
	xstyle=1, xrange=xr1, ystyle=1, yrange=yr, $
	xminor=2, xtickinterval=200, $
	xtitle=strjoin(replicate(' ', 45)) + 'MJD - ' + string(MJD0, format='(I0)'), $
	ytitle='0.3-2 keV Luminosity (10!U39!N erg s!U-1!N)'

xx = xmm_date[2] - MJD0
yy = xmm_lum[2]
yy0 = xmm_lum0[2]
yy1 = xmm_lum1[2]

for i = 0, n_elements(xx) - 1 do begin
	oplot, xx[i]*[1,1], [yy0[i], yy1[i]], color=co[1]
endfor
plotsym, 3, 1.1, /fill
plots, xx, yy, psym=8, color=co[1]
plotsym, 3, 0.8, /fill
plots, xx, yy, psym=8, color=cgcolor('background')

xx = cha_date[0] - MJD0
yy = cha_lum[0]
yy0 = cha_lum0[0]
yy1 = cha_lum1[0]

for i = 0, n_elements(xx) - 1 do begin
	oplot, xx[i]*[1,1], [yy0[i], yy1[i]], color=co[2]
endfor
plotsym, 3, 1.1, /fill
plots, xx, yy, psym=8, color=co[2]
plotsym, 3, 0.8, /fill
plots, xx, yy, psym=8, color=cgcolor('background')

; ######

lx = 150
ly = 14
dy = 1
tx = 50
ty = 0.3
plotsym, 0, 1.1, /fill
plots, lx, ly, psym=8, color=co[0]
plotsym, 0, 0.8, /fill
plots, lx, ly, psym=8, color=cgcolor('background')
xyouts, lx + tx, ly - ty, 'Swift', charsize=cs, font=fnt
plotsym, 3
plots, lx, ly - dy, psym=8, color=co[1]
plotsym, 3, 0.8, /fill
plots, lx, ly - dy, psym=8, color=cgcolor('background')
xyouts, lx + tx, ly - dy - ty, 'XMM-Newton', charsize=cs, font=fnt

plotsym, 4
plots, lx, ly - dy*2, psym=8, color=co[2]
plotsym, 4, 0.8, /fill
plots, lx, ly - dy*2, psym=8, color=cgcolor('background')
xyouts, lx + tx, ly - dy*2 - ty, 'Chandra', charsize=cs, font=fnt


if keyword_set(saveps) then begin
	device, /close
	chepslinewidth, epsname, 40
	set_plot, current_device
endif

stop

end
