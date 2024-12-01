pro pbb, saveps=saveps

distance = 3.4d2 ; distance in unit of 10kpc
files = file_search('./*/fit.txt')

n_obs = n_elements(files)
tbb = dblarr(n_obs, 3)
rbb = dblarr(n_obs, 3)

cstat = dblarr(n_obs)
dof = dblarr(n_obs)

cnet = dblarr(n_obs)
rate = dblarr(n_obs)
rate_err = dblarr(n_obs)
date = dblarr(n_obs)
exposure = dblarr(n_obs)

ebound = [30, 200] ; energy bounds

for i = 0L, n_obs - 1L do begin
	fname = files[i]
	readcol, fname, format='D,D,D', x, e0, e1, /silent
	tbb[i, *] = [x[0], e0[0], e1[0]]
	rbb[i, *] = sqrt([x[1], e0[1], e1[1]]) * distance
	cstat[i] = x[2]
	dof[i] = e0[2]
	
	; directory name
	dir = file_dirname(fname)
	fsrc = dir + '/ulx.pi'
	fbkg = dir + '/ulx_bkg.pi'
	
	dsrc = mrdfits(fsrc, 1, shdr, /silent)
	dbkg = mrdfits(fbkg, 1, bhdr, /silent)
	exptime = sxpar(shdr, 'EXPOSURE')
	
	scales = sxpar(shdr, 'BACKSCAL')
	scaleb = sxpar(bhdr, 'BACKSCAL')
	f = scales / scaleb
	
	; get MJD at the mid-point of the observation
	date0 = sxpar(shdr, 'DATE-OBS')
	date1 = sxpar(shdr, 'DATE-END')
	date[i] = (date_conv(date0, 'MODIFIED') + date_conv(date1, 'MODIFIED')) / 2d
	
	ss = dsrc.counts
	bb = dbkg.counts
	
	csrc = total(ss[ebound[0]:ebound[1] - 1L])
	cbkg = total(bb[ebound[0]:ebound[1] - 1L])
	
	cnet[i] = (csrc - cbkg * f)
	rate[i] = (csrc - cbkg * f) / exptime
	var = csrc + cbkg * f^2d
	rate_err[i] = sqrt(var) / exptime

	exposure[i] = exptime
	
;	print, dir, tbb[i, *], rbb[i, *], $
;		format='(A12, F8.4, "  (", F6.4, " - ", F6.4, ")", E10.2, "  (", E8.2, " - ", E8.2, ")")'	
	
	print, dir, cstat[i], dof[i], tbb[i, 0], format='(A12, F8.2, I5, F8.4)'	
endfor

; print, min(exposure)/1d3, max(exposure)/1d3

if keyword_set(saveps) then begin
	current_device = !d.name ; save current device
	set_plot, 'ps'
	epsname = 'ulx_rt.eps'
	device, filename=epsname, /encapsulated
endif

xr = [0.04, 0.4]
yr = [1d3, 2d5]
cs = 2
fnt = 1
co = [cgcolor('blue'), cgcolor('red'), cgcolor('dark green')]

plot, xr, /xlog, /ylog, /nodata, charsize=cs, font=fnt, $
	xstyle=1, xrange=xr, ystyle=1, yrange=yr, xtickname=' ', $
	xmargin=[5.5, 1], ymargin=[3.5, 1], $
	xtitle='blackbody temperature (keV)', ytitle='blackbody radius (km)'

xtx = [0.05, 0.1, 0.2, 0.3]
xts = ['0.05', '0.1', '0.2', '0.3']
xyouts, xtx, yr[0]/1.7, xts, align=0.5, font=fnt, charsize=cs

xx = tbb[*, 0] ; MJD
yy = rbb[*, 0]
x0 = tbb[*, 1]
x1 = tbb[*, 2]
y0 = rbb[*, 1]
y1 = rbb[*, 2]

plotsym, 0, /fill
oplot, xx, yy, psym=8, color=co[0]
; xyouts, xx, yy, string(cstat, dof, format='(F0.2, "/", I0)'), charsize=1

;#############################################################################

; adding XMM data points
; 07 high state, ulx_agbp
;  11    6   bbodyrad   kT         keV      0.126182     +/-  7.76718E-03  
;  12    6   bbodyrad   norm                1809.24      +/-  997.161      
;  11     0.117901     0.134395    (-0.00824042,0.00825347)
;  12      1006.89      3412.11    (-807.82,1597.4)
; 07 low state, ulx_bp
;   4    4   bbodyrad   kT         keV      5.78170E-02  +/-  9.01660E-03  
;   5    4   bbodyrad   norm                1.60118E+07  +/-  5.13248E+07  
;   4    0.0489621    0.0673099    (-0.00874011,0.00960771)
;   5       765552   8.3238e+08    (-1.59601e+07,8.15655e+08)
; 06, ulx_ebp
;   6    5   bbodyrad   kT         keV      9.44056E-02  +/-  7.01849E-03  
;   7    5   bbodyrad   norm                7448.56      +/-  6524.59      
;   6    0.0872573     0.102021    (-0.00715877,0.00760515)
;   7      3129.14      19464.5    (-4309.82,12025.5)

xmm_tbb = [0.126182, 5.78170E-02, 9.44056E-02]
xmm_tbb0 = [0.117901, 0.0489621, 0.0872573]
xmm_tbb1 = [0.134395, 0.0673099, 0.102021]
xmm_rbb = sqrt([1809.24, 1.60118E+07, 7448.56]) * distance
xmm_rbb0 = sqrt([1006.89, 765552, 3129.14]) * distance
xmm_rbb1 = sqrt([3412.11, 8.3238e+08, 19464.5]) * distance

for i = 0, n_elements(xmm_tbb) - 1 do begin
	oplot, xmm_tbb[i]*[1,1], [xmm_rbb0[i], xmm_rbb1[i]], color=co[1]
	oplot, [xmm_tbb0[i], xmm_tbb1[i]], xmm_rbb[i]*[1,1], color=co[1]
endfor
plotsym, 3
oplot, xmm_tbb, xmm_rbb, psym=8, color=co[1]
plotsym, 3, 0.8, /fill
oplot, xmm_tbb, xmm_rbb, psym=8, color=cgcolor('background')

;#############################################################################

; adding Chandra data
; 12437
;   3    3   bbodyrad   kT         keV      8.85142E-02  +/-  4.31815E-03  
;   4    3   bbodyrad   norm                1.53530E+04  +/-  6213.52      
;     3    0.0844463    0.0930164    (-0.00406794,0.00450221)
;     4      10175.1      22722.8    (-5177.93,7369.79)
; 17547
;   3    3   bbodyrad   kT         keV      0.105386     +/-  5.19391E-03  
;   4    3   bbodyrad   norm                4115.44      +/-  1711.47      
;     3     0.100083     0.111242    (-0.00530276,0.00585563)
;     4      2613.81      6353.58    (-1501.63,2238.14)


cha_tbb = [8.85142E-02, 0.105386]
cha_tbb0 = [0.0844463, 0.100083]
cha_tbb1 = [0.0930164, 0.111242]
cha_rbb = sqrt([1.53530E+04, 4115.44]) * distance
cha_rbb0 = sqrt([10175.1, 2613.81]) * distance
cha_rbb1 = sqrt([22722.8, 6353.58]) * distance

for i = 0, n_elements(cha_tbb) - 1 do begin
	oplot, cha_tbb[i]*[1,1], [cha_rbb0[i], cha_rbb1[i]], color=co[2]
	oplot, [cha_tbb0[i], cha_tbb1[i]], cha_rbb[i]*[1,1], color=co[2]
endfor
plotsym, 4
oplot, cha_tbb, cha_rbb, psym=8, color=co[2]
plotsym, 4, 0.8, /fill
oplot, cha_tbb, cha_rbb, psym=8, color=cgcolor('background')


fx = alog10([xx, xmm_tbb, cha_tbb])
fy = alog10([yy, xmm_rbb, cha_rbb])
a = linfit(fx, fy, sigma=aerr)
oplot, xr, 10d^(a[0] + alog10(xr) * a[1]), linestyle=2

print, 'linear fit:', a, aerr

lx = 0.045
ly = 4d3
dy = 1.5
tx = 1.1
ty = 1.1
plotsym, 0, /fill
plots, lx, ly, psym=8, color=co[0]
xyouts, lx * tx, ly/ty, 'Swift', charsize=cs, font=fnt
plotsym, 3
plots, lx, ly/dy, psym=8, color=co[1]
plotsym, 3, 0.8, /fill
plots, lx, ly/dy, psym=8, color=cgcolor('background')
xyouts, lx * tx, ly/dy/ty, 'XMM-Newton', charsize=cs, font=fnt

plotsym, 4
plots, lx, ly/dy^2, psym=8, color=co[2]
plotsym, 4, 0.8, /fill
plots, lx, ly/dy^2, psym=8, color=cgcolor('background')
xyouts, lx * tx, ly/dy^2/ty, 'Chandra', charsize=cs, font=fnt


if keyword_set(saveps) then begin
	device, /close
	chepslinewidth, epsname, 40
	set_plot, current_device
endif

stop

end
