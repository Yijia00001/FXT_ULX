pro count
; sort and group the observations according to their count rate

files = file_search('./*33469*/ulx.pi')

n_obs = n_elements(files)

count1 = dblarr(n_obs) ; 0.3-1 keV
count2 = dblarr(n_obs) ; 0.3-2 keV
rate1 = dblarr(n_obs)
rate2 = dblarr(n_obs)
dir = strarr(n_obs)

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

ebound = [30, 100, 200] ; energy bounds

for i = 0L, n_obs - 1L do begin
	fsrc = files[i]
	dir[i] = file_dirname(fsrc)
	fbkg = dir[i] + '/ulx_bkg.pi'
	
	dsrc = mrdfits(fsrc, 1, shdr, /silent)
	dbkg = mrdfits(fbkg, 1, bhdr, /silent)
	exptime = sxpar(shdr, 'EXPOSURE')
	
	scales = sxpar(shdr, 'BACKSCAL')
	scaleb = sxpar(bhdr, 'BACKSCAL')
	f = scales / scaleb
	
	ss = dsrc.counts
	bb = dbkg.counts
	
	c1 = total(ss[ebound[0]:ebound[1] - 1L])
	c2 = total(ss[ebound[0]:ebound[2] - 1L])
	b1 = total(bb[ebound[0]:ebound[1] - 1L])
	b2 = total(bb[ebound[0]:ebound[2] - 1L])
	
	count1[i] = c1
	count2[i] = c2
	rate1[i] = (c1 - b1 * f) / exptime
	rate2[i] = (c2 - b2 * f) / exptime
endfor

s = reverse(sort(rate2)) ; descending 
total_count = total(count2)
print, 'Total number of counts:', total_count

c_thres = 120d

ctot = 0d
k = 0L
num = intarr(n_obs)

for i = 0L, n_obs - 1L do begin
	j = s[i]
	ctot = ctot + count2[j]
	num[k]++
	print, dir[j], rate2[j], count2[j]
	if ctot ge c_thres then begin
		print, ctot, num[k], '---------------------'
		ctot = 0d
		k++
	endif
	if i eq n_obs - 1 then print, ctot, num[k], '---------------------'
endfor

num = num[0:k]
num = [3, 4, 5, 7, 9]

print, '================================================'

k = 0L
for i = 0, n_elements(num) - 1 do begin
	print, string(num[i], format='(I0)'), ' ---------------', total(count2[s[k:k+num[i]-1]])
	for j = 0, num[i] - 1 do begin
		print, dir[s[k]], count2[s[k]]
		k++
	endfor
endfor

stop

end
