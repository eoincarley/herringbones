pro plot_burst_stats_raw

;12-Sep-2013 plot burst stats derived from Hough analysis


cd,'~/Data/CALLISTO/20110922'
;--------------------------------------------
;			 First do raw data
;
readcol,'bursts_raw_hough.txt', btall, bfall, biall, format = 'A,D,D'
indices = where(btall eq '-')
indices = [-1, indices]

loadct,39
pos = [0.13, 0.1, 0.95, 0.95]
col_scale=2.5
;-------------------------------------;
;
;		 Intensity v Frequency
;

set_plot,'ps'
!p.font=0
!p.charsize=1.2
device, filename='raw_bfbi.eps', /color, /helvetica, /inches, xsize=7, ysize=7, /encapsulate, yoffset=5

FOR i=0, n_elements(indices)-2 DO BEGIN

	bt = btall[indices[i]+1: indices[i+1]-1]
	bf = bfall[indices[i]+1: indices[i+1]-1]
	bi = biall[indices[i]+1: indices[i+1]-1]
	
	IF i eq 0 THEN BEGIN
		;wset,2
		plot, bf, bi, psym=1, color=0, xr=[45,80], yr = [130, 200], $
		xtitle='Frequency (MHz)', ytitle='Intensity (DN)', position=pos,$
		/normal, title='Raw Herringbones 22-Sep-2011', thick=2, /xs, /ys
		
		oplot, bf, bi, color=(250), psym=1, thick=2
		oplot, bf, bi, color=(250), thick=2
	ENDIF ELSE BEGIN	
		oplot, bf, bi, color = (245 - i*col_scale), psym=1, thick=2
		oplot, bf, bi, color=(245 - i*col_scale), thick=2
	ENDELSE	
ENDFOR	
axis, xaxis=0, xtickname=[' ', ' ', ' ', ' ', ' ', ' ', ' ', ' '], xr=[45,80], xthick=3, /xs
axis, xaxis=1, xtickname=[' ', ' ', ' ', ' ', ' ', ' ', ' ', ' '], xr=[45,80], xthick=3, /xs
axis, yaxis=0, ytickname=[' ', ' ', ' ', ' ', ' ', ' ', ' ', ' '], yr=[130, 200], ythick=3, /ys
axis, yaxis=1, ytickname=[' ', ' ', ' ', ' ', ' ', ' ', ' ', ' '], yr=[130, 200], ythick=3, /ys
device,/close
set_plot,'x'

;-------------------------------------;
;
;		 Frequency v Time
;
set_plot,'ps'
device, filename='raw_btbf.eps', /color, /helvetica, /inches, xsize=7, ysize=7, /encapsulate, yoffset=5	
FOR i=0, n_elements(indices)-2 DO BEGIN

	bt = btall[indices[i]+1: indices[i+1]-1]
	bf = bfall[indices[i]+1: indices[i+1]-1]
	bi = biall[indices[i]+1: indices[i+1]-1]
	
	tsec = anytim(bt[*], /utim) - anytim(bt[0], /utim)
	IF i eq 0 THEN BEGIN
		;wset,3
		plot, tsec, bf, color=0, xr=[0,4], yr=[45,80],position=pos,$
		/normal, title='Raw Herringbones 22-Sep-2011', $
		xtitle='Time (s)', ytitle='Frequency (MHz)', thick=3, /xs, /ys
		
		oplot, tsec, bf, color=(250), thick=3
	ENDIF ELSE BEGIN	
		oplot, tsec, bf, color=(245 - i*col_scale), thick=3
	ENDELSE	

ENDFOR
axis, xaxis=0, xtickname=[' ', ' ', ' ', ' ', ' ', ' '], xr=[0,4], xthick=3, /xs
axis, xaxis=1, xtickname=[' ', ' ', ' ', ' ', ' ', ' '], xr=[0,4], xthick=3, /xs
axis, yaxis=0, ytickname=[' ', ' ', ' ', ' ', ' ', ' ', ' ', ' '], yr=[45,80], ythick=3, /ys
axis, yaxis=1, ytickname=[' ', ' ', ' ', ' ', ' ', ' ', ' ', ' '], yr=[45,80], ythick=3, /ys

device,/close
set_plot,'x'


;-------------------------------------;
;
;		   Time v Intensity
;
set_plot,'ps'
device, filename='raw_btbi.eps', /color, /helvetica, /inches, xsize=7, ysize=7, /encapsulate, yoffset=5	
FOR i=0, n_elements(indices)-2 DO BEGIN

	bt = btall[indices[i]+1: indices[i+1]-1]
	bf = bfall[indices[i]+1: indices[i+1]-1]
	bi = biall[indices[i]+1: indices[i+1]-1]
	tsec = anytim(bt[*], /utim) - anytim(bt[0], /utim)
	IF i eq 0 THEN BEGIN
		;wset,4
		plot, tsec, bi, psym=1, color=0, xr=[0,4], yr=[130, 200], position=pos,$
		/normal, title='Raw Herringbones 22-Sep-2011', $
		xtitle='Time (s)', ytitle='Intensity', thick=2, /ys, /xs
		
		oplot, tsec, bi, color=(250), psym=1, thick=2
		oplot, tsec, bi, color=(250), thick=2
	ENDIF ELSE BEGIN	
		oplot, tsec, bi, color = (245 - i*col_scale), psym=1, thick=2
		oplot, tsec, bi, color=(245 - i*col_scale), thick=2
	ENDELSE		
	
ENDFOR
axis, xaxis=0, xtickname=[' ', ' ', ' ', ' ', ' ', ' ', ' '], xr=[0,4], xthick=3, /xs
axis, xaxis=1, xtickname=[' ', ' ', ' ', ' ', ' ', ' ', ' '], xr=[0,4], xthick=3, /xs
axis, yaxis=0, ytickname=[' ', ' ', ' ', ' ', ' ', ' ', ' ', ' '], yr=[130, 200], ythick=3, /ys
axis, yaxis=1, ytickname=[' ', ' ', ' ', ' ', ' ', ' ', ' ', ' '], yr=[130, 200], ythick=3, /ys
device,/close
set_plot,'x'

;-------------------------------------;
;
;		   Drift v DI/DT
;
alldrift = [0]
allbidrift = [0]
PLOTSYM, 0, /fill               ;Plotting symbol is a circle

set_plot,'ps'
device, filename='raw_dfdt_didt.eps', /color, /helvetica, /inches, xsize=7, ysize=7, /encapsulate, yoffset=5	
FOR i=0, n_elements(indices)-2 DO BEGIN

	bt = btall[indices[i]+1: indices[i+1]-1]
	bf = bfall[indices[i]+1: indices[i+1]-1]
	bi = biall[indices[i]+1: indices[i+1]-1]
	tsec = anytim(bt[*], /utim) - anytim(bt[0], /utim)
	;wset,5
	result = linfit(tsec, bf)
	drift = result[1]
	alldrift = [alldrift, drift]
	
	result = linfit(tsec, bi)
	bidrift = result[1];max(bi)
	allbidrift = [allbidrift, bidrift]
	;print, bidrift
	set_line_color
	IF i eq 0 THEN BEGIN
		plot, [0,drift], [0, bidrift], /noerase, psym=8, color=0, xr=[5,15], yr=[-5, -30], position=pos,$
		/normal, title='Raw Herringbones 22-Sep-2011', $
		xtitle='Drift Rate (MHz s!U-1!N)', ytitle='Rate of change of intensity DI/Dt'
		oplot, [0,drift], [0, bidrift], psym=8, color=5
	ENDIF ELSE BEGIN	
		oplot, 	[0,drift], [0, bidrift], psym=8, color=5
	ENDELSE
	
ENDFOR
axis, xaxis=0, xtickname=[' ', ' ', ' ', ' ', ' ', ' '], xr=[5,15], xthick=3
axis, xaxis=1, xtickname=[' ', ' ', ' ', ' ', ' ', ' '], xr=[5,15], xthick=3
axis, yaxis=0, ytickname=[' ', ' ', ' ', ' ', ' ', ' ', ' '], yr=[-5, -30], ythick=3
axis, yaxis=1, ytickname=[' ', ' ', ' ', ' ', ' ', ' ',' '], yr=[-5, -30], ythick=3
device,/close
set_plot,'x'


alldrift = alldrift[ 1: n_elements(alldrift)-1 ]
remove_nans, alldrift, alldrift
allbidrift = allbidrift[ 1:n_elements(allbidrift)-1 ]
remove_nans, allbidrift, allbidrift

result = correlate(alldrift, allbidrift)
print,'Correlation coefficient: '+string(result)



;---------------------------------------------------;
;
;	 PLOT blank spectrogram with burst detections
;
t1 = anytim(file2time('20110922_105100'),/utim)
t2 = anytim(file2time('20110922_105300'),/utim)

radio_spectro_fits_read,'BIR_20110922_104459_01.fit', data_raw, times, freq
t1_index = closest(times, t1)
t2_index = closest(times, t2)
f1_index = closest(freq,80.0)
f2_index = closest(freq,40.0)
data_bs = constbacksub(data_raw,/auto)
set_plot,'ps'
!p.font=0
!p.charsize=1.5
device, filename='detected_bursts.eps', /color, /helvetica, /inches, xsize=11, ysize=11, /encapsulate, yoffset=7,$
bits_per_pixel=16

loadct,1
spectro_plot, data_bs , times, freq, $
/ys, ytitle='Frequency (MHz)', yr=[80, 40], xr=[t1, t2], /xs, xtitle='Start time:'+anytim(times[t1_index],/yoh),$
position=[0.1, 0.6, 0.9, 0.96], /normal, /noerase

loadct,0

data_raw[*,*] = 220.0
;spectro_plot, bytscl(data_raw, 0, 255), times, freq, yr=[80, 45], xr=[t1, t2], $

tvim, data_raw, position=[0.1, 0.15, 0.9, 0.5], /noaxis
utplot, times, freq, yr=[80, 40], xr=[t1, t2], $
position=[0.1, 0.15, 0.9, 0.5], /normal,thick=3, /noerase, /xs, /ys, ytitle='Frequency (MHz)',$
title='Bursts detected using the Hough trasnform'
col_scale=3.0
stop
loadct,39
FOR i=0, n_elements(indices)-2 DO BEGIN

	bt = btall[indices[i]+1: indices[i+1]-1]
	bf = bfall[indices[i]+1: indices[i+1]-1]
	print, (245 - i*col_scale)
	plots, anytim(bt, /utim), bf, color=(245 - i*col_scale), thick=4
ENDFOR	

device,/close
set_plot,'x'


END