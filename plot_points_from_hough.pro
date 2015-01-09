pro plot_points_from_hough
;
;
;11-Sep-2013 - Code to plot the points detected by the Hough transform
;
;
cd,'~/Data/CALLISTO/20110922'
loadct,1
window,1, xs=1300, ys=600

radio_spectro_fits_read,'BIR_20110922_104459_01.fit', data_raw, times, freq
t1_index = closest(times,anytim(file2time('20110922_105120'),/utim))
t2_index = closest(times,anytim(file2time('20110922_105250'),/utim))
f1_index = closest(freq,80.0)
f2_index = closest(freq,40.0)
data_bs = constbacksub(data_raw,/auto)

spectro_plot, data_bs , times, freq, $
/ys, ytitle='!6Frequency [MHz]', yticks=5, yminor=4, yr = [freq[f1_index],freq[f2_index]], $
xrange=[times[t1_index],times[t2_index]], /xs, xtitle='Start time:'+anytim(times[t1_index],/yoh),$
charsize=1.5
	
set_line_color
restore,'peak_time_freq.sav'
FOR i=0, n_elements(peak_time_freq[*,0])-1 DO BEGIN
    plots,peak_time_freq[i,1:99], peak_time_freq[i,0.0], color=3, psym=1, symsize=1
ENDFOR

;freq = peak_time_freq[*,0]




btimes = 0.0
bf=0.0
drift=0.0
window,5
window,4
window,3
window,2
plot, [0,0], [0,0], xr=[40,80], yr=[140,200]

ft1_index = where(peak_time_freq[80, *] gt 0.0)
ft1 = peak_time_freq[80, ft1_index]


FOR j=1, n_elements(ft1)-1 DO BEGIN
comp_f = ft1[0]
comp_t = ft1[j]

loadct,1
wset,1
spectro_plot, data_bs , times, freq, $
/ys, ytitle='!6Frequency [MHz]', yticks=5, yminor=4, yr = [freq[f1_index],freq[f2_index]], $
xrange=[times[t1_index],times[t2_index]], /xs, xtitle='Start time:'+anytim(times[t1_index],/yoh),$
charsize=1.5


	i=n_elements(peak_time_freq[*, 0])-1
	WHILE i gt 1 DO BEGIN 
	
		ft_index = where(peak_time_freq[i, *] gt 0.0)
		ft = peak_time_freq[i, ft_index]

		f = ft[0]
		index_t = closest(ft, comp_t)
		t = ft[index_t]
		
		
		If abs(comp_t - t) gt 0.5 THEN BEGIN
			i=0
		ENDIF ELSE BEGIN
			set_line_color
			btimes = [btimes, comp_t, t]
			bf = [bf, comp_f, f]
			plots, comp_f, comp_t, color=4, symsize=2, psym=1
			plots, t, f, color=4, symsize=2, psym=1
		ENDELSE
		
		comp_f = f
		comp_t = t
		i=i-1
	
	ENDWHILE
	btimes = btimes[1: n_elements(btimes)-1]
	bf = bf[1: n_elements(bf)-1]
	tindex = btimes
	findex = bf
	
	;---------------------------------------------;
	;				Get profile
	;
	
	FOR i=0, n_elements(btimes)-1 DO BEGIN
		tindex[i] = closest(times, btimes[i])
		findex[i] = closest(freq, bf[i])
	ENDFOR
	prof= interpolate(data_bs, tindex, findex) 
	
	bt = btimes
	bff = bf
	inten = prof
	
	
	write_text, bt, bff, inten
	
	loadct,39
	wset,2
	plot, bf, prof, psym=1,color=(240 - j*3.5), /noerase, xr=[40,80], yr=[140,200],$
	 xtitle='Frequency (MHz)', ytitle='Intensity' 
	oplot, bf, prof, color=(240 - j*3.5)
	
	wset,3
	tsec= btimes[*]-btimes[0]
	plot, tsec, bf,  color=(240 - j*3.5), /noerase, xr=[0,5], yr=[40,80], $
	xtitle='Time (s)', ytitle='Frequency (MHz)'
	
	wset,4
	plot, tsec, prof,  color=(240 - j*3.5), xr=[0,5], yr=[140, 200], /noerase
	
	;---------------------------------------------;
	;				Get drift rate
	;		
	wset,5
	index = where(tsec gt 0.0)
	if index[0] gt -1 then begin
		result = linfit(tsec[index], bf[index])
	ENDIF ELSE BEGIN
		result = linfit(tsec, bf)
	ENDELSE	
	drift = result[1]
	print,'Drift rate: '+string(drift)
	
	
	if index[0] gt -1 then begin
		result = linfit(tsec[index], prof[index])
	ENDIF ELSE BEGIN
		result = linfit(tsec, prof)
	ENDELSE		
	it = result[1]
	maxI = max(prof)
	print, 'Intesnity v time: '+string(it)
	set_line_color
	plot, [0,drift], [0, maxI]  , /noerase, psym=2, xr=[6,14], yr=[160, 180], color=1
	
	
	btimes = 0.0
	bf=0.0
	drift=0.0

ENDFOR
wset,2
set_line_color
axis, xaxis=0, xr=[40,80], color=1, xtitle='Frequency (MHz)'
axis, xaxis=1, xr=[40,80], color=1, xtickname=[' ', ' ', ' ', ' ', ' ']

axis, yaxis=0, xr=[0, 50], color=1, ytitle='Intensity' 
axis, yaxis=1, xr=[0, 50], color=1, ytickname=[' ', ' ', ' ', ' ', ' ',' ']

wset,3
set_line_color
axis, xaxis=0, xr=[0,5], color=1, xtitle='Time (s)'
axis, xaxis=1, xr=[0,5], color=1, xtickname=[' ', ' ', ' ', ' ', ' ']

axis, yaxis=0, yr=[40, 80], color=1, ytitle='Frequency (MHz)'
axis, yaxis=1, yr=[40, 80], color=1, ytickname=[' ', ' ', ' ', ' ', ' ',' ']
;plots, t1, f1, color=4, psym=1, symsize=2
	



END


pro write_text, bt, bff, inten

IF file_test('bursts_bs_hough.txt') eq 1 THEN BEGIN
	readcol,'bursts_bs_hough.txt', btprev, bffprev, intenprev,$
	format = 'A,D,D'
	
	bt = [btprev, '-', anytim(bt, /ccs)]
    bff = [bffprev, !Values.F_NAN, bff]
	inten = [intenprev, !Values.F_NAN, inten]

ENDIF


writecol, 'bursts_bs_hough.txt', anytim(bt, /ccs), bff, inten, fmt='(A,D,D)'


END
