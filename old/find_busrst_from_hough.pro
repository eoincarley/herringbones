pro find_busrst_from_hough
;
;
;11-Sep-2013 - Code to plot the points detected by the Hough transform
;
;
cd,'~/Data/22Sep2011_event/herringbones'
loadct,1
window,1, xs=1300, ys=600

radio_spectro_fits_read,'BIR_20110922_104459_01.fit', data_raw, times, freq
t1_index = closest(times,anytim(file2time('20110922_105100'),/utim))
t2_index = closest(times,anytim(file2time('20110922_105300'),/utim))
f1_index = closest(freq,80.0)
f2_index = closest(freq,40.0)
data_bs = constbacksub(data_raw,/auto)

spectro_plot, data_bs , times, freq, $
/ys, ytitle='!6Frequency [MHz]', yticks=5, yminor=4, yr = [freq[f1_index],freq[f2_index]], $
xrange=[times[t1_index],times[t2_index]], /xs, xtitle='Start time:'+anytim(times[t1_index],/yoh),$
charsize=1.5
	
set_line_color
restore,'peak_time_freq.sav'


;freq = peak_time_freq[*,0]
nfreq = n_elements(peak_time_freq[*,0])

btimes = 0.0
bf=0.0
drift=0.0
;window,5
;window,4;
;window,3
;window,2
plot, [0,0], [0,0], xr=[40,80], yr=[140,200]

ft1_index = where(peak_time_freq[nfreq-1, *] gt 0.0)
ft1 = peak_time_freq[nfreq-1, ft1_index]


FOR j=1, n_elements(ft1)-1 DO BEGIN
comp_f = ft1[0]
comp_t = ft1[j]

loadct,1
wset,1
spectro_plot, data_bs , times, freq, $
/ys, ytitle='!6Frequency [MHz]', yticks=5, yminor=4, yr = [freq[f1_index],freq[f2_index]], $
xrange=[times[t1_index],times[t2_index]], /xs, xtitle='Start time:'+anytim(times[t1_index],/yoh),$
charsize=1.5
FOR i=0, n_elements(peak_time_freq[*,0])-1 DO BEGIN
	set_line_color
    plots,peak_time_freq[i,1:99], peak_time_freq[i,0.0], color=3, psym=1, symsize=1
ENDFOR

	i=n_elements(peak_time_freq[*, 0])-1
	WHILE i gt 1 DO BEGIN 
	
		ft_index = where(peak_time_freq[i, *] gt 0.0)
		ft = peak_time_freq[i, ft_index]

		f = ft[0]
		index_t = closest(ft, comp_t)
		t = ft[index_t]
		
		
		If abs(comp_t - t) gt 0.8 THEN BEGIN
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
		i=i-2
	
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
	prof= interpolate(data_raw, tindex, findex) 
	
	bt = btimes
	bff = bf
	inten = prof
	
	
	write_text, bt, bff, inten
	
	btimes = 0.
	bf = 0.
	prof = 0.
	wait,1.0
ENDFOR
END


pro write_text, bt, bff, inten

IF file_test('bursts_raw_hough.txt') eq 1 THEN BEGIN
	readcol,'bursts_raw_hough.txt', btprev, bffprev, intenprev,$
	format = 'A,D,D'
	
	bt = [btprev, '-', anytim(bt, /ccs)]
    bff = [bffprev, !Values.F_NAN, bff]
	inten = [intenprev, !Values.F_NAN, inten]

ENDIF


writecol, 'bursts_raw_hough.txt', anytim(bt, /ccs), bff, inten, fmt='(A,D,D)'


END
