pro intensity_over_burst

	; Procedure to illustrate what is meant by intensity change over time.

	;---------------------------------------------;
	;			Define variables
	!p.thick = 4
	!p.charsize = 1.5
	loadct, 39
	pos = [0.13, 0.1, 0.95, 0.95]
	col_scale = 1.3
	dimen = 600
	xpos = 500
	ypos = 50 
	folder = '~/Data/2011_sep_22/herringbones'
	cd, folder

	;---------------------------------------------;
	;			Read spectrogram
	loadct, 5
	!p.charsize = 1
	cd, '~/Data/2011_sep_22/herringbones'
	radio_spectro_fits_read, 'BIR_20110922_104459_01.fit', data_raw, times, freq

	t1_index = closest(times, anytim(file2time('20110922_105000'),/utim))
	t2_index = closest(times, anytim(file2time('20110922_105300'),/utim))
	f1_index = closest(freq, 90.0)
	f2_index = closest(freq, 41.0)
	inten0 = -20
	inten1 = 40

	;---------  Chosse intensity clipping and Hough angles  --------;
	data_bs = constbacksub(data_raw, /auto)

	loadct, 0
	window, 1, xs=700, ys=700
	spectro_plot,( bytscl(data_bs, inten0, inten1) ), times, freq, $
		/ys, $
		ytitle = '!6Frequency [MHz]', $
		yticks = 5, $
		yminor = 4, $
		yr = [freq[f1_index], freq[f2_index]], $
		xrange = [times[t1_index], times[t2_index]], $
		/xs, $
		xtitle = 'Start time: '+anytim(times[t1_index], /cc, /trun)+' (UT)', $
		charsize = 2.0
    
    ;---------------------------------------------;
	;			Read burst data
	readcol, 'bursts_ft_second_master_reverse.txt', btall0, bfall0, biall0, format = 'A,D,D'
	plot_bursts, data_bs, times, freq, btall0, bfall0, biall0, color = 170 ; Green



END

pro plot_bursts, data, time, freq, btall, bfall, biall, color = color
  
  indices = where(btall eq '-')
  indices = [-1, indices]
  
  plotsym, 0, /fill
  
	;FOR i=20, n_elements(indices)-2 DO BEGIN 
	    
	    ;i = n_elements(indices)-5
	    i=30
	    bt = anytim(btall[indices[i]+1: indices[i+1]-1], /utim)  ; This selects one burst.
	    bf = bfall[indices[i]+1: indices[i+1]-1]
	    bi = biall[indices[i]+1: indices[i+1]-1]	

	    set_line_color
	    oplot, bt, bf, psym = 8, symsize=0.4, color=0
	    ;oplot, bt, bf, psym = 8, symsize=0.3, color=0
	    loadct, 0, /silent
	    oplot, bt, bf, psym = 8, symsize=0.3, color = color  

	    window, 0, ys = 1000, xs = 400

	    burst_surface = dblarr(n_elements(bf), 100.0)
	    ts_array = dblarr(n_elements(bf), 100.0)
	    fs_array = fltarr(n_elements(bf))
	    FOR j = 0, n_elements(bt)-1 DO BEGIN
	        tpoint = closest(time, bt[0])
	        fpoint = closest(freq, bf[0])

			tbef = bt[j] - 0.8
			taft = bt[j] + 0.8
			ts = (dindgen(100)*(taft - tbef)/99.0) + tbef
			fs = findgen(n_elements(ts))
			fs[*] = closest(freq, bf[j])

			tinds = interpol(findgen(3600), time, ts)
			; Run the interpolation algorithm...
			profile_interp=interpolate(data, tinds, fs, cubic=-0.5)
			if j eq 0 then begin
				window, 0, ys = 1000, xs = 400
				 utplot, ts, profile_interp - j, $
				 		thick=1, $
				 		yr = [-320, 10], $
				 		/ys, $
				 		xr = [time[tinds[0]], time[tinds[n_elements(tinds)-1]+20]], $
				 		title = 'Burst: '+string(i)
				 burst_surface[j, *] = profile_interp[*]
				 ts_array[j, *] = ts
				 fs_array[j] = bf[j]
			endif else begin
				wset, 0
				outplot, ts, profile_interp - j*5.0, thick=1
				print, n_elements(ts)
				burst_surface[j, *] = profile_interp[*]
				ts_array[j, *] = ts
				fs_array[j] = bf[j]
			;print, bf[j]
			endelse	 

		ENDFOR
		print, i

	;ENDFOR

	save, burst_surface, ts_array, fs_array, filename = 'single_burst_suface.sav'

  
END  

