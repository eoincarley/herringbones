function calc_vel, tsec, ftfit, model = model

  ftfit_Hz = 1.0*ftfit*1e6                 ; Hz
  dens = freq2dens(ftfit_Hz)
  rads = density_to_radius(dens, model = model, fold = 1)
  result = linfit(tsec, rads)
  velocity = abs(result[1])*6.955e8  ; m/s
  velocity = [velocity/2.997e8]   ; speed of light units
  kins = list(velocity, rads)
  
  return, kins
  
END

pro max_intensity_check

	; Check that intensity fits are ok.

	!p.thick = 2 
	!p.charsize = 1.5
	!p.charthick = 1.5
	loadct, 39
	pos = [0.13, 0.1, 0.95, 0.95]
	col_scale = 1.5
	dimen = 600
	xpos = 500
	ypos = 50
	folder = '~/Data/2011_sep_22/herringbones'
	cd, folder



	;--------------------------------------;
	;		Read and plot Spectrogram
	;
	radio_spectro_fits_read,'BIR_20110922_104459_01.fit', data_raw, times, freq
	t1_index = closest(times, anytim(file2time('20110922_104700'),/utim))
	t2_index = closest(times, anytim(file2time('20110922_105340'),/utim))

	f1_index = closest(freq, 90.0)
	f2_index = closest(freq, 10.0)
	data_bs = constbacksub(data_raw, /auto)
	xtit = 'Start time: '+anytim(times[t1_index], /cc, /trun)+ ' (UT)'

	;---------------------------------;
	;			 Read data
	;
	readcol, 'bursts_ft_first_master_reverse.txt', btall0, bfall0, biall0, format = 'A,D,D'
	readcol, 'bursts_ft_second_master_reverse.txt', btall1, bfall1, biall1, format = 'A,D,D'
	readcol, 'bursts_ft_second_master_forward.txt', btall2, bfall2, biall2, format = 'A,D,D'
	readcol, 'bursts_ft_first_master_forward.txt', btall3, bfall3, biall3, format = 'A,D,D'

	btall = [btall0, '-', btall1] ;, '-', btall2, '-', btall3]
	bfall = [bfall0, !Values.F_NAN, bfall1] ;, !Values.F_NAN, bfall2, !Values.F_NAN, bfall3]
	biall = [biall0, !Values.F_NAN, biall1]; , !Values.F_NAN, biall2, !Values.F_NAN, biall3]
	rev1f0 = bfall0[0]
	rev2f0 = bfall1[0]
	for2f0 = bfall2[0]
	for1f0 = bfall3[0]

	indices = where(btall eq '-')
	indices = [-1, indices]
	displ = fltarr(n_elements(indices))
	max_bi = fltarr(n_elements(indices))
	
	;setup_ps, 'figures/'+figs_folder+'/inten_time_fits.eps'
	j=0
	window, 0
	window, 1, xs=500, ys=500
	
	FOR i=0, n_elements(indices)-2 DO BEGIN

        bt = btall[indices[i]+1: indices[i+1]-1]
        bf = bfall[indices[i]+1: indices[i+1]-1]
        bi = biall[indices[i]+1: indices[i+1]-1]
        tsec = anytim(bt[*], /utim) - anytim(bt[0], /utim)
    
    	; Do fitting
        result = linfit(tsec, bf, yfit = yfit)
        start = [result[1]]
        
        fit = 'p[0]*x + ' +	string(bf[0], format='(f5.2)')
        result = mpfitexpr(fit, tsec , bf, err, yfit=ftfit, start)
        
        kins = calc_vel(tsec, ftfit, model =  'tcd')
        radii = kins[1]
        displ[i] = abs(radii[0] - radii[n_elements(radii)-1])
        
        ;Plot burst on spectra
        wset, 0
		plot_burst_on_spectra, data_bs, times, freq, f1_index, f2_index, t1_index, t2_index, bt, bf

    	max_bi[i] = max(bi)

    	wset, 1
    	if max_bi[i] gt 34 and max_bi[i] lt 46 and displ[i] gt 0.12 and displ[i] lt 0.16 then begin

	    	plot, [displ[i]], [max_bi[i]], $
		        psym=8, $
		        xr = [0.0, 0.3], $
		        yr = [0.0, 50], $, 
		        /ys, $
		        symsize = 1.1, $
		        xtit = 'Beam displacement (Rsun)', $
		        ytit = 'Maximum intensity (DN)', $
		        color = 200, $
		        /noerase
		    stop    
		endif       

			plot, [displ[i]], [max_bi[i]], $
		        psym=8, $
		        xr = [0.0, 0.3], $
		        yr = [0.0, 50], $, 
		        /ys, $
		        symsize = 1.1, $
		        xtit = 'Beam displacement (Rsun)', $
		        ytit = 'Maximum intensity (DN)', $
		        /noerase
	ENDFOR	
END

pro plot_burst_on_spectra, data_bs, times, freq, f1_index, f2_index, t1_index, t2_index, bt, bf


	t1_index = closest(times, anytim(bt[0],/utim) - 20.0)
	t2_index = closest(times, anytim(bt[n_elements(bt)-1], /utim) + 20.0)
		
	  loadct, 74, /silent
      spectro_plot, bytscl(data_bs, -10, 40) , times, freq, $
          /ys, $
          ytitle = 'Frequency (MHz)', $
          yr = [freq[f1_index], freq[f2_index]], $
          xrange = [times[t1_index],times[t2_index]], $
          /xs, $
          title = 'eCallisto, Birr Castle, Ireland.', $
          xcolor = '100', $
          xticklen = -0.01, $
          yticklen = -0.01

      plotsym, 0, /fill    
      set_line_color
      oplot, anytim(bt, /utim), bf, psym = 8, symsize=0.5, color=4    



END
