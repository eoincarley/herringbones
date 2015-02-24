pro burst_intensity_check

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
	folder = '~/Data/22Sep2011_event/herringbones'
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

	indices = where(btall eq '-')
	indices = [-1, indices]

	restore, 'beam_kins.sav', /verb

	drift = beam_kins[0]
	vels = beam_kins[1]
	displ = beam_kins[2]
	flength =  beam_kins[4]
	figs_folder = beam_kins[5]
	max_bi = fltarr(n_elements(vels))


	;setup_ps, 'figures/'+figs_folder+'/inten_time_fits.eps'
	j=0
	window, 0
	window, 1
	FOR i=0, n_elements(indices)-2 DO BEGIN

        bt = btall[indices[i]+1: indices[i+1]-1]
        bf = bfall[indices[i]+1: indices[i+1]-1]
        bi = biall[indices[i]+1: indices[i+1]-1]
        tsec = anytim(bt[*], /utim) - anytim(bt[0], /utim)
    
        ;Plot burst on spectra
        wset, 0
		plot_burst_on_spectra, data_bs, times, freq, f1_index, f2_index, t1_index, t2_index, bt, bf


		wset, 1
        plot, tsec, bi, $
            psym=1, $
            xr=[0, 6], $
            yr=[0, 45], $
            /xs, $
            /ys, $
            position=pos, $
            /normal, $
            xtitle = 'Frequency (MHz)', $
            ytitle = 'Intensity (DNs)', $
            title = 'Burst number: '+string(j, format='(I02)') 


        ;---------------------------------;
        ;		Check fit
        ;   
        result = linfit(tsec, bi, yfit = yfit)
      
        oplot, tsec, yfit, color=100
     
        stop
    ENDFOR	
	;device, /close

END

pro plot_burst_on_spectra, data_bs, times, freq, f1_index, f2_index, t1_index, t2_index, bt, bf


	t1_index = closest(times, anytim(bt[0],/utim) - 20.0)
	t2_index = closest(times, anytim(bt[n_elements(bt)-1], /utim) + 20.0)
		
	  loadct, 74
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
