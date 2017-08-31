function calc_vel, tsec, ftfit, model = model

   ftfit_Hz = 1.0*ftfit*1e6                 ; Hz
   dens = freq2dens(ftfit_Hz)
   rads = density_to_radius(dens, model = model, fold = 1.0)
   result = linfit(tsec, rads)
   velocity = abs(result[1])*6.955e8  ; m/s
   velocity = [velocity/2.997e8]   ; speed of light units
   kins = list(velocity, rads)

   return, kins
  
END

pro setup_ps, name
  
   set_plot,'ps'
   !p.font=0
   !p.charsize=2.5
   device, filename = name, $
          /color, $
          /helvetica, $
          /inches, $
          xsize=7, $
          ysize=7, $
          /encapsulate, $
          yoffset=5

end

pro burst_column_dens, model

   ; Do a linear fit to (f, t) data.
   ; Convert to density, then height.
   ; Turn into velocity.
   ;!p.thick = 4
   !p.charsize = 1.5
   loadct, 39
   pos = [0.13, 0.1, 0.95, 0.95]
   col_scale = 1.3
   dimen = 600
   xpos = 500
   ypos = 50 
   folder = '~/Data/2011_sep_22/herringbones'
   figs_folder = model
   cd, folder

   ;-------------------------------------;
   ;			    Read data
   ;
   readcol, 'bursts_ft_first_master_reverse.txt', btall0, bfall0, biall0, format = 'A,D,D'
   readcol, 'bursts_ft_second_master_reverse.txt', btall1, bfall1, biall1, format = 'A,D,D'
   readcol, 'bursts_ft_first_master_forward.txt', btall2, bfall2, biall2, format = 'A,D,D'
   readcol, 'bursts_ft_second_master_forward.txt', btall3, bfall3, biall3, format = 'A,D,D'

   btall = [btall0, '-', btall1, '-', btall2, '-', btall3]
   bfall = [bfall0, !Values.F_NAN, bfall1, !Values.F_NAN, bfall2, !Values.F_NAN, bfall3]
   biall = [biall0, !Values.F_NAN, biall1, !Values.F_NAN, biall2, !Values.F_NAN, biall3]
   rev1f0 = bfall0[0]
   rev2f0 = bfall1[0]
   for1f0 = bfall2[0]
   for2f0 = bfall3[0]

   indices = where(btall eq '-')
   indices = [-1, indices]
   drift = fltarr(n_elements(indices))
   vels = fltarr(n_elements(indices))
   displ = fltarr(n_elements(indices))
   flength = fltarr(n_elements(indices))
   start_f = fltarr(n_elements(indices))
   end_f = fltarr(n_elements(indices))
   lifetime =  fltarr(n_elements(indices))

   ;-------------------------------------;
   ;
   ;		Plot frequency v time
   ;
   bsample = [btall2, btall3]    
   ncolors1 = n_elements(where(bsample eq '-')) + 2	; two more bursts than '-'
   colors1 = (findgen(ncolors1)*(80)/(ncolors1-1))
   ncolors2 = n_elements(indices) - ncolors1 - 2
   colors2 = (findgen(ncolors2)*(255 - 150)/(ncolors2-1)) + 150.0
   colors = [colors2, colors1]

   ;setup_ps, 'figures/'+figs_folder+'/freq_time_data.eps'
   window, 0
      FOR i=n_elements(indices)-2, 0, -1 DO BEGIN ;backwards to plot the shortest bursts first.
        bt = btall[indices[i]+1: indices[i+1]-1]  ; This selects one burst.
        bf = bfall[indices[i]+1: indices[i+1]-1]
        bi = biall[indices[i]+1: indices[i+1]-1]	
        tsec = anytim(bt[*], /utim) - anytim(bt[0], /utim)
    
        ;Plot density agains radius here
        dens = freq2dens(bf*1e6)
   		rads = density_to_radius(dens, model = model, fold = 1)

        if i eq n_elements(indices)-2 then begin
        	plot, rads, dens, $
        		;/ylog, $
        		xtitle = 'Heliocentric distance (Rsun)', $
        		ytitle = 'Density (cm!U-3!N)', $
        		/xs, $
        		/ys
	        rads_cm = rads*6.9e10		
	       
	        result = INT_TABULATED( rads_cm, dens )	
	        print, result	
	        column_dens = result
	        displ = abs(rads[0] - rads[n_elements(rads)-1])
	        end_f = bf[n_elements(bf)-1]

			kins_result = linfit(tsec, bf, yfit = yfit)
			start = [kins_result[1]]
			fit = 'p[0]*x + ' +	string(bf[0], format='(f5.2)')
			kins_result = mpfitexpr(fit, tsec , bf, err, yfit=ftfit, start)
			kins = calc_vel(tsec, ftfit, model =  model)
			vels = abs(kins[0]) 

			start_f = bf[0]

        endif else begin
			plot, rads, dens, color = colors[i] ;(245 - i*col_scale)
			rads_cm = rads*6.9e10		
	        result = INT_TABULATED( rads_cm, dens )	
	        column_dens = [column_dens, result]	
	        displ = [displ, abs(rads[0] - rads[n_elements(rads)-1])]
	        end_f = [end_f, bf[n_elements(bf)-1]]


			kins_result = linfit(tsec, bf, yfit = yfit)
			start = [kins_result[1]]
			fit = 'p[0]*x + ' +	string(bf[0], format='(f5.2)')
			kins_result = mpfitexpr(fit, tsec , bf, err, yfit=ftfit, start)
			kins = calc_vel(tsec, ftfit, model =  model)
			vels = [vels, abs(kins[0])]

			start_f = [start_f, bf[0]]


        endelse 
           
      ENDFOR
    
   	order_mag = 1e16
    index_forward = where(round(end_f) lt 35)
    cdf = column_dens[index_forward]

    index_reverse = where(round(end_f) gt 35)
    cdr = column_dens[index_reverse]

  

    good_index = where(column_dens gt 0.0 and column_dens/order_mag lt 1e4)
    column_dens = column_dens[good_index]
    displ = displ[good_index]
    vels = vels[good_index]

    cdr = (cdr)[where(cdr gt 0.0 and cdr/order_mag lt 1e4)]
    cdf = (cdf)[where(cdf gt 0.0 and cdf/order_mag lt 1e4)]


    window, 1
    loadct, 0
	cghistoplot, column_dens/order_mag, binsize=5, /fill, $
		xr = [0.00, 200.0], $
		;yr = [0, 60], $
		polycolor = 240, $
		/xs, $
		xtitle = 'Column Density (1e6 cm!U-2!N)', $
		ytitle='Number of occurences', $
		MININPUT=0.0 

	window, 2
	set_line_color
	plotsym, 0, /fill
	plot, [(displ*6.95e8)/1e6], [column_dens/order_mag], $
		psym=8, $  
		;xr = [10 , 260], $
		;yr=[1, 50], $
		/xs, $
		/ys, $
		;/ylog, $
		xtitle = 'Displacment (Mm)', $
		ytitle = 'Column Density (10!U6!N cm!U-2!N)' 


	nsims = n_elements(displ)
	sym_size = ( dindgen(nsims)*(3.-1.)/(nsims-1) ) + 1.
	vel_sims = ( dindgen(nsims)*(max(vels) - min(vels))/(nsims-1) ) + min(vels)
	syms = interpol(sym_size, vel_sims, vels)



	loadct, 39
    for i =0, n_elements(displ)-1 do begin
		;Green = 170,  Blue = 220,  Orange = 80,  Red = 50

		if start_f[i] eq rev2f0 then color = 240    ;Second reverse (red)
		if start_f[i] eq for2f0 then color = 0     
		if start_f[i] eq for1f0 then color = 80    ;First forward (blue)	;These are cur short because of rfi, hence the two forward drifters are in separate distributions
		if start_f[i] eq rev1f0 then color = 140    ;First reverse (green)

		if start_f[i] ne rev2f0 and start_f[i] ne for2f0 and start_f[i] ne rev1f0 and start_f[i] ne for1f0 then color=200 ;Second forward (yellow)

		oplot, [(displ[i]*6.95e8)/1e6], [(column_dens/order_mag)[i]], $
			psym=8, $
			color = color, $
			symsize = sym_size[i];, $ ;xr = [0.0 , 20], $
			;yr=[0, 20], $
			;/xs, $
			;/ys, $
			;/noerase 
		
    endfor    	

    stop

END