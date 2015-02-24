
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

pro calc_kins, btall, bfall, biall, $
				vels = vels, $
				end_f = end_f, $
				model = model
				
	indices = where(btall eq '-')
	indices = [-1, indices]
	drift = fltarr(n_elements(indices))
	vels = fltarr(n_elements(indices))
	displ = fltarr(n_elements(indices))
	flength = fltarr(n_elements(indices))
	start_f = fltarr(n_elements(indices))
	end_f = fltarr(n_elements(indices))
	lifetime =  fltarr(n_elements(indices))			
				
	loadct, 0
	j = 0  
	FOR i=n_elements(indices)-2, 0, -1 DO BEGIN ;backwards to plot the shortest bursts first.
		bt = btall[indices[i]+1: indices[i+1]-1]  ; This selects one burst.
		bf = bfall[indices[i]+1: indices[i+1]-1]
		bi = biall[indices[i]+1: indices[i+1]-1]	
		tsec = anytim(bt[*], /utim) - anytim(bt[0], /utim)  

		; Do fitting
		result = linfit(tsec, bf, yfit = yfit)
		start = [result[1]]
		if round(bf[0]) eq 43 then intersect = '43.4'   ;second reverse
		if round(bf[0]) eq 32 then intersect = '32.0'   ;first reverse
		if round(bf[0]) eq 42 then intersect = '41.5'   ;second forward
		if round(bf[0]) eq 31 then intersect = '31.2'   ;first forward
		fit = 'p[0]*x + '+	intersect
		result = mpfitexpr(fit, tsec , bf, err, yfit=ftfit, start)

		drift[j] = result[0]
		flength[j] = abs(bf[n_elements(bf)-1] - bf[0])
		start_f[j] = bf[0]
		end_f[j] = bf[n_elements(bf)-1]
		kins = calc_vel(tsec, ftfit, model =  model)
		vels[j] = abs(kins[0]) 
		radii = kins[1]
		displ[j] = abs(radii[0] - radii[n_elements(radii)-1])
		lifetime[j] = tsec[n_elements(tsec)-1]
		j = j+1
	ENDFOR
	
	drift = drift[where(drift ne 0.0)]
	vels = vels[where(vels ne 0.0)]
	flength = flength[where(flength ne 0.0)]
	displ = displ[where(displ ne 0.0)]
	start_f = start_f[where(start_f ne 0.0)]
	end_f = end_f[where(end_f ne 0.0)]
	lifetime = lifetime[where(end_f ne 0.0)]

END
	

pro setup_ps, name
  
  set_plot,'ps'
  !p.font=0
  !p.charsize=1.5
  device, filename = name, $
          /color, $
          /helvetica, $
          /inches, $
          xsize=7, $
          ysize=7, $
          /encapsulate, $
          yoffset=5

end


pro burst_vels_separate, model

	; Do a linear fit to (f, t) data.
	; Convert to density, then height.
	; Turn into velocity.
	!p.thick = 4
	!p.charsize = 1.5
	loadct, 39
	pos = [0.13, 0.1, 0.95, 0.95]
	col_scale = 1.3
	dimen = 600
	xpos = 500
	ypos = 50 
	folder = '~/Data/22Sep2011_event/herringbones'
	figs_folder = model
	cd, folder

	;-------------------------------------;
	;			  Read data
	;
	readcol, 'bursts_ft_first_master_reverse.txt', btall0, bfall0, biall0, format = 'A,D,D'
	readcol, 'bursts_ft_first_master_forward.txt', btall1, bfall1, biall1, format = 'A,D,D'

	readcol, 'bursts_ft_second_master_reverse.txt', btall2, bfall2, biall2, format = 'A,D,D'
	readcol, 'bursts_ft_second_master_forward.txt', btall3, bfall3, biall3, format = 'A,D,D'

	btall_first = [btall0, '-', btall1] 		 
	bfall_first = [bfall0, !Values.F_NAN, bfall1] 
	biall_first = [biall0, !Values.F_NAN, biall1]

	btall_second = [btall2, '-', btall3] 
	bfall_second = [bfall2, !Values.F_NAN, bfall3] 
	biall_second = [biall2, !Values.F_NAN, biall3] 

	calc_kins, btall_first, bfall_first, biall_first, $
			vels = vels_first, $
			end_f = end_f_first, $
			model = model
					
	calc_kins, btall_second, bfall_second, biall_second, $
			vels = vels_second, $
			end_f = end_f_second, $
			model = model
					
	index_forward_first = where( round(end_f_first) lt 35 )
	index_reverse_first = where( round(end_f_first) gt 35 )
	
	index_forward_second = where( round(end_f_second) lt 35 )
	index_reverse_second = where( round(end_f_second) gt 35 )
	
	;--------------------------------;
	;     Velocity histogram
	; 
	; binsize = speed resolution given f and t resolution and model used.
	dens = freq2dens([42.0e6, 42.9e6])
	rads = density_to_radius(dens, model = model, fold = 1)
	binsize = (((rads[0] - rads[1])/(0.5))*6.955e8)/2.9e8 

	setup_ps, 'figures/'+figs_folder+'/hist_vels_first.eps' 

		loadct, 0
		cghistoplot, vels_first, binsize=binsize, /fill, $
			xr = [0.00, 0.35], $
			yr = [0, 35], $
			polycolor = 200, $
			orientation = [45], $
			/xs, $
			xtitle =  '  ', $
			ytitle='Number of occurences', $
			MININPUT=0.0 
		
		loadct, 74
		cghistoplot, vels_first[index_forward_first], binsize=binsize, /line_fill, $
			xr = [0.00, 0.35], $
			yr = [0, 35], $
			polycolor = 50, $
			orientation = [45], $
			/xs, $
			xtitle =  '  ', $
			ytitle='Number of occurences', $
			/noerase, $
			MININPUT=0.0 

		cghistoplot, vels_first[index_reverse_first], binsize=binsize, $
			xr = [0.00, 0.35], $
			/xs, $
			yr = [0, 35], $
			orientation = [-45], $
			polycolor = 170, $
			/line_fill, $
			xtitle = '  ', $
			ytitle=' ', $
			/noerase, $
			MININPUT=0.0 
			
		loadct, 0
		xyouts, 0.55, 0.01, 'Velocity (c)', alignment=0.5, /normal	
		
		loadct, 74
		legend, ['First set forward', 'First set reverse'], psym=[6, 6], color=[50, 170], $
				box=0, /right
		
	device, /close  
	
	;------------------------------------------;
	; Now plot second set of bursts separately
	;------------------------------------------;
	setup_ps, 'figures/'+figs_folder+'/hist_vels_second.eps' 
		loadct, 0
		cghistoplot, vels_second, binsize=binsize, /fill, $
			xr = [0.00, 0.35], $
			yr = [0, 35], $
			polycolor = 200, $
			orientation = [45], $
			/xs, $
			xtitle =  '  ', $
			ytitle='Number of occurences', $
			MININPUT=0.0 
			
		loadct, 74	
		cghistoplot, vels_second[index_reverse_second], binsize=binsize, $
			xr = [0.00, 0.35], $
			/xs, $
			yr = [0, 35], $
			orientation = [-45], $
			polycolor = 220, $
			/line_fill, $
			xtitle = '  ', $
			ytitle=' ', $
			/noerase, $
			MININPUT=0.0 
			
		cghistoplot, vels_second[index_forward_second], binsize=binsize, $
			xr = [0.00, 0.35], $
			/xs, $
			yr = [0, 35], $
			orientation = [40], $
			polycolor = 80, $
			/line_fill, $
			xtitle = '  ', $
			ytitle=' ', $
			/noerase, $
			MININPUT=0.0 	

		loadct, 0
		xyouts, 0.55, 0.01, 'Velocity (c)', alignment=0.5, /normal
		
		loadct, 74
		legend, ['Second set forward', 'Second set reverse'], psym=[6, 6], color=[80, 220], $
				box=0, /right
	device, /close  
	
	;------------------------------------------;
	; 	 Now plot first and second bursts
	;------------------------------------------;
	setup_ps, 'figures/'+figs_folder+'/hist_vels_first_second.eps' 
		loadct, 0
		cghistoplot, [vels_first, vels_second], binsize=binsize, /fill, $
			xr = [0.00, 0.35], $
			yr = [0, 50], $
			polycolor = 200, $
			orientation = [45], $
			/xs, $
			xtitle =  '  ', $
			ytitle='Number of occurences', $
			MININPUT=0.0 
			
		loadct, 74
		cghistoplot, vels_second, binsize=binsize, /fill, $
			xr = [0.00, 0.35], $
			yr = [0, 50], $
			polycolor = 190, $
			/xs, $
			xtitle =  '  ', $
			ytitle='Number of occurences', $
			MININPUT=0.0, $
			/noerase
			
		cghistoplot, vels_first, binsize=binsize, /line_fill, $
			xr = [0.00, 0.35], $
			yr = [0, 50], $
			polycolor = 40, $
			orientation = [90], $
			/xs, $
			xtitle =  '  ', $
			ytitle='Number of occurences', $
			MININPUT=0.0, $
			/noerase
			
		loadct, 0
		xyouts, 0.55, 0.01, 'Velocity (c)', alignment=0.5, /normal	
			
		loadct, 74	
		legend, ['First set', 'Second set'], psym=[6, 6], color=[40, 190], $
				box=0, /right
	device, /close		
 

END  
