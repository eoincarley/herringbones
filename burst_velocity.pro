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


pro burst_velocity, model

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
  ;			        Read data
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
  for2f0 = bfall2[0]
  for1f0 = bfall3[0]
  
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
  ;		   Plot frequency v Time
  ;
  bsample = [btall2, btall3]    
  ncolors1 = n_elements(where(bsample eq '-')) + 2	; two more bursts than '-'
  colors1 = (findgen(ncolors1)*(80)/(ncolors1-1))
  ncolors2 = n_elements(indices) - ncolors1 - 2
  colors2 = (findgen(ncolors2)*(255 - 150)/(ncolors2-1)) + 150.0
  colors = [colors2, colors1]
  
  setup_ps, 'figures/'+figs_folder+'/freq_time_data.eps'
      FOR i=n_elements(indices)-2, 0, -1 DO BEGIN ;backwards to plot the shortest bursts first.
        bt = btall[indices[i]+1: indices[i+1]-1]  ; This selects one burst.
        bf = bfall[indices[i]+1: indices[i+1]-1]
        bi = biall[indices[i]+1: indices[i+1]-1]	
        tsec = anytim(bt[*], /utim) - anytim(bt[0], /utim)
    
        if i eq n_elements(indices)-2 then begin
          ;void = execute(w0)
          plot, tsec, bf, $
              /ys, $
              xr=[0, 6], $
              yr=[10, 90], $
              position=pos, $
              /normal, $
              xtitle='Time (s)', $
              ytitle='Frequency (MHz)'
              loadct, 74, /silent
        endif else begin
          oplot, tsec, bf, color = colors[i] ;(245 - i*col_scale)
        endelse      
      ENDFOR
  device, /close
 
  ;-------------------------------------;
  ;
  ;		  Frequency v Time (FITTING)
  ;
  loadct, 0
  setup_ps, 'figures/'+figs_folder+'/freq_time_fits.eps'  
      j = 0  
      FOR i=n_elements(indices)-2, 0, -1 DO BEGIN ;backwards to plot the shortest bursts first.
        bt = btall[indices[i]+1: indices[i+1]-1]  ; This selects one burst.
        bf = bfall[indices[i]+1: indices[i+1]-1]
        bi = biall[indices[i]+1: indices[i+1]-1]	
        tsec = anytim(bt[*], /utim) - anytim(bt[0], /utim)  
    
        ; Do fitting
        result = linfit(tsec, bf, yfit = yfit)
        start = [result[1]]
        
        fit = 'p[0]*x + ' +	string(bf[0], format='(f5.2)')
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
        
        if i eq n_elements(indices)-2 then begin      
          plot, tsec, ftfit, $
              /ys, $
              xr=[0, 6], $
              yr=[10, 90], $
              position=pos, $
              /normal, $
              xtitle='Time (s)', $
              ytitle='Frequency (MHz)'  
              
           loadct, 74, /silent   
        endif else begin
          oplot, tsec, ftfit, color = colors[i] ;(245 - i*col_scale)
        endelse        
      
        j = j+1
      ENDFOR
  device, /close

  drift = drift[where(drift ne 0.0)]
  vels = vels[where(vels ne 0.0)]
  flength = flength[where(flength ne 0.0)]
  displ = displ[where(displ ne 0.0)]
  start_f = start_f[where(start_f ne 0.0)]
  end_f = end_f[where(end_f ne 0.0)]
  lifetime = lifetime[where(end_f ne 0.0)]

  ;--------------------------------;
  ;     Drift rate histogram
  ;
  index_forward = where(round(end_f) lt 35)
  drift_forward = drift[index_forward]
  
  index_reverse = where(round(end_f) gt 35)
  drift_reverse = drift[index_reverse]
    

  setup_ps, 'figures/'+figs_folder+'/hist_dfdt.eps' 
    ; Callisto has 200 channels over 90 MHz -> ~0.45 MHz per pix. Take two pix as resolution: 0.9 MHz
    ; and two pixels in time = 0.5 seconds -> Drift res is 0.9/0.5 = 1.8 MHz/s
    loadct, 1  
    cghistoplot, drift_reverse, binsize=1.9, /line_fill, $
      xr = [-10, 22], $
      yr = [0, 40], $
      polycolor = 230, $
      ORIENTATION=[-45], $
      /xs, $
      /ys, $
      xtitle = 'Drift rate (MHz s!U-1!N)', $
      ytitle='Number of occurences'
   
     cghistoplot, drift_forward, binsize=1.9, /fill, $
      xr = [-10, 22], $
      yr = [0, 40], $
      /xs, $
      /ys, $
      xtitle = ' ', $
      ytitle=' ', $
      /noerase
    
  device, /close

  ;--------------------------------;
  ;     Velocity histogram
  ; 
  ; binsize = speed resolution given f and t resolution and model used.
  dens = freq2dens([42.0e6, 42.9e6])
  rads = density_to_radius(dens, model = model, fold = 1)
  binsize = (((rads[0] - rads[1])/(0.5))*6.955e8)/2.9e8 
  
  setup_ps, 'figures/'+figs_folder+'/hist_vels_rev_forw.eps' 
    
    loadct, 0
    cghistoplot, vels, binsize=binsize, /fill, $
      xr = [0.00, 0.35], $
      yr = [0, 50], $s
      polycolor = 240, $
      /xs, $
      xtitle = 'Velocity (c)', $
      ytitle='Number of occurences', $
      MININPUT=0.0 
      
    cghistoplot, vels[index_forward], binsize=binsize, $
      xr = [0.00, 0.35], $
      /xs, $
      yr = [0, 50], $
      /fill, $
      xtitle = ' ', $
      ytitle=' ', $
      /noerase, $
      MININPUT=0.0 
      
    loadct, 1  
    cghistoplot, vels[index_reverse], binsize=binsize, $
      xr = [0.00, 0.35], $
      /xs, $
      yr = [0, 50], $
      polycolor = 230, $
      /line_fill, $
      ORIENTATION=[-45], $
      xtitle = ' ', $
      ytitle=' ', $
      /noerase, $
      MININPUT=0.0
      
  device, /close  
  
  ;--------------------------------;
  ;     Lifetime histogram
  ; 
  ; binsize = speed resolution given f and t resolution and model used.
  binsize = 0.5
  
  setup_ps, 'figures/'+figs_folder+'/hist_lifetime.eps' 
   lifetime = round(lifetime*10.0)/10.0
    loadct, 0
    cghistoplot, lifetime, binsize=binsize, /fill, $
      xr = [0.0, 6], $
      yr = [0, 35], $
      polycolor = 240, $
      /xs, $
      xtitle = 'Lifetime (s)', $
      ytitle='Number of occurences', $
      MININPUT=0.0 
      
    cghistoplot, lifetime[index_forward], binsize=binsize, $
      xr = [0.0, 6], $
      /xs, $
      yr = [0, 35], $
      /fill, $
      xtitle = ' ', $
      ytitle=' ', $
      /noerase, $
      MININPUT=0.0 
      
    loadct, 1  
    cghistoplot, lifetime[index_reverse], binsize=binsize, $
      xr = [0.0, 6], $
      /xs, $
      yr = [0, 35], $
      polycolor = 230, $
      /line_fill, $
      ORIENTATION=[-45], $
      xtitle = ' ', $
      ytitle=' ', $
      /noerase, $
      MININPUT=0.0 
      
      lt_forward = round(lifetime[index_forward]*100.0)/100.0
      lt_reverse = round(lifetime[index_reverse]*100.0)/100.0
      
      print, 'Mean lifetime of forward drifters: '+string(mean(lt_forward))
      print, 'Mode lifetime of forward drifters: '+string(mode(lt_forward))
      print,' '
      print, 'Mean lifetime of reverse drifters: '+string(mean(lt_reverse))
      print, 'Mode lifetime of reverse drifters: '+string(mode(lt_reverse))
      
      
  device, /close  
  
  ;--------------------------------;
  ;     Displacement histogram
  ; 
  ; binsize = speed resolution given f and t resolution and model used.
  dens = freq2dens([42.0e6, 44.9e6])
  rads = density_to_radius(dens, model = model, fold = 5)
  binsize = (rads[0] - rads[1])
  ;binsize=0.05


  setup_ps, 'figures/'+figs_folder+'/hist_displ.eps' 
    
    loadct, 0 
    cghistoplot, displ, binsize=binsize, /fill, $
      xr = [0.0, 0.45], $
      yr = [0, 40], $
      polycolor = 240, $
      /xs, $
      xtitle = 'Displacement (Rsun)', $
      ytitle='Number of occurences', $
      MININPUT=0.0
      
    cghistoplot, displ[index_forward], binsize=binsize, $
      xr = [0.0, 0.45], $
      yr = [0, 40], $
      /fill, $
      /xs, $
      xtitle = ' ', $
      ytitle = ' ', $
      /noerase, $
      MININPUT=0.0  
      
    loadct, 1
    cghistoplot, displ[index_reverse], binsize=binsize, $
      xr = [0.0, 0.45], $
      yr = [0, 40], $
      /line_fill, $
      polycolor = 230, $
      orientation = [-45], $
      xtitle = ' ', $
      ytitle=' ', $
      /noerase, $
      MININPUT=0.0

      
  device, /close  

  ;----------------------------------------;
  ;   	Drift and frequency length
  ;
  
  driftabs = abs(drift)
  loadct, 0
  plotsym, 0, /fill
  
  setup_ps, 'figures/'+figs_folder+'/scatter_drift_flen.eps' 
    plot, [driftabs], [flength], $
      psym=8, $  
      xr = [0.0 , 20], $
      ;yr=[0, 0.35], $
      /xs, $
      /ys, $
      xtitle = 'Drift (MHz s!U-1!N)', $
      ytitle = 'Frequency span (MHZ)'
      
      
     loadct, 74
     for i =0, n_elements(flength)-1 do begin
		  ;Green = 170,  Blue = 220,  Orange = 80,  Red = 50
		  if start_f[i] eq rev2f0 then color = 220    ;Second reverse
      if start_f[i] eq for2f0 then color = 80     ;Second forward
      if start_f[i] eq for1f0 then color = 50     ;First forward 
      if start_f[i] eq rev1f0 then color = 170    ;First reverse
	  
		  oplot, [driftabs[i]], [flength[i]], $
			  psym=8, $
			  color = color
     endfor    
  device, /close
  set_plot, 'x'  
  
  ;----------------------------------------;
  ;   	Velocities and displacement
  ;
  
  driftabs = abs(drift)
  loadct, 0
  plotsym, 0, /fill
  
  setup_ps, 'figures/'+figs_folder+'/scatter_vels_displ.eps' 
    plot, [vels], [displ], $
      psym=8, $  
      xr = [0.0 , 0.35], $
      ;yr=[0, 0.35], $
      /xs, $
      /ys, $
      xtitle = 'Velocity (c)', $
      ytitle = 'Displacment (Rsun)'
      
      
     loadct, 74
     for i =0, n_elements(flength)-1 do begin
		  ;Green = 170,  Blue = 220,  Orange = 80,  Red = 50
		  if start_f[i] eq rev2f0 then color = 220    ;Second reverse
      if start_f[i] eq for2f0 then color = 80     ;Second forward
      if start_f[i] eq for1f0 then color = 50     ;First forward 
      if start_f[i] eq rev1f0 then color = 170    ;First reverse
	  
		  oplot, [vels[i]], [displ[i]], $
			  psym=8, $
			  color = color
     endfor    
  device, /close
  set_plot, 'x'  

  ;------------------------------------------------------------------------;
  ;   Lifetime v normalised drift rate (Like Figure 4 Mann et al. 2005)
  ;
  setup_ps, 'figures/'+figs_folder+'/scatter_lifetime_drift.eps' 
    plot, [start_f/driftabs], [lifetime], $
      psym=8, $  
      ;xr = [0.0 , 0.35], $
      ;yr=[0, 0.35], $
      /xs, $
      /ys, $
      xtitle = '(Start frq.)/(Drift rate)', $
      ytitle = 'Lifetime (s)'
      
      
     loadct, 74
     for i =0, n_elements(flength)-1 do begin
      ;Green = 170,  Blue = 220,  Orange = 80,  Red = 50
      if start_f[i] eq rev2f0 then color = 220    ;Second reverse
      if start_f[i] eq for2f0 then color = 80     ;Second forward
      if start_f[i] eq for1f0 then color = 50     ;First forward 
      if start_f[i] eq rev1f0 then color = 170    ;First reverse
    
      oplot, [start_f[i]/driftabs[i]], [lifetime[i]], $
        psym=8, $
        color = color
     endfor    
  device, /close
  set_plot, 'x'  
  

  beam_kins = list(drift, vels, displ, tsec, flength, figs_folder)
  save, beam_kins, filename = 'beam_kins.sav'
             
  dens = freq2dens(end_f*1e6)
  stop_rads = density_to_radius(dens, model = model, fold = 1)
  
  end_45 = end_f(where(round(start_f) eq 43))*1e6                 ; Hz
  dens = freq2dens(end_45)
  rads45 = density_to_radius(dens, model = model, fold = 1)
  
  end_32 = end_f(where(round(start_f) eq 32))*1e6                 ; Hz
  dens = freq2dens(end_32)
  rads32 = density_to_radius(dens, model = model, fold = 1)
  
  
END

;------------ END MAIN PROCEDURE ------------------;