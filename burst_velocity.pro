function calc_vel, tsec, ftfit, model = model

  ftfit_Hz = ftfit*1e6                 ; Hz
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
  
  indices = where(btall eq '-')
  indices = [-1, indices]
  drift = fltarr(n_elements(indices))
  vels = fltarr(n_elements(indices))
  displ = fltarr(n_elements(indices))
  flength = fltarr(n_elements(indices))
  start_f = fltarr(n_elements(indices))
  end_f = fltarr(n_elements(indices))
  
  ;-------------------------------------;
  ;
  ;		   Plot frequency v Time
  ;
  setup_ps, 'figures/'+figs_folder+'/freq_time_data.ps'
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
        endif else begin
          ;wset,0
          oplot, tsec, bf,  color=(245 - i*col_scale)
        endelse      
      ENDFOR
  device, /close
 
  ;-------------------------------------;
  ;
  ;		  Frequency v Time (FITTING)
  ;
  setup_ps, 'figures/'+figs_folder+'/freq_time_fits.ps'  
      j = 0  
      FOR i=n_elements(indices)-2, 0, -1 DO BEGIN ;backwards to plot the shortest bursts first.

        bt = btall[indices[i]+1: indices[i+1]-1]  ; This selects one burst.
        bf = bfall[indices[i]+1: indices[i+1]-1]
        bi = biall[indices[i]+1: indices[i+1]-1]	
        tsec = anytim(bt[*], /utim) - anytim(bt[0], /utim)  
    
        ; Do fitting
        result = linfit(tsec, bf, yfit = yfit)
        start = [result[1]]
        if round(bf[0]) eq 43 then intersect='43.4'   ;second reverse
        if round(bf[0]) eq 32 then intersect='32.0'   ;first reverse
        if round(bf[0]) eq 42 then intersect='41.5'   ;second forward
        if round(bf[0]) eq 31 then intersect='31.2'   ;first forward
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
   
 
        if i eq n_elements(indices)-2 then begin      
          plot, tsec, ftfit, $
              /ys, $
              xr=[0, 6], $
              yr=[10, 90], $
              position=pos, $
              /normal, $
              xtitle='Time (s)', $
              ytitle='Frequency (MHz)'  
        endif else begin
          oplot, tsec, ftfit, color=(245 - i*col_scale)
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

  ;--------------------------------;
  ;     Drift rate histogram
  ;

  index_forward = where(round(end_f) lt 34)
  drift_forward = drift[index_forward]
  
  index_reverse = where(round(end_f) gt 34)
  drift_reverse = drift[index_reverse]

    

  setup_ps, 'figures/'+figs_folder+'/hist_dfdt.ps' 
    ; Callisto has 200 channels over 90 MHz -> ~0.45 MHz per pix. Take two pix as resolution: 0.9 MHz
    ; and two pixels in time = 0.5 seconds -> Drift res is 0.9/0.5 = 1.8 MHz/s
    cghistoplot, drift_reverse, binsize=1.9, /fill, $
      xr = [-10, 20], $
      yr = [0, 40], $
      /xs, $
      /ys, $
      xtitle = 'Drift rate (MHz s!U-1!N)', $
      ytitle='Number of occurences'
   

   loadct, 1  
     cghistoplot, drift_forward, binsize=1.9, /fill, $
      xr = [-10, 20], $
      polycolor = 230, $
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
  ;wset, 3
  
  setup_ps, 'figures/'+figs_folder+'/hist_vels.ps' 
    
    cghistoplot, vels, binsize=binsize, /fill, $
      xr = [0.05, 0.45], $
      yr = [0, 60], $
      /xs, $
      xtitle = 'Velocity (c)', $
      ytitle='Number of occurences'  
      
    loadct, 1  
    cghistoplot, vels[index_forward], binsize=binsize, $
      xr = [0.05, 0.45], $
      /xs, $
      yr = [0, 60], $
      polycolor = 230, $
      /fill, $
      xtitle = ' ', $
      ytitle=' ', $
      /noerase
      
  device, /close  
  stop
  ;--------------------------------;
  ;     Displacement histogram
  ; 
  ; binsize = speed resolution given f and t resolution and model used.
  dens = freq2dens([42.0e6, 44.9e6])
  rads = density_to_radius(dens, model = model, fold = 5)
  binsize = (rads[0] - rads[1])
  binsize=0.05
  
  loadct, 1  
  setup_ps, 'figures/'+figs_folder+'/hist_displ.ps' 
    
    cghistoplot, displ[index_reverse], binsize=binsize, /fill, $
      xr = [0.0, 0.6], $
      yr = [0, 40], $
      polycolor = 230, $
      /xs, $
      xtitle = 'Displacement (Rsun)', $
      ytitle='Number of occurences'
      
    cghistoplot, displ[index_forward], binsize=binsize, $
      xr = [0.0, 0.6], $
      yr = [0, 40], $
      /fill, $
      /xs, $
      xtitle = ' ', $
      ytitle=' ', $
      /noerase
      
  device, /close  

  ;------------------------------;
  ;   Frequency span and drift histogram
  ;
  
  driftabs = abs(drift)
  loadct, 39
  plotsym, 0, /fill
  
  setup_ps, 'figures/'+figs_folder+'/scatter_flen_drift.ps' 
    plot, flength, driftabs, $
      psym=8, $  
      yr = [2 , 15], $
      xr=[0, 45], $
      /xs, $
      /ys, $
      xtitle = 'Frequency span (MHz)', $
      ytitle = 'Drift rate (MHz s!U-1!N)'
      set_line_color
     for i =0, n_elements(flength)-1 do begin
      if round(start_f[i]) eq 43 then color=3 
      if round(start_f[i]) eq 31 then color=5    
      if round(start_f[i]) eq 32 then color=4 
      if round(start_f[i]) eq 42 then color=2 
      
      oplot, [flength[i]], [driftabs[i]], $
          psym=8, $
          color = color
     endfor    
  device, /close
  set_plot, 'x'  
  
  beam_kins = list(drift, vels, displ, tsec, flength, figs_folder)
  save, beam_kins, filename = 'beam_kins.sav'
  
 
END

;------------ END MAIN PROCEDURE ------------------;

