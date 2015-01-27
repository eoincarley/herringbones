function calc_vel, tsec, ftfit

  ftfit_Hz = ftfit*1e6                 ; Hz
  dens = freq2dens(ftfit_Hz)
  rads = density_to_radius(dens, /saito, fold = 1)
  result = linfit(tsec, rads)
  velocity = abs(result[1])*6.955e8  ; m/s
  velocity = [velocity/2.997e8]   ; speed of light units
  kins = list(velocity, rads)
  return, kins
  
END

pro burst_velocity

  ; Do a linear fit to (f, t) data.
  ; Convert to density, then height.
  ; Turn into velocity.
  !p.thick = 4
  !p.charsize = 1.5
  loadct, 39
  pos = [0.13, 0.1, 0.95, 0.95]
  col_scale = 1.5
  dimen = 600
  xpos = 500
  ypos = 50 
  folder = '~/Data/22Sep2011_event/herringbones'
  cd, folder
  
  ;-------------------------------------;
  ;			        Read data
  ;
  readcol, 'bursts_ft_first_master_reverse.txt', btall0, bfall0, biall0, format = 'A,D,D'
  readcol, 'bursts_ft_second_master_reverse.txt', btall1, bfall1, biall1, format = 'A,D,D'
  ;readcol, 'bursts_ft_first_master_positive.txt', btall2, bfall2, biall2, format = 'A,D,D'
  ;readcol, 'bursts_ft_second_master_positive.txt', btall3, bfall3, biall3, format = 'A,D,D'
  
  btall = [btall0, '-', btall1] ;, '-', btall2, '-', btall3]
  bfall = [bfall0, !Values.F_NAN, bfall1] ;, !Values.F_NAN, bfall2, !Values.F_NAN, bfall3]
  biall = [biall0, !Values.F_NAN, biall1] ;, !Values.F_NAN, biall2, !Values.F_NAN, biall3]
  
  
  indices = where(btall eq '-')
  indices = [-1, indices]
  drift = fltarr(n_elements(indices))
  vels = fltarr(n_elements(indices))
  displ = fltarr(n_elements(indices))
  flength = fltarr(n_elements(indices))
  start_f = fltarr(n_elements(indices))
  
  ;-------------------------------------;
  ;
  ;		   Plot frequency v Time
  ;
  setup_ps, 'figures/freq_time_data.ps'
      FOR i=n_elements(indices)-2, 0, -1 DO BEGIN ;backwards to plot the shortest bursts first.

        bt = btall[indices[i]+1: indices[i+1]-1]  ; This selects one burst.
        bf = bfall[indices[i]+1: indices[i+1]-1]
        bi = biall[indices[i]+1: indices[i+1]-1]	
        tsec = anytim(bt[*], /utim) - anytim(bt[0], /utim)
    
        if i eq n_elements(indices)-2 then begin
          ;void = execute(w0)
          plot, tsec, bf, $
              /ys, $
              xr=[0, 5], $
              yr=[20, 80], $
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
  setup_ps, 'figures/freq_time_fits.ps'  
      j = 0  
      FOR i=n_elements(indices)-2, 0, -1 DO BEGIN ;backwards to plot the shortest bursts first.

        bt = btall[indices[i]+1: indices[i+1]-1]  ; This selects one burst.
        bf = bfall[indices[i]+1: indices[i+1]-1]
        bi = biall[indices[i]+1: indices[i+1]-1]	
        tsec = anytim(bt[*], /utim) - anytim(bt[0], /utim)  
    
        ; Do fitting
        result = linfit(tsec, bf, yfit = yfit)
        start = [result[1]]
        if round(bf[0]) gt 40 then intersect='43.4'
        if round(bf[0]) lt 40 then intersect='32.8'
        fit = 'p[0]*x + '+	intersect
        result = mpfitexpr(fit, tsec , bf, err, yfit=ftfit, start)
    
        drift[j] = abs(result[0])
        print, drift[j]
        flength[j] = abs(bf[n_elements(bf)-1] - bf[0])
        start_f[j] = bf[0]
        kins = calc_vel(tsec, ftfit)
        vels[j] = abs(kins[0])
        radii = kins[1]
        displ[j] = abs(radii[0] - radii[n_elements(radii)-1])
   
 
        if i eq n_elements(indices)-2 then begin 
          ;void = execute(w1)      
          plot, tsec, ftfit, $
              /ys, $
              xr=[0, 5], $
              yr=[20, 80], $
              position=pos, $
              /normal, $
              xtitle='Time (s)', $
              ytitle='Frequency (MHz)'  
        endif else begin
          ;wset,1
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
  
  stop
  ;--------------------------------;
  ;     Drift rate histogram
  ;

  index33 = where(round(start_f) lt 40)
  drift_low = drift[index33]

  setup_ps, 'figures/hist_dfdt.ps' 
    ; Callisto has 200 channels over 90 MHz -> ~0.45 MHz per pix. Take two pix as resolution: 0.9 MHz
    ; and two pixels in time = 0.5 seconds -> Drift res is 0.9/0.5 = 1.8 MHz/s
    cghistoplot, drift, binsize=1.8, /fill, $
      xr = [2, 18], $
      yr = [0, 45], $
      /xs, $
      /ys, $
      xtitle = 'Drift rate (MHz s!U-1!N)', $
      ytitle='Number of occurences'
   

   loadct, 1  
     cghistoplot, drift_low, binsize=1.8, /fill, $
      xr = [2, 18], $
      polycolor = 230, $
      yr = [0, 45], $
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
  rads = density_to_radius(dens, /saito, fold = 5)
  binsize = (((rads[0] - rads[1])/(0.5))*6.955e8)/2.9e8 
  ;wset, 3
  
  setup_ps, 'figures/hist_vels.ps' 
    
    cghistoplot, vels, binsize=binsize, /fill, $
      xr = [0.1, 0.6], $
      yr = [0, 60], $
      xtitle = 'Velocity (c)', $
      ytitle='Number of occurences'  
      
    loadct, 1  
    cghistoplot, vels[index33], binsize=binsize, $
      xr = [0.1, 0.6], $
      yr = [0, 60], $
      polycolor = 230, $
      /fill, $
      xtitle = ' ', $
      ytitle=' ', $
      /noerase
      
  device, /close  

  ;------------------------------;
  ;   Displacement histogram
  ;
  loadct, 39
  plotsym, 0, /fill
  
  setup_ps, 'figures/scatter_flen_drift.ps' 
    plot, flength, drift, $
      psym=8, $  
      yr = [3 , 15], $
      xr=[0, 45], $
      /xs, $
      /ys, $
      xtitle = 'Frequency span (MHz)', $
      ytitle = 'Drift rate (MHz s!U-1!N)'
      set_line_color
     for i =0, n_elements(flength)-1 do begin
      if round(start_f[i]) eq 43 then color=3 
      if round(start_f[i]) eq 33 then color=5    
      if round(start_f[i]) eq 31 then color=4 
      if round(start_f[i]) eq 42 then color=2 
      
      oplot, [0, flength[i]], [0, drift[i]], $
          psym=8, $
          color = color
     endfor    
  device, /close
  set_plot, 'x'  
  
  beam_kins = list(drift, vels, displ, tsec, flength)
  save, beam_kins, filename = 'beam_kins.sav'
  

END

;------------ END MAIN PROCEDURE ------------------;

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

