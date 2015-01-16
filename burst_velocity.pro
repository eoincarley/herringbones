function calc_vel, tsec, yfit

  yfit = yfit*1e6                 ; Hz
  dens = freq2dens(yfit)
  rads = density_to_radius(dens)
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
  !p.thick=1
  !p.charsize=1.5
  loadct, 39
  pos = [0.13, 0.1, 0.95, 0.95]
  col_scale = 2.0
  dimen = 600
  xpos = 500
  ypos = 50
  window, 0, xs=dimen, ys=dimen, xpos = xpos, ypos = ypos
  window, 1, xs=dimen, ys=dimen, xpos = xpos + 1.0*dimen, ypos = ypos
  window, 2, xs=dimen, ys=dimen, xpos = xpos, ypos = ypos + 1.0*dimen
  window, 3, xs=dimen, ys=dimen, xpos = xpos + 1.0*dimen, ypos = ypos + 1.0*dimen
  window, 4, xs=dimen, ys=dimen, xpos = xpos + 2.0*dimen, ypos = ypos 
  folder = '/Users/ecarley/Data/hough_hbones/20110922'
  cd, folder
  ;-------------------------------------;
  ;			 Read data
  ;
  readcol,'bursts_bs_hough.txt', btall, bfall, biall, format = 'A,D,D'
  indices = where(btall eq '-')
  indices = [-1, indices]
  drift = fltarr(n_elements(indices))
  vels = fltarr(n_elements(indices))
  displ = fltarr(n_elements(indices))
  
  ;-------------------------------------;
  ;
  ;		 Frequency v Time
  ;
  j = 0
  FOR i=n_elements(indices)-2, 0, -1 DO BEGIN ;backwards to plot the shortest bursts first.

    bt = btall[indices[i]+1: indices[i+1]-1]
    bf = bfall[indices[i]+1: indices[i+1]-1]
    bi = biall[indices[i]+1: indices[i+1]-1]	
    tsec = anytim(bt[*], /utim) - anytim(bt[0], /utim)
    
    tsec = anytim(bt[*], /utim) - anytim(bt[0], /utim)
    result = linfit(tsec, bf, yfit = yfit)
    
    start = [result[1]]
	  fit = 'p[0]*x + 46.0'			
	  result = mpfitexpr(fit, tsec , bf, err, yfit=yfit, start)
    
    drift[j] = result[0]
    
    IF i eq n_elements(indices)-2 THEN BEGIN
    
      wset, 0
      plot, tsec, bf, $
          /ys, $
          xr=[0, 5], $
          yr=[45, 80], $
          position=pos, $
          /normal, $
          title='Back-sub Herringbones 22-Sep-2011', $
          xtitle='Time (s)', $
          ytitle='Frequency (MHz)'
          
      oplot, tsec, bf, color=(250)    
          
      wset, 1    
      plot, tsec, yfit, $
          /ys, $
          xr=[0, 5], $
          yr=[45, 80], $
          position=pos, $
          /normal, $
          title='Back-sub Herringbones 22-Sep-2011', $
          xtitle='Time (s)', $
          ytitle='Frequency (MHz)'    
      
      oplot, tsec, yfit, color=250 
      
      kins = calc_vel(tsec, yfit)
      vels[j] = kins[0]
      radii = kins[1]
      displ[j] = abs(radii[n_elements(radii)-1] - radii[0])
      
      
    ENDIF ELSE BEGIN
      wset, 0
       oplot, tsec, bf,  color=(245 - i*col_scale)
      
      wset, 1	
      oplot, tsec, yfit, color=(245 - i*col_scale)
    
      kins = calc_vel(tsec, yfit)
      vels[j] = kins[0]
      radii = kins[1]
      displ[j] = abs(radii[n_elements(radii)-1] - radii[0])
      
    ENDELSE	

    j = j+1
  ENDFOR

  wset, 2
  pdf = histogram(drift, locations = xbin)
  plot, xbin, pdf, $
    xr=[4, 20], $
    psym=10, $
    /xs, $
    xtitle = 'Drift rate (MHz/s)', $
    ytitle='Number of occurences'
    
   
  wset, 3
  pdf = histogram(vels, binsize=0.05, locations = xbin)
  plot, xbin, pdf, $
    psym=10, $
    /xs, $
    xtitle = 'Velocity (c)', $
    ytitle='Number of occurences'  
    
  wset, 4
  plot, displ, vels, $
    psym=4, $  
    yr = [0.0 , 1.0], $
    /xs, $
    xtitle = 'Beam distance (Rsun)', $
    ytitle = 'Beam velocity (c)'
  
  beam_kins = list(drift, vels, displ, tsec)
  save, beam_kins, filename = 'beam_kins.sav'
  
END

;------------ END MAIN PROCEDURE ------------------;

