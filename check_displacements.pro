function calc_vel, tsec, ftfit, model = model

  ftfit_Hz = 2.0*ftfit*1e6                 ; Hz
  dens = freq2dens(ftfit_Hz)
  rads = density_to_radius(dens, model = model, fold = 1)
  result = linfit(tsec, rads)
  velocity = abs(result[1])*6.955e8  ; m/s
  velocity = [velocity/2.997e8]   ; speed of light units
  kins = list(velocity, rads)
  if ftfit[n_elements(ftfit)-1] gt 88.0 then stop
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


pro check_displacements
	model = 'tcd'

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
  folder = '~/Data/2011_sep_22/herringbones'
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
  ;setup_ps, 'figures/'+figs_folder+'/freq_time_fits.eps'  
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

   		if end_f[j] ge 88 then stop  
        
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
  ;device, /close


END