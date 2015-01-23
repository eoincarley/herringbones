pro burst_intensity

  !p.thick = 2 
  !p.charsize = 1.5
  !p.charthick = 1.5
  loadct, 39
  pos = [0.13, 0.1, 0.95, 0.95]
  col_scale = 2.5
  dimen = 600
  xpos = 500
  ypos = 50
  ;window, 0, xs=dimen, ys=dimen, xpos = xpos, ypos = ypos
  ;window, 1, xs=dimen, ys=dimen, xpos = xpos + 1.0*dimen, ypos = ypos
  ;window, 2, xs=dimen, ys=dimen, xpos = xpos, ypos = ypos + 1.0*dimen
  ;window, 3, xs=dimen, ys=dimen, xpos = xpos + 1.0*dimen, ypos = ypos + 1.0*dimen
  ;window, 4, xs=dimen, ys=dimen, xpos = xpos + 2.0*dimen, ypos = ypos 
  ;window, 5, xs=dimen, ys=dimen, xpos = xpos + 2.0*dimen, ypos = ypos + 1.0*dimen
  folder = '~/Data/22Sep2011_event/herringbones'
  cd, folder
  
  ;------------------------;
  ;			 Read data
  ;
  readcol,'bursts_ft_first_master_reverse.txt', btall1, bfall1, biall1, format = 'A,D,D'
  readcol,'bursts_ft_second_master_reverse.txt', btall2, bfall2, biall2, format = 'A,D,D'
  ;readcol,'bursts_bs_hough_second.txt', btall3, bfall3, biall3, format = 'A,D,D'
  
  btall = [btall1, '-', btall2];, '-', btall3]
  bfall = [bfall1, !Values.F_NAN, bfall2];, !Values.F_NAN, bfall3]
  biall = [biall1, !Values.F_NAN, biall2];, !Values.F_NAN, biall3]
  
  indices = where(btall eq '-')
  indices = [-1, indices]
  
  restore, 'beam_kins.sav', /verb
  drift = beam_kins[0]
  vels = beam_kins[1]
  displ = beam_kins[2]
  max_bi = fltarr(n_elements(vels))
  ;-----------------------------------;
  ;			 First intensity v time
  ;
  setup_ps, 'figures/freq_inten_data.ps'
      FOR i=n_elements(indices)-2, 0, -1 DO BEGIN

        bt = btall[indices[i]+1: indices[i+1]-1]
        bf = bfall[indices[i]+1: indices[i+1]-1]
        bi = biall[indices[i]+1: indices[i+1]-1]
        tsec = anytim(bt[*], /utim) - anytim(bt[0], /utim)
    
        IF i eq n_elements(indices)-2 THEN BEGIN
    
          ;wset, 5
          plot, bf, bi, $
            psym=1, $
            xr=[32, 80], $
            yr=[0, 45], $
            /xs, $
            /ys, $
            position=pos, $
            /normal, $
            xtitle='Frequency (MHz)', $
            ytitle='Intensity (DNs)'       
     
          oplot, bf, bi, color = i*col_scale, psym=1
          oplot, bf, bi, color = i*col_scale 
      
        ENDIF ELSE BEGIN	
  
          oplot, bf, bi, color = i*col_scale, psym=1
          oplot, bf, bi, color = i*col_scale
      
        ENDELSE		

      ENDFOR	
  device, /close
  
  set_plot, 'x'
  ;-----------------------------------;
  ;			 First intensity v time
  ;
  alldrift = fltarr(n_elements(indices))
  bidrift = fltarr(n_elements(indices))
  start_f = fltarr(n_elements(indices))
  
  setup_ps, 'figures/inten_time_data.ps'
      j=0
      FOR i=n_elements(indices)-2, 0, -1 DO BEGIN
    
        bt = btall[indices[i]+1: indices[i+1]-1]
        bf = bfall[indices[i]+1: indices[i+1]-1]
        bi = biall[indices[i]+1: indices[i+1]-1]
        tsec = anytim(bt[*], /utim) - anytim(bt[0], /utim)
    
          IF i eq n_elements(indices)-2 THEN BEGIN
    
       
            plot, tsec, bi, $
                psym=1, $
                xr=[0, 3], $
                yr=[0, 45], $
                /ys, $
                position=pos, $
                /normal, $
                xtitle='Time (s)', $
                ytitle='Intensity (DNs)'
        
            max_bi[j] = max(bi)  
            oplot, tsec, bi, color = i*col_scale, psym=1
            oplot, tsec, bi, color = i*col_scale
        
      
          ENDIF ELSE BEGIN
      
            max_bi[j] = max(bi) 
            oplot, tsec, bi, color = i*col_scale, psym=1
            oplot, tsec, bi, color = i*col_scale
          ENDELSE
          j = j+1
      ENDFOR  
  device, /close

  ;set_plot, 'x'
  
  setup_ps, 'figures/inten_time_fits.ps'
    j=0
    FOR i=n_elements(indices)-2, 0, -1 DO BEGIN  
    
        bt = btall[indices[i]+1: indices[i+1]-1]
        bf = bfall[indices[i]+1: indices[i+1]-1]
        bi = biall[indices[i]+1: indices[i+1]-1]
        tsec = anytim(bt[*], /utim) - anytim(bt[0], /utim)
    
        result = linfit(tsec, bi, yfit = yfit)
        start = [result[1]]
        if bf[0] lt 41.0 then intersect='31.0'
        if bf[0] gt 41.0 then intersect='42.0'
        fit = 'p[0]*x + '+	intersect
        result = mpfitexpr(fit, tsec , bi, err, yfit=yfit, start)
        alldrift[j] = result[0]
        result = linfit(tsec, bi, yfit = yfit)
        bidrift[j] = result[1];max(bi)
        start_f[j] = bf[0]
    
        IF i eq n_elements(indices)-2 THEN BEGIN
          plot, tsec, yfit, $
            xr=[0, 3], $
            yr=[0, 50], $
            position=pos, $
            /normal, $
            xtitle='Time (s)', $
            ytitle='Intensity (DNs)'
      
        ENDIF ELSE BEGIN	
    
          oplot, tsec, yfit, color = i*col_scale
      
        ENDELSE
  
      j = j+1
    ENDFOR	
  device, /close
 
  
  bidrift = bidrift[where(bidrift ne 0.0)]
  alldrift = alldrift[where(alldrift ne 0.0)]
  start_f = start_f[where(start_f ne 0.0)]
  
  
  loadct, 39
  ;-----------------------------------;
  ;			 Intensity v Velocity
  ;
  pos_k = where(bidrift ne 0.0)    
  pos_j = where(alldrift ne 0.0)    
  pos_i = [pos_j, pos_k]
  
  ;wset, 2
  setup_ps, 'figures/scatter_vels_maxi.ps'
    plotsym, 0, /fill 
    plot, [vels], [max_bi], $
        psym=8, $
        xr = [0, 0.6], $
        xtit = 'Beam velocity (c)', $
        ytit = 'Maximum intensity (DN)'
      
    set_line_color    
    for i =0, n_elements(vels)-1 do begin
      if start_f[i] lt 41 then color=3 else color=5    
      oplot, [vels[i]], [max_bi[i]], $
          psym=8, $
          color = color
    endfor    
  device, /close    
    
  loadct, 39
  ;-----------------------------------;
  ;			 Intensity v Displacement
  ;    
  
  setup_ps, 'figures/scatter_displ_maxi.ps'
    plot, [displ], [max_bi], $
        psym=8, $
        xr = [0.0, 0.35], $
        yr = [5.0, 50], $, 
        /ys, $
        xtit = 'Beam displacement (Rsun)', $
        ytit = 'Maximum intensity (DN)'
  
    set_line_color     
    for i =0, n_elements(displ)-1 do begin
      if start_f[i] lt 41 then color=3 else color=5    
        oplot, [displ[i]], [max_bi[i]], $
          psym=8, $
          color = color
    endfor
  device, /close
  
  loadct, 39
  ;-----------------------------------------;      
  ;     Frequency drift v Intensity drift
  ;    
  setup_ps, 'figures/scatter_fdrift_idrift.ps'
    plot, [alldrift[0]], [bidrift[0]], $
          psym=8, $
          xr = [0, 20], $
          yr = [-25, 5], $
          xtit = 'Drift rate (MHz s!U-1!N)', $
          ytit = 'dI/dt (DN s!U-1!N)'
        
    set_line_color         
    for i =0, n_elements(alldrift)-1 do begin
      if start_f[i] lt 41 then color=3 else color=5
      oplot, [drift[i]], [bidrift[i]], $
          psym=8, $
          color=color
          
    endfor      
   device, /close   
   set_plot, 'x'

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