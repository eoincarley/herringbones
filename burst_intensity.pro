pro burst_intensity

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
  
  ;------------------------;
  ;			 Read data
  ;
  readcol, 'bursts_ft_first_master_reverse.txt', btall0, bfall0, biall0, format = 'A,D,D'
  readcol, 'bursts_ft_second_master_reverse.txt', btall1, bfall1, biall1, format = 'A,D,D'
  readcol, 'bursts_ft_second_master_forward.txt', btall2, bfall2, biall2, format = 'A,D,D'
  readcol, 'bursts_ft_first_master_forward.txt', btall3, bfall3, biall3, format = 'A,D,D'
 
  
  btall = [btall0, '-', btall1, '-', btall2, '-', btall3]
  bfall = [bfall0, !Values.F_NAN, bfall1, !Values.F_NAN, bfall2, !Values.F_NAN, bfall3]
  biall = [biall0, !Values.F_NAN, biall1, !Values.F_NAN, biall2, !Values.F_NAN, biall3]
  
  indices = where(btall eq '-')
  indices = [-1, indices]
  
  restore, 'beam_kins.sav', /verb
  
  drift = beam_kins[0]
  vels = beam_kins[1]
  displ = beam_kins[2]
  flength =  beam_kins[4]
  figs_folder = beam_kins[5]
  max_bi = fltarr(n_elements(vels))
  ;-----------------------------------;
  ;			 First intensity v time
  ;
  setup_ps, 'figures/'+figs_folder+'/freq_inten_data.ps'
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
  
  ;-----------------------------------;
  ;			 First intensity v time
  ;
  alldrift = fltarr(n_elements(indices))
  bidrift = fltarr(n_elements(indices))
  start_f = fltarr(n_elements(indices))
  end_f = fltarr(n_elements(indices))
  
  setup_ps, 'figures/'+figs_folder+'/inten_time_data.ps'
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
        
            max_bi[j] = max(bi)  - min(bi)

            oplot, tsec, bi, color = i*col_scale, psym=1
            oplot, tsec, bi, color = i*col_scale
        
      
          ENDIF ELSE BEGIN
      
            max_bi[j] = max(bi) - min(bi) 

            oplot, tsec, bi, color = i*col_scale, psym=1
            oplot, tsec, bi, color = i*col_scale
          ENDELSE
          j = j+1
      ENDFOR  
  device, /close

  ;set_plot, 'x'
  
  setup_ps, 'figures/'+figs_folder+'/inten_time_fits.ps'
    j=0
    FOR i=n_elements(indices)-2, 0, -1 DO BEGIN  
    
        bt = btall[indices[i]+1: indices[i+1]-1]
        bf = bfall[indices[i]+1: indices[i+1]-1]
        bi = biall[indices[i]+1: indices[i+1]-1]
        tsec = anytim(bt[*], /utim) - anytim(bt[0], /utim)
        
        result = linfit(tsec, bi, yfit = yfit)
        bidrift[j] = result[1]  
        start_f[j] = bf[0]
        end_f[j] = bf[n_elements(bf)-1]
       
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
  start_f = start_f[where(start_f ne 0.0)]
  end_f = end_f[where(end_f ne 0.0)]
  
  ipos = where(drift gt 0.0)
  
  vels = vels[ipos]
  max_bi = max_bi[ipos]
  displ = displ[ipos]
  drift = drift[ipos]
  bidrift = bidrift[ipos]
  start_f = start_f[ipos]
  end_f = end_f[ipos]
  flength = flength[ipos]
  
  loadct, 39
  ;-----------------------------------;
  ;			 Intensity v Velocity
  ;
  
  ;wset, 2
  setup_ps, 'figures/'+figs_folder+'/scatter_vels_maxi.ps'
    plotsym, 0, /fill 
    plot, [vels], [max_bi], $
        psym=8, $
        xr = [0, 0.6], $
        xtit = 'Beam velocity (c)', $
        ytit = 'Maximum intensity (DN)'
      
    set_line_color  

    for i =0, n_elements(vels)-1 do begin
      if round(start_f[i]) eq 43 then color=3 
      if round(start_f[i]) eq 31 then color=5    
      if round(start_f[i]) eq 32 then color=4 
      if round(start_f[i]) eq 42 then color=2    

      oplot, [vels[i]], [max_bi[i]], $
          psym=8, $
          color = color
    endfor    
  device, /close    
    
  loadct, 39
  
  ;-----------------------------------;
  ;			 Intensity v Displacement
  ;    
  
  setup_ps, 'figures/'+figs_folder+'/scatter_displ_maxi.ps'
    plot, [flength[0]], [max_bi[0]], $
        psym=8, $
        xr = [0.0, 40], $
        yr = [0.0, 45], $, 
        /ys, $
        xtit = 'Frequency range (MHz)', $
        ytit = 'Maximum intensity (DN)'
  
    set_line_color     
    for i =0, n_elements(displ)-1 do begin
      if round(start_f[i]) eq 43 then color=3 
      if round(start_f[i]) eq 31 then color=5    
      if round(start_f[i]) eq 32 then color=4 
      if round(start_f[i]) eq 42 then color=2    
        oplot, [flength[i]], [max_bi[i]], $
          psym=8, $
          color = color
    endfor
    PRINT, 'Intensity and displacement correlation coefficient: '+ string(CORRELATE(displ, max_bi))

  device, /close
  
  loadct, 39
  ;-----------------------------------------;      
  ;    Frequency drift v Intensity drift
  ;    
  setup_ps, 'figures/'+figs_folder+'/scatter_fdrift_idrift.ps'
    plot, [drift[0]], [bidrift[0]], $
          psym=8, $
          xr = [-10, 20], $
          yr = [-25, 5], $
          /ys, $
          xtit = 'Drift rate (MHz s!U-1!N)', $
          ytit = 'dI/dt (DN s!U-1!N)'
        
    set_line_color         
    for i =0, n_elements(drift)-1 do begin
      if round(start_f[i]) eq 43 then color=3 
      if round(start_f[i]) eq 31 then color=5    
      if round(start_f[i]) eq 32 then color=4 
      if round(start_f[i]) eq 42 then color=2 
      oplot, [drift[i]], [bidrift[i]], $
          psym=8, $
          color=color
          
    endfor

    PRINT, 'Intensity drift and frequency drift correlation coefficient: '+ string(CORRELATE(drift, bidrift))      

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