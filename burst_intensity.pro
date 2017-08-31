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
  folder = '~/Data/2011_sep_22/herringbones'
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
  rev1f0 = bfall0[0]
  rev2f0 = bfall1[0]
  for2f0 = bfall2[0]
  for1f0 = bfall3[0]
  
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
  ;			 First intensity v frequency
  ;
  setup_ps, 'figures/'+figs_folder+'/freq_inten_data.eps'
      FOR i=n_elements(indices)-2, 0, -1 DO BEGIN

        bt = btall[indices[i]+1: indices[i+1]-1]
        bf = bfall[indices[i]+1: indices[i+1]-1]
        bi = biall[indices[i]+1: indices[i+1]-1]
        tsec = anytim(bt[*], /utim) - anytim(bt[0], /utim)
    
        IF i eq n_elements(indices)-2 THEN BEGIN
    
          ;wset, 5
          plot, bf, bi, $
            psym=1, $
            xr=[10, 80], $
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
  
  setup_ps, 'figures/'+figs_folder+'/inten_time_data.eps'
      j=0
      FOR i=n_elements(indices)-2, 0, -1 DO BEGIN
    
        bt = btall[indices[i]+1: indices[i+1]-1]
        bf = bfall[indices[i]+1: indices[i+1]-1]
        bi = biall[indices[i]+1: indices[i+1]-1]
        tsec = anytim(bt[*], /utim) - anytim(bt[0], /utim)
    
          IF i eq n_elements(indices)-2 THEN BEGIN
    
       
            plot, tsec, bi, $
                psym=1, $
                xr=[0, 6], $
                yr=[0, 45], $
                /ys, $
                position=pos, $
                /normal, $
                xtitle='Time (s)', $
                ytitle='Intensity (DNs)'
        
            max_bi[j] = max(bi) - min(bi)

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
  
  setup_ps, 'figures/'+figs_folder+'/inten_time_fits.eps'
    j=0
    FOR i=n_elements(indices)-2, 0, -1 DO BEGIN  
    
        bt = btall[indices[i]+1: indices[i+1]-1]
        bf = bfall[indices[i]+1: indices[i+1]-1]
        bi = biall[indices[i]+1: indices[i+1]-1]
        tsec = anytim(bt[*], /utim) - anytim(bt[0], /utim)
        
        if bf[0] eq rev2f0 then begin         ;Second reverse
          index_max = where(bi eq max(bi))
          index_min = where(bi eq min(bi))
          ;index_gt2 = where(tsec gt 2.2)
          bi = bi[index_max : n_elements(bi)-1]
          tsec = tsec[index_max : n_elements(tsec)-1]
        endif  

        result = linfit(tsec, bi, yfit = yfit)
        bidrift[j] = result[1]  
        start_f[j] = bf[0]
        end_f[j] = bf[n_elements(bf)-1]
       
        IF i eq n_elements(indices)-2 THEN BEGIN
          plot, tsec, yfit, $
            xr=[0, 6], $
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
  
  ipos = where(drift gt 0.0)	; Indices of reverse drift only
  
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
  setup_ps, 'figures/'+figs_folder+'/scatter_vels_maxi.eps'
    plotsym, 0, /fill 
    plot, [vels], [max_bi], $
        psym=8, $
        symsize=1.1, $
        xr = [0.05, 0.25], $
        /ys, $
        /xs, $
        xtit = 'Beam velocity (c)', $
        ytit = 'Maximum intensity (DN)'
      
    loadct, 74  
    for i =0, n_elements(vels)-1 do begin
      ;Green = 170,  Blue = 220,  Orange = 80,  Red = 50
      if start_f[i] eq rev2f0 then color = 220 		;Second reverse
      if start_f[i] eq for2f0 then color = 80    	;Second forward
      if start_f[i] eq for1f0 then color = 50    	;First forward 
      if start_f[i] eq rev1f0 then color = 170 		;First reverse

      oplot, [vels[i]], [max_bi[i]], $
          psym=8, $
          color = color
    endfor 
    set_line_color   
    idft = linfit(vels, max_bi, yfit = yfit)
    oplot, vels, yfit, color=0   
    cc = string(CORRELATE(vels, max_bi), format='(f5.2)')
    PRINT, 'Velocity and intensity correlation coefficient:'+ cc     
    legend, ['Total CC: '+cc], box=0, pos = [0.19, 0.92], /normal, $
            charsize=1.2, thick=3
    
  device, /close    

  ;-----------------------------------;
  ;			 Intensity v Displacement
  ;    
  loadct, 0
  setup_ps, 'figures/'+figs_folder+'/scatter_displ_maxi.eps'
    plot, [displ], [max_bi], $
        psym=8, $
        ;xr = [0.0, 40], $
        yr = [0.0, 45], $, 
        /ys, $
        color = 0, $
        symsize = 1.1, $
        xtit = 'Beam displacement (R  )', $
        ytit = 'Maximum intensity (DN)'
  
    loadct, 74    
    for i =0, n_elements(displ)-1 do begin
      ;Green = 170,  Blue = 220,  Orange = 80,  Red = 50
      if start_f[i] eq rev2f0 then color = 220    ;Second reverse
      if start_f[i] eq for2f0 then color = 80     ;Second forward
      if start_f[i] eq for1f0 then color = 50     ;First forward 
      if start_f[i] eq rev1f0 then color = 170    ;First reverse
      
      if end_f[i] gt 88.0 then begin
      	void = execute('PLOTSYM, 7, 2.5, thick=3, /FILL') 
      endif else begin
      	void = execute('PLOTSYM, 0, /FILL')
      endelse
          
        oplot, [displ[i]], [max_bi[i]], $
          psym = 8, $
          color = color
    endfor
    
    set_line_color   
    idft = linfit(displ, max_bi, yfit = yfit)
    oplot, displ, yfit, color=0
    cc_tot = string(CORRELATE(displ, max_bi), format='(f5.2)')
    ;PRINT, 'Intensity and displacement correlation coefficient: '+ cc     

    
     ; Now correlated separately
    indeces43 = where(start_f eq rev2f0)
    indeces32 = where(start_f eq rev1f0)
    loadct, 74
    idft = linfit(displ[indeces43], max_bi[indeces43], yfit = yfit)
    oplot, displ[indeces43], yfit, color=220   
    cc43 = string(CORRELATE(displ[indeces43], max_bi[indeces43]), format='(f5.2)')
  
    idft = linfit(displ[indeces32], max_bi[indeces32], yfit = yfit)
    set_line_color
    oplot, displ[indeces32], yfit, color=0, thick=3.0
    loadct, 74
    oplot, displ[indeces32], yfit, color=140   
    cc32 = string(CORRELATE(displ[indeces32], max_bi[indeces32]), format='(f5.2)')
    
    loadct, 39
    legend, ['Total CC: '+cc_tot, 'Reverse 43 MHz CC: '+cc43, 'Reverse 32 MHz CC: '+cc32], $
            color = [0, 80, 180], linestyle = [0,0,0], box=0, pos = [0.18, 0.92], /normal, $
            charsize=1.2, thick=3
    
  device, /close
  

  ;-----------------------------------------;      
  ;    Frequency drift v Intensity drift
  ;    
  loadct,39
  setup_ps, 'figures/'+figs_folder+'/scatter_fdrift_idrift.eps'
    plot, [abs(drift)], [bidrift], $
          psym=8, $
          xr = [0, 15], $
          yr = [-25, 5], $
          symsize=1.1, $
          /ys, $
          xtit = 'Drift rate (MHz s!U-1!N)', $
          ytit = 'dI/dt (DN s!U-1!N)'
        
    loadct, 74        
    for i =0, n_elements(drift)-1 do begin
      ;Green = 170,  Blue = 220,  Orange = 80,  Red = 50
      if start_f[i] eq rev2f0 then color = 220    ;Second reverse
      if start_f[i] eq for2f0 then color = 80     ;Second forward
      if start_f[i] eq for1f0 then color = 50     ;First forward 
      if start_f[i] eq rev1f0 then color = 170    ;First reverse
          
      oplot, [abs(drift[i])], [bidrift[i]], $
          psym=8, $
          color=color
          
    endfor
    zero = drift
    zero[*] = 0.0
    oplot, findgen(21), zero, color=0, linestyle = 3
    
    set_line_color   
    idft = linfit(drift, bidrift, yfit = yfit)
    oplot, drift, yfit, color=0
    cc = string(CORRELATE(drift, bidrift), format='(f5.2)')      

    PRINT, 'Intensity drift and frequency drift correlation coefficient: '+ cc
    legend, ['Correlation Coefficient: '+cc], box=0, pos = [0.45, 0.92], /normal
    
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