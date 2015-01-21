pro burst_intensity

  !p.thick = 2 
  !p.charsize = 1.5
  !p.charthick = 1.5
  loadct, 39
  reverse_ct
  pos = [0.13, 0.1, 0.95, 0.95]
  col_scale = 3.0
  dimen = 600
  xpos = 500
  ypos = 50
  window, 0, xs=dimen, ys=dimen, xpos = xpos, ypos = ypos
  window, 1, xs=dimen, ys=dimen, xpos = xpos + 1.0*dimen, ypos = ypos
  window, 2, xs=dimen, ys=dimen, xpos = xpos, ypos = ypos + 1.0*dimen
  window, 3, xs=dimen, ys=dimen, xpos = xpos + 1.0*dimen, ypos = ypos + 1.0*dimen
  window, 4, xs=dimen, ys=dimen, xpos = xpos + 2.0*dimen, ypos = ypos 
  window, 5, xs=dimen, ys=dimen, xpos = xpos + 2.0*dimen, ypos = ypos + 1.0*dimen
  folder = '~/Data/22Sep2011_event/herringbones'
  cd, folder
  
  ;------------------------;
  ;			 Read data
  ;
  readcol,'bursts_bs_hough_first1.txt', btall1, bfall1, biall1, format = 'A,D,D'
  readcol,'bursts_bs_hough_first2.txt', btall2, bfall2, biall2, format = 'A,D,D'
  readcol,'bursts_bs_hough_second.txt', btall3, bfall3, biall3, format = 'A,D,D'
  
  btall = [btall1, '-', btall2, '-', btall3]
  bfall = [bfall1, !Values.F_NAN, bfall2, !Values.F_NAN, bfall3]
  biall = [biall1, !Values.F_NAN, biall2, !Values.F_NAN, biall3]
  
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

  FOR i=n_elements(indices)-2, 0, -1 DO BEGIN

    bt = btall[indices[i]+1: indices[i+1]-1]
    bf = bfall[indices[i]+1: indices[i+1]-1]
    bi = biall[indices[i]+1: indices[i+1]-1]
    tsec = anytim(bt[*], /utim) - anytim(bt[0], /utim)
    
    IF i eq n_elements(indices)-2 THEN BEGIN
    
      wset, 5
      plot, bf, bi, $
        psym=1, $
        xr=[45, 80], $
        yr=[0, 50], $
        position=pos, $
        /normal, $
        title='Back-sub Herringbones 22-Sep-2011', $
        xtitle='Frequency (MHz)', $
        ytitle='Intensity'       
     
      oplot, bf, bi, color=(250), psym=1
      oplot, bf, bi, color=(250) 
      
    ENDIF ELSE BEGIN	
    
      wset,5
  
      oplot, bf, bi, color = (245 - i*col_scale), psym=1
      oplot, bf, bi, color = (245 - i*col_scale)
      
    ENDELSE		

  ENDFOR	
  ;-----------------------------------;
  ;			 First intensity v time
  ;
  alldrift = fltarr(n_elements(indices))
  bidrift = fltarr(n_elements(indices))
  start_f = fltarr(n_elements(indices))
  
  j=0
  FOR i=n_elements(indices)-2, 0, -1 DO BEGIN

    bt = btall[indices[i]+1: indices[i+1]-1]
    bf = bfall[indices[i]+1: indices[i+1]-1]
    bi = biall[indices[i]+1: indices[i+1]-1]
    tsec = anytim(bt[*], /utim) - anytim(bt[0], /utim)
    
    result = linfit(tsec, bf, yfit = yfit)
    
    start = [result[1]]
    if bf[0] lt 46.00 then intersect='33.25'
    if bf[0] eq 46.25 then intersect='46.25'
    
    
	  fit = 'p[0]*x + '+	intersect
	  result = mpfitexpr(fit, tsec , bf, err, yfit=yfit, start)
	  
    alldrift[j] = result[0]
    
    result = linfit(tsec, bi, yfit = yfit)
    bidrift[j] = result[1];max(bi)
    start_f[j] = bf[0]
    
    IF i eq n_elements(indices)-2 THEN BEGIN
    
      wset, 0
      plot, tsec, bi, $
        psym=1, $
        xr=[0, 3], $
        yr=[0, 50], $
        position=pos, $
        /normal, $
        title='Back-sub Herringbones 22-Sep-2011', $
        xtitle='Time (s)', $
        ytitle='Intensity'
        
      max_bi[j] = max(bi)  
      oplot, tsec, bi, color=(250), psym=1
      oplot, tsec, bi, color=(250)
      
      wset, 1
      plot, tsec, yfit, $
        xr=[0, 3], $
        yr=[0, 50], $
        position=pos, $
        /normal, $
        title='Back-sub Herringbones 22-Sep-2011', $
        xtitle='Time (s)', $
        ytitle='Intensity'
      
    ENDIF ELSE BEGIN	
    
      wset,0
      max_bi[j] = max(bi) 
      oplot, tsec, bi, color = (245 - i*col_scale), psym=1
      oplot, tsec, bi, color = (245 - i*col_scale)
      
      wset, 1
      oplot, tsec, yfit, color = (245 - i*col_scale)
    ENDELSE		
    j = j+1
  ENDFOR	
  
  loadct, 39
  reverse_ct
  ;-----------------------------------;
  ;			 Intensity v Velocity
  ;
  pos_k = where(bidrift ne 0.0)    
  pos_j = where(alldrift ne 0.0)    
  pos_i = [pos_j, pos_k]
  
  wset, 2
  plotsym, 0, /fill 
  plot, [0, vels], [0, max_bi], $
      psym=8, $
      xr = [0, 0.6], $
      xtit = 'Beam velocity (c)', $
      ytit = 'Maximum intensity (DN)'
      
  set_line_color    
  for i =0, n_elements(vels)-1 do begin
    if start_f[i] lt 45 then color=3 else color=5    
      oplot, [0, vels[i]], [0, max_bi[i]], $
        psym=8, $
        color = color
  endfor    
      
    
  loadct, 39
  reverse_ct
  ;-----------------------------------;
  ;			 Intensity v Displacement
  ;    
  wset, 3
  plot, [0, displ], [0, max_bi], $
      psym=8, $
      xr = [0.0, 0.2], $
      xtit = 'Beam displacement (Rsun)', $
      ytit = 'Maximum intensity (DN)', $
      color = 255
  
  set_line_color     
  for i =0, n_elements(displ)-1 do begin
    if start_f[i] lt 45 then color=3 else color=5    
      oplot, [0, displ[i]], [0, max_bi[i]], $
        psym=8, $
        color = color
  endfor
  
  loadct, 39
  reverse_ct
  ;-----------------------------------------;      
  ;     Frequency drift v Intensity drift
  ;    
  wset, 4
  ;pos_i = where(bidrift
  plot, [0, alldrift[0]], [0, bidrift[0]], $
        psym=8, $
        yr = [-40, 10], $
        xr = [0, 20], $
        xtit = 'Drift (MHz/s)', $
        ytit = 'di/dt (DN)', $
        color = 255
        
  set_line_color         
  for i =0, n_elements(alldrift)-1 do begin
    if start_f[i] lt 45 then color=3 else color=5
    oplot, [0, alldrift[i]], [0, bidrift[i]], $
        psym=8, $
        color = color
  endfor      
       
stop

END