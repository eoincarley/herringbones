pro burst_intensity

  !p.thick=1
  !p.charsize=1.5
  loadct, 39
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
  readcol,'bursts_bs_hough_second.txt', btall2, bfall2, biall2, format = 'A,D,D'
  
  btall = btall2 ;[btall1, '-', btall2]
  bfall = bfall2 ;[bfall1, !Values.F_NAN, bfall2]
  biall = biall2 ;[biall1, !Values.F_NAN, biall2]
  
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
  
  ;-----------------------------------;
  ;			 Intensity v velocity
  ;
  pos_k = where(bidrift ne 0.0)    
  pos_j = where(alldrift ne 0.0)    
  pos_i = [pos_j, pos_k]
  
  wset, 2
  plot, vels[pos_i], max_bi[pos_i], $
      psym=4, $
      xr = [0, 0.6], $
      xtit = 'Beam velocity (c)', $
      ytit = 'Maximum intensity (DN)'
      
  wset, 3
  plot, displ[pos_i], max_bi[pos_i], $
      psym=4, $
      xr = [0.0, 0.4], $
      xtit = 'Beam displacement (Rsun)', $
      ytit = 'Maximum intensity (DN)'
      
  wset, 4
  ;pos_i = where(bidrift
  plot, alldrift[pos_i], bidrift[pos_i], $
      psym=4, $
      yr = [20, -60], $
      xr = [0, 25], $
      xtit = 'Drift (MHz/s)', $
      ytit = 'di/dt (DN)'   
      
       
stop

END