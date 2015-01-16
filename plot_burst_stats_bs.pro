pro plot_burst_stats_bs

  ;12-Sep-2013 plot burst stats derived from Hough analysis

  ;Last time this was run need to fix ranges on plot.
  !p.thick=1
  !p.charsize=1.5
  loadct, 39
  pos = [0.13, 0.1, 0.95, 0.95]
  col_scale = 2.0
  dimen = 600
  xpos = 900
  ypos = 50
  window, 0, xs=dimen, ys=dimen, xpos = xpos, ypos = ypos
  window, 1, xs=dimen, ys=dimen, xpos = xpos + 1.0*dimen, ypos = ypos
  window, 2, xs=dimen, ys=dimen, xpos = xpos, ypos = ypos + 1.0*dimen
  window, 3, xs=dimen, ys=dimen, xpos = xpos + 1.0*dimen, ypos = ypos + 1.0*dimen
  ;window, 4, xs=dimen, ys=dimen, xpos = 10, ypos = ypos + 4.0*dimen
  
  folder = '/Users/ecarley/Data/hough_hbones/20110922'
  cd, folder
  ;-------------------------------------
  ;			 First do raw data
  ;
  readcol,'bursts_bs_hough.txt', btall, bfall, biall, format = 'A,D,D'
  indices = where(btall eq '-')
  indices = [-1, indices]


  ;-------------------------------------;
  ;
  ;		 Intensity v Frequency
  ;

      ;set_plot,'ps'
      ;!p.font=0
      ;!p.charsize=1.5
      ;device, filename='bs_bfbi.eps', $
              ;/color, $
              ;/helvetica, $
              ;/inches, $
              ;xsize=7, $
              ;ysize=7, $
              ;/encapsulate, $
              ;yoffset=5

  FOR i=0, n_elements(indices)-2 DO BEGIN

    bt = btall[indices[i]+1: indices[i+1]-1]
    bf = bfall[indices[i]+1: indices[i+1]-1]
    bi = biall[indices[i]+1: indices[i+1]-1]
    
    IF i eq 0 THEN BEGIN
      wset, 0
      plot, bf, bi, $
          /xs, $
          /ys, $
          xr=[45, 65], $
          yr=[0, 55], $
          psym=1, $
          xtitle='Frequency (MHz)', $
          ytitle='Intensity (DN)', $
          position=pos, $
          /normal, $
          title='Back-sub Herringbones 22-Sep-2011'
          
      oplot, bf, bi, color=(250), psym=1
      oplot, bf, bi, color=(250)
    ENDIF ELSE BEGIN	
      oplot, bf, bi, color = (245 - i*col_scale), psym=1
      oplot, bf, bi, color=(245 - i*col_scale)
    ENDELSE	
    
  ENDFOR	
  


  ;-------------------------------------;
  ;
  ;		 Frequency v Time
  ;
  FOR i=n_elements(indices)-2, 0, -1 DO BEGIN ;backwards to plot the shortest bursts first.
  

    bt = btall[indices[i]+1: indices[i+1]-1]
    bf = bfall[indices[i]+1: indices[i+1]-1]
    bi = biall[indices[i]+1: indices[i+1]-1]	
    tsec = anytim(bt[*], /utim) - anytim(bt[0], /utim)
    IF i eq n_elements(indices)-2 THEN BEGIN
    
      wset, 1
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
    ENDIF ELSE BEGIN	
      oplot, tsec, bf, color=(245 - i*col_scale)
    ENDELSE	

  ENDFOR



  ;-------------------------------------;
  ;
  ;		   Time v Intensity
  ;
  FOR i=n_elements(indices)-2, 0, -1 DO BEGIN

    bt = btall[indices[i]+1: indices[i+1]-1]
    bf = bfall[indices[i]+1: indices[i+1]-1]
    bi = biall[indices[i]+1: indices[i+1]-1]
    tsec = anytim(bt[*], /utim) - anytim(bt[0], /utim)
    
    IF i eq n_elements(indices)-2 THEN BEGIN
    
      wset, 2
      plot, tsec, bi, $
        psym=1, $
        xr=[0, 3], $
        yr=[0, 50], $
        position=pos, $
        /normal, $
        title='Back-sub Herringbones 22-Sep-2011', $
        xtitle='Time (s)', $
        ytitle='Intensity'
      
      oplot, tsec, bi, color=(250), psym=1
      oplot, tsec, bi, color=(250)
    ENDIF ELSE BEGIN	
      oplot, tsec, bi, color = (245 - i*col_scale), psym=1
      oplot, tsec, bi, color = (245 - i*col_scale)
    ENDELSE		
  
  ENDFOR	



  ;-------------------------------------;
  ;
  ;		   Drift v DI/DT
  ;
  alldrift = [0]
  allbidrift = [0]
  PLOTSYM, 0, /fill               ;Plotting symbol is a circle

  FOR i=0, n_elements(indices)-2 DO BEGIN

    bt = btall[indices[i]+1: indices[i+1]-1]
    bf = bfall[indices[i]+1: indices[i+1]-1]
    bi = biall[indices[i]+1: indices[i+1]-1]
    tsec = anytim(bt[*], /utim) - anytim(bt[0], /utim)
   
    result = linfit(tsec, bf)
    drift = result[1]
    alldrift = [alldrift, drift]

    result = linfit(tsec, bi)
    bidrift = result[1];max(bi)
    allbidrift = [allbidrift, bidrift]
    loadct, 39
    print, bidrift
    
    IF i eq 0 THEN BEGIN
      wset, 3
      plot, [0,drift], [0, bidrift], $
        /noerase, $
        psym=8, $
        xr=[5, 15], $
        yr=[0, -15], $
        position=pos, $
        /normal, $
        title='Back-sub Herringbones 22-Sep-2011', $
        xtitle='Drift Rate (MHz s!U-1!N)', $
        ytitle='Rate of change of intensity DI/Dt'
      
      oplot, [0,drift], [0, bidrift], psym=8, color=5
    ENDIF ELSE BEGIN	
      oplot, 	[0,drift], [0, bidrift], psym=8, color = (245 - i*col_scale)
    ENDELSE
  
  ENDFOR	

stop

  alldrift = alldrift[ 1: n_elements(alldrift)-1 ]
  remove_nans, alldrift, alldrift
  allbidrift = allbidrift[ 1:n_elements(allbidrift)-1 ]
  remove_nans, allbidrift, allbidrift

  result = correlate(alldrift, allbidrift)
  print,'Correlation coefficient: '+string(result)

  ;-------------------------------------;
  ;
  ;		   Drift v Max I
  ;
  PLOTSYM, 0, /fill               ;Plotting symbol is a circle

  FOR i=0, n_elements(indices)-2 DO BEGIN

    bt = btall[indices[i]+1: indices[i+1]-1]
    bf = bfall[indices[i]+1: indices[i+1]-1]
    bi = biall[indices[i]+1: indices[i+1]-1]
    tsec = anytim(bt[*], /utim) - anytim(bt[0], /utim)
    result = linfit(tsec, bf)
    drift = result[1]
  
    bimax = max(bi)
    set_line_color
    print, bidrift
    IF i eq 0 THEN BEGIN
      plot, [0,drift], [0, bimax], /noerase, psym=8, color=0, xr=[5,15], yr=[20, 50], position=pos,$
      /normal, title='Back-sub Herringbones 22-Sep-2011', $
      xtitle='Drift Rate (MHz s!U-1!N)', ytitle='Max Intensity'
      oplot, [0,drift], [0, bimax], psym=8, color=5
    ENDIF ELSE BEGIN	
      oplot, 	[0,drift], [0, bimax], psym=8, color=5
    ENDELSE
  
  ENDFOR	
  axis, xaxis=0, xtickname=[' ', ' ', ' ', ' ', ' ', ' '], xr=[5, 15], xthick=3
  axis, xaxis=1, xtickname=[' ', ' ', ' ', ' ', ' ', ' '], xr=[5, 15], xthick=3
  axis, yaxis=0, ytickname=[' ', ' ', ' ', ' ', ' ', ' '], yr=[20, 50], ythick=3
  axis, yaxis=1, ytickname=[' ', ' ', ' ', ' ', ' ', ' '], yr=[20, 50], ythick=3



END