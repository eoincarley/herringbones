pro plot_points_figure_paper

  cd, '~/Data/22Sep2011_event/herringbones'
  !p.charsize = 1.5
  loadct, 1
  reverse_ct
  window, 3, xs=1300, ys=800


  ;-------------------------------------;
  ;			 Read and plot Spectrogram
  ;

  radio_spectro_fits_read,'BIR_20110922_104459_01.fit', data_raw, times, freq
  t1_index = closest(times, anytim(file2time('20110922_104700'),/utim))
  t2_index = closest(times, anytim(file2time('20110922_105400'),/utim))

  f1_index = closest(freq, 90.0)
  f2_index = closest(freq, 10.0)
  data_bs = constbacksub(data_raw, /auto)
  xtit = 'Start time: '+anytim(times[t1_index], /cc)+ '(UT)'
  
  
  spectro_plot, bytscl(data_bs, -10, 40) , times, freq, $
      /ys, $
      ytitle = '!6Frequency [MHz]', $
      yr = [freq[f1_index], freq[f2_index]], $
      xrange = [times[t1_index],times[t2_index]], $
      /xs, $
      xtitle = xtit

 
  ;-------------------------------------;
  ;			     Read burst data
  ;
  readcol, 'bursts_ft_first_master_reverse.txt', btall0, bfall0, biall0, format = 'A,D,D'
  readcol, 'bursts_ft_second_master_reverse.txt', btall1, bfall1, biall1, format = 'A,D,D'
  readcol, 'bursts_ft_second_master_positive.txt', btall2, bfall2, biall2, format = 'A,D,D'
  readcol, 'bursts_ft_first_master_positive.txt', btall3, bfall3, biall3, format = 'A,D,D'
  
  btall = [btall0, '-', btall1, '-', btall2, '-', btall3]
  bfall = [bfall0, !Values.F_NAN, bfall1, !Values.F_NAN, bfall2, !Values.F_NAN, bfall3]
  biall = [biall0, !Values.F_NAN, biall1, !Values.F_NAN, biall2, !Values.F_NAN, biall3]

  
  indices = where(btall eq '-')
  indices = [-1, indices]
  
  ;-------------------------------------;
  ;			    Plot burst data
  ;
  plotsym, 0, /fill
  set_line_color
  FOR i=0, n_elements(indices)-2 DO BEGIN 
        bt = anytim(btall[indices[i]+1: indices[i+1]-1], /utim)  ; This selects one burst.
        bf = bfall[indices[i]+1: indices[i+1]-1]
        bi = biall[indices[i]+1: indices[i+1]-1]	
    
        oplot, bt, bf, psym = 8, symsize=0.6, color=0
           
  ENDFOR
  
  FOR i=0, n_elements(indices)-2 DO BEGIN 
        bt = anytim(btall[indices[i]+1: indices[i+1]-1], /utim)  ; This selects one burst.
        bf = bfall[indices[i]+1: indices[i+1]-1]
        bi = biall[indices[i]+1: indices[i+1]-1]	
    
        oplot, bt, bf, psym = 8, symsize=0.4, color=10
           
  ENDFOR
  
  btall = [btall0];, '-', btall1, '-'] ;, btall3]
  bfall = [bfall0];, !Values.F_NAN, bfall1] ;, !Values.F_NAN, bfall3]
  biall = [biall0];, !Values.F_NAN, biall1] ;, !Values.F_NAN, biall3]

  
  indices = where(btall eq '-')
  indices = [-1, indices]
  
  ;-------------------------------------;
  ;			    Plot burst data
  ;
  plotsym, 0, /fill
  set_line_color
  FOR i=0, n_elements(indices)-2 DO BEGIN 
        bt = anytim(btall[indices[i]+1: indices[i+1]-1], /utim)  ; This selects one burst.
        bf = bfall[indices[i]+1: indices[i+1]-1]
        bi = biall[indices[i]+1: indices[i+1]-1]	
    
        oplot, bt, bf, psym = 8, symsize=0.6, color=0
           
  ENDFOR
  
  FOR i=0, n_elements(indices)-2 DO BEGIN 
        bt = anytim(btall[indices[i]+1: indices[i+1]-1], /utim)  ; This selects one burst.
        bf = bfall[indices[i]+1: indices[i+1]-1]
        bi = biall[indices[i]+1: indices[i+1]-1]	
    
        oplot, bt, bf, psym = 8, symsize=0.4, color=4
           
  ENDFOR
  
  

END