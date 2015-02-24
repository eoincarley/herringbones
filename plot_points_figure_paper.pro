pro setup_ps, name
  
  set_plot,'ps'
  !p.font=0
  !p.charsize=1.5
  device, filename = name, $
          /color, $
          /helvetica, $
          /inches, $
          xsize=13, $
          ysize=6, $
          /encapsulate, $
          yoffset=5

end


pro plot_points_figure_paper

  cd, '~/Data/22Sep2011_event/herringbones'
  !p.charsize = 1.5
  loadct, 0
  ;-------------------------------------;
  ;			 Read and plot Spectrogram
  ;

  radio_spectro_fits_read,'BIR_20110922_104459_01.fit', data_raw, times, freq
  t1_index = closest(times, anytim(file2time('20110922_104700'),/utim))
  t2_index = closest(times, anytim(file2time('20110922_105340'),/utim))

  f1_index = closest(freq, 90.0)
  f2_index = closest(freq, 10.0)
  data_bs = constbacksub(data_raw, /auto)
  xtit = 'Start time: '+anytim(times[t1_index], /cc, /trun)+ ' (UT)'
  
  setup_ps, 'figures/burst_detections.eps'
      spectro_plot, bytscl(data_bs, -10, 40) , times, freq, $
          /ys, $
          ytitle = 'Frequency (MHz)', $
          yr = [freq[f1_index], freq[f2_index]], $
          xrange = [times[t1_index],times[t2_index]], $
          /xs, $
          xtitle = xtit, $
          title = 'eCallisto, Birr Castle, Ireland.', $
          xcolor = '100', $
          xticklen = -0.01, $
          yticklen = -0.01
  
      ; Plot twice because IDL's shitty fucking reverse_ct routine will NOT print the fucking 
      ; axes in black. Christ.
      loadct, 0    
      reverse_ct    
      spectro_plot, bytscl(data_bs, -10, 40) , times, freq, $
          /ys, $
          ytitle = '  ', $
          yr = [freq[f1_index], freq[f2_index]], $
          xrange = [times[t1_index],times[t2_index]], $
          /xs, $
          xtitle = ' ', $
          xticks=1, $
          yticks=1, $
          xtickname = [' ', ' ', ' ', ' ', ' ', ' ', ' ', ' '], $
          /noerase

      ;-------------------------------------;
      ;			     Read burst data
      ;
      readcol, 'bursts_ft_first_master_reverse.txt', btall0, bfall0, biall0, format = 'A,D,D'
      readcol, 'bursts_ft_second_master_reverse.txt', btall1, bfall1, biall1, format = 'A,D,D'
      readcol, 'bursts_ft_second_master_forward.txt', btall2, bfall2, biall2, format = 'A,D,D'
      readcol, 'bursts_ft_first_master_forward.txt', btall3, bfall3, biall3, format = 'A,D,D'
  
      ;-------------------------------------;
      ;			     Plot burst data
      ;
      plot_bursts, btall0, bfall0, biall0, color = 170 ; Green
      plot_bursts, btall1, bfall1, biall1, color = 220 ; Blue
      plot_bursts, btall2, bfall2, biall2, color = 80  ; Orange
      plot_bursts, btall3, bfall3, biall3, color = 50  ; Red

  device, /close
  set_plot, 'x'  

END


pro plot_bursts, btall, bfall, biall, color = color
  
  indices = where(btall eq '-')
  indices = [-1, indices]
  
  plotsym, 0, /fill
  
  FOR i=0, n_elements(indices)-2 DO BEGIN 
        bt = anytim(btall[indices[i]+1: indices[i+1]-1], /utim)  ; This selects one burst.
        bf = bfall[indices[i]+1: indices[i+1]-1]
        bi = biall[indices[i]+1: indices[i+1]-1]	
    
        set_line_color
        oplot, bt, bf, psym = 8, symsize=0.4, color=0
        ;oplot, bt, bf, psym = 8, symsize=0.3, color=0
        loadct, 74, /silent
        oplot, bt, bf, psym = 8, symsize=0.3, color = color         
  ENDFOR
  
END  
