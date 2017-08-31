pro read_goes, goes_array
  
  file = '~/Data/2011_sep_22/20110922_Gp_xr_1m.txt'

 readcol, file, y, m, d, hhmm, mjd, sod, short_channel, long_channel;, $
          format = 'A, A, A, A, A, A, L, L'
 

  ;-------- Time in correct format --------
  time  = strarr(n_elements(y))
  time[*] = string(y[*], format='(I04)') + string(m[*], format='(I02)') $
  		  + string(d[*], format='(I02)') + '_' + string(hhmm[*], format='(I04)')
  time = anytim(file2time(time), /utim)

  ;------- Build data array --------------
 
  goes_array = dblarr(3, n_elements(y))
  goes_array[0,*] = time
  goes_array[1,*] = long_channel
  goes_array[2,*] = short_channel
  
  
END  

pro plot_callisto_goes

  ;------------------------------;
  ;			   READ THE DATA			   ;
  ;


  cd, '~/Data/2011_sep_22/herringbones'
  
  xstart = anytim(file2time('20110922_100000'),/utim)
  xend = anytim(file2time('20110922_120000'),/utim)
  xstart_display = anytim(file2time('20110922_100000'),/yoh, /trun)

  files = findfile('*01.fit')
  radio_spectro_fits_read, files[0], low1data, l1time, low_freq
  radio_spectro_fits_read, files[1], low2data, l2time, low_freq


  low_data = [low1data, low2data]
  low_times = [l1time, l2time]

  ;------------------------------------------;
  ;			   PLOT THE DATA			   ;
  ;
  x0 = 0.12
  x1 = 0.95
  ;angstrom = '!6!sA!r!u!9 %!6!n' 
  !p.charsize = 1
  !p.font = 0
  set_plot,'ps'
  device, filename = 'goes_callisto_20110922.eps', $
          /color, $
          bits=8, $
          /inches, $
          /helvetica, $
          /encapsulate, $
          xsize=7, $
          ysize=10

      read_goes, goes_data
      
      set_line_color
      utplot, goes_data[0,*], goes_data[1,*], $
          /ylog, $
          xr = [xstart,xend], $
          /xs, $
          thick = 4, $
          yrange = [1e-7,1e-3], $
          xtitle='Start Time (2011-Sep-22 10:00:00 UT)', $
          position = [x0, 0.7, x1, 0.99], $
          color = 3
      
     
      xyouts, 0.03, 0.82, 'Watts m!U-2!N', /normal, orientation=90
      axis,yaxis=1,ytickname=[' ','C','M','X',' ']
      
      oplot, goes_data[0,*], goes_data[2,*], $
            color=5, $
            thick=4

      loadct,0
      i_gen = dindgen(1001)*((1.0e-3 - 1.0e-7)/1000)+1.0e-7
      tcon = anytim(file2time('20110922_104700'),/utim)
      plots, tcon, i_gen, linestyle=2
      i_gen = dindgen(1001)*((1.0e-3 - 1.0e-7)/1000)+1.0e-7
      tcon = anytim(file2time('20110922_105330'),/utim)
      set_line_color
      plots, tcon, i_gen, linestyle=2

      plots, [x0, 0.445], [0.63, 0.7], /normal, linestyle=2
      plots, [0.95, 0.485], [0.63, 0.7], /normal, linestyle=2
      
      t1 = anytim(file2time('20110922_10000'), /utim)
      t2 = anytim(file2time('20110922_12000'), /utim)
      indices = where(goes_data[0,*] gt t1 and goes_data[0,*] lt t2)
      plots, goes_data[0, indices], 1e-6, thick=1, color=0
      plots, goes_data[0, indices], 1e-5, thick=1, color=0
      plots, goes_data[0, indices], 1e-4, thick=1, color=0
      
      legend, ['GOES15 0.1-0.8 nm', 'GOES15 0.05-0.4 nm'],$
            linestyle=[0,0], $
            color=[3, 5], $
            box=0, $
            pos=[0.13, 0.99], $
            /normal


      ;---------------------------------------;
      ;		      CALLISTO PLOTS
      ;
      loadct, 0
      low_data_bg = constbacksub(low_data, /auto)
      spectro_plot, low_data_bg > (-5) < 40 , $
          low_times, $
          low_freq, $
          /xs, $
          /ys, $
          xrange='2011-sep-22 '+['10:50:30','10:53:10'], $
          yr=[100, 10], $
          xticklen=-0.01, $
          yticklen=-0.01, $
          ytitle=' ', $
          yminor=2, $
          position=[x0, 0.38, 0.95, 0.63], $
          /noerase

      loadct, 74
      low_data_bg = constbacksub(low_data, /auto)
      spectro_plot, low_data_bg > (-5) < 40 , $
          low_times, $
          low_freq, $
          /xs, $
          /ys, $
          xtitle = ' ', $
          tick_unit = 60.0*20.0, $
          xtickname = [' ', ' ', ' '], $
          yticks = 2, $
          ytickname = [' ', ' ', ' '], $
          xrange='2011-sep-22 '+['10:50:30','10:53:10'], $
          yr=[100, 10], $
          xticklen=0, $
          yticklen=0, $
          ytitle=' ', $s
          position=[x0, 0.38, 0.95, 0.63], $
          /noerase
          
      ;---------------------------------------;
      ;	 Plot back projected Hough transform
      ;   
      loadct, 0
      restore, 'back_proj_hough.sav'
      spectro_plot, normal_back > 0.6 <1.0, time_set, freq_set, $
          /xs, $
          /ys, $
          yr = [90, 45], $
          xrange='2011-sep-22 '+['10:50:30','10:53:10'], $
          ytitle='Frequency (MHz)', $
          title='Backprojected Hough Transform', $
          position=[x0, 0.06, 0.95, 0.30], $
          xticklen=-0.01, $
          yticklen=-0.01, $
          yminor=2, $
          /normal, $
          /noerase
          
      loadct, 74    
      spectro_plot, normal_back > 0.6 <1.0, time_set, freq_set, $
          /xs, $
          /ys, $
          yr = [90, 45], $
          xtitle = ' ', $
          tick_unit = 60.0*20.0, $
          xtickname = [' ', ' ', ' '], $
          yticks = 2, $
          ytickname = [' ', ' ', ' '], $
          xrange='2011-sep-22 '+['10:50:30','10:53:10'], $
          xticklen=0, $
          yticklen=0, $
          ytitle=' ', $
          position=[x0, 0.06, 0.95, 0.30], $
          /normal, $
          /noerase


  device,/close
  set_plot,'x'


END