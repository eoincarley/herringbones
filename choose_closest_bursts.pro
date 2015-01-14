pro choose_closest_bursts
  
  ; Procedure to chose the start and end of the spikes in herringbones.
  ; Then finds the closest points in peak_time_freq.sav.

  cd,'~/Data/22Sep2011_event/herringbones'
  !p.charsize = 1.5
  loadct, 1
  window, 1, xs=2300, ys=600

  radio_spectro_fits_read, 'BIR_20110922_104459_01.fit', data_raw, times, freq
  t1_index = closest(times,anytim(file2time('20110922_105000'),/utim))
  t2_index = closest(times,anytim(file2time('20110922_105300'),/utim))
  f1_index = closest(freq,80.0)
  f2_index = closest(freq,45.0)
  data = constbacksub(data_raw,/auto)

  spectro_plot, data , times, freq, $
      /ys, $
      ytitle='!6Frequency [MHz]', $
      yticks=5, yminor=4, $
      yr = [freq[f1_index], freq[f2_index]], $
      xrange=[times[t1_index], times[t2_index]], $
      /xs, $
      xtitle='Start time:'+anytim(times[t1_index],/yoh)+' (UT)'
  
  set_line_color
  restore,'peak_time_freq.sav'
  FOR i=0, n_elements(peak_time_freq[*,0])-1 DO BEGIN
      plots, peak_time_freq[i,1:99], peak_time_freq[i,0.0], color=3, psym=1, symsize=1
  ENDFOR




END