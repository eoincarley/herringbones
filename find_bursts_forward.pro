pro find_bursts_forward, save_bursts = save_bursts, $
          FIRST = first, SECOND = second
  ;
  ;
  ;11-Sep-2013 - Code to plot the points detected by the Hough transform
  ;
  ;
  cd, '~/Data/2011_sep_22/herringbones'
  !p.charsize = 1.5

  radio_spectro_fits_read,'BIR_20110922_104459_01.fit', data_raw, times, freq
  
  if keyword_set(first) then begin
      frequency1 = 33.0
      frequency2 = 12.0
      t1_index = closest(times,anytim(file2time('20110922_104730'),/utim))
      t2_index = closest(times,anytim(file2time('20110922_105100'),/utim))
      f1_index = closest(freq, frequency1)
      f2_index = closest(freq, frequency2)
      inten0 = -40
      inten1 = 10
      peak_file = 'peak_ft_first_master_forward.sav'
      burst_file = 'bursts_ft_first_master_forward.txt'
      zero_burst = -5.0
  endif
  
  if keyword_set(second) then begin
      frequency1 = 42.0
      frequency2 = 17.0
      t1_index = closest(times,anytim(file2time('20110922_105110'),/utim))
      t2_index = closest(times,anytim(file2time('20110922_105500'),/utim))
      f1_index = closest(freq, frequency1)
      f2_index = closest(freq, frequency2)
      peak_file = 'peak_ft_second_master_forward.sav'
      burst_file = 'bursts_ft_second_master_forward.txt'
      inten0 = -20
      inten1 = 50
      zero_burst = 0.0
  endif

  data_bs = smooth(constbacksub(data_raw, /auto), 3)
  xtit = 'Start time: '+anytim(times[t1_index], /cc)+ '(UT)'
  
  in = data_bs
  get_bg, in, times, freq, frequency1, frequency2, $
          out_bg, sig
          
  finuse = freq[f1_index:f2_index]
  data_displ = bytscl(data_bs, inten0, inten1)
  plot_herb = 'spectro_plot, data_displ , times, freq, /ys, ytitle = "Frequency [MHz]", ticks = 5, yminor = 4, yr = [freq[f1_index], freq[f2_index]], xrange = [times[t1_index],times[t2_index]], /xs, xtitle = xtit'
  
  restore, peak_file
  nfreqs = (size(peak_time_freq))[1]
  nburst = (size(peak_time_freq))[2]
  
  ;FOR i=0, nfreqs-1 DO BEGIN
  ;    plots, peak_time_freq[i, 1:(nburst-1)], peak_time_freq[i, 0.0], color=3, psym=1, symsize=1
  ;ENDFOR

  btimes = 0.0
  bf = 0.0
  drift = 0.0
  
  ft1_index = where(peak_time_freq[1, *] gt 0.0) ;start at second frequency column. First does not go through all points.
  ft1 = peak_time_freq[1, ft1_index]


  FOR j=2, n_elements(ft1)-1 DO BEGIN
    loadct, 5, /silent
    reverse_ct
    window, 1, xs=2300, ys=800          
    void = execute(plot_herb)
  
    comp_f = ft1[0]
    comp_t = ft1[j]
 
      i = 2
      WHILE i lt n_elements(peak_time_freq[*, 0])-1 DO BEGIN 
  
          ft_index = where(peak_time_freq[i, *] gt 0.0)
          ft = peak_time_freq[i, ft_index]

          f = ft[0]
          index_t = closest(ft, comp_t)
          t = ft[index_t]
          
          ;also put a check in here. If intensity of found point is within 1 stdev of nominal back ground then i=0
          
          sample0 = closest(times, t)    
          sample1 = closest(freq, f) 
          intensity = data_bs[sample0, sample1]
          samplef_i = closest(finuse, f)
          backg_sample = out_bg[samplef_i]
          sig_sample = sig[samplef_i]
          
          ;print, intensity, backg_sample, sig_sample
          ;if intensity le (backg_sample + 0.3*sig_sample) then i=0

          
          If abs(comp_t - t) gt 0.8 THEN BEGIN
            i = n_elements(peak_time_freq[*, 0])
          ENDIF ELSE BEGIN
            btimes = [btimes,  t]
            bf = [bf, f]
            wset,1
            plotsym, 0, /fill
            ;plots, comp_f, comp_t, color=4, symsize=0.5, psym=4
            set_line_color
            plots, t, f, color=4, symsize=0.5, psym=8
          ENDELSE
    
          comp_f = f
          comp_t = t
          i=i+1
  
      ENDWHILE
     
      negative_check = btimes[n_elements(btimes)-1] - btimes[0]

      if n_elements(btimes) gt 5 then begin
          btimes = btimes[1: n_elements(btimes)-1]
          bf = bf[1: n_elements(bf)-1]
          tindex = btimes
          findex = bf
          ;---------------------------------------------;
          ;				        Get profile
          ;
          FOR i=0, n_elements(btimes)-1 DO BEGIN
            tindex[i] = closest(times, btimes[i])
            findex[i] = closest(freq, bf[i])
          ENDFOR
          inten = interpolate(data_bs, tindex, findex) 

          window, 3, ysize=500, xpos = 400, ypos=800
          loadct, 0
          plot, bf, smooth(inten, 5), /xs, /ys, xr = [frequency1, frequency2], xtit='Frequency (MHz)', ytit='Intensity (DNs)'
          zero = fltarr(n_elements(bf))
          zero[*] = zero_burst
          oplot, bf, zero, linestyle=2
  
          ;---------------------------------------------;
          ;       Manually choose stopping frequency
          ;
          print, 'Click below 25 to store burst.'
          cursor, x, y, /data
          if x lt 25 then begin

              ;manual = 'n'
              ;READ, manual, prompt = 'Enter stop frequency manually? (y/n)'
              ;if manual eq 'y' then begin
              ;
              
              y=0
              wait,0.5
              print,'Choose stop frequency: '
              cursor, user_f, user_I, /data
              print,'Manual stop frequency: '+string(user_f)

              user_stop_f = closest(bf, user_f)
              bf = bf[ 0:user_stop_f ]
              btimes = btimes[ 0:user_stop_f ]
              loadct, 5, /silent
              reverse_ct
              window, 1, xs=2300, ys=800          
              void = execute(plot_herb)
              plotsym, 0, /fill
              set_line_color
              plots, btimes, bf, color=4, symsize=0.4, psym=8
              wait,2.0
              if keyword_set(save_bursts) then write_text, btimes, bf, inten, burst_file
        
          endif                
      endif
        
      btimes = 0.0
      bf=0.0
      drift=0.0
    
  ENDFOR

END


pro write_text, bt, bff, inten, burst_filename

  IF file_test(burst_filename) eq 1 THEN BEGIN
    readcol,burst_filename, btprev, bffprev, intenprev,$
    format = 'A,D,D'
  
    bt = [btprev, '-', anytim(bt, /ccs)]
      bff = [bffprev, !Values.F_NAN, bff]
    inten = [intenprev, !Values.F_NAN, inten]

  ENDIF

  writecol, burst_filename, anytim(bt, /ccs), bff, inten, fmt='(A,D,D)'

END


pro get_bg, in, times, freq, frequency1, frequency2, $
          out, sig

  index0 = closest(times, anytim(file2time('20110922_104730'),/utim))
  index1 = closest(times, anytim(file2time('20110922_104800'),/utim))
  index2 = closest(freq, frequency1)
  index3 = closest(freq, frequency2)

  in = in[index0:index1, index2:index3]
  out = average(in, 2)
  sig = stddev(in, dimension = 1)


END

