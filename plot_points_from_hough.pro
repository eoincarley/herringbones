pro plot_points_from_hough
  ;
  ;
  ;11-Sep-2013 - Code to plot the points detected by the Hough transform
  ;
  ;
  cd, '~/Data/22Sep2011_event/herringbones'
  !p.charsize = 1.5
  loadct, 1
  reverse_ct
  window, 1, xs=1300, ys=600
  ;window, 2, xs=1300, ys=600

  radio_spectro_fits_read,'BIR_20110922_104459_01.fit', data_raw, times, freq
  ;t1_index = closest(times, anytim(file2time('20110922_105000'),/utim))
  ;t2_index = closest(times, anytim(file2time('20110922_105300'),/utim))
  ;f1_index = closest(freq, 80.0)
  ;f2_index = closest(freq, 45.0)
  
  ; Region of first set of herringbones. Choose angles 190 and 210.
       t1_index = closest(times,anytim(file2time('20110922_104730'),/utim))
       t2_index = closest(times,anytim(file2time('20110922_105000'),/utim))
       f1_index = closest(freq, 60.0)
       f2_index = closest(freq, 32.0)

  
  
  data_bs = constbacksub(data_raw, /auto)
  xtit = 'Start time: '+anytim(times[t1_index], /cc)+ '(UT)'
  
  in = data_bs
  
          
  finuse = freq[f1_index:f2_index]
  
  spectro_plot, bytscl(data_bs, -20, 50) , times, freq, $
      /ys, $
      ytitle = '!6Frequency [MHz]', $
      yticks = 5, $
      yminor = 4, $
      yr = [freq[f1_index], freq[f2_index]], $
      xrange = [times[t1_index],times[t2_index]], $
      /xs, $
      xtitle = xtit
     
    
  
  set_line_color
  restore,'peak_time_freq_0.sav'
  nfreqs = (size(peak_time_freq))[1]
  nburst = (size(peak_time_freq))[2]
  
  ;FOR i=0, nfreqs-1 DO BEGIN
  ;    plots, peak_time_freq[i, 1:(nburst-1)], peak_time_freq[i, 0.0], color=3, psym=1, symsize=1
  ;ENDFOR


  btimes = 0.0
  bf = 0.0
  drift = 0.0
  inten_sample = fltarr(5)
  
  ;window,5
  ;window,4
  ;window,3
  ;window,2
  ;plot, [0,0], [0,0], xr=[40,80], yr=[140,200]
  
  
  ft1_index = where(peak_time_freq[nfreqs - 3, *] gt 0.0) ;start at second frequency column. First does not go through all points.
  ft1 = peak_time_freq[nfreqs - 3, ft1_index]

  FOR j=1, n_elements(ft1)-1 DO BEGIN
    comp_f = ft1[0]
    comp_t = ft1[j]

    in = data_bs
    get_bg, in, times, freq, comp_t, $
          outback, sig
          
        
 
      i = n_elements(peak_time_freq[*, 0]) - 2
      WHILE i gt 1 DO BEGIN 
  
          ft_index = where(peak_time_freq[i, *] gt 0.0)
          ft = peak_time_freq[i, ft_index]

          f = ft[0]
          index_t = closest(ft, comp_t)
          t = ft[index_t]
          
          ;also put a check in here. If intensity of found point is within 1 stdev of nominal back ground then i=0
          
          sample0 = closest(times, t)    
          sample1 = closest(freq, f) 
          intensity = data_bs[sample0, sample1]
                    ;samplef_i = closest(finuse, f)
                    ;backg_sample = outback[samplef_i]
                    ;sig_sample = sig[samplef_i]
          
          inten_sample = [inten_sample, intensity]
          inten_diff = inten_sample[n_elements(inten_sample)-1] - inten_sample[n_elements(inten_sample)-3]
          ;print, inten_sample, 
           ;print, 'Diff: '+string(inten_diff) 
           ;print, 'STDEV: '+string(stdev(inten_sample))
          if inten_diff lt (-5.0) then begin
            
            i=0
          endif  
         
                    ;print, intensity, backg_sample, sig_sample
                    ;if intensity le (backg_sample + 0.0*sig_sample) then i=0
       
          wset,1
          If abs(comp_t - t) gt 1.0 THEN BEGIN
            i=0
          ENDIF ELSE BEGIN
            set_line_color
            btimes = [btimes, comp_t, t]
            bf = [bf, comp_f, f]
            plots, comp_f, comp_t, color=4, symsize=0.5, psym=1
            plots, t, f, color=4, symsize=0.5, psym=1
          ENDELSE
    
          comp_f = f
          comp_t = t
          i=i-1
  
      ENDWHILE
      ;print, 'STDEV: '+string(stdev(inten_sample))
     
      
      
      if n_elements(btimes) gt 2 then begin
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
        prof= interpolate(data_bs, tindex, findex) 
  
        bt = btimes
        bff = bf
        inten = prof
  
        ;write_text, bt, bff, inten
       
      endif
        
      btimes = 0.0
      bf=0.0
      drift=0.0

  ENDFOR

END


pro write_text, bt, bff, inten

IF file_test('bursts_bs_hough.txt') eq 1 THEN BEGIN
	readcol,'bursts_bs_hough.txt', btprev, bffprev, intenprev,$
	format = 'A,D,D'
	
	bt = [btprev, '-', anytim(bt, /ccs)]
    bff = [bffprev, !Values.F_NAN, bff]
	inten = [intenprev, !Values.F_NAN, inten]

ENDIF


writecol, 'bursts_bs_hough.txt', anytim(bt, /ccs), bff, inten, fmt='(A,D,D)'


END


pro get_bg, in, times, freq, compt,  $
          out, sig

  t1 = anytim(compt, /utim) - 10.0 > 0
  t2 = anytim(compt, /utim) + 10.0
  index0 = closest(times, t1)
  index1 = closest(times, t2)
  index2 = closest(freq, 60.0)
  index3 = closest(freq, 32.0)
  
  in = in[index0:index1, index2:index3]

  out = fltarr(n_elements(in[0, *]))
  for i=0, n_elements(in[0, *])-1 DO BEGIN
    out[i] = mean(in[*, i])
  endfor   
  sig = stddev(in, dimension = 1)


END