pro hough_herringbone_20110922_v6, angle1, angle2, normal_back, PLOT_HOUGH=plot_hough

  ;v5 first finds all bursts for a particular frequency slice and then succesively moves onto the
  ;next frequencies. Need to do it the other way around e.g., get the first peak and trace this peak 
  ;all the way trough so I can get an array indicating individual bursts in columns

  ;This version (v6) is now fucntionalised. 

  ; best performance is angle1 = 180, angle2 = 200
  loadct, 1
  !p.charsize = 1
  cd,'~/Data/22Sep2011_event/herringbones'
  radio_spectro_fits_read, 'BIR_20110922_104459_01.fit', data_raw, times,freq

  ;----------  Define Time and Frequency Interval  ----------- 
  t1_index = closest(times,anytim(file2time('20110922_105000'),/utim))
  t2_index = closest(times,anytim(file2time('20110922_105300'),/utim))
  f1_index = closest(freq, 80.0)
  f2_index = closest(freq, 45.0)


  ;---------  Chosse intensity clipping and Hough angles  --------;

  data_section = data_raw[t1_index:t2_index, f1_index:f2_index]
  data_clip =  gradient(bytscl(constbacksub( data_section, /auto), -20.0, 50.0))

  theta1 = angle1*!dtor
  theta2 = angle2*!dtor
  n_points = 1000.0
  theta = (dindgen(n_points+1)*(theta2-theta1)/(n_points))+theta1

  result = HOUGH(data_clip, RHO=rho, THETA=theta, /gray) 

  ;--------  Plot the Hough transform  ----------
  IF keyword_set(plot_hough) THEN BEGIN
    window,3
    result_n=result/max(result)
    plot_image, result_n > (-2.0) < 2.0, $
        xticks=1, $
        xtickname=[' ',' '], $
        xtitle='Theta (degrees)', $
        yticks=1, $
        ytickname=[' ',' '], $
        title='Radius-Angle parameter space of Hough Transform', $
        position=[0.1, 0.1, 0.9, 0.9], $
        /normal
  
    axis, xaxis=0, xr=[ theta1*!radeg, theta2*!radeg ], /xs
    axis, yaxis=0, yticks=10,yr=[rho[0],rho[n_elements(rho)-1]], $
    ytitle='Radius length from centre (pixel units)'
  ENDIF

  ;---------- Plot the Hough transform backprojection  -----------

  window, 4, xs=1200, ys=500
  t_range = ((t2_index+1.0)-t1_index)
  f_range = (f2_index-f1_index) + 1.0
  backproject = HOUGH(result, /BACKPROJECT, RHO=rho, THETA=theta, nx = t_range, ny = f_range) 
  normal_back = smooth( backproject/max(backproject) , 3)

  freq_set = freq[f1_index:f2_index]
  time_set = times[t1_index:t2_index]
  spectro_plot, bytscl(normal_back, 0.5, 1.0), time_set, freq_set, $
      /xs, $
      /ys, $
      charsize=1.5, $
      ytitle='Frequency (MHz)', $
      title='Backprojected Hough Transform'


  ;         			  END OF HOUGH IMAGE PROCESSING
  ;
  ;************************************************************************;




  ;************************************************************************;
  ;			INTERPOLATION SEEMS TO HAVE DONE THE JOP OF 'JOINING THE DOTS'
  ;     IN A SMOOTH WAY
  ;
  burst_times = 0.0
  index = 0.0
  peak_time_freq = dblarr(n_elements(freq_set), 250)  ;pad it with zeros
  FOR k=n_elements(freq_set)-1, 0, -1 DO BEGIN
  
      f_slice = k ;closest(freq_set,44) ;choose frequency slice and plot
      get_interp_profile, normal_back, freq_set, time_set, f_slice, $
                  time_4profile, profile_interp
              
      intensity = normal_back[*,f_slice]
      

      ; PEAK FINDING:
      ; In the back projected Hough transform choose the left edge of a dark peak. 
      ; This point corresponds to a peak in the original image e.g., center of a 
      ; burst in the original dynamic spectrum. In a light curve of the transform 
      ; at a particular frequency, the left edge is half way between a peak and a trough.
      ; Prodecure: Build a code that finds the peaks and troughs and finds the time half 
      ; way in between this
    
    
      find_peaks, time_set, intensity, $
              peak_times, peak_intensity ;Get peak and trough points
    


      get_half_intensity_time, peak_intensity, peak_times, time_4profile, profile_interp, $
                     burst_time, half_inten, times_from_interp   
                  ; Get half way point burst_time is just from the half-way point in between troughs
                  ; times_from_interp is the half-way point taken from the interpolated
                  ; profile.
    

      loadct, 9
      wset, 4
      spectro_plot, bytscl(normal_back, 0.5, 1.0), time_set, freq_set, $
          /xs, $
          /ys, $
          charsize=1.5, $
          ytitle='Frequency (MHz)', $
          title='Backprojected Hough Transform'
    
      set_line_color
      plots, times_from_interp, freq_set[f_slice], psym=1, color=3, symsize=2

      left_over = abs(n_elements(times_from_interp) - 249.0)
      padding = dblarr(left_over)
      concat = [times_from_interp,padding]
      peak_time_freq[k,0.0] = freq_set[f_slice] ;particular frequency
      peak_time_freq[k,1:249] = concat  ;<----All times of peaks for a particular freqeuncy.
                        ;This array is columns of times at partiuclar freq
                        ;with the frequency in the first element of the column
 
  ENDFOR  

  plots, burst_times, freq_set, psym=1, color=3, symsize=2

  ;-------------- PLOT ALL DATA POINTS FOUND -----------------
  
  loadct,1
  reverse_ct
  window, 1, xs=2400, ys=700
  spectro_plot,( bytscl(constbacksub(data_raw,/auto), -20, 50) ), times, freq, $
      /ys, $
      ytitle = '!6Frequency [MHz]', $
      yticks = 5, $
      yminor = 4, $
      yr = [freq[f1_index],freq[f2_index]], $
      xrange = [times[t1_index],times[t2_index]], $
      /xs, $
      xtitle = 'Start time: '+anytim(times[t1_index], /cc, /trun)+' (UT)', $
      charsize = 2.0
    
  set_line_color
  save, peak_time_freq, filename='peak_time_freq.sav'
  FOR i=0,n_elements(freq_set)-1 do begin
      plots,peak_time_freq[i,1:99], peak_time_freq[i,0.0], color=3, psym=1, symsize=1
  ENDFOR


END
;
;						END OF MAIN CODE						   ;
;																   
;******************************************;



;******************************************;
;																   
;			Various Peak Finding Procedure				   

pro get_interp_profile, normal_back, freq_set, time_set, f_slice,  $
                        time_4profile, profile_interp

	npoints=100000.0 
	xloc = (n_elements(time_set)-1.0)*dindgen(npoints+1.0)/(npoints) ;as was before
	freq_test=dblarr(npoints+1.0) ;as was before
	freq_test[*]= f_slice
	index = dindgen(n_elements(freq_set))+1.0
	yloc = interpol(index, freq_set, freq_test)  
	yloc = dblarr(npoints+1.0)
	yloc[*] = f_slice
	; Run the interpolation algorithm...
	profile_interp=interpolate(normal_back,xloc,yloc,cubic=-0.5)
	time_4profile = dindgen(npoints+1.0)*( (time_set[n_elements(time_set)-1.0]- time_set[0])/(npoints)  )$
			+time_set[0]

END

pro find_peaks, time_set, intensity, $
                peak_times, peak_intensity

; Procedure to find peak and trough points
	peak_times=0.0
	FOR i=1.0, n_elements(time_set)-2 DO BEGIN
   		IF intensity[i] gt intensity[i-1] and intensity[i] gt intensity[i+1] THEN BEGIN
   			IF peak_times[0] eq 0.0 THEN BEGIN
   				peak_times = time_set[i]
   				peak_intensity = intensity[i]
   			ENDIF ELSE BEGIN	
   		    	peak_times = [peak_times,time_set[i]]
   		    	peak_intensity = [peak_intensity,intensity[i]]
   			ENDELSE    
  		 ENDIF
   
   		IF intensity[i] lt intensity[i-1] and intensity[i] lt intensity[i+1] THEN BEGIN
   			IF peak_times[0] eq 0.0 THEN BEGIN
   				peak_times = time_set[i]
   				peak_intensity = intensity[i]
   			ENDIF ELSE BEGIN	
   		    	peak_times = [peak_times,time_set[i]]
   		    	peak_intensity = [peak_intensity,intensity[i]]
   			ENDELSE    
   		ENDIF
	ENDFOR 
	
END	

pro get_half_intensity_time, peak_intensity, peak_times, time_4profile, profile_interp,$
					         burst_time, half_inten, times_from_interp
	
  burst_time = 0.0d
  half_inten = 0.0d
  burst_time_test = 0.0d
  times_from_interp = 0.0d
	FOR i=1, n_elements(peak_times)-1 DO BEGIN
  		IF peak_intensity[i] gt peak_intensity[i-1] THEN BEGIN
  			; Find intensity value half way between peak intensities
  			; Find the time value half way between peak times
    		burst_time = [burst_time, peak_times[i-1] +  (peak_times[i] - peak_times[i-1])/2.0]
    		half_inten = [half_inten, peak_intensity[i-1] +  (peak_intensity[i] - peak_intensity[i-1])/2.0]
    
        indeces = where(time_4profile gt peak_times[i-1] and time_4profile lt peak_times[i])
        time_profile = time_4profile[indeces]
        profile_slice = profile_interp[indeces]
        index_I = closest(profile_slice, half_inten[n_elements(half_inten)-1])

        plots, time_profile[index_I], profile_slice[index_I], psym=6, color=200 
        times_from_interp = [times_from_interp, time_profile[index_I]]
	
  		ENDIF  
	ENDFOR
	
END	


		