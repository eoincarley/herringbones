pro determine_detection_threshold, angle1, angle2, normal_back, $
          PLOT_HOUGH=plot_hough, SAVE_POINTS = save_points, $
          FIRST = first, SECOND = second

  ;To detect a burst, the Hough transform must identify a potential straight-line feature in 
  ;the background subtracted image that has been passed through a gradient function. It is 
  ;difficult to calculate the theoretical detection limit, given all of these processing steps. 
  ;However, we can perform a simple simulation to determine a detection threshold: 

  ;We firstly find a particularly clean area of our spectrogram in the same frequency range 
  ;as we detected the herringbones. In a 30x30 area of the spectrogram (about the size of 
  ;a herringbone) we get the mean value (noise floor) and standard deviation (estimate 
  ;of the Gaussian noise in our system). We next assign a line of pixels in our 30x30 
  ;area a very low DN e.g., noise floor + 0.01xStandard deviation. This simulates a 
  ;very weak single herringbone, barely above the noise floor. To make the simulation 
  ;more realistic we make the line broken and patchy along its length. We also apply 
  ;the rate of intensity fall-off that we found for our observed bursts. We pass this 
  ;simulated herringbone through the same algorithm (or set of algorithms) we used on 
  ;the real herringbones. We incrementally increase the line intensity in steps of 0.01 
  ;standard deviations until we get a positive detection. In this way we determined 
  ;the detection threshold for a single weak herringbone burst.

  
  loadct, 5
  !p.charsize = 1
  cd, '~/Data/2011_sep_22/herringbones'
  radio_spectro_fits_read, 'BIR_20110922_104459_01.fit', data_raw, times, freq

  ; Region of first second bursts.  Best performance is angle1 = 182, angle2 = 195
  t1 = anytim(file2time('20110922_105810'),/utim)
  t2 = anytim(file2time('20110922_105950'),/utim)
  t1_index = closest(times, t1)
  t2_index = closest(times, t2)
  f1_index = closest(freq, 60.0)
  f2_index = closest(freq, 40.0)
  inten0 = -20
  inten1 = 40
  filename_save = 'peak_ft_second_master_reverse.sav'


  ;Sample background
  t1_bg = anytim(file2time('20110922_105900'),/utim)
  t2_bg = anytim(file2time('20110922_105926'),/utim)
  t1_bg_index = closest(times, t1_bg)
  t2_bg_index = closest(times, t2_bg)
  f1_bg_index = closest(freq, 60.0)
  f2_bg_index= closest(freq, 40.0)
  data_bg_sample = data_raw[t1_bg_index:t2_bg_index, f1_bg_index:f2_bg_index]

  ;window, 1
  ;plot_image, data_bg_sample > 100 < 240
  noise_floor = mean(data_bg_sample)
  noise = stdev(data_bg_sample)

  ;***********Simulate burst here**************
  t1_sim_burst = anytim(file2time('20110922_105900'),/utim)
  t2_sim_burst = anytim(file2time('20110922_105905'),/utim)
  t1_sim_burst_index = closest(times, t1_sim_burst)
  t2_sim_burst_index = closest(times, t2_sim_burst)
  f1_sim_burst_index = closest(freq, 60.0)
  f2_sim_burst_index = closest(freq, 40.0)
  fsims_index = reverse((dindgen(100)*(f2_sim_burst_index - f1_sim_burst_index)/99.) + f1_sim_burst_index)
  tsims_index = (dindgen(100)*(t2_sim_burst_index - t1_sim_burst_index)/99.) + t1_sim_burst_index


  t1_sim_burst = anytim(file2time('20110922_105900'),/utim)
  t2_sim_burst = anytim(file2time('20110922_105905'),/utim)
  t1_sim_burst_index = closest(times, t1_sim_burst)
  t2_sim_burst_index = closest(times, t2_sim_burst)
  f1_sim_burst_index = closest(freq, 60.0)
  f2_sim_burst_index = closest(freq, 45.0)
  region_stdev = stdev( data_raw[t1_sim_burst_index:t2_sim_burst_index , f1_sim_burst_index:f2_sim_burst_index] )

  local_mean = dindgen(100)
  local_stdev = dindgen(100)
  local_size = reverse(dindgen(100)*(5 - 0)/99.) 
  randsize = RANDOMN(1, 100)*reverse( (dindgen(100)*(2.0)/99.) )
  local_size = round(local_size + randsize)
  randnoise = RANDOMN(1, 100)*reverse((dindgen(100)*(2.0)/99.))
  
  ;sim_inten = dblarr(10, 100)
  ;for i =0, 99 do sim_inten[*, i] = 
  ;gaussian( dindgen(10), [local_mean[i], 5, 1+local_width*2])
  enhance = dindgen(100)*(2.0 - 0.0)/99.

  burst_life = dindgen(100)*(5)/99.0


  window, 0
  for j = 0, n_elements(enhance)-1 do begin

    for i=0, 99 do begin
      if i eq 0 then x_index=closest(dindgen(n_elements(times)), tsims_index[i]) else x_index=[x_index, closest(dindgen(n_elements(times)), tsims_index[i])]
      if i eq 0 then y_index=closest(dindgen(n_elements(freq)), fsims_index[i]) else y_index=[y_index, closest(dindgen(n_elements(freq)), fsims_index[i])] 
      
      local_mean[i] = mean( data_raw[x_index[i]-5:x_index[i]+5 , y_index[i]-5:y_index[i]+5] )
      local_stdev[i] = stdev( data_raw[x_index[i]-5:x_index[i]+5 , y_index[i]-5:y_index[i]+5] )

      local_noise = local_stdev[i];*randnoise[i]
      local_inten = local_mean[i] ;+ local_noise
      local_width = abs(local_size[i])
      local_profile = local_inten - RANDOMN(1, 1+local_width*2)*(local_noise*0.5)

      burst_inten = region_stdev*enhance[j] - 2.0*burst_life > 0.0

      local_gauss = local_profile + gaussian(dindgen(n_elements(local_profile)), [burst_inten[i], n_elements(local_profile)/2, local_width*(enhance[j]/10.)])
     
      data_raw[ (x_index[i]-local_width) : (x_index[i]+local_width), y_index[i]] = local_gauss;*local_inten;*local_gauss
    endfor

    wset, 0
    spectro_plot, data_raw > 100 , times, freq, $
        /xs, $
        /ys, $
        charsize=1.5, $
        ytitle='Frequency (MHz)', $
        title='Backprojected Hough Transform', $
        xr = [t1, t2], $
        yr=[60, 40]

    print, enhance[j]
    wait, 0.1    
  endfor    
  
stop
  ;---------  Chosse intensity clipping and Hough angles  --------;
  data_bs = smooth(constbacksub(data_raw, /auto), 1)
  data_section = data_bs[t1_index:t2_index, f1_index:f2_index]
  data_clip =  gradient(bytscl(data_section, inten0, inten1))

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
  normal_back = smooth( backproject/max(backproject), 3)

  freq_set = freq[f1_index:f2_index]
  time_set = times[t1_index:t2_index]
  spectro_plot, bytscl(normal_back, 0.5, 1.0), time_set, freq_set, $
      /xs, $
      /ys, $
      charsize=1.5, $
      ytitle='Frequency (MHz)', $
      title='Backprojected Hough Transform'
  ;save, normal_back, time_set, freq_set, filename = 'back_proj_hough.sav'

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
    

      loadct, 9, /silent
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
      peak_time_freq[k, 0] = freq_set[f_slice] ;particular frequency
      peak_time_freq[k, 1:249] = concat  ;<----All times of peaks for a particular freqeuncy.
                        ;This array is columns of times at partiuclar freq
                        ;with the frequency in the first element of the column
 
  ENDFOR  

  plots, burst_times, freq_set, psym=1, color=3, symsize=2

  ;-------------- PLOT ALL DATA POINTS FOUND -----------------
  
  loadct, 74
  window, 1, xs=2400, ys=700
  spectro_plot,( bytscl(data_bs, inten0, inten1) ), times, freq, $
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
  plotsym, 0, /fill
  cd,'~/Data/2011_sep_22/herringbones'
  ;if keyword_set(save_points) then save, peak_time_freq, filename=filename_save

  stop
  
  FOR i=0, n_elements(freq_set)-1 do begin
      plots, peak_time_freq[i,1:99], peak_time_freq[i,0.0], color=4, psym=8, symsize=0.5
      plots, peak_time_freq[i,1:99], peak_time_freq[i,0.0], color=0, psym=8, symsize=0.4    
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
  ;utplot, time_4profile, profile_interp, /xs, /ys

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

        ;plots, time_profile[index_I], profile_slice[index_I], psym=6, color=200 
        ;times_from_interp = [times_from_interp, time_profile[index_I]]
	
  		ENDIF  
	ENDFOR
	
END	


		