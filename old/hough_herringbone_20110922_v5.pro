pro hough_herringbone_20110922_v5,angle1,angle2,normal_back, PLOT_HOUGH=plot_hough

;This version is new and imporved v3. It finds the half intensity time using the 
;result of the interpolate function


cd,'/Users/eoincarley/Data/CALLISTO/20110922'
radio_spectro_fits_read,'BIR_20110922_104459_01.fit',data_raw,times,freq

; ************* Define Time and Frequency Interval *****************
t1_index = closest(times,anytim(file2time('20110922_105000'),/utim))
t2_index = closest(times,anytim(file2time('20110922_105300'),/utim))
f1_index = closest(freq,75.0)
f2_index = closest(freq,45.0)

;* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
loadct,9
!p.color=255
!p.background=0
stretch,255,60


;***************Chosse clipping and Hough Angles********************

data_section = data_raw[t1_index:t2_index, f1_index:f2_index]
data_clip =  gradient(bytscl(constbacksub( data_section ,/auto),-20.0,100.0))


theta1 = (angle1*!Pi)/180.0
theta2 = (angle2*!Pi)/180.0
n_points = 1000.0
theta = (dindgen(n_points+1)*(theta2-theta1)/(n_points))+theta1
		;theta_1st = (dindgen(n_points+1)*(175.0-160.0)/n_points)+160.0
		;theta_2nd = (dindgen(n_points+1)*(190.0-185.0)/n_points)+185.0
		;theta = [theta_1st,theta_2nd]

result = HOUGH(data_clip, RHO=rho, THETA=theta, /gray) 


;***************Plot the Hough transform****************
IF keyword_set(plot_hough) THEN BEGIN
	window,3;,xs=1000,ys=600
	result_n=result/max(result)
	plot_image,(bytscl(result_n,-0.9,0.9)),xticks=1,xtickname=[' ',' '],xtitle='Theta (degrees)',$
	charsize=1.5,yticks=1,ytickname=[' ',' '],title='Radius-Angle parameter space of Hough Transform',$
	position=[0.1,0.1,0.9,0.9],/normal
	
	axis,xaxis=0,xr=[(theta1*180.0)/!Pi,(theta2*180.0)/!Pi],/xs,charsize=1.5
	axis,yaxis=0,yticks=10,yr=[rho[0],rho[n_elements(rho)-1]],charsize=1.5,$
	ytitle='Radius length from centre (pixel units)'
ENDIF

;***************Plot the Hough transform backprojection*************

window,4,xs=600,ys=500
t_range = ((t2_index+1.0)-t1_index)
f_range = (f2_index-f1_index) + 1.0
backproject = HOUGH(result, /BACKPROJECT, RHO=rho, THETA=theta, nx = t_range, ny = f_range) 
normal_back = smooth(backproject/max(backproject),3)

freq_set = freq[f1_index:f2_index]
time_set = times[t1_index:t2_index]
spectro_plot,bytscl(normal_back,0.5,1.0),time_set,freq_set,$
/xs,/ys,charsize=1.5,$
ytitle='Frequency (MHz)',title='Backprojected Hough Transform'

;************************************************************************
;         			  END OF HOUGH IMAGE PROCESSING
;************************************************************************




;************************************************************************************************
;				INTERPOLATION SEEMS TO HAVE DONE THE JOP OF 'JOINING THE DOTS'
;        		IN A SMOOTH WAY
;
peak_time_freq = dblarr(n_elements(freq_set), 250)  ;pad it with zeros
FOR k=0,n_elements(freq_set)-1 DO BEGIN
	f_slice = k;closest(freq_set,42.0) ;choose frequency slice and plot
	npoints=100000.0 
	xloc = (n_elements(time_set)-1.0)*dindgen(npoints+1.0)/(npoints) ;as was before
	freq_test=dblarr(npoints+1.0) ;as was before
	freq_test[*]= f_slice

	index = dindgen(n_elements(freq_set))+1.0
	yloc = interpol(index, freq_set, freq_test)  


	yloc = dblarr(npoints+1.0)
	yloc[*] = f_slice
	;Run the interpolation algorithm...
	profile_interp=interpolate(normal_back,xloc,yloc,cubic=-0.5)

	time_4profile = dindgen(npoints+1.0)*( (time_set[n_elements(time_set)-1.0]- time_set[0])/(npoints)  )$
			+time_set[0]

	wset,5;,xs=1200
	intensity = normal_back[*,f_slice]
	utplot,time_set[*],intensity,/xs,/ys,psym=4,ytitle='Intesnity',charsize=1.5
	oplot,time_set[*],intensity,color=240.0
	oplot,time_4profile,profile_interp


	;*******************PEAK FINDING******************

	;In the back projected Hough transform choose the right edge of a dark peak. This point corresponds 
	;to a peak in the original image e.g., center of a burst in the original dynamic spectrum
	;In a light curve of the transform at a particular frequency, the right edge is half 
	;way between a peak and a trough.
	;Prodecure: Build a code that finds the peaks and troughs and finds the time half way in between this
    
    ;*****************Get peak and trough points*****************
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
	loadct,39
	plots,peak_times,peak_intensity,psym=1,color=240,symsize=2
    
    
    ;*******************Get half way point***************
	burst_time = 0.0d
	half_inten = 0.0d
	burst_time_test = 0.0d
	times_from_interp = 0.0d
	FOR i=1, n_elements(peak_times)-1 DO BEGIN
  		IF peak_intensity[i] gt peak_intensity[i-1] THEN BEGIN
  				;Find intensity value half way between peak intensities
  				;Find the time value half way between peak times
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
	plots,burst_time,half_inten,psym=1,color=100,symsize=2

	wset,4
	plots,burst_time[4],freq_set[f_slice],psym=1,color=240,symsize=1

	loadct,9
	!p.color=255
	!p.background=0
	stretch,255,60
	wset,1
		spectro_plot,(bytscl(constbacksub(data_raw,/auto),-20,200)), times, freq, $
		/ys, ytitle='!6Frequency [MHz]', yticks=5, yminor=4, yr = [freq[f1_index],freq[f2_index]], $
		xrange=[times[t1_index],times[t2_index]], /xs, xtitle='Start time:'+anytim(times[t1_index],/yoh),$
		charsize=1.5

	loadct,39

	;plots,burst_time_test,freq_set[f_slice],psym=1,color=240,symsize=1  ;points closest to half intensity
	;plots,burst_time,freq_set[f_slice],psym=1,color=100,symsize=1       ;points closest to half intensity
	
	plots,times_from_interp, freq_set[f_slice], psym=3, color=230, symsize=2
     

	left_over = abs(n_elements(times_from_interp) - 249.0)
	padding = dblarr(left_over)
	concat = [times_from_interp,padding]
	peak_time_freq[k,0.0] = freq_set[f_slice] ;particular frequency
	peak_time_freq[k,1:249] = concat  ;<----All times of peaks for a particular freqeuncy.
									  ;This array is columns of times at partiuclar freq
									  ;with the frequency in the first element of the column

ENDFOR
	
;**************PLOT ALL DATA POINTS FOUND************
	loadct,9
	!p.color=255
	!p.background=0
	stretch,255,60
        wset,1
		spectro_plot,(bytscl(constbacksub(data_raw,/auto),-20,200)), times, freq, $
		/ys, ytitle='!6Frequency [MHz]', yticks=5, yminor=4, yr = [freq[f1_index],freq[f2_index]], $
		xrange=[times[t1_index],times[t2_index]], /xs, xtitle='Start time:'+anytim(times[t1_index],/yoh),$
		charsize=1.5
	loadct,39
	
	FOR i=0,n_elements(freq_set)-1 do begin
	    plots,peak_time_freq[i,1:99], peak_time_freq[i,0.0], color=230,psym=1,symsize=0.3
	ENDFOR
;******************************************************







;*******************************WTD*******************************
loadct,9
!p.color=255
!p.background=0
stretch,255,60
window,6,xs=500,ys=500
wtd = dblarr(n_elements(times_from_interp))
FOR i=2,n_elements(times_from_interp)-1 DO BEGIN
 wtd[i-2] = times_from_interp[i] - times_from_interp[i-1]
ENDFOR 

cumu_t = dindgen(1001)*(7.0)/1000.0
p = dblarr(n_elements(cumu_t))

FOR j=0.,n_elements(cumu_t)-1 DO BEGIN
time = cumu_t[j]
	FOR i=0.,n_elements(wtd)-1 DO BEGIN
  		IF wtd[i] gt time THEN p[j] = p[j]+1.0
	ENDFOR
	
ENDFOR
index_gt = where( p gt 0.0)

plot, cumu_t[index_gt], p[index_gt], $
charsize=1.5,ytitle='No. Burst > '+Greek('delta', /CAPITAL)+' t',$
xtitle=Greek('delta', /CAPITAL)+' t',xstyle=1,yr=[0.0,150.0],ystyle=1,$
psym=10


sim_cumu = dblarr(n_elements(cumu_t))
sim_cumu[*] = igamma(floor(cumu_t[*])+1.0,2.5)/factorial(floor(cumu_t[*]))
;oplot,cumu_t,sim_cumu 



hist = HISTOGRAM(wtd,binsize=0.3) 
bins = FINDGEN(N_ELEMENTS(hist)) + MIN(wtd) 
loadct,9
!p.color=0
!p.background=255
window,7,xs=500,ys=500
plot,bins,hist,psym=10,charsize=1.5,xtitle='Waiting time [s]',ytitle='No. of burts'

x = dindgen(101)*(10.0-0.0)/100.0
pois = dblarr(n_elements(x))
aver = 3.5
pois[*] = (aver^x[*])*exp(-aver)/factorial(x[*])
loadct,39
!p.color=0
!p.background=255
oplot,x,pois*n_elements(times_from_interp),color=240,thick=2;*n_elements(times_from_interp)
stop

END