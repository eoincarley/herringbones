pro hough_herringbone_20110922_v2,angle1,angle2,normal_back, PLOT_HOUGH=plot_hough


cd,'/Users/eoincarley/Data/CALLISTO/20110922'
radio_spectro_fits_read,'BIR_20110922_104459_01.fit',data_raw,times,freq

; ************* Define Time and Frequency Interval *****************
t1_index = closest(times,anytim(file2time('20110922_105130'),/utim))
t2_index = closest(times,anytim(file2time('20110922_105200'),/utim))
f1_index = closest(freq,80.0)
f2_index = closest(freq,40.0)

;* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

;***********************Plot Dynamic Spectrum***********************
loadct,9
!p.color=255
!p.background=0
stretch,255,60
wset,1

;spectro_plot,(bytscl(constbacksub(data_raw,/auto),-20,200)), times, freq, $
;/ys, ytitle='!6Frequency [MHz]', yticks=5, yminor=4, yr = [80.0,40.0], $
;xrange='2011-sep-22 '+['10:51:30','10:52:20'], /xs, xtitle='!6Start Time (2011-Sep-22 10:51:00 UT)',$
;charsize=1.5


;***************Choose clipping and Hough Angles********************

data_section = data_raw[t1_index:t2_index ,f1_index:f2_index]
;window,3
data_clip =  gradient(bytscl(constbacksub( data_section ,/auto),-100,100))
;plot_image,congrid(data_clip,1000,500)

theta1 = (angle1*!Pi)/180.0
theta2 = (angle2*!Pi)/180.0
n_points = 1000.0
theta = (dindgen(n_points+1)*(theta2-theta1)/(n_points))+theta1

		;theta_1st = (dindgen(n_points+1)*(175.0-160.0)/n_points)+160.0
		;theta_2nd = (dindgen(n_points+1)*(190.0-185.0)/n_points)+185.0
		;theta = [theta_1st,theta_2nd]

result = HOUGH(data_clip, RHO=rho, THETA=theta, /gray) 
 ;result = (result - 50) > 0

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

;***************Plot the Hough transform backprojection***************
window,4
t_range = ((t2_index+1)-t1_index)
f_range = (f2_index-f1_index)
backproject = HOUGH(result, /BACKPROJECT, RHO=rho, THETA=theta, nx = t_range, ny = f_range) 
	;backproject = gradient(backproject)
normal_back = backproject/max(backproject)
	;normal_back=filter_image(normal_back,/median)
answer = 'y'


;WHILE answer eq 'y' DO BEGIN
wset,4

set_freq = freq[f1_index:f2_index]
set_time = times[t1_index:t2_index]
spectro_plot,bytscl(normal_back,0.4,1.0),set_time,set_freq,$
/xs,/ys,charsize=1.5,$
ytitle='Frequency (MHz)',title='Backprojected Hough Transform'
	;point,x,y
	;wset,1
	;plots,x,y,color=0
;ENDWHILE
;************Take a slice through the backprojection*****************




;********************************************************************

f_slice = closest(set_freq,50.0) ;choose frequency slice and plot


wset,5;,xs=1200
intensity = smooth(normal_back[*,f_slice],1)
utplot,set_time[*],intensity,/xs,/ys,psym=4,ytitle='Intesnity',charsize=1.5
oplot,set_time[*],intensity


	;In the back projected Hough transform choose the right edge of a dark peak. This point corresponds 
	;to a peak in the original image e.g., center of a burst in the original dynamic spectrum
	;In a light curve of the transform at a particular frequency, the right edge is half 
	;way between a peak and a trough.
	;Prodecure: Build a code that finds the peaks and troughs and finds the time half way in between this
    
    ;*****************Get peak and trough points*****************
	peak_times=0.0
	FOR i=1.0, n_elements(set_time)-2 DO BEGIN
   		IF intensity[i] gt intensity[i-1] and intensity[i] gt intensity[i+1] THEN BEGIN
   			IF peak_times[0] eq 0.0 THEN BEGIN
   				peak_times = set_time[i]
   				peak_intensity = intensity[i]
   			ENDIF ELSE BEGIN	
   		    	peak_times = [peak_times,set_time[i]]
   		    	peak_intensity = [peak_intensity,intensity[i]]
   			ENDELSE    
  		 ENDIF
   
   		IF intensity[i] lt intensity[i-1] and intensity[i] lt intensity[i+1] THEN BEGIN
   			IF peak_times[0] eq 0.0 THEN BEGIN
   				peak_times = set_time[i]
   				peak_intensity = intensity[i]
   			ENDIF ELSE BEGIN	
   		    	peak_times = [peak_times,set_time[i]]
   		    	peak_intensity = [peak_intensity,intensity[i]]
   			ENDELSE    
   		ENDIF
	ENDFOR   
	loadct,39
	plots,peak_times,peak_intensity,psym=1,color=240,symsize=2
    
    ;*******************Get half way point***************
	burst_time = 0.0d
	half_inten = 0.0d
	FOR i=1, n_elements(peak_times)-1 DO BEGIN
  		IF peak_intensity[i] gt peak_intensity[i-1] THEN BEGIN
  
    		burst_time = [burst_time, peak_times[i-1] +  (peak_times[i] - peak_times[i-1])/2.0]
    		half_inten = [half_inten, peak_intensity[i-1] +  (peak_intensity[i] - peak_intensity[i-1])/2.0]
    	
  		ENDIF  
	ENDFOR
	plots,burst_time,half_inten,psym=1,color=100,symsize=2
	
	
	
	


wset,4

plots,burst_time[4],set_freq[f_slice],psym=1,color=240,symsize=1


loadct,9
!p.color=255
!p.background=0
stretch,255,60
wset,1
spectro_plot,(bytscl(constbacksub(data_raw,/auto),-20,200)), times, freq, $
/ys, ytitle='!6Frequency [MHz]', yticks=5, yminor=4, yr = [freq[f1_index],freq[f2_index]], $
xrange=[times[t1_index],times[t2_index]], /xs, xtitle='!6Start Time (2011-Sep-22 10:51:00 UT)',$
charsize=1.5
loadct,39
;stop
plots,burst_time,set_freq[f_slice],psym=1,color=240,symsize=1


END