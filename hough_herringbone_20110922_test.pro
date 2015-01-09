pro hough_herringbone_20110922_test,angle1,angle2,normal_back, PLOT_HOUGH=plot_hough


cd,'/Users/eoincarley/Data/CALLISTO/20110922'
radio_spectro_fits_read,'BIR_20110922_104459_01.fit',data_raw,times,freq

; ************* Define Time and Frequency Interval *****************
t1_index = closest(times,anytim(file2time('20110922_104500'),/utim))
t2_index = closest(times,anytim(file2time('20110922_105200'),/utim))
f1_index = closest(freq,100.0)
f2_index = closest(freq,10.0)

;* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

;***********************Plot Dynamic Spectrum***********************
loadct,1
!p.color=255
!p.background=0
stretch,255,60
;window,1

;spectro_plot,(bytscl(constbacksub(data_raw,/auto),-20,200)), times, freq, $
;/ys, ytitle='!6Frequency [MHz]', yticks=5, yminor=4, yr = [80.0,40.0], $
;;xrange='2011-sep-22 '+['10:51:00','10:52:00'], /xs, xtitle='!6Start Time (2011-Sep-22 10:51:00 UT)',$
;charsize=1.5

set_plot,'ps'
device,filename = 'Hough_herringbones.ps',/color,bits=8,/inches,/encapsulate,$
xsize=10,ysize=8,xoffset=9.0

;***************Chosse clipping and Hough Angles********************

data_section = data_raw[t1_index:t2_index ,f1_index:f2_index]
;window,3
data_clip = gradient(bytscl(data_section,100,200)); gradient(bytscl( data_section,-100,100))
;plot_image,congrid(data_clip,1000,500)

theta1 = (angle1*!Pi)/180.0
theta2 = (angle2*!Pi)/180.0
n_points = 1000.0
theta = (dindgen(n_points)*(theta2-theta1)/(n_points-1.0))+theta1

		;theta_1st = (dindgen(n_points)*(175-110)/n_points)+110
		;theta_2nd = (dindgen(n_points)*(260-185)/n_points)+185
		;theta = [theta_1st,theta_2nd]

result = HOUGH(data_clip, RHO=rho, THETA=theta, /gray) 
 ;result = (result - 50) > 0

;***************Plot the Hough transform****************
IF keyword_set(plot_hough) THEN BEGIN
	;window,3;,xs=1000,ys=600
	result_n=result/max(result)
	plot_image,alog(bytscl(result_n,-0.1,0.7)),xticks=1,xtickname=[' ',' '],xtitle='Theta (degrees)',$
	charsize=1.5,yticks=1,ytickname=[' ',' '],title='Radius-Angle parameter space of Hough Transform',$
	position=[0.15,0.1,0.95,0.9],/normal,color=254
	
	axis,xaxis=0,xr=[(theta1*180.0)/!Pi,(theta2*180.0)/!Pi],/xs,charsize=1.5,color=254
	axis,yaxis=0,yticks=10,yr=[rho[0],rho[n_elements(rho)-1]],charsize=1.5,$
	ytitle='Radius length from centre (pixel units)',color=254
ENDIF


;***************Plot the Hough transform backprojection***************
;window,4
t_range = (t2_index-t1_index)
f_range = (f2_index-f1_index)

backproject = HOUGH(result, /BACKPROJECT, RHO=rho, THETA=theta, nx = t_range, ny = f_range) 
	;backproject = gradient(backproject)
normal_back = backproject/max(backproject)
	;normal_back=filter_image(normal_back,/median)
answer = 'y'


;WHILE answer eq 'y' DO BEGIN


;spectro_plot,bytscl(normal_back,0.5,1.0),times[t1_index:t2_index],freq[f1_index:f2_index],$
;/xs,/ys,charsize=1.5,$
;ytitle='Frequency (MHz)',title='Backprojected Hough Transform'
;	point,x,y
	;wset,1
	;plots,x,y,color=0
;ENDWHILE
;************Take a slice through the backprojection*****************

device,/close
set_plot,'x'
END