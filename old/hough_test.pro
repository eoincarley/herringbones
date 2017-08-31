pro hough_test,angle1,angle2,light_curve_hough

cd,'/Users/eoincarley/Data/CALLISTO/20110607/Fits'
radio_spectro_fits_read,'BIR_20110607_062400_10.fit',z,x,y

x1=450
x2=550
;y1=closest(y,40)
;y2=closest(y,70)
y1=0
y2=150


;=======================Plot normal spectrum===========================
;================Choose harsh clipping, best is 180-200================
window,1
spectro_plot,(bytscl(z[x1:x2,y1:y2],150,160)),x[x1:x2],y[y1:y2],/xs,/ys,$
ytitle='Frequency (Mhz)',title='Normal dynamic spectrum',charsize=1.5

;=======================Plot spectrum gradient=========================
window,2
spectro_plot,sobel(bytscl(z[x1:x2,y1:y2],150,160)),x[x1:x2],y[y1:y2],/xs,/ys,$
ytitle='Frequency (Mhz)',title='Dynamic spectrum gradient',charsize=1.5

;==========================Hough Transform=============================
;=================Scale the test spectrum harshly also=================
test = (bytscl(z[x1:x2,y1:y2],150,160))
theta1 = (angle1*!pi)/180
theta2 = (angle2*!pi)/180
n_points = 1000
theta = (dindgen(n_points)*(theta2-theta1)/n_points)+theta1

result = HOUGH(test, RHO=rho, THETA=theta, /gray) 

;======================Plot the Hough transform==========================
window,3,xs=1000,ys=600
result_n=result/max(result)
plot_image,bytscl(result_n,0,0.7),xticks=1,xtickname=[' ',' '],xtitle='Theta (degrees)',$
charsize=1.5,yticks=1,ytickname=[' ',' '],title='Radius-Angle parameter space of Hough Transform'
axis,xaxis=0,xr=[(theta1*180)/!pi,(theta2*180)/!pi],/xs,charsize=1.5
axis,yaxis=0,yticks=10,yr=[rho[0],rho[n_elements(rho)-1]],charsize=1.5,$
ytitle='Radius length from centre (pixel units)'

;===================Plot the Hough transform backprojection==============
window,4
backproject = HOUGH(result, /BACKPROJECT, RHO=rho, THETA=theta,nx=x2-x1,ny=y2-y1) 
;backproject = gradient(backproject)
normal_back = backproject/max(backproject)
;normal_back=filter_image(normal_back,/median)
a = anytim(file2time('20110607_062800'),/utim)
b = anytim(file2time('20110607_062900'),/utim)
;spectro_plot,normal_back^3,x[x1:x2],y[y1:y2],/xs,/ys,charsize=1.5,$
spectro_plot,bytscl((normal_back),0.4,0.8), x[x1:x2],y[y1:y2],/xs,/ys,charsize=1.5,$

ytitle='Frequency (MHz)',title='Backprojected Hough transform'

;plot_image,congrid((normal_back),800,600),$
;title='Back-projected'
;window,5
;utplot,x[x1:x2],normal_back[*,50],charsize=1.5

light_curve_hough = dblarr(2,n_elements(normal_back[*,50]))
light_curve_hough(0,*) = x[x1:x2-1]
light_curve_hough(1,*) = normal_back[*,50]


end