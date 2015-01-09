pro hough_herringbone,angle1,angle2,normal_back

!p.multi=[0,1,1]


cd,'/Users/eoincarley/Data/CALLISTO/20110607/Fits'
radio_spectro_fits_read,'BIR_20110607_062400_10.fit',z,x,y

x1=600
x2=1200
y1=closest(y,75)
y2=closest(y,35)

;=======================Plot normal spectrum===========================
;================Choose harsh clipping, best is 100-200================
window,1
spectro_plot,(bytscl(z[x1:x2,y1:y2],100,200)),x[x1:x2],y[y1:y2],/xs,/ys,$
ytitle='Frequency (Mhz)',title='Normal dynamic spectrum',charsize=1.5

;=======================Plot spectrum gradient=========================
window,2
spectro_plot,gradient(z[x1:x2,y1:y2]),x[x1:x2],y[y1:y2],/xs,/ys,$
ytitle='Frequency (Mhz)',title='Dynamic spectrum gradient',charsize=1.5

;==========================Hough Transform=============================
;=====================Choose appropriate angles========================
test = gradient(bytscl(z[x1:x2,y1:y2],0,200))
theta1 = (angle1*!pi)/180.0
theta2 = (angle2*!pi)/180.0
n_points = 1000
theta = (dindgen(n_points)*(theta2-theta1)/n_points)+theta1

;theta_1st = (dindgen(n_points)*(175-110)/n_points)+110
;theta_2nd = (dindgen(n_points)*(260-185)/n_points)+185
;theta = [theta_1st];,theta_2nd]

result = HOUGH(test, RHO=rho, THETA=theta, /gray) 

;======================Plot the Hough transform==========================
window,3,xs=1000,ys=600
result_n=result/max(result)
plot_image,bytscl(result_n,-1.0,1.0),xticks=1,xtickname=[' ',' '],xtitle='Theta (degrees)',$
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


;plot_image,congrid((smooth(normal_back,2))^3,700,400)  

spectro_plot,bytscl(normal_back,0,1),x[x1:x2],y[y1:y2],/xs,/ys,charsize=1.5,$
ytitle='Frequency (MHz)',title='Backprojected Hough transform'

END