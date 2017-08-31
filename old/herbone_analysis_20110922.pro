pro herbone_analysis_20110922

!p.multi=[0,1,1]

cd,'/Users/eoincarley/Data/CALLISTO/20110922'
radio_spectro_fits_read,'BIR_20110922_104459_01.fit',data,times,freq

t1=closest(times,anytim(file2time('20110922_105200'),/utim))
t2=closest(times,anytim(file2time('20110922_105300'),/utim))
y1=closest(y,80)
y2=closest(y,40)

loadct,9
!p.color=255
!p.background=0
stretch,255,60

!p.multi=[0,1,2]

a = closest(freq,45.0)
b  = closest(freq,70.0)

IF keyword_set(do_plot) THEN BEGIN
	FOR i=a,b,-1.0 DO BEGIN

		loadct,9
		!p.color=255
		!p.background=0
		stretch,255,60

		;*************Plot Spectra****************
		spectro_plot,bytscl(constbacksub(data,/auto),0,40),times,freq,/xs,/ys,ytitle='!6Frequency [MHz]',$
		charsize=1.5,xtitle='!6Start Time (2011-Sep-22 10:51:00 UT)',yr=[100,10],yticks=9,yminor=2,$
		xrange='2011-sep-22 '+['10:51:00','10:53:00']


		;*************Plot Lightcurve*************
		loadct,39
		freq_line = dblarr(n_elements(x))
		freq_line[*] = freq[i]
	
		plots,times[*],freq_line[*],color=240,/data
		;prof = closest(y,freq_line[0])
		utplot,times,data[*,i],/ys,xr=[times[t1], times[t2]],charsize=1.5,color=0,$
		title='Frequency: '+string(freq[i])
  	;****Best frequency 'slice' is at 53.938 MHz, y[104]
	ENDFOR
ENDIF

j_junk=0.0
answer = 'y'
data = constbacksub(data,/auto)
full_data = constbacksub(data,/auto)
full_freq = freq
WHILE answer ne 'n' DO BEGIN
;*****************Plot Overall lightcurve and select a peak***************
t1=closest(times,anytim(file2time('20110922_105100'),/utim))
t2=closest(times,anytim(file2time('20110922_105300'),/utim))
loadct,9
!p.color=0
!p.background=255
!p.multi=[0,1,1]
window,1;,xs=1500,ys=500
;utplot,x,z[*,104],/ys,xr=[x[x1],x[x2]],charsize=1.5,color=0,$
;title='Frequency: '+string(y[104]),yr=[145,190]
spectro_plot,bytscl(full_data,0,40),times,full_freq,/xs,/ys,ytitle='!6Frequency [MHz]',$
		charsize=1.5,xtitle='!6Start Time (2011-Sep-22 10:51:00 UT)',yr=[100,10],yticks=9,yminor=2,$
		xrange=[times[t1],times[t2]]


print,'Choose peak' 
cursor,peak_time,junk,/data
jump1:
t1_orig = closest(times,peak_time-10)
t2_orig = closest(times,peak_time+10)

spectro_plot,bytscl(constbacksub(full_data,/auto),0,40),times,full_freq,/xs,/ys,ytitle='!6Frequency [MHz]',$
		charsize=1.5,xtitle='!6Start Time (2011-Sep-22 10:51:00 UT)',yr=[100,10],yticks=9,yminor=2,$
		xr=[times[t1_orig],times[t2_orig]]
	;**** place marker on dyn spectra where selection was made
	loadct,39		
	plots,peak_time,53.938,psym=2,color=250
	loadct,9
	!p.color=255
	!p.background=0
	stretch,255,60
    ;**********************************************************


;x1=closest(x,anytim(file2time('20110922_105040'),/utim))
;x2=closest(x,anytim(file2time('20110922_105043'),/utim))
;x3=closest(x,anytim(file2time('20110922_105100'),/utim))


;******************Plot Zoom of selected Peak***********************

window,2  
utplot,times,data[*,104],/ys,xr=[times[t1_orig],times[t2_orig]],charsize=1.5,$
title='Frequency: '+string(freq[104]);,yr=[155,170]
		;**** place marker on lightcurve where selection was made
        line = dindgen(100.0)*(200.0 - 0.0)/99.0 
        peak_time_line = dblarr(n_elements(line))
        peak_time_line[*] = peak_time
		oplot,peak_time_line,line,linestyle=1
        ;******
        
print,'*********************'
print,'Choose start time:' 
print,'*********************'
cursor,tim1,junk,/data
wait,1
print,'*********************'
print,'Choose end time:' 
print,'*********************'
cursor,tim2,junk,/data
t1 = closest(times,tim1)
t2 = closest(times,tim2)
utplot,times,data[*,104],ystyle=1,xr=[times[t1],times[t2]],charsize=1.5,$
title='Frequency: '+string(freq[104])



;************Gaussfit**************
data_f = data[t1:t2,104]
data_t =times[t1:t2] - times[t1]

res = gaussfit(data_t,data_f,gaus_params,nterms=4)
    ;        Produce the Gauss
	tau0 =  dindgen(1001)*(data_t[n_elements(data_t)-1] - data_t[0])/1000.0
	tau1 = (tau0 - gaus_params[1])/gaus_params[2]
	y=dblarr(n_elements(tau0))
	y[*] = gaus_params[0]*exp( (-tau1[*]^2.0)/2.0 )+gaus_params[3]
	

;**********Oplot the fit**************
loadct,39
oplot,tau0+times[t1],y,color=60
freqs = ( dindgen(1001)*(170.0 - 120.0)/1000.0 )+150.0
tim = dblarr(n_elements(freqs))
time = gaus_params[1]+times[t1]
plots,time,freqs,color=240
print,'*********************'
print,'Try again?'
print,'*********************'
 yn=' '
 read,yn
  IF yn eq 'y' THEN GOTO, jump1
  k_junk=0
IF k_junk eq 0 THEN restore,'peak_times.sav'
  k_junk=1
  
peak_times = [peak_times,time]
save,peak_times,filename='peak_times.sav'
IF j_junk eq 0.0 THEN BEGIN 
	allthe_ys = fltarr(n_elements(y))
    allthe_ts = fltarr(n_elements(y))
ENDIF    
j_junk=1.0
allthe_ys = [allthe_ys,y]
allthe_ts = [allthe_ts,tau0+times[t1]]


test1=closest(times,anytim(file2time('20110922_105100'),/utim))
test2=closest(times,anytim(file2time('20110922_105300'),/utim))
loadct,39
!p.color=0
!p.background=255
!p.multi=[0,1,1]
wset,3
utplot,times,data[*,104],/ys,xr=[times[test1],times[test2]],charsize=1.5,color=0,$
title='Frequency: '+string(data[104]),yr=[0.0,40.0],/noerase
;oplot,t0+x[x1],y,color=60
oplot,allthe_ts,allthe_ys,color=60,psym=3
stop
oplot,times,freq,color=60
print,'*********************'
print,'Keep going?'
print,'*********************'
 answer=' '
 read,answer

ENDWHILE
END