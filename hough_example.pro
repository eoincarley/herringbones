PRO hough_example  
  
  ;Create an image with a random set of pixels  
  seed = 12345   ; remove this line to get different random images  
  array = RANDOMU(seed,128,128) GT 0.95  
  
  ;Draw three lines in the image  
  x = FINDGEN(32)*4  
  array[x,0.5*x+20] = 1b  
  array[x,0.5*x+30] = 1b  
  array[-0.5*x+100,x] = 1b  
  
  ;Create display window, set graphics properties  
  WINDOW, XSIZE=330,YSIZE=630, TITLE='Hough Example'  
  !P.BACKGROUND = 255  ; white  
  !P.COLOR = 0 ; black  
  !P.FONT=2  
  ERASE  
  
  XYOUTS, .1, .94, 'Noise and Lines', /NORMAL  
  ;Display the image. 255b changes black values to white:  
  TVSCL, 255b - array, .1, .72, /NORMAL    
  stop
  ;Calculate and display the Hough transform  
  result = HOUGH(array, RHO=rho, THETA=theta)  
  XYOUTS, .1, .66, 'Hough Transform', /NORMAL  
  TVSCL, 255b - result, .1, .36, /NORMAL  
  
  ;Keep only lines that contain more than 20 points:  
  result = (result - 20) > 0  
  
  ;Find the Hough backprojection and display the output  
  backproject = HOUGH(result, /BACKPROJECT, RHO=rho, THETA=theta)  
  XYOUTS, .1, .30, 'Hough Backprojection', /NORMAL  
  TVSCL, 255b - backproject, .1, .08, /NORMAL  
  stop
END  