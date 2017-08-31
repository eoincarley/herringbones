pro setup_ps, name
  
  set_plot,'ps'
  !p.font=0
  !p.charsize=1.5
  !p.charthick=0.5
  !p.thick = 0.5
  device, filename = name, $
          /color, $
          /helvetica, $
          /inches, $
          xsize=7, $
          ysize=7, $
          /encapsulate, $
          yoffset=5

end

pro plot_burst_surface

	; Followinf save file produced using intensity_over_burst.pro
	cd, '~/Data/2011_sep_22/herringbones/'
	restore, 'single_burst_suface.sav'
	white_space = (dindgen(20)*(5 - 0.7)/19. ) + 0.7
	set_line_color

	FOR i = 0, n_elements(fs_array)-1 DO BEGIN
       	
       	times = ts_array[i, *]
		inten = burst_surface[i, *]

		if i eq 0 then begin
			setup_ps, 'burst_surface.eps'
			tim_range = [ times[0], times[n_elements(times)-1]+4.0 ]

			utplot, times, inten - i, $
					thick=1, $
					yr = [-140, 20], $
					/ys, $
					/xs, $
					ytickname = [' ', ' ', ' ', ' '], $
					xr = tim_range, $
					xtickname = ['2', '4', '6'], $
					tick_unit = 2.0, $
					position = [0.1, 0.1, 0.9, 0.9]

			
			 
			 
		endif else begin
			utplot, times, inten - i*2.0, $
				thick=1, $
			 	yr = [-140, 20], $
			 	xr = tim_range, $
			 	/ys, $
			 	/xs, $
			 	ytickname = [' ', ' ', ' ', ' '], $
			 	/noerase, $
			 	xtickname = ['2', '4', '6'], $
				tick_unit = 2.0, $
				position = [0.1, 0.1, 0.9, 0.9]
			
		endelse	 

		if i ne n_elements(fs_array)-1 then begin
			for j = 0, n_elements(white_space)-1 do begin

				utplot, times, inten - i*2.0 - white_space[j], $
						thick=1, $
						yr = [-140, 20], $
						/ys, $
						/xs, $
						xr = tim_range, $
						ytickname = [' ', ' ', ' ', ' '], $
						color = 1, $
						/noerase, $
						xtickname = ['2', '4', '6'], $
						tick_unit = 2.0, $
						position = [0.1, 0.1, 0.9, 0.9]


			endfor
		endif
		
		if i mod 4 eq 0.0 and fs_array[i] gt 45.0 then xyouts, $
			times[n_elements(times)-1] + 0.05, inten[n_elements(inten)-1] - i*2.0, $
			string(fs_array[i], format = '(f4.1)')+' MHz', $
			charsize=1.3

		;plot intensity scale	
		if i eq 6 then begin
			outplot, [times[0], times[0]], [-130, -120], thick=4		
			outplot, [times[0], times[0]], [-130, -120], psym=1, thick=4
		endif
	ENDFOR

	device, /close
	set_plot, 'x'
END