pro determine_detection_theshold_all_sd
	
	; Simple procedure to imput multipe standard deviations to see when the Hough algorithm 
	; reliably detects a simulated burst.

	sd = dindgen(30)*(3.5 - 0.0)/29

	for i=0, n_elements(sd)-1 do begin 
		determine_detection_threshold2, 222, 225, sd[i]
		wait, 5
	endfor	

END