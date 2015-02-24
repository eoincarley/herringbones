pro tcd_ne_model_20110922

  ;--------------------------------;
  ;				Read the data
  ;--------------------------------;
  cd,'~/Data/22sep2011_event/density_mag'
  restore,'cartesian_density_map_22_sep.sav', /verb
  files_c2 = findfile('*.fts')
  c2data = lasco_readfits(files_c2[0], hdr)
  rsun_arcsec = get_solar_radius(hdr)

  ;--------------------------------;
  ;				Plot the data
  ;--------------------------------;

  window,1
  plot_image, alog10(CAR_DEN_ALL.data)
  xcen = 300.0
  ycen = 300.0
  rhos = dindgen(250)

  ;---------------------------------------------------;
  ;	Get line profiles between two angles degrees lat
  ;---------------------------------------------------;
  start_angle= 180.0
  stop_angle = 225.0
  n_profiles = stop_angle - start_angle
  angles = ( findgen(n_profiles)*(stop_angle - start_angle)/(n_profiles-1) ) + start_angle
  prof_array = dblarr(n_profiles, 250)

  FOR i=0, n_profiles-1 DO BEGIN
    angle = angles[i]*!DTOR 
    xline = (COS(angle) * rhos + xcen)*1.0 ;working in pixel units ;hdr.cdelt1
    yline = (SIN(angle) * rhos + ycen)*1.0 ;working in pixel units ;hdr.cdelt1
    plots, xline, yline, color=255, thick=1
    line_profile = interpolate(car_den_all.data, xline, yline)
    prof_array[i, *] = line_profile
  ENDFOR
  freq = 75e6
  n = freq2dens(freq)

  window,3
  plot_image, alog10(prof_array)


  ;Remove the nans from the array
  remove_nans, prof_array, junk, junk, nan_pos
  prof_array[nan_pos[*]] = 0.0
   
  tcd_ne = average(prof_array, 1) 
  
  xarcsec = abs(xline - 300.0)*CAR_DEN_ALL.dx
  yarcsec = abs(yline - 300.0)*CAR_DEN_ALL.dy
  rsun_asec = get_solar_radius(hdr)
  xrsun = xarcsec/rsun_asec
  yrsun = yarcsec/rsun_asec
  tcd_rads = sqrt(xrsun^2.0 + yrsun^2.0)
  
  save, tcd_rads, tcd_ne, filename = 'tcd_model_20110922.sav'
  
  stop
END