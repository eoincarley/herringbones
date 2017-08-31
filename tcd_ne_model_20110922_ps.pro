pro setup_ps, name
  
  set_plot,'ps'
  !p.font=0
  !p.charsize=1.5
  device, filename = name, $
          /color, $
          /helvetica, $
          /inches, $
          bits_per_pixel = 16, $
          xsize=10, $
          ysize=10, $
          /encapsulate, $
          yoffset=5

end

pro tcd_ne_model_20110922_ps

  ;--------------------------------;
  ;				Read the data
  ;--------------------------------;
  cd,'~/Data/2011_sep_22/density_mag'
  restore,'cartesian_density_map_22_sep.sav', /verb
  files_c2 = findfile('*.fts')
  c2data = lasco_readfits(files_c2[0], hdr)
  rsun_arcsec = get_solar_radius(hdr)

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
   ; plots, xline, yline, color=255, thick=1
    line_profile = interpolate(car_den_all.data, xline, yline)
    prof_array[i, *] = line_profile
  ENDFOR
  freq = 75e6
  n = freq2dens(freq)


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
  
  ;save, tcd_rads, tcd_ne, filename = 'tcd_model_20110922.sav'
  
 
  ;--------------------------------;
  ;     Plot the in postscript
  ;--------------------------------;

  ;setup_ps, 'tcd_density_map_aa.eps'
    loadct, 5
    FOV = [120.0, 120.0]
    car_den_all.data = alog10(car_den_all.data)
    tit = '2011-Sep-22 00:00:00 UT';strsplit(anytim(hdr.date_obs, /cc, /trun), 'T', /extract)+' UT'
    plot_map, car_den_all, $
        fov = FOV, $
        title = tit, $
        xticklen = -0.02, $
        yticklen = -0.02, $
        dmin = 4, $
        dmax = 9

    loadct, 74
    plot_map, car_den_all, $
        fov = FOV, $
        title = ' ', $
        xticklen = -0.02, $
        yticklen = -0.02, $
        dmin = 4, $
        dmax = 9, $
        /noerase, $
        /noaxes

    set_line_color
    rhos = dindgen(150)
    ;suncenter = get_sun_center(hdr)
    xline = (COS(180*!DTOR) * rhos + 0)*hdr.cdelt1 ;working in pixel units ;hdr.cdelt1
    yline = (SIN(180*!DTOR) * rhos + 0)*hdr.cdelt2 ;working in pixel units ;hdr.cdelt1
    ;plots, xline, yline, color = 0, thick = 3

    xline = (COS(225*!DTOR) * rhos + 0)*hdr.cdelt1 ;working in pixel units ;hdr.cdelt1
    yline = (SIN(225*!DTOR) * rhos + 0)*hdr.cdelt2 ;working in pixel units ;hdr.cdelt1
    ;plots, xline, yline, color = 0, thick = 3

    tvcircle, 960.0, 0.0, 0.0, color=0, /data, /fill

    restore,'~/Data/2011_sep_22/density_mag/Density_map_disk_Zucca_et_al_20110922.sav', /verb
    plot_map, density_map, $
          /composite, $
          /average



    set_line_color
    plot_helio, hdr.date_obs, /over, gstyle=0, gthick=3.0, gcolor=1, grid_spacing=15.0

    data = (car_den_all.data)
    remove_nans, data, data_out;, /zero
    data_round = round(data_out*10.0)/10.0
    dens1 = alog10(freq2dens(32e6))
    dens1 = round(dens1*10.0)/10.0

    dens2 = alog10(freq2dens(43e6))
    dens2 = round(dens2*10.0)/10.0

    ind32 = where(data_round eq dens1)
    xy32 = array_indices(data, ind32)
    xy32 = (xy32 - 300.0)*car_den_all.dx 

    ind43 = where(data_round eq dens2)
    xy43 = array_indices(data, ind43)
    xy43 = (xy43 - 300.0)*car_den_all.dx 
    plotsym, 0, /fill
    plots, xy32[0, *], xy32[1, *], color=0, thick=1, psym=8, symsize=0.2
    plots, xy43[0, *], xy43[1, *], color=0, thick=1, psym=8, symsize=0.2
    

    loadct, 74
    cgColorbar, Range=[4.0, 9.0], $
        /vertical, $
        /right, $
        position = [0.88, 0.15, 0.89, 0.85], $
        ticklen=-0.35

    loadct, 0    
    xyouts, 0.94, 0.5, 'Electron Density log!L10!N(n!Le!N) (cm!U-3!N)', orient = 270, alignment = 0.5, /normal

  ;device, /close
  ;set_plot, 'x'


END