pro diffuse_emission_check

    ; Look at the nature of the diffuse emission near herringbones
    ; on 2011 sep 22.

    loadct, 74, /silent
    reverse_ct
    window, 1
    !p.charsize = 1

    cd, '~/Data/2011_sep_22/herringbones'
    radio_spectro_fits_read, 'BIR_20110922_104459_01.fit', data_raw, times, freq
    data_bs = constbacksub(data_raw, /auto)
    data_bs = (data_bs - mean(data_bs))/stdev(data_bs)

    t1_index = closest(times,anytim(file2time('20110922_104800'),/utim))
    t2_index = closest(times,anytim(file2time('20110922_105000'),/utim))
    f1_index = closest(freq, 80.0)
    f2_index = closest(freq, 30.0)

    spectro_plot, smooth(data_bs, 1) > (-1) < 3, times, freq, $
        /ys, $
        ytitle = '!6Frequency [MHz]', $
        yticks = 5, $
        yminor = 4, $
        yr = [freq[f1_index],freq[f2_index]], $
        xrange = [times[t1_index],times[t2_index]], $
        /xs, $
        xtitle = 'Start time: '+anytim(times[t1_index], /cc, /trun)+' (UT)'

END        