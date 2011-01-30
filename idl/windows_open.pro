PRO WINDOWS_OPEN,  slice_array, ps, cs, display

  n_slices = N_ELEMENTS(slice_array) - 1

  FOR i_slice=1, n_slices DO BEGIN
    s = slice_array(i_slice)
    IF(s.sw EQ 'on' AND s.status EQ 1) THEN BEGIN
        IF(s.type EQ cs.plane3) THEN resol=[ps.win_hsize_xy,ps.win_vsize_xy]
        IF(s.type EQ cs.plane2) THEN resol=[ps.win_hsize_xz,ps.win_vsize_xz]
        IF(s.type EQ cs.plane1) THEN resol=[ps.win_hsize_yz,ps.win_vsize_yz]
        IF(display EQ 'y') THEN BEGIN
          WINDOW,i_slice,xsize=resol[0],ysize=resol[1], $
                  title = s.name
        ENDIF
     ENDIF
  ENDFOR

END
