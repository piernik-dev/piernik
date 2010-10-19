PRO WINDOWS_OPEN,  slice_array, display

COMMON coordinates,  ni, nj, nk, x1a, x2a, x3a, x1b, x2b, x3b, dims

COMMON coordsys    , coordsys, coordnames $
                   , coord1,coord2,coord3,plane1,plane2,plane3
COMMON plot_coords $
      , win_vsize, plot_h0, plot_v0, plot_h1, plot_v1 $
      , win_hsize_xy, win_vsize_xy,  plot_hsize_xy, plot_vsize_xy, position_xy $
      , win_hsize_xz, win_vsize_xz,  plot_hsize_xz, plot_vsize_xz, position_xz $
      , win_hsize_yz, win_vsize_yz,  plot_hsize_yz, plot_vsize_yz, position_yz $
      , win_hsize_b1, win_vsize_b1,  plot_hsize_b1, plot_vsize_b1, position_b1


  n_slices = N_ELEMENTS(slice_array) - 1

  FOR i_slice=1, n_slices DO BEGIN
    s = slice_array(i_slice)
    IF(s.sw EQ 'on') THEN BEGIN


      IF(dims EQ '2d') THEN BEGIN
        IF(s.type EQ plane3) THEN BEGIN
          resol=[win_hsize_xy,win_vsize_xy]
         IF(display EQ 'y') THEN BEGIN
            WINDOW,i_slice,xsize=resol[0],ysize=resol[1], $
                             title = s.name
         ENDIF
        ENDIF


      ENDIF ELSE IF(dims EQ '3d') THEN BEGIN
        IF(s.type EQ plane3) THEN resol=[win_hsize_xy,win_vsize_xy]
        IF(s.type EQ plane2) THEN resol=[win_hsize_xz,win_vsize_xz]
        IF(s.type EQ plane1) THEN resol=[win_hsize_yz,win_vsize_yz]
        IF(display EQ 'y') THEN BEGIN
          WINDOW,i_slice,xsize=resol[0],ysize=resol[1], $
                  title = s.name
        ENDIF
      ENDIF

     ENDIF

  END

END
