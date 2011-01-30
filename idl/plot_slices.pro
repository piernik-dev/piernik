PRO PLOT_SLICES,  data, slice_array, ps, cs, plp, disp_par, units, n_vectors, $
                  display_frames, step, png_output, png_dir, win_open, first_call

  time = data.attr.time

  nx  = data.attr.nxd
  ny  = data.attr.nyd
  nz  = data.attr.nzd

  vars = data.vars

; CHANGING UNITS (should be done some place else)

  vars_mag = WHERE(STRCMP(strmid(vars,0,3),'mag',/FOLD_CASE))
  IF(N_ELEMENTS(vars_mag) EQ 3) THEN BEGIN
     data.magx.arr=0.459*data.magx.arr
     data.magy.arr=0.459*data.magy.arr
     data.magz.arr=0.459*data.magz.arr
  ENDIF

  n_slices = N_ELEMENTS(slice_array)-1

  FOR i_slice=1, n_slices DO BEGIN
    s = slice_array(i_slice)

    ivar_valid = WHERE(s.scal_disp EQ data.vars)
    IF(ivar_valid(0) GT -1) THEN BEGIN
       slice_array(i_slice).status=1
    ENDIF ELSE BEGIN
       PRINT, 'WARNING: "'+s.scal_disp+'" in panel ('+s.panel_name+') is not a valid variable to display'
    ENDELSE

    IF(s.type EQ'xy' AND (nx EQ 1 OR ny EQ 1)) THEN BEGIN
       slice_array(i_slice).status=0
       PRINT, 'WARNING: "'+s.type+'"   in panel ('+s.panel_name+') is not a valid slice type to display'
    ENDIF

    IF(s.type EQ'xz' AND (nx EQ 1 OR nz EQ 1)) THEN BEGIN
       slice_array(i_slice).status=0
       PRINT, 'WARNING: "'+s.type+'"   in panel ('+s.panel_name+') is not a valid slice type to display'
    ENDIF

    IF(s.type EQ'yz' AND (ny EQ 1 OR nz EQ 1)) THEN BEGIN
       slice_array(i_slice).status=0
       PRINT, 'WARNING: "'+s.type+'"   in panel ('+s.panel_name+') is not a valid slice type to display'
    ENDIF
  END
  PRINT, ' '

  IF(first_call EQ 'y') THEN BEGIN
    WINDOWS_OPEN,  slice_array,  ps, cs, display_frames
    first_call = 'n'
  ENDIF


  FOR i_slice=1, n_slices DO BEGIN
    s = slice_array(i_slice)

    IF(s.sw EQ 'on' AND s.status EQ 1) THEN BEGIN

      IF(s.type EQ cs.plane3) THEN resol=[ps.win_hsize_xy,ps.win_vsize_xy]
      IF(s.type EQ cs.plane2) THEN resol=[ps.win_hsize_xz,ps.win_vsize_xz]
      IF(s.type EQ cs.plane1) THEN resol=[ps.win_hsize_yz,ps.win_vsize_yz]


      IF(display_frames EQ 'y') THEN BEGIN
        WSHOW, i_slice
        WSET, i_slice
      ENDIF

      MAKE_SLICE, data, cs, slice_array, i_slice, n_vectors, sldata

      SET_PLOT, 'Z'
      DEVICE, SET_RESOLUTION=resol,SET_COLORS=ncolors
     !P.BACKGROUND = disp_par.backg
     !P.COLOR = disp_par.color



      SHOW_SLICE, data, cs, plp, disp_par,  units, slice_array, ps, i_slice, time, sldata

      frame = TVRD()

      IF(display_frames EQ 'y') THEN BEGIN
        SET_PLOT,'X'
        TV, frame
      ENDIF

      IF(png_output EQ 'y') THEN BEGIN
        TVLCT,red,green,blue,/GET
        pngfile = png_dir+'/'+s.name+'_'+string(step,format='(I4.4)')+'.png'
        WRITE_PNG,  pngfile,frame,red,green,blue, /verbose
      ENDIF
     ENDIF
  END

SKIP:

END
