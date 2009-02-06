PRO PLOT_SLICES,  slice_array, n_vectors, $
                  display_frames, hdf_num, png_output, png_dir, time,step

COMMON coordsys    , coordsys, coordnames $
                   , coord1,coord2,coord3,plane1,plane2,plane3
COMMON coordinates,  ni, nj, nk, x1a, x2a, x3a, x1b, x2b, x3b, dims
COMMON plot_coords $
      , win_vsize, plot_h0, plot_v0, plot_h1, plot_v1 $
      , win_hsize_xy, win_vsize_xy,  plot_hsize_xy, plot_vsize_xy, position_xy $
      , win_hsize_xz, win_vsize_xz,  plot_hsize_xz, plot_vsize_xz, position_xz $
      , win_hsize_yz, win_vsize_yz,  plot_hsize_yz, plot_vsize_yz, position_yz $
      , win_hsize_b1, win_vsize_b1,  plot_hsize_b1, plot_vsize_b1, position_b1  
COMMON frames, first_frame, i_frame

  n_slices = N_ELEMENTS(slice_array)-1

  FOR i_slice=1, n_slices DO BEGIN
    s = slice_array(i_slice)
    IF(s.sw EQ 'on') THEN BEGIN

      IF(s.type EQ plane3) THEN resol=[win_hsize_xy,win_vsize_xy]

      IF((s.type EQ plane1) OR (s.type EQ plane2) AND (dims EQ '2d')) THEN GOTO, SKIP

      IF(s.type EQ plane2) THEN resol=[win_hsize_xz,win_vsize_xz]
      IF(s.type EQ plane1) THEN resol=[win_hsize_yz,win_vsize_yz]


      IF(display_frames EQ 'y') THEN BEGIN
        WSHOW, i_slice
        WSET, i_slice
      ENDIF
    
      MAKE_SLICE, slice_array, i_slice, n_vectors

      SET_PLOT, 'Z' 
      DEVICE, SET_RESOLUTION=resol,SET_COLORS=ncolors
      SHOW_SLICE, slice_array, i_slice, time

      frame = TVRD() 

      IF(display_frames EQ 'y') THEN BEGIN
        SET_PLOT,'X'
        TV, frame
      ENDIF

      IF(png_output EQ 'y') THEN BEGIN
        TVLCT,red,green,blue,/GET
        pngfile = png_dir+'/'+hdf_num+'_'+s.name+'_'+string(step,format='(I3.3)')+'.png'
        WRITE_PNG,  pngfile,frame,red,green,blue, /verbose
      ENDIF
     ENDIF
  END

SKIP:

END
