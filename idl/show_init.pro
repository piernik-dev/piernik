PRO SHOW_INIT, s, slice_array, n_slices

COMMON parameters,   nhy, time $
                   , gamma, eta1, eta2, jcrit

COMMON coordsys    , coordsys, coordnames $
                   , coord1,coord2,coord3,plane1,plane2,plane3

COMMON coordinates,  ni, nj, nk, x1a, x2a, x3a, x1b, x2b, x3b, dims
COMMON quantites,    dd, ee, gp, v1, v2, v3, b1, b2, b3, ecr
COMMON labels,       x1alab, x2alab, x3alab, x1blab, x2blab, x3blab  $
                     , ddlab, eelab, gplab, v1lab                      $
                     , v2lab, v3lab, b1lab, b2lab, b3lab,ecrlab

COMMON quantites2,    dd1, ee1, ecr1, jj

COMMON files, fname,fext

COMMON whatisshown, den, ene, grav, velo, magf, crene, cur

COMMON data_min_max, min_dd, max_dd, min_dd0, max_dd0, min_dd1, max_dd1 $
                   , min_ee, max_ee, min_ee0, max_ee0, min_ee1, max_ee1 $
                   , min_ecr, max_ecr, min_ecr0, max_ecr0, min_ecr1, max_ecr1 $
                   , min_jj, max_jj

COMMON mag_field, b1_mean, b2_mean, bh_mean, b1_m

COMMON sldata  $
       , xdis, ydis, udis, vdis, sdis  $
       , xrange_, yrange_ $
       , time_str, slice_str $
       , lunit, tunit, bunit, vunit, dunit, eunit, ecrunit


COMMON plot_params, chars_, xmarg_, ymarg_, color_bar_

COMMON plot_coords $
      , win_vsize, plot_h0, plot_v0, plot_h1, plot_v1 $
      , win_hsize_xy, win_vsize_xy,  plot_hsize_xy, plot_vsize_xy, position_xy $
      , win_hsize_xz, win_vsize_xz,  plot_hsize_xz, plot_vsize_xz, position_xz $
      , win_hsize_yz, win_vsize_yz,  plot_hsize_yz, plot_vsize_yz, position_yz $
      , win_hsize_b1, win_vsize_b1,  plot_hsize_b1, plot_vsize_b1, position_b1


COMMON disp, ncolors, min_data, max_data, mincolor, maxcolor

COMMON anim, animation_name,nframes

COMMON endian, swap_endian

COMMON mode, batchmode
COMMON dumptype, dump

;COMMON first, iframe
;iframe = 0

;------------------------------------------------------------------------------
s = {name:      ' ', sw:       ' ',$
     type:      ' ', coord:        0.0, panel_name: ' ',     $
     vect_disp: ' ', vect_scaling: ' ', vect_scale: 1.0,     $
     scal_disp: ' ', scal_pert:    ' ',scal_scaling:' ',     $
     scal_log:  ' ', scal_scale:fltarr(2)}
slice_array = [s]
n_slices = 0
;------------------------------------------------------------------------------


swap_endian   = 'n'
data_type = 'hdf1'   ; 'usr', 'hdf1' (ZeusMP) or 'hdf2' (Zeus3D)


                     ; Specify the directory with data files
batchmode = 0
dump = 'hdfdump'


IF(!D.NAME NE 'X') THEN BEGIN
  device,/close
  set_plot,'X'
ENDIF

;WINDOW, /free, /pixmap, colors=-10
;WDELETE, !d.window
;DEVICE, get_visual_depth=depth
;print, 'Display depth:    ', depth
;print, 'Color table size: ', !d.table_size

color_table = 39
ncolors =!d.table_size-1
LOADCT, color_table
tvlct,red,green,blue,/get

chars_ = 0.8
xmarg_ = [16,3]
ymarg_ = [4,3]
color_bar_ = 'y'


lunit=''          ;'[pc]'
tunit=''          ;'[Myr]'
vunit=''          ;'kms!U-1!N'
bunit=''          ;'[!4l!3G]'





END
