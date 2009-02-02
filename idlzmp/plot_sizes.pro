PRO PLOT_SIZES, file, win_vert_size, coord_sys

;ON_ERROR, 1

COMMON coordinates,  ni, nj, nk, x1a, x2a, x3a, x1b, x2b, x3b, dims


COMMON plot_params, chars_, xmarg_, ymarg_, color_bar_

COMMON plot_coords $
      , win_vsize, plot_h0, plot_v0, plot_h1, plot_v1 $
      , win_hsize_xy, win_vsize_xy,  plot_hsize_xy, plot_vsize_xy, position_xy $
      , win_hsize_xz, win_vsize_xz,  plot_hsize_xz, plot_vsize_xz, position_xz $
      , win_hsize_yz, win_vsize_yz,  plot_hsize_yz, plot_vsize_yz, position_yz 

;IF(color_bar_ EQ 'y') THEN BEGIN
  colbar = 1 
;ENDIF ELSE BEGIN
;  colbar=0
;ENDELSE

;READ_DIMS_AXES, file, x1b, x2b, x3b, ni, nj, nk
LOAD_DIMS_HDF, file, pdims=pdims, pcoords=pcoords, dims=dims, $
       nxd=nxd,nyd=nyd,nzd=nzd, nxb=nxb,nyb=nyb,nzb=nzb, nb=nb,$
       xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, zmin=zmin, zmax=zmax

ni = nxd
nj = nyd
nk = nzd
dx = (xmax - xmin) / nxd
dy = (ymax - ymin) / nyd
dz = (zmax - zmin) / nzd
print, dx,dy,dz
x1b = findgen(ni)
x2b = findgen(nj)
x3b = findgen(nk)

for i = 0,ni-1 do  x1b[i] = xmin + 0.5*dx + dx*x1b[i]
for i = 0,nj-1 do  x2b[i] = ymin + 0.5*dy + dy*x2b[i]
for i = 0,nk-1 do  x3b[i] = zmin + 0.5*dz + dz*x3b[i]

plot,x1b

IF(nj EQ 1 AND nk EQ 1) THEN dims= '1d' ELSE IF(nk EQ 1) THEN dims= '2d' ELSE dims= '3d' 


IF(dims EQ '2d') THEN BEGIN

  plot_h0 = FIX(chars_*!D.X_CH_SIZE*xmarg_(0))
  plot_v0 = FIX(chars_*!D.Y_CH_SIZE*ymarg_(0))
  plot_h1 = FIX(chars_*!D.X_CH_SIZE*xmarg_(1))
  plot_v1 = FIX(chars_*!D.Y_CH_SIZE*ymarg_(1))

  win_vsize_xy = win_vert_size 
  plot_vsize_xy = FIX(win_vsize_xy - chars_*!D.Y_CH_SIZE*(ymarg_(0)+ymarg_(1)))
  plot_hsize_xy = FIX(plot_vsize_xy*(MAX(x1b)-MIN(x1b))/(MAX(x2b)-MIN(x2b)))
  win_hsize_xy  = FIX(plot_hsize_xy + chars_*!D.X_CH_SIZE*(xmarg_(0)+xmarg_(1)))
  win_hsize_xy = win_hsize_xy + colbar*(plot_h0 + 20)

  position_xy = $
    [ plot_h0/FLOAT(win_hsize_xy) $
    , plot_v0/FLOAT(win_vsize_xy) $
    , plot_h0/FLOAT(win_hsize_xy) + plot_hsize_xy/FLOAT(win_hsize_xy) $
    , plot_v0/FLOAT(win_vsize_xy) + plot_vsize_xy/FLOAT(win_vsize_xy) ] 

ENDIF ELSE IF(dims EQ '3d') THEN BEGIN

  win_vsize = win_vert_size

  win_vsize_xz = win_vsize 
  win_vsize_yz = win_vsize 


  plot_vsize_xz = FIX(win_vsize_xz - chars_*!D.Y_CH_SIZE*(ymarg_(0)+ymarg_(1)))

  plot_hsize_xz = FIX(plot_vsize_xz*(MAX(x1b)-MIN(x1b))/(MAX(x3b)-MIN(x3b)))
  win_hsize_xz  = FIX(plot_hsize_xz + chars_*!D.X_CH_SIZE*(xmarg_(0)+xmarg_(1)))

  plot_vsize_yz = plot_vsize_xz  
  plot_hsize_yz = FIX(plot_vsize_yz*(MAX(x2b)-MIN(x2b))/(MAX(x3b)-MIN(x3b)))
  win_hsize_yz  = FIX(plot_hsize_yz + chars_*!D.X_CH_SIZE*(xmarg_(0)+xmarg_(1)))

  IF(coord_sys EQ 'zrp') THEN BEGIN
    plot_hsize_xy = FIX(plot_vsize_yz*(MAX(x2b)-MIN(x2b))/(MAX(x3b)-MIN(x3b)))
    plot_vsize_xy = FIX(plot_hsize_xy*(MAX(x1b)-MIN(x1b))/(MAX(x2b)-MIN(x2b)))
  ENDIF ELSE BEGIN  
    plot_hsize_xy = plot_hsize_xz
    plot_vsize_xy = FIX(plot_hsize_xy*(MAX(x2b)-MIN(x2b))/(MAX(x1b)-MIN(x1b)))
  ENDELSE
  
  win_hsize_xy  = FIX(plot_hsize_xy + chars_*!D.X_CH_SIZE*(xmarg_(0)+xmarg_(1)))
  win_vsize_xy  = FIX(plot_vsize_xy + chars_*!D.Y_CH_SIZE*(ymarg_(0)+ymarg_(1)))

  plot_h0 = FIX(chars_*!D.X_CH_SIZE*xmarg_(0))
  plot_v0 = FIX(chars_*!D.Y_CH_SIZE*ymarg_(0))
  plot_h1 = FIX(chars_*!D.X_CH_SIZE*xmarg_(1))
  plot_v1 = FIX(chars_*!D.Y_CH_SIZE*ymarg_(1))

  win_hsize_xz = win_hsize_xz + colbar*(plot_h0 + 20)
  win_hsize_yz = win_hsize_yz + colbar*(plot_h0 + 20)
  win_hsize_xy = win_hsize_xy + colbar*(plot_h0 + 20)

 
  position_xz = $
    [ plot_h0/FLOAT(win_hsize_xz) $
    , plot_v0/FLOAT(win_vsize_xz) $
    , plot_h0/FLOAT(win_hsize_xz) + plot_hsize_xz/FLOAT(win_hsize_xz) $
    , plot_v0/FLOAT(win_vsize_xz) + plot_vsize_xz/FLOAT(win_vsize_xz) ] 

  position_yz = $
    [ plot_h0/FLOAT(win_hsize_yz) $
    , plot_v0/FLOAT(win_vsize_yz) $
    , plot_h0/FLOAT(win_hsize_yz) + plot_hsize_yz/FLOAT(win_hsize_yz) $
    , plot_v0/FLOAT(win_vsize_yz) + plot_vsize_yz/FLOAT(win_vsize_yz) ] 

  position_xy = $
    [ plot_h0/FLOAT(win_hsize_xy) $
    , plot_v0/FLOAT(win_vsize_xy) $
    , plot_h0/FLOAT(win_hsize_xy) + plot_hsize_xy/FLOAT(win_hsize_xy) $
    , plot_v0/FLOAT(win_vsize_xy) + plot_vsize_xy/FLOAT(win_vsize_xy) ] 
    

ENDIF

END
