PRO PLOT_SIZES, data, win_vert_size, crdsys, psiz, plp

;ON_ERROR, 1

  attr = data.attr

  nx = attr.nxd
  ny = attr.nyd
  nz = attr.nzd

  x = data.axes.x
  y = data.axes.y
  z = data.axes.z

  psiz = {win_vsize:1, $
  	  plot_h0:1, $
	  plot_v0:1, $
	  plot_h1:1, $
	  plot_v1:1, $
      	  win_hsize_xy:1, $
	  win_vsize_xy:1, $
	  plot_hsize_xy:1, $
	  plot_vsize_xy:1, $
	  position_xy:FLTARR(4), $
      	  win_hsize_xz:1, $
	  win_vsize_xz:1, $
	  plot_hsize_xz:1, $
	  plot_vsize_xz:1, $
	  position_xz:FLTARR(4), $
      	  win_hsize_yz:1, $
	  win_vsize_yz:1, $
	  plot_hsize_yz:1, $
	  plot_vsize_yz:1, $
	  position_yz:FLTARR(4)}

  psiz.win_vsize = win_vert_size

  psiz.plot_h0 = FIX(plp.chars*!D.X_CH_SIZE*plp.xmarg(0))
  psiz.plot_v0 = FIX(plp.chars*!D.Y_CH_SIZE*plp.ymarg(0))
  psiz.plot_h1 = FIX(plp.chars*!D.X_CH_SIZE*plp.xmarg(1))
  psiz.plot_v1 = FIX(plp.chars*!D.Y_CH_SIZE*plp.ymarg(1))

IF(nx EQ 1) THEN BEGIN

  psiz.win_vsize_yz = win_vert_size
  psiz.plot_vsize_yz = FIX(psiz.win_vsize_yz - plp.chars*!D.Y_CH_SIZE*(plp.ymarg(0)+plp.ymarg(1)))
  psiz.plot_hsize_yz = FIX(psiz.plot_vsize_yz*(MAX(y)-MIN(y))/(MAX(z)-MIN(z)))
  psiz.win_hsize_yz  = FIX(psiz.plot_hsize_yz + plp.chars*!D.X_CH_SIZE*(plp.xmarg(0)+plp.xmarg(1)))
  psiz.win_hsize_yz = psiz.win_hsize_yz + plp.colbar*(psiz.plot_h0 + 20)

  psiz.position_yz = $
    [ psiz.plot_h0/FLOAT(psiz.win_hsize_yz) $
    , psiz.plot_v0/FLOAT(psiz.win_vsize_yz) $
    , psiz.plot_h0/FLOAT(psiz.win_hsize_yz) + psiz.plot_hsize_yz/FLOAT(psiz.win_hsize_yz) $
    , psiz.plot_v0/FLOAT(psiz.win_vsize_yz) + psiz.plot_vsize_yz/FLOAT(psiz.win_vsize_yz) ]

ENDIF ELSE IF(ny EQ 1) THEN BEGIN

  psiz.win_vsize_xz = win_vert_size
  psiz.plot_vsize_xz = FIX(psiz.win_vsize_xz - plp.chars*!D.Y_CH_SIZE*(plp.ymarg(0)+plp.ymarg(1)))
  psiz.plot_hsize_xz = FIX(psiz.plot_vsize_xz*(MAX(x)-MIN(x))/(MAX(z)-MIN(z)))
  psiz.win_hsize_xz  = FIX(psiz.plot_hsize_xz + plp.chars*!D.X_CH_SIZE*(plp.xmarg(0)+plp.xmarg(1)))
  psiz.win_hsize_xz = psiz.win_hsize_xz + plp.colbar*(psiz.plot_h0 + 20)

  psiz.position_xz = $
    [ psiz.plot_h0/FLOAT(psiz.win_hsize_xz) $
    , psiz.plot_v0/FLOAT(psiz.win_vsize_xz) $
    , psiz.plot_h0/FLOAT(psiz.win_hsize_xz) + psiz.plot_hsize_xz/FLOAT(psiz.win_hsize_xz) $
    , psiz.plot_v0/FLOAT(psiz.win_vsize_xz) + psiz.plot_vsize_xz/FLOAT(psiz.win_vsize_xz) ]

ENDIF ELSE IF(nz EQ 1) THEN BEGIN

  psiz.win_vsize_xy = win_vert_size
  psiz.plot_vsize_xy = FIX(psiz.win_vsize_xy - plp.chars*!D.Y_CH_SIZE*(plp.ymarg(0)+plp.ymarg(1)))
  psiz.plot_hsize_xy = FIX(psiz.plot_vsize_xy*(MAX(x)-MIN(x))/(MAX(y)-MIN(y)))
  psiz.win_hsize_xy  = FIX(psiz.plot_hsize_xy + plp.chars*!D.X_CH_SIZE*(plp.xmarg(0)+plp.xmarg(1)))
  psiz.win_hsize_xy = psiz.win_hsize_xy + plp.colbar*(psiz.plot_h0 + 20)

  psiz.position_xy = $
    [ psiz.plot_h0/FLOAT(psiz.win_hsize_xy) $
    , psiz.plot_v0/FLOAT(psiz.win_vsize_xy) $
    , psiz.plot_h0/FLOAT(psiz.win_hsize_xy) + psiz.plot_hsize_xy/FLOAT(psiz.win_hsize_xy) $
    , psiz.plot_v0/FLOAT(psiz.win_vsize_xy) + psiz.plot_vsize_xy/FLOAT(psiz.win_vsize_xy) ]

ENDIF ELSE IF(nx NE 1 AND ny NE 1 AND nz NE 1) THEN BEGIN

  psiz.win_vsize_xz = psiz.win_vsize
  psiz.win_vsize_yz = psiz.win_vsize


  psiz.plot_vsize_xz = FIX(psiz.win_vsize_xz - plp.chars*!D.Y_CH_SIZE*(plp.ymarg(0)+plp.ymarg(1)))

  psiz.plot_hsize_xz = FIX(psiz.plot_vsize_xz*(MAX(x)-MIN(x))/(MAX(z)-MIN(z)))
  psiz.win_hsize_xz  = FIX(psiz.plot_hsize_xz + plp.chars*!D.X_CH_SIZE*(plp.xmarg(0)+plp.xmarg(1)))

  psiz.plot_vsize_yz = psiz.plot_vsize_xz
  psiz.plot_hsize_yz = FIX(psiz.plot_vsize_yz*(MAX(y)-MIN(y))/(MAX(z)-MIN(z)))
  psiz.win_hsize_yz  = FIX(psiz.plot_hsize_yz + plp.chars*!D.X_CH_SIZE*(plp.xmarg(0)+plp.xmarg(1)))

  IF(crdsys EQ 'zrp') THEN BEGIN
    psiz.plot_hsize_xy = FIX(psiz.plot_vsize_yz*(MAX(y)-MIN(y))/(MAX(z)-MIN(z)))
    psiz.plot_vsize_xy = FIX(psiz.plot_hsize_xy*(MAX(x)-MIN(x))/(MAX(y)-MIN(y)))
  ENDIF ELSE BEGIN
    psiz.plot_hsize_xy = psiz.plot_hsize_xz
    psiz.plot_vsize_xy = FIX(psiz.plot_hsize_xy*(MAX(y)-MIN(y))/(MAX(x)-MIN(x)))
  ENDELSE

  psiz.win_hsize_xy  = FIX(psiz.plot_hsize_xy + plp.chars*!D.X_CH_SIZE*(plp.xmarg(0)+plp.xmarg(1)))
  psiz.win_vsize_xy  = FIX(psiz.plot_vsize_xy + plp.chars*!D.Y_CH_SIZE*(plp.ymarg(0)+plp.ymarg(1)))

  psiz.win_hsize_xz = psiz.win_hsize_xz + plp.colbar*(psiz.plot_h0 + 20)
  psiz.win_hsize_yz = psiz.win_hsize_yz + plp.colbar*(psiz.plot_h0 + 20)
  psiz.win_hsize_xy = psiz.win_hsize_xy + plp.colbar*(psiz.plot_h0 + 20)


  psiz.position_xz = $
    [ psiz.plot_h0/FLOAT(psiz.win_hsize_xz) $
    , psiz.plot_v0/FLOAT(psiz.win_vsize_xz) $
    , psiz.plot_h0/FLOAT(psiz.win_hsize_xz) + psiz.plot_hsize_xz/FLOAT(psiz.win_hsize_xz) $
    , psiz.plot_v0/FLOAT(psiz.win_vsize_xz) + psiz.plot_vsize_xz/FLOAT(psiz.win_vsize_xz) ]

  psiz.position_yz = $
    [ psiz.plot_h0/FLOAT(psiz.win_hsize_yz) $
    , psiz.plot_v0/FLOAT(psiz.win_vsize_yz) $
    , psiz.plot_h0/FLOAT(psiz.win_hsize_yz) + psiz.plot_hsize_yz/FLOAT(psiz.win_hsize_yz) $
    , psiz.plot_v0/FLOAT(psiz.win_vsize_yz) + psiz.plot_vsize_yz/FLOAT(psiz.win_vsize_yz) ]

  psiz.position_xy = $
    [ psiz.plot_h0/FLOAT(psiz.win_hsize_xy) $
    , psiz.plot_v0/FLOAT(psiz.win_vsize_xy) $
    , psiz.plot_h0/FLOAT(psiz.win_hsize_xy) + psiz.plot_hsize_xy/FLOAT(psiz.win_hsize_xy) $
    , psiz.plot_v0/FLOAT(psiz.win_vsize_xy) + psiz.plot_vsize_xy/FLOAT(psiz.win_vsize_xy) ]

ENDIF ELSE BEGIN

    PRINT, 'Plots available only in 2D or 3D'

ENDELSE

END
