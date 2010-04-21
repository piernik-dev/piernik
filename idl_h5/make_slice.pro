PRO MAKE_SLICE, data, cs, slice_array, i_slice, n_vectors, sldata

;ON_ERROR, 1

s = slice_array(i_slice)

x = data.axes.x
y = data.axes.y
z = data.axes.z

nx = data.attr.nxd
ny = data.attr.nyd
nz = data.attr.nzd

IF(s.type EQ cs.plane3) THEN BEGIN
  mz = min((z-s.coord)^2,iz)

  IF(s.vect_disp EQ 'mag') THEN BEGIN
    u = reform(data.magx.arr(*,*,iz), nx,ny)
    v = reform(data.magy.arr(*,*,iz), nx,ny)
  ENDIF ELSE IF(STRMID(s.vect_disp,0,3) EQ 'vel') THEN BEGIN
    flind = STRMID(s.vect_disp,3,1)
    velx = 'vlx'+flind
    vely = 'vly'+flind
    ok=EXECUTE('u = REFORM(data.'+velx+'.arr(*,*,iz), nx,ny)')
    ok=EXECUTE('v = REFORM(data.'+vely+'.arr(*,*,iz), nx,ny)')
  ENDIF


    var_scd = s.scal_disp

    ok = Execute('scd = data.'+var_scd+'.arr') &
    ok = Execute('min_scd = data.'+var_scd+'.min')
    ok = Execute('max_scd = data.'+var_scd+'.max')

    sdis = reform(scd(*,*,iz), nx,ny)


    xdis = CONGRID(x,n_vectors(0),/interp)
    ydis = CONGRID(y,n_vectors(1),/interp)
    udis = congrid(u, n_vectors(0),n_vectors(1),/interp)
    vdis = congrid(v, n_vectors(0),n_vectors(1),/interp)

    xrange = [data.attr.xmin,data.attr.xmax]
    yrange = [data.attr.ymin,data.attr.ymax]

ENDIF ELSE IF(s.type EQ cs.plane1) THEN BEGIN
  mx = min((x-s.coord)^2,ix)

  IF(s.vect_disp EQ 'mag') THEN BEGIN
    u = reform(data.magy.arr(ix,*,*), ny,nz)
    v = reform(data.magz.arr(ix,*,*), ny,nz)
  ENDIF ELSE IF(STRMID(s.vect_disp,0,3) EQ 'vel') THEN BEGIN
    flind = STRMID(s.vect_disp,3,1)
    vely = 'vly'+flind
    velz = 'vlz'+flind
    ok=EXECUTE('u = REFORM(data.'+vely+'.arr(ix,*,*), ny,nz)')
    ok=EXECUTE('v = REFORM(data.'+velz+'.arr(ix,*,*), ny,nz)')
  ENDIF

    var_scd = s.scal_disp

    ok = Execute('scd = data.'+var_scd+'.arr') &
    ok = Execute('min_scd = data.'+var_scd+'.min')
    ok = Execute('max_scd = data.'+var_scd+'.max')

    sdis = reform(scd(ix,*,*), ny,nz)

  xdis = CONGRID(y,n_vectors(1),/interp)
  ydis = CONGRID(z,n_vectors(2),/interp)
  udis = congrid(u, n_vectors(1),n_vectors(2),/interp)
  vdis = congrid(v, n_vectors(1),n_vectors(2),/interp)

  xrange = [data.attr.ymin,data.attr.ymax]
  yrange = [data.attr.zmin,data.attr.zmax]

ENDIF ELSE IF(s.type EQ cs.plane2) THEN BEGIN
  my = min((y-s.coord)^2,iy)

  IF(s.vect_disp EQ 'mag') THEN BEGIN
    u = reform(data.magx.arr(*,iy,*), nx,nz)
    v = reform(data.magz.arr(*,iy,*), nx,nz)
  ENDIF ELSE IF(STRMID(s.vect_disp,0,3) EQ 'vel') THEN BEGIN
    flind = STRMID(s.vect_disp,3,1)
    velx = 'vlx'+flind
    velz = 'vlz'+flind
    ok=EXECUTE('u = REFORM(data.'+velx+'.arr(*,iy,*), nx,nz)')
    ok=EXECUTE('v = REFORM(data.'+velz+'.arr(*,iy,*), nx,nz)')
  ENDIF

    var_scd = s.scal_disp

    ok = Execute('scd = data.'+var_scd+'.arr') &
    ok = Execute('min_scd = data.'+var_scd+'.min')
    ok = Execute('max_scd = data.'+var_scd+'.max')

    sdis = reform(scd(*,iy,*), nx,nz)

  xdis = CONGRID(x,n_vectors(0),/interp)
  ydis = CONGRID(z,n_vectors(2),/interp)
  udis = congrid(u, n_vectors(0),n_vectors(2),/interp)
  vdis = congrid(v, n_vectors(0),n_vectors(2),/interp)

  xrange = [data.attr.xmin,data.attr.xmax]
  yrange = [data.attr.zmin,data.attr.zmax]


ENDIF

IF(ABS(xrange(0)) LT 0.01*ABS(xrange(1))) THEN xrange(0) = 0.0
IF(ABS(yrange(0)) LT 0.01*ABS(yrange(1))) THEN yrange(0) = 0.0
IF(ABS(xrange(1)) LT 0.01*ABS(xrange(0))) THEN xrange(1) = 0.0
IF(ABS(yrange(1)) LT 0.01*ABS(yrange(0))) THEN yrange(1) = 0.0


sldata = {xdis:xdis, $
          ydis:ydis, $
          udis:udis, $
          vdis:vdis, $
          sdis:sdis, $
          xrange:xrange, $
          yrange:yrange}

END
