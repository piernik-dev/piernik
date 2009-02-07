PRO MAKE_SLICE, slice_array, i_slice, n_vectors

;ON_ERROR, 1

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

COMMON data_min_max, min_dd0, max_dd0, min_ee0, max_ee0 $
                   , min_dd1, max_dd1, min_ee1, max_ee1, min_jj, max_jj $
                   , min_ecr0, max_ecr0, min_ecr1, max_ecr1  

COMMON mag_field, b1_mean, b2_mean, bh_mean, b1_m 

;COMMON slices, slice_array, i_slice

COMMON sldata  $
       , xdis, ydis, udis, vdis, sdis  $
       , xrange_, yrange_ $ 
       , time_str, slice_str  


s = slice_array(i_slice)
IF(s.type EQ plane1) THEN begin
   hdx = (x2b(1)-x2b(0))/2.
   hdy = (x3b(1)-x3b(0))/2.
endif else begin
   hdx = (x1b(1)-x1b(0))/2.
   hdy = (x2b(1)-x2b(0))/2.
endelse
IF(dims EQ '3d') THEN hdz = (x3b(1)-x3b(0))/2.


IF(s.type EQ plane3) THEN BEGIN
  mz = min((x3b-s.coord)^2,iz)

  IF(s.vect_disp EQ 'b') THEN BEGIN
    u = reform(b1(*,*,iz), ni,nj)
    v = reform(b2(*,*,iz), ni,nj)
  ENDIF ELSE IF(s.vect_disp EQ 'v') THEN BEGIN
    u = reform(v1(*,*,iz), ni,nj) 
    v = reform(v2(*,*,iz), ni,nj)
  ENDIF

  IF(s.scal_disp EQ 'd') THEN BEGIN
    IF(s.scal_pert EQ '') THEN sdis = reform(dd(*,*,iz), ni,nj) 
    IF(s.scal_pert EQ 'p') THEN sdis = reform(dd1(*,*,iz), ni,nj) 
  ENDIF ELSE IF(s.scal_disp EQ 'e') THEN BEGIN
    IF(s.scal_pert EQ '') THEN sdis = reform(ee(*,*,iz), ni,nj) 
    IF(s.scal_pert EQ 'p') THEN sdis = reform(ee1(*,*,iz), ni,nj) 
  ENDIF ELSE IF(s.scal_disp EQ 'cr') THEN BEGIN
    IF(s.scal_pert EQ '') THEN sdis = reform(ecr(*,*,iz), ni,nj) 
    IF(s.scal_pert EQ 'p') THEN sdis = reform(ecr1(*,*,iz), ni,nj) 
  ENDIF ELSE IF(s.scal_disp EQ 'j') THEN BEGIN
;    sdis = reform(jj(*,*,iz), ni,nj) 
    ASCURL,x1b,x2b,x3b,b1,b2,b3,cur,s.type,iz
    sdis=jj
    min_jj = 0.0
    max_jj = MAX(jj(5:ni-6,5:nj-6))
  ENDIF
  
  
    xdis = CONGRID(x1b,n_vectors(0),/interp)
    ydis = CONGRID(x2b,n_vectors(1),/interp) 
    udis = congrid(u, n_vectors(0),n_vectors(1),/interp)
    vdis = congrid(v, n_vectors(0),n_vectors(1),/interp)

    xrange_ = [MIN(x1b)-hdx,MAX(x1b)+hdx]
    yrange_ = [MIN(x2b)-hdy,MAX(x2b)+hdy]

  IF((s.type EQ 'zr') AND (dims EQ '3d')) THEN BEGIN  
    xrange_ = [MIN(x2b)-hdy,MAX(x2b)+hdy]
    yrange_ = [MIN(x1b)-hdx,MAX(x1b)+hdx]
    sdis = ROTATE(sdis,3)
    xdis = CONGRID(x2b,n_vectors(1),/interp)
    ydis = CONGRID(x1b,n_vectors(0),/interp) 
    ydis = ROTATE(ydis,2)
    udis = congrid(v, n_vectors(0),n_vectors(1),/interp)
    vdis = congrid(u , n_vectors(0),n_vectors(1),/interp)
    udis = ROTATE(udis,3)
    vdis = ROTATE(vdis,3)
  ENDIF

ENDIF ELSE IF(s.type EQ plane1) THEN BEGIN
  mx = min((x1b-s.coord)^2,ix)

  IF(s.vect_disp EQ 'b') THEN BEGIN
    u = reform(b2(ix,*,*), nj,nk)
    v = reform(b3(ix,*,*), nj,nk)
  ENDIF ELSE IF(s.vect_disp EQ 'v') THEN BEGIN
    u = reform(v2(ix,*,*), nj,nk) 
    v = reform(v3(ix,*,*), nj,nk)
  ENDIF

  IF(s.scal_disp EQ 'd') THEN BEGIN
    IF(s.scal_pert EQ '') THEN sdis = reform(dd(ix,*,*), nj,nk) 
    IF(s.scal_pert EQ 'p') THEN sdis = reform(dd1(ix,*,*), nj,nk) 
  ENDIF ELSE IF(s.scal_disp EQ 'e') THEN BEGIN
    IF(s.scal_pert EQ '') THEN sdis = reform(ee(ix,*,*), nj,nk) 
    IF(s.scal_pert EQ 'p') THEN sdis = reform(ee1(ix,*,*), nj,nk) 
  ENDIF ELSE IF(s.scal_disp EQ 'cr') THEN BEGIN
    IF(s.scal_pert EQ '') THEN sdis = reform(ecr(ix,*,*), nj,nk) 
    IF(s.scal_pert EQ 'p') THEN sdis = reform(ecr1(ix,*,*), nj,nk) 
  ENDIF ELSE IF(s.scal_disp EQ 'j') THEN BEGIN
    ;sdis = reform(jj(ix,*,*), nj,nk) 
    ASCURL,x1b,x2b,x3b,b1,b2,b3,cur,s.type,ix
    sdis=jj
    min_jj = 0.0
    max_jj = MAX(jj(5:ni-6,5:nk-6))
ENDIF
  
  
  xdis = CONGRID(x2b,n_vectors(1),/interp)
  ydis = CONGRID(x3b,n_vectors(2),/interp) 
  udis = congrid(u, n_vectors(1),n_vectors(2),/interp)
  vdis = congrid(v, n_vectors(1),n_vectors(2),/interp)

  xrange_ = [MIN(x2b)-hdy,MAX(x2b)+hdy]
  yrange_ = [MIN(x3b)-hdz,MAX(x3b)+hdz]

ENDIF ELSE IF(s.type EQ plane2) THEN BEGIN
  my = min((x2b-s.coord)^2,iy)

  IF(s.vect_disp EQ 'b') THEN BEGIN
    u = reform(b1(*,iy,*), ni,nk)
    v = reform(b3(*,iy,*), ni,nk)
  ENDIF ELSE IF(s.vect_disp EQ 'v') THEN BEGIN
    u = reform(v1(*,iy,*), ni,nk) 
    v = reform(v3(*,iy,*), ni,nk)
  ENDIF

  IF(s.scal_disp EQ 'd') THEN BEGIN
    IF(s.scal_pert EQ '') THEN sdis = reform(dd(*,iy,*), ni,nk) 
    IF(s.scal_pert EQ 'p') THEN sdis = reform(dd1(*,iy,*), ni,nk) 
  ENDIF ELSE IF(s.scal_disp EQ 'e') THEN BEGIN
    IF(s.scal_pert EQ '') THEN sdis = reform(ee(*,iy,*), ni,nk) 
    IF(s.scal_pert EQ 'p') THEN sdis = reform(ee1(*,iy,*), ni,nk) 
  ENDIF ELSE IF(s.scal_disp EQ 'cr') THEN BEGIN
    IF(s.scal_pert EQ '') THEN sdis = reform(ecr(*,iy,*), ni,nk) 
    IF(s.scal_pert EQ 'p') THEN sdis = reform(ecr1(*,iy,*), ni,nk) 
  ENDIF ELSE IF(s.scal_disp EQ 'j') THEN BEGIN
    ;sdis = reform(jj(*,iy,*), ni,nk) 
    ASCURL,x1b,x2b,x3b,b1,b2,b3,jj,s.type,iy
    sdis=jj
    min_jj = 0.0
    max_jj = MAX(jj(5:ni-6,5:nk-6))
  ENDIF
  xdis = CONGRID(x1b,n_vectors(0),/interp)
  ydis = CONGRID(x3b,n_vectors(2),/interp) 
  udis = congrid(u, n_vectors(0),n_vectors(2),/interp)
  vdis = congrid(v, n_vectors(0),n_vectors(2),/interp)

  xrange_ = [MIN(x1b)-hdx,MAX(x1b)+hdx]
  yrange_ = [MIN(x3b)-hdz,MAX(x3b)+hdz]


;  IF(MAX(udis)^2+MAX(vdis)^2 LT 1.0e-6) THEN BEGIN
;    vdis(0,0) = -1.0e-6 
;  ENDIF

ENDIF 

IF(ABS(xrange_(0)) LT 0.01*ABS(xrange_(1))) THEN xrange_(0) = 0.0
IF(ABS(yrange_(0)) LT 0.01*ABS(yrange_(1))) THEN yrange_(0) = 0.0
IF(ABS(xrange_(1)) LT 0.01*ABS(xrange_(0))) THEN xrange_(1) = 0.0
IF(ABS(yrange_(1)) LT 0.01*ABS(yrange_(0))) THEN yrange_(1) = 0.0

END
