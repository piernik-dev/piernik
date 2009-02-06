PRO READ_DATA, dir,prefix,step, vars, time

; Reads data from user text file and HDF file specified by parameters

; COMMON parameters,   nhy, time $
;                   , gamm, eta1, eta2, jcrit


  COMMON coordinates,  nx1, nx2, nx3, x1a, x2a, x3a, x1b, x2b, x3b, dims
  COMMON quantites,    dd, ee, gp, v1, v2, v3, b1, b2, b3, ecr
  COMMON labels,       x1alab, x2alab, x3alab, x1blab, x2blab, x3blab  $
                     , ddlab, eelab, gplab, v1lab                      $
                     , v2lab, v3lab, b1lab, b2lab, b3lab, ecrlab

COMMON data_min_max, min_dd, max_dd, min_dd0, max_dd0, min_dd1, max_dd1 $
                   , min_ee, max_ee, min_ee0, max_ee0, min_ee1, max_ee1 $
                   , min_ecr, max_ecr, min_ecr0, max_ecr0, min_ecr1, max_ecr1 $ 
                   , min_jj, max_jj
 

  
  COMMON whatisshown, den, ene, grav, velo, magf, crene, cur

  nvar = n_elements(vars)

;print, 'file= ', file

FOR ivar=0, nvar-1 DO BEGIN
  var = vars(ivar)
  data = LOAD_DATA_HDF(dir,prefix, step, var, $
            xcoord=x, ycoord = y, zcoord = z, nxa=nxa,nya=nya,nza=nza, $
            time = time,allblocks='y')
  status = 0
;  READ_VARIABLE, file, var, data, time, status
  IF(status NE 0) THEN GOTO, SKIP

  IF(var EQ 'den1') THEN BEGIN
    dd = data
    min_dd = MIN(dd)
    max_dd = MAX(dd)
  ENDIF
  IF(var EQ 'den2') THEN BEGIN
    dd = data
    min_dd = MIN(dd)
    max_dd = MAX(dd)
  ENDIF
  IF(var EQ 'ene1') THEN BEGIN
    ee = data
    min_ee = MIN(ee)
    max_ee = MAX(ee)
  ENDIF
  IF(var EQ 'ene2') THEN BEGIN
    ee = data
    min_ee = MIN(ee)
    max_ee = MAX(ee)
  ENDIF
  IF(var EQ 'er') THEN BEGIN
    er = data
    min_er = MIN(er)
    max_er = MAX(er)
  ENDIF
  IF(var EQ 'encr') THEN BEGIN
    ee = data
    min_ee = MIN(ee)
    max_ee = MAX(ee)
  ENDIF
  
  IF(var EQ 'curz') THEN BEGIN
     dd = data
     min_dd = MIN(dd)
     max_dd = MAX(dd)
  endif
  
  IF(var EQ 'vlx1') THEN v1 = data
  IF(var EQ 'vly1') THEN v2 = data
  IF(var EQ 'vlz1') THEN v3 = data
  IF(var EQ 'vlx2') THEN v1 = data
  IF(var EQ 'vly2') THEN v2 = data
  IF(var EQ 'vlz2') THEN v3 = data
  IF(var EQ 'magx') THEN b1 = data
  IF(var EQ 'magy') THEN b2 = data
  IF(var EQ 'magz') THEN b3 = data

  SKIP:
ENDFOR

 
END
