PRO READ_DATA, file, vars, time

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


FOR ivar=0, nvar-1 DO BEGIN
  var = vars(ivar)
  READ_VARIABLE, file, var, data, time, status
  IF(status NE 0) THEN GOTO, SKIP

  IF(var EQ 'dd') THEN BEGIN
    dd = data
    min_dd =MIN(dd)
    max_dd = MAX(dd)
  ENDIF
  IF(var EQ 'ee') THEN BEGIN
    ee = data
    min_ee = MIN(ee)
    max_ee = MAX(ee)
  ENDIF
  IF(var EQ 'er') THEN BEGIN
    er = data
    min_er = MIN(er)
    max_er = MAX(er)
  ENDIF
  IF(var EQ 'ec') THEN BEGIN
    ec = data
    min_ec = MIN(ec)
    max_ec = MAX(ec)
  ENDIF

  IF(var EQ 'v1') THEN v1 = data
  IF(var EQ 'v2') THEN v2 = data
  IF(var EQ 'v3') THEN v3 = data
  IF(var EQ 'b1') THEN b1 = data
  IF(var EQ 'b2') THEN b2 = data
  IF(var EQ 'b3') THEN b3 = data

  SKIP:
ENDFOR

 
END
