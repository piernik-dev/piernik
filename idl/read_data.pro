PRO READ_DATA, dir, prefix, iframe, data, noarrays=noarrays

; Reads data from  HDF5 file
; The filename is constructed on the base of input variables: dir, prefix
; and step
;
; The output 'data' structure contains all data attributes, x,y,z-grids,
; variable names and data arrays together with  minima and and maxima
; of all variables.
;

  file = dir + '/' + prefix + '_'+STRING(FORMAT = '((I4.4))', iframe) + '.h5'

  nvar = n_elements(vars)

;READING ATTRIBUTES

  attr = read_attr_h5(file)
  data = CREATE_STRUCT("attr",attr)

  time = attr.time
  nxd  = attr.nxd
  nyd  = attr.nyd
  nzd  = attr.nzd

;CREATING AXES

  dx = (attr.xmax-attr.xmin)/float(attr.nxd)
  dy = (attr.ymax-attr.ymin)/float(attr.nyd)
  dz = (attr.zmax-attr.zmin)/float(attr.nzd)

  x = attr.xmin + dx*(0.5+FINDGEN(attr.nxd))
  y = attr.ymin + dy*(0.5+FINDGEN(attr.nyd))
  z = attr.zmin + dz*(0.5+FINDGEN(attr.nzd))

  axes = {x:x,y:y,z:z,dx:dx,dy:dy,dz:dz}

  data = CREATE_STRUCT("axes",axes,data)

; READING 3D VARIABLES

  start = [0,0,0]
  count = [nxd,nyd,nzd]

  h5vars = h5datasets(file)
  n_h5vars=n_elements(h5vars)

  data = CREATE_STRUCT("vars",h5vars,data)

  IF KEYWORD_SET(noarrays) THEN BEGIN
    PRINT, ''
    PRINT, 'READING FILE: ',file, ', time = ', strtrim(string(time),1)
    PRINT, 'DATA ATTRIBUTES, AXES AND VARNAMES'
    RETURN
  ENDIF

  PRINT, ''
  PRINT, 'READING FILE: ',file, ', time = ', strtrim(string(time),1)
  PRINT,  strtrim(n_h5vars-1,1), ' VARIABLES:',h5vars
  PRINT, ''

  IF KEYWORD_SET(noarrays) THEN RETURN

  FOR ivar=1, n_h5vars-1 DO BEGIN

     data_var  = H5VARS(ivar)
     data_arr  = READ_VAR_H5(file=file, var=data_var, start=start, count=count)

; COMPUTING MINIMA AND MAXIMA OF DATA ARRAYS

     data_min  = MIN(data_arr)
     data_max  = MAX(data_arr)

; STORING EVERYTHING IN THE STRUCTURE DATA

     data_item = {arr:data_arr,min:data_min,max:data_max}
     data = CREATE_STRUCT(data_var,data_item,data)

  ENDFOR


END
