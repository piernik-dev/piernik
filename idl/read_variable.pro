PRO READ_VARIABLE, filename, varname, data, time, status

; Written by M.Hanasz

 status =-1

 FileID   = HDF_SD_START(filename, /READ)
 HDF_SD_FILEINFO, FileID, n_datasets, n_attributes
 FOR sds_index = 0, n_datasets - 1 DO BEGIN
   sds_id = HDF_SD_SELECT( FileID, sds_index )
      IF ( HDF_SD_ISCOORDVAR(sds_id) NE 1 ) THEN BEGIN

        HDF_SD_GETINFO, sds_id, LABEL = label, DIMS = dims
        lab = str_sep(label,' AT ')
        tmp = str_sep(label,'=')
        time = float(tmp(1))

        IF(varname EQ 'dd' AND lab(0) EQ 'DENSITY') THEN BEGIN
                HDF_SD_GETDATA, sds_id, data 
                status = 0
        ENDIF ELSE IF (varname EQ 'ee' AND lab(0) EQ 'GAS ENERGY')  THEN  BEGIN
                HDF_SD_GETDATA, sds_id, data 
                status = 0
        ENDIF ELSE IF (varname EQ 'er' AND lab(0) EQ 'RADIATION T') THEN  BEGIN
                HDF_SD_GETDATA, sds_id, data 
                status = 0
        ENDIF ELSE IF (varname EQ 'v1' AND lab(0) EQ '1-VELOCITY')  THEN  BEGIN
                HDF_SD_GETDATA, sds_id, data 
                status = 0
        ENDIF ELSE IF (varname EQ 'v2' AND lab(0) EQ '2-VELOCITY')  THEN  BEGIN
                HDF_SD_GETDATA, sds_id, data  
                status = 0
        ENDIF ELSE IF (varname EQ 'v3' AND lab(0) EQ '3-VELOCITY')  THEN  BEGIN
                HDF_SD_GETDATA, sds_id, data 
                status = 0
        ENDIF ELSE IF (varname EQ 'b1' AND lab(0) EQ '1-MAG FIELD') THEN  BEGIN
                HDF_SD_GETDATA, sds_id, data 
                status = 0
        ENDIF ELSE IF (varname EQ 'b2' AND lab(0) EQ '2-MAG FIELD') THEN  BEGIN
                HDF_SD_GETDATA, sds_id, data  
                status = 0
        ENDIF ELSE IF (varname EQ 'b3' AND lab(0) EQ '3-MAG FIELD') THEN  BEGIN
                HDF_SD_GETDATA, sds_id, data 
                status = 0
        ENDIF ELSE BEGIN
;
        ENDELSE
      ENDIF

      HDF_SD_ENDACCESS, sds_id
    ENDFOR

  HDF_SD_END, FileID



END
