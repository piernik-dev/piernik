
PRO READ_DIMS_AXES, filename, x, y, z, nx, ny, nz

; Written by M.Hanasz

ON_ERROR, 1

 FileID   = HDF_SD_START(filename, /READ)
 HDF_SD_FILEINFO, FileID, n_datasets, n_attributes
 FOR sds_index = 0, n_datasets - 1 DO BEGIN
   sds_id = HDF_SD_SELECT( FileID, sds_index )
      IF ( HDF_SD_ISCOORDVAR(sds_id) EQ 1 ) THEN BEGIN
        HDF_SD_GETINFO, sds_id, NAME = name, LABEL = label, DIMS = dim
        CASE name OF
          'fakeDim1': BEGIN
                   HDF_SD_GETDATA, sds_id, x
                   nx = dim
                 END
          'fakeDim2': BEGIN
                   HDF_SD_GETDATA, sds_id, y
                   ny = dim
                 END
          'fakeDim3': BEGIN
                   HDF_SD_GETDATA, sds_id, z
                   nz = dim
                 END
        ELSE:
        ENDCASE
      ENDIF
      HDF_SD_ENDACCESS, sds_id
 ENDFOR

 HDF_SD_END, FileID

;IF(nz EQ 1) THEN BEGIN
  x_ = x & nx_ = nx
  x  = y & nx  = ny
  y  = x_ & ny = nx_
;ENDIF

END
