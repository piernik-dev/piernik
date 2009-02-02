PRO DIRS, data_dir, frame_dir, png_output

ds = STR_SEP(data_dir,'/')
ns = N_ELEMENTS(ds)
subsub_dir = ds(ns-2)

IF(png_output EQ 'y') THEN BEGIN
  dlength = STRLEN( FINDFILE(frame_dir ) )
  IF ( dlength[0] EQ 0 ) THEN BEGIN
    FILE_MKDIR, frame_dir+'/'
  ENDIF
    png_frame_dir = frame_dir +'/'
    dlength = STRLEN( FINDFILE(png_frame_dir ) )
  IF ( dlength[0] EQ 0 ) THEN BEGIN
    FILE_MKDIR, png_frame_dir
  ENDIF
ENDIF

END
