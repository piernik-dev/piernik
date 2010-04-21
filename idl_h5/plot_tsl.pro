PRO PLOT_TSL, dir

; Written by M. Hanasz, 2003, 2007

batchmode = 0
IF (n_params(0) ge 1) THEN BEGIN
  batchmode = 1
  set_plot,'Z'
  PRINT, 'IDL:  Entering TSL in the batch mode'
ENDIF

batchmode_output = 'png'

colortable=0 ; 39
LOADCT, colortable

ncolors =!d.table_size-1
white = ncolors
!P.BACKGROUND = white
!P.COLOR = 0
!P.POSITION = [0.17,0.22,0.95,0.92]
!P.CHARSIZE  = 1


resolution = [400,300]
ps_xsize = 12
ps_ysize = 9

START:

;===========================================================================

ans=''
var=''
output=''
line = ''


; determine current working directory

CD, './', CURRENT=dir

NEXTFILE:

CASE batchmode OF
  0: BEGIN
       work_dir = ''
READWD:
       READ, PROMPT='WORKING DIRECTORY: ', work_dir
       PRINT, ''
       IF(work_dir EQ '') THEN BEGIN
          work_dir=last_dir
         GOTO, displayvar
       ENDIF ELSE BEGIN
         wdir_exist= FILE_TEST(work_dir,/directory)
         IF(wdir_exist NE 1) THEN GOTO, readwd
         last_dir=work_dir
       ENDELSE

     END
  1: work_dir = dir
ENDCASE

; make list of all '*.tsl' filenames

CD, work_dir

tsl_files = FILE_SEARCH('*.tsl')
ntsl = N_ELEMENTS(tsl_files)
start_times = fltarr(ntsl)

; determine start times  of all tsl data

FOR itsl=0, ntsl-1 do begin
  tsl_file = tsl_files(itsl)
  PRINT, 'tsl file #',itsl+1,':      ', tsl_file
  file_name = work_dir + '/'+ tsl_file
  OPENR, tsl_lun, file_name, /GET_LUN,ERROR=err
    READF,tsl_lun,line
    namearr = STRSPLIT(line, ' ',/extract,count=nvar)
    itime = WHERE(namearr EQ 'time') -1
    READF,tsl_lun,line      ; empty line

    data_line=FLTARR(nvar)
    READF,tsl_lun, data_line
    start_times(itsl) = data_line(itime)
  FREE_LUN, tsl_lun
ENDFOR


; determine the number of variables in tsl files

if(namearr(0) EQ '#') then namearr = namearr(1:nvar-1)
nvar = n_elements(namearr)
IF(nvar MOD 2 EQ 0) THEN BEGIN
  list_nvar = nvar
  list_namearr = namearr
ENDIF ELSE BEGIN
  list_nvar = nvar + 1
  list_namearr = [namearr,'none']
ENDELSE

firststep = 'y'

message, /reset_error_state

; read data from tsl files

FOR itsl = 0, ntsl-1 DO BEGIN

  file_name = tsl_files(itsl)

  OPENR, tsl_lun, file_name, /GET_LUN,ERROR=err

  READF,tsl_lun,line
  READF,tsl_lun,line
  darrayn=fltarr(nvar)

  WHILE (NOT EOF(tsl_lun)) DO BEGIN
    READF,tsl_lun, darrayn

     IF(itsl LT ntsl-1) THEN IF(darrayn(itime) GT start_times(itsl+1)) THEN GOTO, SKIP

     IF(firststep EQ 'y') THEN BEGIN
       darray = darrayn
       firststep = 'n'
       npoint = 1
     ENDIF ELSE BEGIN
       darray = [[darray],[darrayn]]
       npoint = npoint + 1
     ENDELSE

  ENDWHILE

  SKIP:

  FREE_LUN,tsl_lun

ENDFOR

GOTO, displayvar

;=============================================================================

DISPLAYVAR:

CASE batchmode OF
  0: BEGIN
       VARIABLES:
       PRINT, ''
       READ, PROMPT='VARIABLE #/(l)ist/(n)ew data/e(x)it: ', var
       IF(var EQ  'x') THEN GOTO, finish
       IF(var EQ  'n') THEN GOTO, start
       IF(var EQ  'l') THEN BEGIN
         ih = list_nvar/2
         FOR i = 0, ih-1 DO BEGIN
             PRINT, FORMAT='(2(i4,": ",a12,10x))', i, STRTRIM(list_namearr(i),1), i+ih, STRTRIM(list_namearr(i+ih),1)
         ENDFOR
       GOTO, variables
       ENDIF
       ivar = FIX(var)
       IF(ivar GE nvar) THEN BEGIN
         PRINT, '# out of range'
         GOTO, DISPLAYVAR
       ENDIF
       dispind = [ivar]
       disparr = [namearr(ivar)]
       ndisp = 1
     END
  1: BEGIN
       file_name = 'plot_tsl.var'
       OPENR, var_lun, file_name, /GET_LUN,ERROR=err
         READF, var_lun, line
       FREE_LUN, var_lun
       disparr = STRSPLIT(line, ' ',/extract,count=ndisp)
       dispind = INTARR(ndisp)
       FOR idisp = 0, ndisp-1 DO BEGIN
         dispind(idisp) = WHERE(namearr EQ disparr(idisp))
       ENDFOR
     END
ENDCASE

PLOTTING:

FOR i=0, ndisp-1 DO BEGIN

  CASE batchmode OF
    0: BEGIN
         PRINT,''
         READ, PROMPT='(s)crean, (ps), (png), (dat) or (c)ancel ? ', output
       END
    1: BEGIN
         output = batchmode_output
         ivar = dispind(i)
       END
  ENDCASE

IF (output EQ 'c') THEN BEGIN
  GOTO, displayvar
ENDIF ELSE IF (output EQ 'x') THEN BEGIN
  set_plot,'x'
  LOADCT, colortable
  ncolors =!d.table_size-1
  white = ncolors
  !P.BACKGROUND = white
  !P.COLOR = 0
ENDIF ELSE IF (output EQ 'dat') THEN BEGIN
  plotname0=strtrim(namearr(ivar))
  plotname1=STRSPLIT(plotname0 ,/EXTRACT)
  plotname = STRJOIN(plotname1,'_')

  dat_dir= work_dir+'/dat'
  dat_exist= FILE_TEST(dat_dir,/directory)
  IF(dat_exist NE 1) THEN SPAWN, 'mkdir ' + dat_dir

  file_name = work_dir+ '/dat/'+ plotname +'.dat'
  PRINT, 'DATA FILE: ', file_name
  OPENW,dat_lun,file_name,/GET_LUN

  FOR np=0, npoint-1 DO BEGIN
    PRINTF, FORMAT='(f12.8,1x,e12.4)', dat_lun, darray(itime,np),darray(ivar,np)
  ENDFOR
  FREE_LUN,dat_lun
ENDIF ELSE IF (output EQ 's' OR output EQ 'ps' OR output EQ 'png') THEN BEGIN
  plotname0=strtrim(namearr(ivar))
  plotname1=STRSPLIT(plotname0 ,/EXTRACT)
  plotname = STRJOIN(plotname1,'_')

  plot_dir= work_dir+'/plots'
  plot_exist= FILE_TEST(plot_dir,/directory)
  IF(plot_exist NE 1) THEN SPAWN, 'mkdir ' + plot_dir



  IF(output EQ 'ps') THEN BEGIN
    set_plot,'ps'
    psfile=plot_dir+ '/' + plotname +'.ps'
    PRINT, 'PS FILE: ', psfile

    device, file=psfile $
          , xsize=ps_xsize,ysize=ps_ysize $
          , xoffset=2,yoffset=2 $
          ,  bits=8,/color
  ENDIF
  IF(output EQ 'png') THEN BEGIN
    SET_PLOT, 'Z'
    DEVICE, SET_RESOLUTION=resolution
    !P.BACKGROUND = white
    !P.COLOR = 0
  ENDIF
  IF(output EQ 's') THEN BEGIN
    SET_PLOT, 'X'
    !P.BACKGROUND = white
    !P.COLOR = 0
    WINDOW, ivar, xsize= resolution(0), ysize=resolution(1)
  ENDIF

  xrange = [min(darray(itime,*)), max(darray(itime,*))]
  range  = [min(darray(ivar,*)),  max(darray(ivar,*))]


  CASE namearr(ivar) OF

; special wishes

    'emag':  plot, darray(itime,*),darray(ivar,*) $
                 , title='MAGNETIC ENERGY'   $
                 , xrange = xrange, xstyle=1 $
                 , xtitle = 'time [Myr]' $
                 , yrange = yrange $
                 , ystyle=1 $
                 , /ylog

  ELSE: BEGIN

; default plot format

             plot, darray(itime,*),darray(ivar,*),title=namearr(ivar) $
                 , xtitle = 'time [Myr]' $
                 , xrange = xrange, xstyle=1
        END
  ENDCASE

  plottime = SYSTIME()
  XYOUTS, 0.03,0.05, plottime, /norm

  dirsep = STR_SEP(work_dir,'/')
  ns = N_ELEMENTS(dirsep)
  sub_dir = dirsep(ns-1)
  XYOUTS, 0.03,0.01, sub_dir, /norm

  IF(output EQ 'png') THEN BEGIN
    frame = TVRD()
    pngfile = plot_dir+ '/'+ plotname +'.png'
    WRITE_PNG,  pngfile,frame,red,green,blue, /verbose
  ENDIF
  IF(output EQ 'ps' OR output EQ 'png') THEN BEGIN
    device, /close
    set_plot,'x'
    !P.BACKGROUND = white
    !P.COLOR = 0
  ENDIF

ENDIF

ENDFOR

CASE batchmode OF
  0: GOTO, displayvar
  1: GOTO, finish
ENDCASE

FINISH:
IF(batchmode EQ 1) THEN BEGIN
  PRINT, 'IDL: Leaving TSL '
ENDIF

END
