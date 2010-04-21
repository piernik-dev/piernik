PRO SHOW_INIT, s, slice_array, n_slices, plot_par, display_par, units, batchmode


;------------------------------------------------------------------------------
s = {name:      ' ', sw:       ' ',     status: 0,           $
     type:      ' ', coord:        0.0, panel_name: ' ',     $
     vect_disp: ' ', vect_scaling: ' ', vect_scale: 1.0,     $
     scal_disp: ' ', scal_pert:    ' ',scal_scaling:' ',     $
     scal_log:  ' ', scal_scale:fltarr(2)}
slice_array = [s]
n_slices = 0
;------------------------------------------------------------------------------


batchmode = 0
colbar=1


IF(!D.NAME NE 'X') THEN BEGIN
  device,/close
  set_plot,'X'
ENDIF


color_table = 39
ncolors =!d.table_size-1
LOADCT, color_table,file='./colors1.tbl',/SILENT
tvlct,red,green,blue,/get

IF(color_table EQ 42) THEN BEGIN
   backg = ncolors/2
   color = ncolors
   vec_color = color
ENDIF ELSE BEGIN
   backg = ncolors
   color = 0
   vec_color = backg
ENDELSE

!P.BACKGROUND = backg
!P.COLOR = color

; Default values of plot params

plot_par = {chars:1.0, $
            xmarg:[12,3],  $
            ymarg:[4,3],  $
            colbar:1, $
            colbar_tick_form:'(E9.2)', $
            slice_str_form: '(E9.2)', $
            time_str_form: '(E9.2)', $
            xticks:4, $
            yticks:4, $
            zticks:4, $
            xminor:6, $
            yminor:6, $
            zminor:6, $
            xtickf:'(E9.2)', $
            ytickf:'(E9.2)', $
            ztickf:'(E9.2)'}

; User defined values of plot params

plot_par.chars  =  1.0
plot_par.xmarg  = [12,3]
plot_par.ymarg  = [4,3]
plot_par.colbar =  colbar
plot_par.colbar_tick_form ='(F5.1)'
plot_par.slice_str_form = '(F5.1)'
plot_par.time_str_form = '(I5)'
plot_par.xticks = 4
plot_par.yticks = 4
plot_par.zticks = 4
plot_par.xminor = 6
plot_par.yminor = 6
plot_par.zminor = 6
plot_par.xtickf = '(F5.1)'
plot_par.ytickf = '(F5.1)'
plot_par.ztickf = '(F5.1)'


; Default units
lunit={name:'',value:1}
tunit={name:'',value:1}

; Galactic disk units
;lunit={name:'kpc',value:1000}
;tunit={name:'Myr',value:1}

vunit='kms!U-1!N'
bunit='!4l!3G'

units = {len:lunit, $
         tim:tunit, $
         vel:vunit, $
         mag:bunit}

display_par = {ncolors:ncolors,     $
               color:color,         $
               vec_color:vec_color, $
               backg:backg}










END
