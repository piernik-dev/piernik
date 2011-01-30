PRO SHOW_SLICE, data, cs, plp, disp_par, units, slice_array, ps, i_slice, time, sldata


;ON_ERROR, 1

s = slice_array(i_slice)

xdis=sldata.xdis
ydis=sldata.ydis
udis=sldata.udis
vdis=sldata.vdis
sdis=sldata.sdis
xrange=sldata.xrange
yrange=sldata.yrange

lunit=units.len.name
tunit=units.tim.name
vunit=units.vel
bunit=units.mag

ncolors=disp_par.ncolors
color=disp_par.color
background  = disp_par.backg
vec_color = disp_par.vec_color

ok = Execute('mini = data.'+s.scal_disp+'.min')
ok = Execute('maxi = data.'+s.scal_disp+'.max')

IF((maxi - mini) LT 2.e-5) THEN BEGIN
  mini = mini - 1.e-5
  maxi = maxi + 1.e-5
ENDIF


IF(s.type EQ cs.plane3) THEN BEGIN
  plot_hsize = ps.plot_hsize_xy
  plot_vsize = ps.plot_vsize_xy
  win_hsize = ps.win_hsize_xy
  win_vsize = ps.win_vsize_xy
  position_ = ps.position_xy

  xtitle = cs.coord1
  ytitle = cs.coord2
  IF(lunit NE '') THEN BEGIN
    xtitle = xtitle+' ['+lunit+']'
    ytitle = ytitle+' ['+lunit+']'
  ENDIF

  slice_str = '  '+cs.coord3+'='+ STRING(s.coord, FORMAT=plp.slice_str_form)+' '+lunit
  !x.ticks = plp.xticks
  !y.ticks = plp.yticks
  !x.minor = plp.xminor
  !y.minor = plp.yminor
   xtickf  = plp.xtickf
   ytickf  = plp.ytickf
ENDIF


IF(s.type EQ cs.plane2) THEN BEGIN
  plot_hsize = ps.plot_hsize_xz
  plot_vsize = ps.plot_vsize_xz
  win_hsize = ps.win_hsize_xz
  win_vsize = ps.win_vsize_xz
  position_ = ps.position_xz

  xtitle = cs.coord1
  ytitle = cs.coord3
  IF(lunit NE '') THEN BEGIN
    xtitle = xtitle+' ['+lunit+']'
    ytitle = ytitle+' ['+lunit+']'
  ENDIF

  slice_str = '  '+cs.coord2+'='+ STRING(s.coord, FORMAT=plp.slice_str_form)+' '+lunit
  !x.ticks = plp.xticks
  !y.ticks = plp.yticks
  !x.minor = plp.xminor
  !y.minor = plp.yminor
   xtickf  = plp.xtickf
   ytickf  = plp.ztickf
ENDIF

IF(s.type EQ cs.plane1) THEN BEGIN
  plot_hsize = ps.plot_hsize_yz
  plot_vsize = ps.plot_vsize_yz
  win_hsize = ps.win_hsize_xy
  win_vsize = ps.win_vsize_xy
  position_ = ps.position_yz

  xtitle = cs.coord2
  ytitle = cs.coord3
  IF(lunit NE '') THEN BEGIN
    xtitle = xtitle+' ['+lunit+']'
    ytitle = ytitle+' ['+lunit+']'
  ENDIF

  slice_str = '  '+cs.coord1+'='+ STRING(s.coord, FORMAT=plp.slice_str_form)+' '+lunit
  !x.ticks = plp.yticks
  !y.ticks = plp.zticks
  !x.minor = plp.yminor
  !y.minor = plp.zminor
   xtickf  = plp.ytickf
   ytickf  = plp.ztickf
ENDIF


time_str = 't='+ STRING(time, FORMAT=plp.time_str_form)+' '+tunit



image = BYTSCL(sdis,MIN=mini,MAX=maxi,TOP=ncolors-2)


IF(s.vect_scaling EQ 'fix') THEN BEGIN
  vect_length = MAX(SQRT(udis^2+vdis^2))*s.vect_scale
ENDIF ELSE IF(s.vect_scaling EQ 'free') THEN BEGIN
  vect_length = s.vect_scale
ENDIF

;print,  'vect_length =', vect_length

uvdismx=MAX(ABS(udis)+ABS(vdis))

IF(uvdismx LT 1.e-10) THEN BEGIN
  udis(0,0) = udis(0,0)-1.e-10
  vdis(0,0) = udis(0,0)-1.e-10
ENDIF

vsiz = SIZE(vdis)
nv1 = vsiz(1)
nv2 = vsiz(2)
vnul = FLTARR(nv1,nv2)

xdis = xdis/units.len.value
ydis = ydis/units.len.value
xrange = xrange/units.len.value
yrange = yrange/units.len.value

velovect1,vnul,vnul,xdis,ydis,length=vect_length ,charsize=plp.chars $
        , position = position_ $
        , xrange = xrange, yrange=yrange, xstyle=1,ystyle=1 $
        , xtickformat=xtickf,ytickformat=ytickf $
        , color=color,background=background,dots=1,/nodata


SX = plot_hsize
SY = plot_vsize
TV,CONGRID(image,plot_hsize,plot_vsize),ps.plot_h0,ps.plot_v0

; Plotting vectors
velovect1,udis,vdis,xdis,ydis,length=vect_length, charsize=plp.chars $
        , position = position_ $
        , xrange = xrange, yrange=yrange, xstyle=1,ystyle=1 $
        , xtickformat=xtickf,ytickformat=ytickf $
        , color=vec_color , background =background , maxvector=maxvector, dots=1,/noerase

IF(s.panel_name NE '-') THEN BEGIN
  panel ='('+s.panel_name+')'
ENDIF ELSE BEGIN
  panel = '   '
ENDELSE

; Plotting frame
velovect1, vnul,vnul,xdis,ydis,length=vect_length $
        , title= +time_str  +' '+ slice_str $
        , position = position_ $
        , /noerase ,charsize=plp.chars $
        , xrange = xrange, yrange=yrange, xstyle=1,ystyle=1 $
        , xtickformat=xtickf,ytickformat=ytickf $
        , color=color,bacground=background, /nodata

XYOUTS, ps.plot_h0/6, ps.plot_v0+plot_vsize+ps.plot_v1/3, panel, color=color,align=0.0 $
       ,/device, charsize =1.2*plp.chars

XYOUTS, ps.plot_h0+plot_hsize/2, 5,xtitle, color=color,align=0.5 $
       ,/device, charsize =plp.chars

XYOUTS, ps.plot_h0/4, ps.plot_v0 + plot_vsize/2, ytitle, color=color,align=0.5 $
       ,orientation=90,/device, charsize =plp.chars


; THINGS TO CLARIFY !!!
max_arrow_device = maxvector.arrow/(data.attr.ymax-data.attr.ymin)*plot_vsize*units.len.value
max_step_device = maxvector.step/(data.attr.ymax-data.attr.ymin)*plot_vsize*units.len.value
;max_arrow_device = maxvector.arrow*plot_vsize
;max_step_device = maxvector.step*plot_vsize

c_unit = FIX(maxvector.arrow/maxvector.step)
unit_arrow_device = maxvector.arrow/c_unit
unit_norm         = maxvector.norm/c_unit

arrow_device = max_arrow_device
;arrow_device = unit_arrow_device

norm = maxvector.norm
;norm = unit_norm

IF(plp.colbar EQ 1) THEN BEGIN
ARROW, ps.plot_h0+plot_hsize+0.5*plp.chars*!D.X_CH_SIZE  $
     , ps.plot_v0+plot_vsize+ps.plot_v1-0.75*plp.chars*!D.Y_CH_SIZE-arrow_device   $
     , ps.plot_h0+plot_hsize+0.5*plp.chars*!D.X_CH_SIZE   $
     , ps.plot_v0+plot_vsize+ps.plot_v1-0.75*plp.chars*!D.Y_CH_SIZE $
     , color=color,hsize = 0.75*!D.X_CH_SIZE


IF(s.vect_disp EQ 'mag') THEN BEGIN
  title_vec ='B:'
  unit_vec = bunit
  norm_disp = norm
ENDIF ELSE IF(strmid(s.vect_disp,0,3) EQ 'vel') THEN BEGIN
  flind = STRMID(s.vect_disp,3,1)
  title_vec ='V'+flind+':'
  unit_vec = vunit
  norm_disp = norm
ENDIF

IF(norm_disp GE 10.0) THEN BEGIN
  norm_str = STRING(norm_disp, format='(I4)')
ENDIF ELSE  IF(norm_disp LE 0.1 AND norm_disp GT 0.01 ) THEN BEGIN
  norm_str = STRING(norm_disp, format='(F5.2)')
ENDIF ELSE  IF(norm_disp LE 0.01 AND norm_disp GT 0.001 ) THEN BEGIN
  norm_str = STRING(norm_disp, format='(F6.3)')
ENDIF ELSE  IF(norm_disp LE 0.001 AND norm_disp GT 0.0001 ) THEN BEGIN
  norm_str = STRING(norm_disp, format='(F7.4)')
ENDIF ELSE  IF(norm_disp LE 0.0001 AND norm_disp GT 0.00001 ) THEN BEGIN
  norm_str = STRING(norm_disp, format='(F7.5)')
ENDIF ELSE BEGIN
  norm_str = STRING(norm_disp, format='(f4.2)')
ENDELSE

XYOUTS, ps.plot_h0+plot_hsize - ps.plot_h1/2. , ps.plot_v0+plot_vsize+ps.plot_v1/4. $
        , title_vec $
        , color=color,align=0.5,/device, charsize =1.2*plp.chars


XYOUTS, ps.plot_h0+plot_hsize + 0.6*ps.plot_h1, ps.plot_v0+plot_vsize+0.25*ps.plot_v1  $
      , color = color, align=0.0,/device, charsize=1.2*plp.chars $
      , '='+norm_str + unit_vec

ENDIF


; PLOTTING THE COLOR BAR
; ----------------------

color_bar = FLTARR(20,plot_vsize)

FOR i=0,plot_vsize-1 DO color_bar(*,i)=mini+(maxi-mini)*i/FLOAT(plot_vsize-1)

pos_cbar = [2*ps.plot_h0 + plot_hsize , ps.plot_v0, $
               2*ps.plot_h0 + plot_hsize+20, ps.plot_v0 + plot_vsize]

IF(plp.colbar EQ 1) THEN BEGIN
TV, BYTSCL(color_bar,MIN=mini,MAX=maxi,TOP=ncolors-2) $
           , pos_cbar(0),pos_cbar(1) ,/device

title_cbar = s.scal_disp

IF(s.scal_disp EQ 'd' AND s.scal_pert EQ 'p') THEN title_cbar ='!4Dq/q!3!D0!N'
IF(s.scal_disp EQ 'd' AND s.scal_pert EQ '') THEN title_cbar ='!3log!D10!N(!4q!3)'
IF(s.scal_disp EQ 'e' AND s.scal_pert EQ 'p') THEN title_cbar ='!4D!3e/e'
IF(s.scal_disp EQ 'e' AND s.scal_pert EQ '') THEN title_cbar ='e'
IF(s.scal_disp EQ 'cr' AND s.scal_pert EQ 'p') THEN title_cbar ='!4D!3e!Dcr!N/e!Dcr!N!Dcr!N'
IF(s.scal_disp EQ 'cr' AND s.scal_pert EQ '') THEN title_cbar ='!3log!D10!N(e!Dcr!N)'
;IF(s.scal_disp EQ 'cr' AND s.scal_pert EQ '') THEN title_cbar ='V!Dphi!N'
;IF(s.scal_disp EQ 'cr' AND s.scal_pert EQ '') THEN title_cbar ='B!D!4u!3!N'
IF(s.scal_disp EQ 'j')                        THEN title_cbar ='j'


IF(ABS(maxi-mini) LE 1.e-5 ) THEN BEGIN
  ytickf='(F10.6)' ;'(E13.6)'
ENDIF ELSE  IF(ABS(maxi-mini) LE 1.e-4 ) THEN BEGIN
  ytickf='(F10.6)' ;'(E12.5)'
ENDIF ELSE  IF(ABS(maxi-mini) LE 1.e-3 ) THEN BEGIN
  ytickf='(F9.5)' ;'(E11.4)'
ENDIF ELSE IF(ABS(maxi-mini) LE 1.e-2 ) THEN BEGIN
  ytickf='(F8.4)' ;'(E10.3)'
ENDIF ELSE IF(ABS(maxi-mini) LE 1.e-1 ) THEN BEGIN
  ytickf='(F7.3)' ;'(E9.2)'
ENDIF ELSE IF(ABS(maxi-mini) LE 1.e+3 ) THEN BEGIN
  ytickf='(F7.2)' ;'(E9.2)'
;  ytickf='(F7.1)' ;'(E9.2)'
ENDIF ELSE BEGIN
  ytickf='(E9.2)'
ENDELSE

plot,color_bar $
     , position=pos_cbar , color = color , charsize = plp.chars $
     , xticks = 1, xminor=1, xtickname=[' ',' '] $
     , yticklen=0.08,yticks=4,ytickformat=ytickf, ystyle =1 $
     ,/noerase,/nodata,/device

XYOUTS, 2*ps.plot_h0+plot_hsize + 10 , ps.plot_v0/2.  $
      , color = color, align=0.66,/device, charsize=1.4*plp.chars $
      ,  title_cbar

ENDIF

SKIP:


END
