PRO SHOW_SLICE, slice_array, i_slice, time


COMMON coordsys    , coordsys, coordnames $
                   , coord1,coord2,coord3,plane1,plane2,plane3

COMMON coordinates,  ni, nj, nk, x1a, x2a, x3a, x1b, x2b, x3b, dims
COMMON quantites,    dd, ee, gp, v1, v2, v3, b1, b2, b3, ecr
COMMON labels,       x1alab, x2alab, x3alab, x1blab, x2blab, x3blab  $
                     , ddlab, eelab, gplab, v1lab                      $
                     , v2lab, v3lab, b1lab, b2lab, b3lab,ecrlab

COMMON quantites2,    dd1, ee1, ecr1, jj


COMMON data_min_max, min_dd, max_dd, min_dd0, max_dd0, min_dd1, max_dd1 $
                   , min_ee, max_ee, min_ee0, max_ee0, min_ee1, max_ee1 $
                   , min_ecr, max_ecr, min_ecr0, max_ecr0, min_ecr1, max_ecr1 $ 
                   , min_jj, max_jj
 
COMMON mag_field, b1_mean, b2_mean, bh_mean, b1_m 


COMMON sldata  $
       , xdis, ydis, udis, vdis, sdis  $
       , xrange_, yrange_ $ 
       , time_str, slice_str $
       , lunit, tunit, bunit, vunit, dunit, eunit, ecrunit  

COMMON plot_params, chars_, xmarg_, ymarg_, color_bar_

COMMON plot_coords $
      , win_vsize, plot_h0, plot_v0, plot_h1, plot_v1 $
      , win_hsize_xy, win_vsize_xy,  plot_hsize_xy, plot_vsize_xy, position_xy $
      , win_hsize_xz, win_vsize_xz,  plot_hsize_xz, plot_vsize_xz, position_xz $
      , win_hsize_yz, win_vsize_yz,  plot_hsize_yz, plot_vsize_yz, position_yz $
      , win_hsize_bm, win_vsize_bm,  plot_hsize_bm, plot_vsize_bm, position_bm 


COMMON disp, ncolors, min_data, max_data, mincolor, maxcolor

COMMON max_vector, max_arrow, max_norm, max_step

COMMON anim, animation_name,nframes

;on_error, 1

s = slice_array(i_slice)

IF (s.scal_disp EQ 'j') THEN BEGIN
  IF(s.scal_scaling EQ 'free') THEN BEGIN
    mini = min_jj 
    maxi = max_jj 
  ENDIF ELSE IF(s.scal_scaling EQ 'fix') THEN BEGIN 
    mini = s.scal_scale(0)
    maxi = s.scal_scale(1) 
  ENDIF
ENDIF ELSE IF (s.scal_disp EQ 'e') THEN BEGIN
  IF(s.scal_scaling EQ 'free') THEN BEGIN
    IF(s.scal_pert EQ 'p') THEN BEGIN
      mini = min_ee1 
      maxi = max_ee1 
      amax = MAX([ABS(mini),ABS(maxi)])
      mini = -amax
      maxi = amax
    ENDIF ELSE IF(s.scal_pert EQ '') THEN BEGIN
      mini = min_ee 
      maxi = max_ee 
    ENDIF
  ENDIF ELSE IF(s.scal_scaling EQ 'fix0') THEN BEGIN 
      mini = min_ee0 
      maxi = max_ee0 
      amax = MAX([ABS(mini),ABS(maxi)])
      mini = -amax
      maxi = amax
  ENDIF ELSE IF(s.scal_scaling EQ 'fix') THEN BEGIN 
    mini = s.scal_scale(0)
    maxi = s.scal_scale(1) 
  ENDIF
ENDIF ELSE IF (s.scal_disp EQ 'd') THEN BEGIN
  IF(s.scal_scaling EQ 'free') THEN BEGIN
    IF(s.scal_pert EQ 'p') THEN BEGIN
      mini = min_dd1 
      maxi = max_dd1 
      amax = MAX([ABS(mini),ABS(maxi)])
      mini = -amax
      maxi = amax
    ENDIF ELSE IF(s.scal_pert EQ '') THEN BEGIN
      mini = min_dd 
      maxi = max_dd 
    ENDIF
  ENDIF ELSE IF(s.scal_scaling EQ 'fix0') THEN BEGIN 
      mini = min_dd0 
      maxi = max_dd0 
      amax = MAX([ABS(mini),ABS(maxi)])
      mini = -amax
      maxi = amax
  ENDIF ELSE IF(s.scal_scaling EQ 'fix') THEN BEGIN 
      mini = s.scal_scale(0)
      maxi = s.scal_scale(1) 
  ENDIF
ENDIF ELSE IF (s.scal_disp EQ 'cr') THEN BEGIN
  IF(s.scal_scaling EQ 'free') THEN BEGIN
    IF(s.scal_pert EQ 'p') THEN BEGIN
      mini = min_ecr1 
      maxi = max_ecr1 
      amax = MAX([ABS(mini),ABS(maxi)])
      mini = -amax
      maxi = amax
    ENDIF ELSE IF(s.scal_pert EQ '') THEN BEGIN
      mini = min_ecr 
      maxi = max_ecr
    ENDIF
  ENDIF ELSE IF(s.scal_scaling EQ 'fix0') THEN BEGIN 
      mini = min_ecr0 
      maxi = max_ecr0 
      amax = MAX([ABS(mini),ABS(maxi)])
      mini = -amax
      maxi = amax
  ENDIF ELSE IF(s.scal_scaling EQ 'fix') THEN BEGIN 
    mini = s.scal_scale(0)
    maxi = s.scal_scale(1) 
  ENDIF
ENDIF

IF((maxi - mini) LT 2.e-5) THEN BEGIN
  mini = mini - 1.e-5
  maxi = maxi + 1.e-5
ENDIF


IF(s.type EQ plane3) THEN BEGIN
  plot_hsize = plot_hsize_xy
  plot_vsize = plot_vsize_xy
  win_hsize = win_hsize_xy
  win_vsize = win_vsize_xy
  position_ = position_xy
  IF(dims EQ '3d') THEN BEGIN
    IF(s.type EQ 'zr') THEN BEGIN
       xtit_ = coord2+lunit
       ytit_ = coord1+lunit
    ENDIF ELSE BEGIN 
       xtit_ = coord1+lunit
       ytit_ = coord2+lunit
    ENDELSE

  ENDIF ELSE IF(dims EQ '2d') THEN BEGIN
    xtit_ = coord1+lunit
    ytit_ = coord2+lunit
  ENDIF

  slice_str = ' '+coord3+'='+ STRING(s.coord, FORMAT='(E9.2)')+lunit
;  slice_str = '  '+coord3+'='+ STRING(s.coord, FORMAT='(E6.0)')+lunit
;  slice_str = '  '+coord3+'='+ STRING(s.coord, FORMAT='(I4)')+lunit
 !x.ticks = 4
  !y.ticks = 4
  !x.minor = 6;2 
  !y.minor = 6;2
  xtickf_='(E9.2)'
  ytickf_='(E9.2)'
;  xtickf_='(I5)'
;  ytickf_='(I5)'
ENDIF  


IF(s.type EQ plane2) THEN BEGIN 
  plot_hsize = plot_hsize_xz
  plot_vsize = plot_vsize_xz
  win_hsize = win_hsize_xz
  win_vsize = win_vsize_xz
  position_ = position_xz
  xtit_ = coord1+lunit
  ytit_ = coord3+lunit
;s.coord=0  
  slice_str = ' '+coord2+'='+ STRING(s.coord, FORMAT='(E9.2)')+lunit
;  slice_str = '  '+coord2+'='+ STRING(s.coord, FORMAT='(E6.0)')+lunit
;  slice_str = '  '+coord2+'='+ STRING(s.coord, FORMAT='(I4)')+lunit
  !x.ticks = 4
  !y.ticks = 4
  !x.minor = 2 
  !y.minor = 2
  xtickf_='(E9.2)'
  ytickf_='(E9.2)'
;  xtickf_='(I5)'
;  ytickf_='(I5)'
ENDIF

IF(s.type EQ plane1) THEN BEGIN
  plot_hsize = plot_hsize_yz
  plot_vsize = plot_vsize_yz
  win_hsize = win_hsize_xy
  win_vsize = win_vsize_xy
  position_ = position_yz
  xtit_ = coord2+lunit
  ytit_ = coord3+lunit
  slice_str = ' '+coord1+'='+ STRING(s.coord, FORMAT='(E9.2)')+lunit
  ;slice_str = '  '+coord1+'='+ STRING(s.coord, FORMAT='(E6.0)')+lunit
  ;slice_str = '  '+coord1+'='+ STRING(s.coord, FORMAT='(I4)')+lunit
  !x.ticks = 4
  !y.ticks = 4
  !x.minor = 2 
  !y.minor = 2
  xtickf_='(E9.2)'
  ytickf_='(E9.2)'
;  xtickf_='(I5)'
;  ytickf_='(I5)'
ENDIF 


time_str = ' t='+ STRING(time, FORMAT='(E8.2)')+tunit
;time_str = ' t='+ STRING(time, FORMAT='(E6.0)')+tunit
;time_str = 't='+ STRING(time, FORMAT='(I3)')+tunit



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


velovect1,vnul,vnul,xdis,ydis,length=vect_length ,charsize=chars_ $
        , position = position_ $
        , xrange = xrange_, yrange=yrange_, xstyle=1,ystyle=1 $
        , xtickformat=xtickf_,ytickformat=ytickf_ $
        , color=0,background=ncolors,dots=1,/nodata


SX = plot_hsize
SY = plot_vsize 
TV,CONGRID(image,plot_hsize,plot_vsize),plot_h0,plot_v0

velovect1,udis,vdis,xdis,ydis,length=vect_length, charsize=chars_ $
        , position = position_ $
        , xrange = xrange_, yrange=yrange_, xstyle=1,ystyle=1 $
        , xtickformat=xtickf_,ytickformat=ytickf_ $
        , color=ncolors , background = ncolors,dots=1,/noerase

IF(s.panel_name NE '-') THEN BEGIN
  panel ='('+s.panel_name+')'
ENDIF ELSE BEGIN
  panel = '   '
ENDELSE


velovect1, vnul,vnul,xdis,ydis,length=vect_length $
        , title= +time_str  +' '+ slice_str $
        , position = position_ $
        , /noerase ,charsize=chars_ $
        , xrange = xrange_, yrange=yrange_, xstyle=1,ystyle=1 $
        , xtickformat=xtickf_,ytickformat=ytickf_ $
        , color=0,bacground=0, /nodata

XYOUTS, plot_h0/6, plot_v0+plot_vsize+plot_v1/3, panel, color=0,align=0.0 $
       ,/device, charsize =1.2*chars_ 

XYOUTS, plot_h0+plot_hsize/2, 5, xtit_, color=0,align=0.5 $
       ,/device, charsize =chars_ 

XYOUTS, plot_h0/8, plot_v0 + plot_vsize/2, ytit_, color=0,align=0.5 $
       ,orientation=90,/device, charsize =chars_ 


max_arrow_device = max_arrow/(MAX(x2b)-MIN(x2b))*plot_vsize_xy
max_step_device = max_step/(MAX(x2b)-MIN(x2b))*plot_vsize_xy
c_unit = FIX(max_arrow/max_step)
unit_arrow_device = max_arrow_device/c_unit
unit_norm         = max_norm/c_unit 

arrow_device = max_arrow_device 
;arrow_device = unit_arrow_device 

norm = max_norm
;norm = unit_norm

IF(color_bar_ EQ 'y') THEN BEGIN
ARROW, plot_h0+plot_hsize+0.5*chars_*!D.X_CH_SIZE  $
     , plot_v0+plot_vsize+plot_v1-0.75*chars_*!D.Y_CH_SIZE-arrow_device   $
     , plot_h0+plot_hsize+0.5*chars_*!D.X_CH_SIZE   $
     , plot_v0+plot_vsize+plot_v1-0.75*chars_*!D.Y_CH_SIZE $
     ,color=0,hsize = 0.75*!D.X_CH_SIZE;,/device


IF(s.vect_disp EQ 'b') THEN BEGIN
  title_vec ='B:' 
  unit_vec = bunit 
  norm_disp = norm       
ENDIF ELSE IF(s.vect_disp EQ 'v') THEN BEGIN
  title_vec ='V:' 
  unit_vec = vunit
  norm_disp = norm
ENDIF

IF(norm_disp GE 10.0) THEN BEGIN
  norm_str = STRING(norm_disp, format='(E7.1)')
ENDIF ELSE  IF(norm_disp LT 0.1) THEN BEGIN
  norm_str = STRING(norm_disp, format='(E7.1)')
ENDIF ELSE BEGIN
  norm_str = STRING(norm_disp, format='(f4.2)')
ENDELSE

XYOUTS, plot_h0+plot_hsize - plot_h1/2. , plot_v0+plot_vsize+plot_v1/4. $
        , title_vec $
        , color=0,align=0.5,/device, charsize =1.2*chars_ 


XYOUTS, plot_h0+plot_hsize + 0.6*plot_h1, plot_v0+plot_vsize+0.25*plot_v1  $
      , color = 0, align=0.0,/device, charsize=1.2*chars_ $
      , '='+norm_str + unit_vec

ENDIF

;xsurface,sdis , $
;  shade =b2TSCL(sdis,MIN=color_scale*(max_imag-min_data),MAX=max_imag,TOP=ncolors-2)


;print, 'ph =',plot_hsize,'   pv =',plot_vsize


; PLOTTING THE COLOR BAR
; ----------------------

color_bar = FLTARR(20,plot_vsize)

FOR i=0,plot_vsize-1 DO color_bar(*,i)=mini+(maxi-mini)*i/FLOAT(plot_vsize-1)

pos_cbar = [2*plot_h0 + plot_hsize , plot_v0, $
               2*plot_h0 + plot_hsize+20, plot_v0 + plot_vsize]

IF(color_bar_ EQ 'y') THEN BEGIN
TV, BYTSCL(color_bar,MIN=mini,MAX=maxi,TOP=ncolors-2) $
           , pos_cbar(0),pos_cbar(1) ,/device


IF(s.scal_disp EQ 'd' AND s.scal_pert EQ 'p') THEN title_cbar ='!4Dq/q!3!D0!N' 
IF(s.scal_disp EQ 'd' AND s.scal_pert EQ '') THEN title_cbar ='!4q!3' 
IF(s.scal_disp EQ 'e' AND s.scal_pert EQ 'p') THEN title_cbar ='!4D!3e/e' 
IF(s.scal_disp EQ 'e' AND s.scal_pert EQ '') THEN title_cbar ='e' 
IF(s.scal_disp EQ 'cr' AND s.scal_pert EQ 'p') THEN title_cbar ='!4D!3e!Dcr!N/e!Dcr!N!Dcr!N' 
IF(s.scal_disp EQ 'cr' AND s.scal_pert EQ '') THEN title_cbar ='e!Dcr!N' 
IF(s.scal_disp EQ 'j')                        THEN title_cbar ='j' 


IF(ABS(maxi-mini) LE 1.e-5 ) THEN BEGIN
  ytickf_='(F10.6)' ;'(E13.6)'  
ENDIF ELSE  IF(ABS(maxi-mini) LE 1.e-4 ) THEN BEGIN
  ytickf_='(F10.6)' ;'(E12.5)'  
ENDIF ELSE  IF(ABS(maxi-mini) LE 1.e-3 ) THEN BEGIN
  ytickf_='(F9.5)' ;'(E11.4)'  
ENDIF ELSE IF(ABS(maxi-mini) LE 1.e-2 ) THEN BEGIN
  ytickf_='(F8.4)' ;'(E10.3)'  
ENDIF ELSE IF(ABS(maxi-mini) LE 1.e-1 ) THEN BEGIN
  ytickf_='(F7.3)' ;'(E9.2)'
ENDIF ELSE IF(ABS(maxi-mini) LE 1.e+3 ) THEN BEGIN
  ytickf_='(F7.2)' ;'(E9.2)'
ENDIF ELSE BEGIN 
  ytickf_='(E9.2)'
ENDELSE

plot,color_bar $
     , position=pos_cbar , color = 0 , charsize = chars_ $ 
     , xticks = 1, xminor=1, xtickname=[' ',' '] $
     , yticklen=0.08,yticks=4,ytickformat=ytickf_, ystyle =1 $
     ,/noerase,/nodata,/device

XYOUTS, 2*plot_h0+plot_hsize + 10 , plot_v0/2.  $
      , color = 0, align=0.5,/device, charsize=1.4*chars_ $
      ,  title_cbar 

ENDIF

SKIP:


END
