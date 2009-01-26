; $Id: vectors.pro 2 2007-10-29 13:14:38Z kowalik $
;
; Copyright (c) 1983-2001, Research Systems, Inc.  All rights reserved.
;	Unauthorized reproduction prohibited.

;
;+
; NAME:
;	VELOVECT
;
; PURPOSE:
;	Produce a two-dimensional velocity field plot.
;
;	A directed arrow is drawn at each point showing the direction and
;	magnitude of the field.
;
; CATEGORY:
;	Plotting, two-dimensional.
;
; CALLING SEQUENCE:
;	VELOVECT, U, V [, X, Y]
;
; INPUTS:
;	U:	The X component of the two-dimensional field.
;		U must be a two-dimensional array.
;
;	V:	The Y component of the two dimensional field.  Y must have
;		the same dimensions as X.  The vector at point [i,j] has a
;		magnitude of:
;
;			(U[i,j]^2 + V[i,j]^2)^0.5
;
;		and a direction of:
;
;			ATAN2(V[i,j],U[i,j]).
;
; OPTIONAL INPUT PARAMETERS:
; 	X:	Optional abcissae values.  X must be a vector with a length
;		equal to the first dimension of U and V.
;
;	Y:	Optional ordinate values.  Y must be a vector with a length
;		equal to the first dimension of U and V.
;
; KEYWORD INPUT PARAMETERS:
;	COLOR:	The color index used for the plot.
;
;	DOTS:	Set this keyword to 1 to place a dot at each missing point.
;		Set this keyword to 0 or omit it to draw nothing for missing
;		points.  Has effect only if MISSING is specified.
;
;	LENGTH:	Length factor.  The default of 1.0 makes the longest (U,V)
;		vector the length of a cell.
;
;       MISSING: Missing data value.  Vectors with a LENGTH greater
;		than MISSING are ignored.
;
;       OVERPLOT: Set this keyword to make VELOVECT "overplot".  That is, the
;               current graphics screen is not erased, no axes are drawn, and
;               the previously established scaling remains in effect.
;
;
;	Note:   All other keywords are passed directly to the PLOT procedure
;		and may be used to set option such as TITLE, POSITION,
;		NOERASE, etc.
; OUTPUTS:
;	None.
;
; COMMON BLOCKS:
;	None.
;
; SIDE EFFECTS:
;	Plotting on the selected device is performed.  System
;	variables concerning plotting are changed.
;
; RESTRICTIONS:
;	None.
;
; PROCEDURE:
;	Straightforward.  Unrecognized keywords are passed to the PLOT
;	procedure.
;
; MODIFICATION HISTORY:
;	DMS, RSI, Oct., 1983.
;	For Sun, DMS, RSI, April, 1989.
;	Added TITLE, Oct, 1990.
;	Added POSITION, NOERASE, COLOR, Feb 91, RES.
;	August, 1993.  Vince Patrick, Adv. Visualization Lab, U. of Maryland,
;		fixed errors in math.
;	August, 1993. DMS, Added _EXTRA keyword inheritance.
;	January, 1994, KDB. Fixed integer math which produced 0 and caused
;		            divide by zero errors.
;	December, 1994, MWR. Added _EXTRA inheritance for PLOTS and OPLOT.
;	June, 1995, MWR. Removed _EXTRA inheritance for PLOTS and changed
;			 OPLOT to PLOTS.
;       September, 1996, GGS. Changed denominator of x_step and y_step vars.
;       February, 1998, DLD.  Add support for CLIP and NO_CLIP keywords.
;       June, 1998, DLD.  Add support for OVERPLOT keyword.
;-
;
PRO VECTORS,UI,VI,XI,YI, Missing = Missing, Length = length, Dots = dots, NVECTORS=nvectors  $
       , Color=color, CLIP=clip, NOCLIP=noclip, OVERPLOT=overplot, Vcolor=vcolor, Vthick=vthick $
       , MAXVAL = maxval,_EXTRA=extra

COMPILE_OPT strictarr

;        on_error,2                      ;Return to caller if an error occurs


;  IF N_ELEMENTS(stride) LE 0 THEN stride = [1,1]
;  IF N_ELEMENTS(stride) NE 2 THEN BEGIN
;    MESSAGE, 'STRIDE must be 2-elements vector.'
;  ENDIF
;  ni = N_ELEMENTS( xi ) / stride[0]
;  nj = N_ELEMENTS( yi ) / stride[1]
 
  IF N_ELEMENTS(nvectors) LE 0 THEN BEGIN
    u = ui & v = vi & x = xi & y = yi
  ENDIF ELSE BEGIN
    IF N_ELEMENTS(nvectors) EQ 2 THEN BEGIN
      ni = nvectors[0] & nj = nvectors[1]

      uo = CONGRID( ui, ni, nj, /interp, /minus_one )
      vo = CONGRID( vi, ni, nj, /interp, /minus_one )
      xo = CONGRID( xi, ni, /interp, /minus_one )
      yo = CONGRID( yi, nj, /interp, /minus_one )

  xmin = MIN(xi) & xmax = MAX(xi)
  ymin = MIN(yi) & ymax = MAX(yi)
  IF ( N_ELEMENTS(extra.xrange) EQ 2 ) THEN BEGIN
    xmin = extra.xrange[0]
    xmax = extra.xrange[1]
  ENDIF
  IF ( N_ELEMENTS(extra.yrange) EQ 2 ) THEN BEGIN
    ymin = extra.yrange[0]
    ymax = extra.yrange[1]
  ENDIF

  sizeuv = SIZE(ui)
  sx = sizeuv[1]
  sy = sizeuv[2]

  di = sx/ni
  dj = sy/nj

  dx = (xmax-xmin)/(ni)
  dy = (ymax-ymin)/(nj)


     x = xmin+dx*(FINDGEN(ni) + 0.5)
     y = ymin+dy*(FINDGEN(nj) + 0.5)

     u = INTERPOLATE( ui, di*(FINDGEN(ni) + 0.5), dj*(FINDGEN(nj) + 0.5), /GRID)
     v = INTERPOLATE( vi, di*(FINDGEN(ni) + 0.5), dj*(FINDGEN(nj) + 0.5), /GRID)
     
;     x=xo
;     y=yo
;     u=uo
;     v=vo
      
    ENDIF ELSE BEGIN
      MESSAGE, 'NVECTORS must be 2-elements vector.'
    ENDELSE
  ENDELSE

  xmin = MIN(x) & xmax = MAX(x)
  ymin = MIN(y) & ymax = MAX(y)
  IF ( N_ELEMENTS(extra.xrange) EQ 2 ) THEN BEGIN
    xmin = extra.xrange[0]
    xmax = extra.xrange[1]
  ENDIF
  IF ( N_ELEMENTS(extra.yrange) EQ 2 ) THEN BEGIN
    ymin = extra.yrange[0]
    ymax = extra.yrange[1]
  ENDIF

  s = SIZE(u)
  t = SIZE(v)

        if s[0] ne 2 then begin
baduv:   message, 'U and V parameters must be 2D and same size.'
                endif
        if total(abs(s[0:2]-t[0:2])) ne 0 then goto,baduv
;
        if n_params(0) lt 3 then x = findgen(s[1]) else $
                if n_elements(x) ne s[1] then begin
badxy:                  message, 'X and Y arrays have incorrect size.'
                        endif
        if n_params(1) lt 4 then y = findgen(s[2]) else $
                if n_elements(y) ne s[2] then goto,badxy
;
        if n_elements(missing) le 0 then missing = 1.0e30
        if n_elements(length) le 0 then length = 1.0

        mag = sqrt(u^2.+v^2.)             ;magnitude.
                ;Subscripts of good elements
        nbad = 0                        ;# of missing points
        if n_elements(missing) gt 0 then begin
                good = where(mag lt missing)
                if keyword_set(dots) then bad = where(mag ge missing, nbad)
        endif else begin
                good = lindgen(n_elements(mag))
        endelse

    ugood = u[good]
    vgood = v[good]

    x0 = MIN(x)                     ;get scaling
    x1 = MAX(x)
    y0 = MIN(y)
    y1 = MAX(y)

	x_step = ( x1 - x0 ) / ( s[1] - 1.0 )   ; Convert to float. Integer math
	y_step = ( y1 - y0 ) / ( s[2] - 1.0 )   ; could result in divide by 0

	maxmag = MAX( [MAX( ABS( ugood / x_step ) ), MAX( ABS( vgood / y_step ) )] )

;    print, 'vec: ', x1, x0, y1, y0, x_step, y_step, maxmag

    IF ( N_ELEMENTS( maxval ) EQ 1 ) THEN maxmag = maxval

	sina = length * (ugood/maxmag)
	cosa = length * (vgood/maxmag)

;mh+
        COMMON max_vector, max_arrow, max_norm, max_step
        IF(maxmag NE 0.0) THEN BEGIN
          max_arrow = MAX(SQRT(sina^2 + cosa^2))
          max_norm  = MAX(SQRT(ugood^2 + vgood^2))
          max_step  = MAX([x_step,y_step])
        ENDIF
;mh-
;
        if n_elements(title) le 0 then title = ''
        ;--------------  plot to get axes  ---------------
        if n_elements(color) eq 0 then color = !p.color
        if n_elements(noclip) eq 0 then noclip = 1
        x_b0=x0-x_step
	x_b1=x1+x_step
	y_b0=y0-y_step
	y_b1=y1+y_step
        if (not keyword_set(overplot)) then begin
          if n_elements(position) eq 0 then begin
            plot,[x_b0,x_b1],[y_b1,y_b0],/nodata,  $
              color=color, _EXTRA = extra, xstyle=21, ystyle=21
          endif else begin
;            plot,[x_b0,x_b1],[y_b1,y_b0],/nodata,/xst,/yst, $
;              color=color, _EXTRA = extra
          endelse
        endif
        if n_elements(clip) eq 0 then $
            clip = [!x.crange[0],!y.crange[0],!x.crange[1],!y.crange[1]]
;
        r = .3                          ;len of arrow head
        angle = 22.5 * !dtor            ;Angle of arrowhead
        st = r * sin(angle)             ;sin 22.5 degs * length of head
        ct = r * cos(angle)
;
        for i=0L,n_elements(good)-1 do begin     ;Each point
                x0 = x[good[i] mod s[1]]        ;get coords of start & end
                dx = sina[i]
                x1 = x0 + dx
                y0 = y[good[i] / s[1]]
                dy = cosa[i]
                y1 = y0 + dy
;                k = good[i] MOD s[1]
;                l = good[i] / s[1]
;               IF (k MOD stride[0] EQ 0) AND (l MOD stride[1] EQ 0) THEN BEGIN
                xd=x_step
                yd=y_step
;                IF ( x1 LE xmax AND x1 GE xmin AND y1 LE ymax AND y1 GE ymin ) THEN BEGIN
                plots,[x0,x1,x1-(ct*dx/xd-st*dy/yd)*xd, $
                      x1,x1-(ct*dx/xd+st*dy/yd)*xd], $
                      [y0,y1,y1-(ct*dy/yd+st*dx/xd)*yd, $
                      y1,y1-(ct*dy/yd-st*dx/xd)*yd], $
                      color=color,thick=1.0*vthick,clip=clip,noclip=noclip
;                plots,[x0,x1,x1-(ct*dx/xd-st*dy/yd)*xd, $
;                      x1,x1-(ct*dx/xd+st*dy/yd)*xd], $
;                      [y0,y1,y1-(ct*dy/yd+st*dx/xd)*yd, $
;                      y1,y1-(ct*dy/yd-st*dx/xd)*yd], $
;                      color=vcolor,thick=vthick,clip=clip,noclip=noclip
;                ENDIF ELSE BEGIN
;                  IF ( x1 LT xmin ) THEN x1 = xmin 
;                  IF ( x1 GT xmax ) THEN x1 = xmax 
;                  IF ( y1 LT ymin ) THEN y1 = ymin 
;                  IF ( y1 GT ymax ) THEN y1 = ymax 

;                  PLOTS,[x0,x1], [y0,y1], $
;                      color=color,thick=1.0*vthick,clip=clip,noclip=noclip
;                  PLOTS,[x0,x1], [y0,y1], $
;                      color=vcolor,thick=vthick,clip=clip,noclip=noclip
;                ENDELSE
;               endif
        endfor
        if nbad gt 0 then $             ;Dots for missing?
                PLOTS, x[bad mod s[1]], y[bad / s[1]], psym=3, color=color, $
                       clip=clip,noclip=noclip
end
