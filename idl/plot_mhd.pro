PRO PLOT_MHD

; Written by: M. Hanasz, December 2005 - January 2006

  COMMON dims, nv,nx,ny,nz,nxd,nyd,nzd,nb
  COMMON mesh, x,y,z
  COMMON state,iter,t
  COMMON params, Eexpl,d0,e0,cx,cy,cz
  COMMON vars, u,b,wa
  ncolors = !D.N_COLORS -1
  white = 255
  amin = 0. ;0.197 ;2.9 ;-1.0e-10
  amax = 1.2/20;0.203 ;3.5 ;e10
  extrema = 'f'  ; l-local, f-fixed, n-normalized

  plots  ='n'		; direction for plotting
  oplot = 'n'		; 'y' or 'n' or 'c'=yes and use previous plot
  plotwin = 7		; number of window for plotting
  slices ='y'		; 'y' - yes, 'n' - no slices
  slicwin = 4		; number of window for drawing slices
  linie = 'n'		; lines standing for borders of blocks/phys_domains
  compare = 'n'		; n-null, d,e-differ, s,t-sum, v,x-divide, o-oplot
  scalefaq=1		; scale factor for comparing data
  oplotfaq=1		; scale factor for oplotted data
  colourbar='y'		; colorbar under the slice
  pngwrite= 'n'		; 'y' -yes, 'n' - no for writing into png file
  nsim = 1		; number of data into one plot with oplots
  sf = 1		; scale factor to magnifying slices
  plotwait= 0.0001	; time to wait between plots
  dirstep = 200000	; stepnumber for changing dir into dir2

;DUSTKEPLER i pochodne
;   dir = '/raid16/wolt/dustkepler/test17/'
;   dir = '/scrh1/wolt/runkeplerorg/smout20so-iso-none/'
   dir = '/scrh1/wolt/dustkepler/smout21ss-gp-simx-colm-10/'
;   dir = '/raid_balaur/wolt/runkeplerorg/smout17db-smd-0/dbg/'
;   prefix = 'keplsupp_tst'
   prefix = 'kepler_tst'
;   prefix = 'dustkepler_tst'
;STREAMING INSTABILITY
;   dir = '/scrh1/wolt/stream/linb3/'
;   prefix = 'streaming_lnB'
;   dir = '/scrh1/wolt/slab2/'
;   prefix = 'slab_uns'

 dir1=dir 
 dir2 = dir ;'/scrh/wolt/multipiernik/runkeplerorg/supp1j/'
 dir3 = '/scrh/wolt/multipiernik/runkeplerorg/supp2j/'
 dir4 = '/scrh/wolt/multipiernik/runkeplerorg/supp3j/'
; dir2 = '/scrh/wolt/multipiernik/runkepler8suppnon/'
; dir2 = '/scrh/wolt/multipiernik/runkeplsupp1/simpledens/adiab/'
; dir2 = '/scrh1/wolt/keplsupptest/test1/simplevel/subsound/'
; dir3 = '/scrh1/wolt/keplsupptest/test1/simplemom/subsound/'
; dir3 = '/scrh/wolt/multipiernik/runkeplsupp1/simpledens2/adiab/'
; dir3 = '/scrh/wolt/multipiernik/runkeplsupp1/none/adiab/'
 prefix1=prefix
 prefix2=prefix
 prefix3=prefix
 prefix4=prefix

  first=3
  last =3
  freq =10

  var='den2'
  var1=var
  var2='vrt1'
  var3=var
  var4=var
  step2=0

  unitfactor=1.0 ;9.1828e19	

  ylog=0
;  if(var EQ 'dens' OR var EQ 'ener') then ylog=1


  frame   = string(0, format = '(i2.2)')+'_' $
           +string(0, format = '(i2.2)')+'_' $
           +string(0, format = '(i2.2)')+'_' $
           +string(first, format = '(i4.4)')
 
  filename = data_files(0)

  LOAD_DIMS_HDF, filename, pdims=pdims, pcoords=pcoords, dims=dims, $
                           nxd,nyd,nzd, nxb,nyb,nzb, nb, $
                           xmin, xmax, ymin, ymax, zmin, zmax 
  nproc = pdims(0)*pdims(1)*pdims(2)

  if(nxb eq dims(0)) then begin
    nx=nxb
    ny=nyb
    nz=nzb
    nb = 0
  endif else begin
    nx = dims(0)
    ny = dims(1)
    nz = dims(2)
  endelse

  ix = nx/2;+nb  ;nxd/2 ; nx/2 ; 3*nx/4 
  ix = ix*pdims(0)
  iy = ny/2;+nb ; nyd/2 ; ny/2 ; 3*ny/4
  iy = iy*pdims(1)
  iy = nyd/2  


  if(nz eq 1) then begin 
  iz = 1
  endif else begin
  iz = nz/2;+nb ;nz/2 ;nzd/2 ; nz/2 ; 3*nz/4 
  iz = iz*pdims(2)
  endelse
;  iz = 1
  ix = 1
  iy = 1


nxds = sf*nx*pdims(0)
nyds = sf*ny*pdims(1)
nzds = sf*nz*pdims(2)
nzbs = sf*nzb
nbs = sf*nb

lxb = nbs
rxb = nxds-nbs
lyb = nbs
ryb = nyds-nbs
lzb = nbs
rzb = nzds-nbs

if(nxb eq dims(0)) then begin
  lzb = nbs
  rzb = nbs+nzds
endif else begin
  lzb = nbs
  rzb = nbs+nzds
endelse


  if((plots  NE 'n') AND (oplot NE 'c'))  then WINDOW, plotwin, title='plot of '+var+' - window '+string(plotwin,format='(i2.2)')
  if(nz ne 1) then begin
  if(slices EQ 'y') then WINDOW, slicwin, xsize=sf*(nx*pdims(0)+ny*pdims(1)+1)+2, $
                                    ysize=sf*(nz*pdims(2)+ny*pdims(1)+1)+2 $
        ,title=var+' - window '+string(slicwin,format='(i2.2)')  ;+timestr
  endif else begin
  if(colourbar EQ 'y') then begin
  if(slices EQ 'y') then WINDOW, slicwin, ysize=sf*(ny*pdims(1)*1.25+1)+2, $
                                    xsize=sf*(nx*pdims(0)+1)+2, $
                                    title=var+' - window '+string(slicwin,format='(i2.2)')  ;+timestr
  endif else begin
  if(slices EQ 'y') then WINDOW, slicwin, ysize=sf*(ny*pdims(1)+1)+2, $
                                    xsize=sf*(nx*pdims(0)+1)+2, $
                                    title=var+' - window '+string(slicwin,format='(i2.2)')  ;+timestr

  endelse
  endelse


  for step=first,last,freq do begin
   if(step GE dirstep) then dir=dir2
     for istep=1,nsim,1 do begin
       if(istep EQ 1) then begin
         dir=dir1
         prefix=prefix1
         var=var1
       endif
       if(istep EQ 2) then begin
         dir=dir2
         prefix=prefix2
         var=var2
         oplot='c'
       endif
       if(istep EQ 3) then begin
         dir=dir3
         prefix=prefix3
         var=var3
       endif
       if(istep EQ 4) then begin
         dir=dir4
         prefix=prefix4
         var=var4
       endif
  

        filepref= dir+prefix  
        a = LOAD_DATA_HDF(filepref, step, var, $
                          xcoord = x, ycoord = y, zcoord = z, $
                          nxa=nxa,nya=nya,nza=nza, $
                          time = t)

  if(compare NE 'n') then begin
  filepref2 = dir2+prefix2
  if(compare EQ 'd') then begin
  a = a -LOAD_DATA_HDF(filepref2, step2, var2, $
                       xcoord = x, ycoord =y, zcoord = z, $
                       nxa=nxa,nya=nya,nza=nza, $
		       time = t)
  endif
  if(compare EQ 's') then begin
  a = a +LOAD_DATA_HDF(filepref2, step2, var2, $
                       xcoord = x, ycoord =u, zcoord = z, $
                       nxa=nxa,nya=nya,nza=nza, $
                       time = t)
  endif
  if(compare EQ 'e') then begin
  a = a*scalefaq -LOAD_DATA_HDF(filepref2, step, var2, $
                       xcoord = x, ycoord =y, zcoord = z, $
                       nxa=nxa,nya=nya,nza=nza, $
		       time = t)
  endif
  if(compare EQ 't') then begin
  a = a*scalefaq +LOAD_DATA_HDF(filepref2, step, var2, $
                       xcoord = x, ycoord =u, zcoord = z, $
                       nxa=nxa,nya=nya,nza=nza, $
                       time = t)
  endif
  if(compare EQ 'v') then begin
  a = a /LOAD_DATA_HDF(filepref2, step, var2, $
                       xcoord = x, ycoord =y, zcoord = z, $
                       nxa=nxa,nya=nya,nza=nza, $
		       time = t)
  endif
  if(compare EQ 'x') then begin
  a = a /LOAD_DATA_HDF(filepref2, step2, var2, $
                       xcoord = x, ycoord =u, zcoord = z, $
                       nxa=nxa,nya=nya,nza=nza, $
                       time = t)
  endif
  if(compare EQ 'o') then begin
  b = oplotfaq * LOAD_DATA_HDF(filepref2, step, var2, $
                       xcoord = x, ycoord =y, zcoord = z, $
                       nxa=nxa,nya=nya,nza=nza, $
		       time = t)
  endif
  endif
  timestr = '    t='+strtrim(string(t,format='(f7.2)'),0)+'  step='+strtrim(string(step,format='(i4)'),0)


;locmax = 0.0
;locmin = 0.0
;;;locmax = max(a((ix-nx/4):(ix+nx/4),(iy-ny/4):(iy+ny/4),(iz-nz/4):(iz+nz/4)))
;;;locmin = min(a((ix-nx/4):(ix+nx/4),(iy-ny/4):(iy+ny/4),(iz-nz/4):(iz+nz/4)))
;;;print, 'maximum = ',locmax,' minimum = ',locmin
;;;print, 'amplitude = ',locmax-locmin
;;;print, step, locmax-locmin, locmax, locmin, t
a=unitfactor*a
if(extrema EQ 'l') then begin
amin = min(a)
amax = max(a)
endif
if(extrema EQ 'f') then begin
;do nothing
endif
if(extrema EQ 'n') then begin
amin = min(a)
amax = max(a)
extrema = 'd'
endif
;  if(amin EQ 0.0)  then amin = min(a)
;  if(amax EQ 0.0)  then amax = max(a)
;amin = min(a)
;amax = 2.0
yrange=[amin,amax]

IF(plots NE 'n') THEN BEGIN
IF(((oplot NE 'y') OR (step EQ first)) AND (oplot NE 'c')) THEN BEGIN
;IF(oplot NE 'y') THEN BEGIN
  WSET, plotwin
;  WSHOW, plotwin

zero=replicate(0.,nx)

  IF(plots EQ 'x') THEN BEGIN
    px = a(*,iy-1,iz-1)
;    print, px
    PLOT,  x, px, line=0, title=plots+' plot of '+var+timestr, $
;              xstyle=1, xtitle=plots, ytitle='density: gas & '+strtrim(string(oplotfaq,format='(i2)'),0)+' x dust', $
              xstyle=1, ytitle=var, xtitle=plots, $
              yrange=yrange,ystyle=1 ,ylog=ylog 
;    OPLOT,x,zero,line=1
  IF(compare EQ 'o') THEN BEGIN
    px = b(*,iy-1,iz-1)
    OPLOT,  x, px, line=2
  ENDIF
  ENDIF ELSE IF (plots EQ 'y') THEN BEGIN
    py = a(ix-1,*,iz-1)
    PLOT,  y, py, line=0, title=plots+' plot of '+var+timestr, $
            xstyle=1, ytitle=var, xtitle=plots, $
            yrange=yrange,ystyle=1 ,ylog=ylog 
  IF(compare EQ 'o') THEN BEGIN
    py = b(ix-1,*,iz-1)
    OPLOT,  y, py, line=1
  ENDIF
  ENDIF ELSE IF (plots EQ 'z') THEN BEGIN
    pz = a(ix-1,iy-1,*)
;    print, pz
    PLOT,  z, pz, line=0, title=plots+' plot of '+var+timestr, $
            xstyle=1, ytitle=var, xtitle=plots, $
            yrange=yrange,ystyle=1 ,ylog=ylog 
  IF(compare EQ 'o') THEN BEGIN
    pz = b(ix-1,iy-1,*)
    OPLOT,  z, pz, line=1
  ENDIF
  ENDIF
ENDIF ELSE BEGIN
  WSET, plotwin
;  WSHOW, plotwin
  IF(plots EQ 'x') THEN BEGIN
    px = oplotfaq*a(*,iy-1,iz-1)
    OPLOT,  x, px, line=istep ;step-1
  ENDIF ELSE IF (plots EQ 'y') THEN BEGIN
    py = oplotfaq*a(ix-1,*,iz-1)
    OPLOT,  y, py, line=istep ;step-1
  ENDIF ELSE IF (plots EQ 'z') THEN BEGIN
    pz = oplotfaq*a(ix-1,iy-1,*)
    OPLOT,  z, pz, line=istep ;step-1
  ENDIF
ENDELSE
ENDIF

IF(slices EQ 'y') THEN BEGIN
  WSET, slicwin
;  WSHOW, slicwin

  if(compare NE 'n') then begin
  if(compare NE 'o') then begin
  if((compare EQ 'd') OR (compare EQ 'e')) then begin
  znak='-'
  endif
  if((compare EQ 's') OR (compare EQ 't')) then begin
  znak='+'
  endif
  if((compare EQ 'v') OR (compare EQ 'x')) then begin
  znak='/'
  endif
  print, 'MIN(',var,znak,var2,')=',amin,'(',min(a),') Max(',var,znak,var2,')=',amax,'(',max(a),')'
  endif else begin
  print, 'MIN(',var,') =',amin,'(',min(a),') Max(',var,') =',amax,'(',max(a),')'
  print, 'MIN(',var2,') =',amin,'(',min(b),') Max(',var2,') =',amax,'(',max(b),')'
  endelse
  endif else begin
  print, 'MIN(',var,') =',amin,'(',min(a),') Max(',var,') =',amax,'(',max(a),')'
  endelse
;;;  print,''

  if(nza NE 1)  then begin

  image_yz = REFORM(a(ix,*,*),nya,nza)
  image_yz = REBIN(image_yz,sf*(nya),sf*(nza),/sample)
  image_yz = BYTSCL(image_yz,min=amin,max=amax,top=254)


if(linie EQ 'y') then begin
; Linie pionowe
;  for pc=0,pdims(1)-1 do begin
;    image_yz(0+pc*ny*sf,*)       = white    ; skrajna lewa
;    image_yz((pc+1)*ny*sf-1,*)   = white    ; skrajna prawa
;    if(nb NE 0) then begin
;      image_yz(nbs-1+pc*ny*sf,*)   = white    ; lewe  brzegi bloku 
;      image_yz((pc+1)*ny*sf-nbs,*) = white    ; prawe brzegi bloku
;    endif
;  end 

; Linie poziome
;  for pc=0,pdims(2)-1 do begin
;    image_yz(*,0+pc*nz*sf)       = white    ; skrajna lewa
;    image_yz(*,(pc+1)*nz*sf-1)   = white    ; skrajna prawa
;    if(nb NE 0) then begin
;      image_yz(*,nbs-1+pc*nz*sf)   = white    ; lewe  brzegi bloku
;      image_yz(*,(pc+1)*nz*sf-nbs) = white    ; prawe brzegi bloku
;    endif
;  end 
endif

  TV, image_yz,1,sf*nya+sf+1
if(pngwrite EQ 'y') then begin
filenm = prefix+'yz_'+var+string(step, format = '(i4.4)')+'.png'
tvlct, red,green,blue,/get
write_png,filenm,image_yz, red,green,blue
endif

  image_xz = REFORM(a(*,iy,*),nxa,nza)
  image_xz = REBIN(image_xz,sf*nxa,sf*nza,/sample)
  image_xz = BYTSCL(image_xz,min=amin,max=amax,top=254)

if(linie EQ 'y') then begin
; Linie pionowe
;  for pc=0,pdims(0)-1 do begin
;    image_xz(0+pc*nx*sf,*)       = white    ; skrajna lewa
;    image_xz((pc+1)*nx*sf-1,*)   = white    ; skrajna prawa
;    if(nb NE 0) then begin
;      image_xz(nbs-1+pc*nx*sf,*)   = white    ; lewe  brzegi bloku
;      image_xz((pc+1)*nx*sf-nbs,*) = white    ; prawe brzegi bloku
;    endif
;  end 

; Linie poziome
;  for pc=0,pdims(2)-1 do begin
;    image_xz(*,0+pc*nz*sf)       = white    ; skrajna lewa
;    image_xz(*,(pc+1)*nz*sf-1)   = white    ; skrajna prawa
;    if(nb NE 0) then begin
;      image_xz(*,nbs-1+pc*nz*sf)   = white    ; lewe  brzegi bloku
;      image_xz(*,(pc+1)*nz*sf-nbs) = white    ; prawe brzegi bloku
;    endif
;  end 
endif

  TV, image_xz, sf*nya+sf+1,sf*nya+sf+1

  endif

  image_xy = REFORM(a(*,*,iz-1),nxa,nya)
  image_xy = REBIN(image_xy,sf*nxa,sf*nya,/sample)
  image_xy = BYTSCL(image_xy,min=amin,max=amax,top=254)

if(linie EQ 'y') then begin
; Linie pionowe
;  for pc=0,pdims(0)-1 do begin
;    image_xy(0+pc*nx*sf,*)       = white    ; skrajna lewa
;    image_xy((pc+1)*nx*sf-1,*)   = white    ; skrajna prawa
;    if(nb NE 0) then begin
;      image_xy(nbs-1+pc*nx*sf,*)   = white    ; lewe  brzegi bloku
;      image_xy((pc+1)*nx*sf-nbs,*) = white    ; prawe brzegi bloku
;    endif
;  end 

; Linie poziome
;  for pc=0,pdims(1)-1 do begin
;    image_xy(*,0+pc*ny*sf)       = white    ; skrajna lewa
;    image_xy(*,(pc+1)*ny*sf-1)   = white    ; skrajna prawa
;    if(nb NE 0) then begin
;      image_xy(*,nbs-1+pc*ny*sf)   = white    ; lewe  brzegi bloku
;      image_xy(*,(pc+1)*ny*sf-nbs) = white    ; prawe brzegi bloku
;    endif
;  end 
endif
;  image_xy = rotate(image_xy,-45)


if(colourbar EQ 'y') then begin
  TV, image_xy,1,sf*nxa*0.25+1+sf
if(nz EQ 1) then begin
cbposition=[0.1, 0.1, 0.9, 0.15]
;print, pc*ny*sf
;print, (pc+1)*ny+sf-1
bgrnd=image_xy(*,pdims(1)*ny*sf*0.25:((pdims(1)+1)*ny*sf-1)*0.25+sf-1)
bgrnd(*,*)=white
TV, bgrnd,1,1
endif else begin
cbposition=[0.55, 0.18, 0.95, 0.23]
bgrnd=image_xy
bgrnd(*,*)=white
TV, bgrnd,sf*nya+sf+1,1
endelse
COLORBAR, BOTTOM=0, CHARSIZE=0.8, COLOR=0, DIVISIONS=3, $
   FORMAT='(f4.2)', NCOLORS=255, FONT=0, POSITION=cbposition, $
   MINRANGE=amin, MAXRANGE=amax
endif else begin
  TV, image_xy,1,1
endelse
if(pngwrite EQ 'y') then begin
filenm = prefix+'xy_'+var+string(step, format = '(i4.4)')+'.png'
tvlct, red,green,blue,/get
write_png,filenm,image_xy ,red,green,blue
endif

  nxa=20
  nya=20

endif
if((nsim NE 1) AND (istep EQ nsim)) then begin
oplot='n'
endif
endfor

wait, plotwait

endfor

SKIP:

end
