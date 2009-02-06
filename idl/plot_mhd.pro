PRO PLOT_MHD ;,first,last,var

; Written by: M. Hanasz, December 2005 - January 2006

  COMMON dims, nv,nx,ny,nz,nxd,nyd,nzd,nb
  COMMON mesh, x,y,z
  COMMON state,iter,t
  COMMON params, Eexpl,d0,e0,cx,cy,cz
  COMMON vars, u,b,wa
  ncolors = !D.N_COLORS -1
  white = 255
  plots  ='y'
  slices ='y'
  dump   ='n'

  compare = 'n'

  sf = 2

;  dir0 = '/home/mhanasz/work/piernik/trunk/obj-sedov1d/'
;  var0 = 

   DEVICE,DECOMPOSED=0,RETAIN=2

   dir1 = '../obj/'
   prefix1 = 'sedov_tst'
   var1='den1'

   dir0    = dir1
   prefix0 = prefix1
   var0    = var1


;===============================================================================
data_files = FINDFILE(dir1+'/'+prefix1+'_00_00_00_*.hdf', Count=n_files)
IF(n_files EQ 0) THEN GOTO, SKIP
;===============================================================================


  first = 0
  last  = n_files -1
  freq = 1

  var = var0


; if(var EQ 'dens' OR var EQ 'ener' ) then ylog=1


    amin= 0.0
    amax= 0.0
  if(var EQ 'den1') then begin
    amin= 0.0
    amax= 0.0
  endif
  if(var EQ 'ener') then begin
    amin= 0.0 
    amax= 0.0 
  endif
  if(var EQ 'encr') then begin
    amin= 0.0 
    amax= 0.0 
  endif
  if(var EQ 'eint') then begin
    amin= 0.0
    amax= 0.0 
  endif
  if(var EQ 'temp') then begin
    amin= 0.0
    amax= 0.0
  endif
  if(var EQ 'velx') then begin
    amin= 0.0 
    amax= 0.0 
  endif
  if(var EQ 'vely') then begin
    amin= 0.0
    amax= 0.0
  endif
  if(var EQ 'velz') then begin
    amin= 0.0
    amax= 0.0
  endif
  if(var EQ 'magx') then begin
    amin= 0.0 
    amax= 0.0 
  endif
  if(var EQ 'magy') then begin
    amin= 0.0
    amax= 0.0
  endif
  if(var EQ 'magz') then begin
    amin= 0.0
    amax= 0.0
  endif
  if(var EQ 'gpot') then begin
    amin= 0.0
    amax= 0.0
  endif

  frame   = string(0, format = '(i2.2)')+'_' $
           +string(0, format = '(i2.2)')+'_' $
           +string(0, format = '(i2.2)')+'_' $
           +string(first, format = '(i4.4)')
 
  filename = data_files(0)
  
  TVLCT, red, green, blue, /GET

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

  ix = nxd/2

  iy = nyd/2


  if(nz eq 1) then begin 
    iz = 1
  endif else begin
    iz = nzd/2;+nb ;nz/2 ;nzd/2 ; nz/2 ; 3*nz/4 
  endelse


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


  if(plots  NE 'n') then begin
    WINDOW, 1, XSIZE=600, YSIZE=500
  endif

  if(nz ne 1) then begin
  if(slices EQ 'y' AND dump NE 'y') then WINDOW, 2, xsize=sf*(nx*pdims(0)+ny*pdims(1)+1)+2, $
                                    ysize=sf*(nz*pdims(2)+ny*pdims(1)+1)+2, $
                                    title=var  ;+timestr
  endif else begin
  if(slices EQ 'y' AND dump NE 'y') then WINDOW, 2, ysize=sf*(ny*pdims(1)+1)+2, $
                                    xsize=sf*(nx*pdims(0)+1)+2, $
                                    title=var  ;+timestr

  endelse


  for step=first,last,freq do begin

        filepref= dir1+prefix1  
        a0 = LOAD_DATA_HDF(dir0,prefix0, step, var0, xcoord = x, ycoord = y, zcoord = z, $
                          nxa=nxa,nya=nya,nza=nza, $
                          time = t)
        a1 = LOAD_DATA_HDF(dir1,prefix1, step, var1, xcoord = x, ycoord = y, zcoord = z, $
                          nxa=nxa,nya=nya,nza=nza, $
                          time = t)

        IF(compare EQ 'y') THEN a = a1-a0 ELSE a = a0 
        
  timestr = '    t='+strtrim(string(t,format='(e10.4)'),0)

  if(amin EQ 0.0)  then amin = min(a)
  if(amax EQ 0.0)  then amax = max(a)
;  a = alog(a+1.e-15)
  yrange = [amin,amax]
  amin = min(a)
  amax = max(a)


IF(plots NE 'n') THEN BEGIN
  WSET, 1
  

  IF(plots EQ 'x') THEN BEGIN
    px  =  a(*,iy-1,iz-1)
    px = reform(px)
    pxr= rotate(px,2) 
    PLOT,  x, px, line=0, title=var+timestr, $
              xstyle=1;, $
    oplot, x, pxr, line =1,color=254

  ENDIF ELSE IF (plots EQ 'y') THEN BEGIN
    py  =  a(ix-1,*,iz-1)
    py = reform(py);, nyd)
    pyr= rotate(py,2) 
    PLOT,  y, py, line=0, title=var+timestr, $
            xstyle=1, $
            yrange=yrange,ystyle=1 ,ylog=ylog
    oplot, y, pyr, line =1,color=254

  ENDIF ELSE IF (plots EQ 'z') THEN BEGIN
    pz  =  a(ix-1,iy-1,*)
    help, pz
    pz = reform(pz)
    pzr= rotate(pz,2) 

    PLOT,  z, pz, line=0, title=var+timestr, $
            xstyle=1, $
            yrange=yrange,ystyle=1 ,ylog=ylog 
    oplot, z,pzr, line =1,color=254
   

  ENDIF

ENDIF

IF(slices EQ 'y') THEN BEGIN
   IF(dump NE 'y') THEN BEGIN
     WSET, 2
   ENDIF ELSE BEGIN
    xy_name = 'xy'+var+'_'+ prefix +'_'+ string(step,FORMAT='(I4.4)') +'.png'
    xz_name = 'xz'+var+'_'+ prefix +'_'+ string(step,FORMAT='(I4.4)') +'.png'
    yz_name = 'yz'+var+'_'+ prefix +'_'+ string(step,FORMAT='(I4.4)') +'.png'
   ENDELSE

   print, 'MIN(',var,') =',min(a), '    Max(',var,') =',max(a),  $
                                   '    Tot(',var,') =',total(a) 
  print,''

  if(nza NE 1)  then begin

  image_yz = REFORM(a(ix,*,*),nya,nza)
  image_yz = REBIN(image_yz,sf*(nya),sf*(nza),/sample)
  image_yz = BYTSCL(image_yz,min=amin,max=amax,top=254)


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

  IF(dump NE 'y') THEN BEGIN
     TV, image_yz,1,sf*nya+sf+1
  ENDIF ELSE BEGIN
     WRITE_PNG,yz_name, image_yz, red,green,blue
  ENDELSE

  image_xz = REFORM(a(*,iy,*),nxa,nza)
  image_xz = REBIN(image_xz,sf*nxa,sf*nza,/sample)
  image_xz = BYTSCL(image_xz,min=amin,max=amax,top=254)

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

   IF(dump NE 'y') THEN BEGIN
     TV, image_xz, sf*nya+sf+1,sf*nya+sf+1
   ENDIF ELSE BEGIN
     WRITE_PNG,xz_name, image_xz, red,green,blue
   ENDELSE

  endif

  image_xy = REFORM(a(*,*,iz-1),nxa,nya)
  image_xy = REBIN(image_xy,sf*nxa,sf*nya,/sample)
  image_xy = BYTSCL(image_xy,min=amin,max=amax,top=254)

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
;  image_xy = rotate(image_xy,-45)


   IF(dump NE 'y') THEN BEGIN
     TV, image_xy,1,1
   ENDIF ELSE BEGIN
     WRITE_PNG,xy_name, image_xy, red,green,blue
   ENDELSE


endif

endfor

SKIP:

end
