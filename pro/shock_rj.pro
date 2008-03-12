
PRO PLOTMHD,ff=first,lf=last,var=var,zoom=sf,log=llog,min=amin,max=amax,fix=run

; Written by: M. Hanasz, December 2005 - January 2006

  COMMON dims, nv,nx,ny,nz,nxd,nyd,nzd,nb
  COMMON mesh, x,y,z
  COMMON state,iter,t
  COMMON params, Eexpl,d0,e0,cx,cy,cz
  COMMON vars, u,b,wa

  DEVICE,DECOMPOSED=0,RETAIN=2
;  LOADCT,3

  ncolors = !D.N_COLORS -1
  white = 255
  plots  ='x'
  slices ='n'
  blocks ='n'
  dump   ='n'

;  sf = 8

   dir = '/raid20/kowalik/jeans_in/'
   dir = '/raid20/kowalik/gammie/'
   dir = '../run/'
;   dir = '/raid20/mhanasz/runs/crdynamo-d-r100/'
;  dir = '/scr/symulacje/cr-tests/crshear-res/'
    prefix = 'cr_shear_r10'
    prefix = 'keplerian_disk_n4c'
    prefix = 'jeans_in_dsk'
    prefix = '2d_hydro_tst'
    prefix = 'shock1d_r' + run
;    prefix = 'slab_uns'
;    prefix = 'shock2d_tst'
;    prefix = 'otvortex_rs6'
; first = 0
;  last  = 0
  freq = 1

;  var='dens'

;  ylog=1
 if(var EQ 'dens' OR var EQ 'ener') then ylog=1

  if(var EQ 'dns1') then begin
    amin = 0.0
    amax = 0.0
  endif  
   if(var EQ 'dns2') then begin
    amin = 0.01
    amax = 2.2
  endif
  if(var EQ 'vy_1') then begin
    amin = -0.05
    amax=  0.05
  endif
  if(var EQ 'dens') then begin
    amin= 0.0; 0.499280 
    amax= 0.0; 0.51 ; 1.3
  endif
  if(var EQ 'ener') then begin
    amin= 0.0 ;0.01
    amax= 0.0 ;75.  
  endif
  if(var EQ 'eint') then begin
    amin= 0.0 ;0.1
    amax= 0.0 ;1.8
  endif
  if(var EQ 'temp') then begin
    amin= 0.
    amax= 0.
  endif
  if(var EQ 'velx') then begin
    amin= 0.0 ;-1.8 
    amax= 0.0 ;0.8
  endif
  if(var EQ 'vely') then begin
    amin= 0.0
    amax= 0.0
  endif
  if(var EQ 'velz') then begin
    amin= -1.0
    amax= 1.0
  endif
  if(var EQ 'magx') then begin
    amin= 3.0 ;-0.5
    amax= 33.0 ;0.5
  endif
  if(var EQ 'magy') then begin
    amin= 0.0;
    amax= 12.0
  endif
  if(var EQ 'magz') then begin
    amin= 0.
    amax= 0.
  endif
  if(var EQ 'esrc') then begin
    amin=  -0.0
    amax=  0.0
  endif
  if(var EQ 'gpot') then begin
    amin=  0.0
    amax= 0.0
  endif
  if(var EQ 'ccur') then begin
    amin = 0.0
    amax = 0.05
  endif

  frame   = string(0, format = '(i2.2)')+'_' $
           +string(0, format = '(i2.2)')+'_' $
           +string(0, format = '(i2.2)')+'_' $
           +string(first, format = '(i4.4)')

  filename= dir+prefix+'_'+frame+'.hdf'  
  
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

  ix = nxd/2;+nb  ;nxd/2 ; nx/2 ; 3*nx/4 
;  ix = ix*pdims(0)
  iy = nyd/2 ; nyd/2 ; ny/2 ; 3*ny/4
;  iy = iy*pdims(1)

  if(nz eq 1) then begin 
    iz = 1
  endif else begin
    iz = nzd/2;+nb ;nz/2 ;nzd/2 ; nz/2 ; 3*nz/4 
  endelse


;  iz = iz*pdims(2)
;  iz = 1
;  ix = 1
;  iy = 1


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
;    WINDOW, 1, XSIZE=600, YSIZE=500
;    WINDOW, 3, XSIZE=600, YSIZE=500
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

        filepref= dir+prefix  
        a = LOAD_DATA_HDF(dir,prefix, step, var, $
                          xcoord = x, ycoord = y, zcoord = z, $
                          nxa=nxa,nya=nya,nza=nza, $
                          time = t)
    
;        filepref2 = dir2+prefix
;        a2 = LOAD_DATA_HDF(filepref2, step, var, $
;                          xcoord = x, ycoord = y, zcoord = z, $
;                          nxa=nxa,nya=nya,nza=nza, $
;                          time = t)
    
  timestr = '    t='+strtrim(string(t,format='(f7.2)'),0)
  if(llog EQ 1) then  a = alog10(a)
  if(amin EQ 0.0)  then amin = min(a)
  if(amax EQ 0.0)  then amax = max(a)

;  a = (a - shift(a,0,1))/(y(1)-y(0))
;  a(*,0) = a(*,1)
  yrange = [amin,amax]
  amin = min(a)
  amax = max(a)


IF(plots NE 'n') THEN BEGIN
  WSET, 1
  WSHOW, 1
  

  IF(plots EQ 'x') THEN BEGIN
    px  =  a(*,iy-1,iz-1)
    px = reform(px)
    pxr= rotate(px,2) 
;    px2 = a2(*,iy-1,iz-1)
;    PLOT,  x, (px*x)^2, line=0, title=var+timestr, $
    PLOT,  x, px, line=0, title=var+timestr, $
              xstyle=1,yrange=yrange,PSYM=4, SYMSIZE=0.5;, $
;              yrange=yrange,ystyle=1 ,ylog=ylog 
;plot,x,px, xrange=[-1,101], yrange=[0.98*amin,1.02*amax], PSYM=4, SYMSIZE=0.5, title=var

;    OPLOT, x, px2, thick=2
;    oplot, x, pxr, line =1,color=254

;    WSET, 3
;    WSHOW, 3
;    plot, x, px-px2

;    print, max(px-px2)

  ENDIF ELSE IF (plots EQ 'y') THEN BEGIN
    py  =  a(ix-1,*,iz-1)
    py = reform(py);, nyd)
    pyr= rotate(py,2) 
;    py2 = a2(ix-1,*,iz-1)
    PLOT,  y, py, line=0, title=var+timestr, $
            xstyle=1, $
            yrange=yrange,ystyle=1 ,ylog=ylog
;    OPLOT, y, py2, thick=2
;    oplot, y, pyr, line =1,color=254

;    WSET, 3
;    WSHOW, 3
;    plot, y, py-py2

;    print, max(py-py2)

  ENDIF ELSE IF (plots EQ 'z') THEN BEGIN
    pz  =  a(ix-1,iy-1,*)
    pz = reform(pz, nzd)
    pzr= rotate(pz,2) 
;    pz2 = a2(ix-1,iy-1,*)
    PLOT,  z, pz, line=0, title=var+timestr, $
            xstyle=1, $
            yrange=yrange,ystyle=1 ,ylog=ylog 
;    OPLOT, z, pz2, thick=2
;    oplot, z,pzr, line =1,color=254
   
;    WSET, 3
;    WSHOW, 3
;    plot, z, pz-pz2

;    print, max(pz-pz2)

  ENDIF

ENDIF

IF(slices EQ 'y') THEN BEGIN
   IF(dump NE 'y') THEN BEGIN
    WSET, 2
;    WSHOW, 2
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


 if(blocks EQ 'y') then begin
; Linie pionowe
  for pc=0,pdims(1)-1 do begin
    image_yz(0+pc*ny*sf,*)       = white    ; skrajna lewa
    image_yz((pc+1)*ny*sf-1,*)   = white    ; skrajna prawa
    if(nb NE 0) then begin
      image_yz(nbs-1+pc*ny*sf,*)   = white    ; lewe  brzegi bloku 
      image_yz((pc+1)*ny*sf-nbs,*) = white    ; prawe brzegi bloku
    endif
  end 

; Linie poziome
  for pc=0,pdims(2)-1 do begin
    image_yz(*,0+pc*nz*sf)       = white    ; skrajna lewa
    image_yz(*,(pc+1)*nz*sf-1)   = white    ; skrajna prawa
    if(nb NE 0) then begin
      image_yz(*,nbs-1+pc*nz*sf)   = white    ; lewe  brzegi bloku
      image_yz(*,(pc+1)*nz*sf-nbs) = white    ; prawe brzegi bloku
    endif
  end 
  endif
  IF(dump NE 'y') THEN BEGIN
     TV, image_yz,1,sf*nya+sf+1
  ENDIF ELSE BEGIN
     WRITE_PNG,yz_name, image_yz, red,green,blue
  ENDELSE

  image_xz = REFORM(a(*,iy,*),nxa,nza)
  image_xz = REBIN(image_xz,sf*nxa,sf*nza,/sample)
  image_xz = BYTSCL(image_xz,min=amin,max=amax,top=254)

 if(blocks EQ 'y') then begin
; Linie pionowe
  for pc=0,pdims(0)-1 do begin
    image_xz(0+pc*nx*sf,*)       = white    ; skrajna lewa
    image_xz((pc+1)*nx*sf-1,*)   = white    ; skrajna prawa
    if(nb NE 0) then begin
      image_xz(nbs-1+pc*nx*sf,*)   = white    ; lewe  brzegi bloku
      image_xz((pc+1)*nx*sf-nbs,*) = white    ; prawe brzegi bloku
    endif
  end 

; Linie poziome
  for pc=0,pdims(2)-1 do begin
    image_xz(*,0+pc*nz*sf)       = white    ; skrajna lewa
    image_xz(*,(pc+1)*nz*sf-1)   = white    ; skrajna prawa
    if(nb NE 0) then begin
      image_xz(*,nbs-1+pc*nz*sf)   = white    ; lewe  brzegi bloku
      image_xz(*,(pc+1)*nz*sf-nbs) = white    ; prawe brzegi bloku
    endif
  end 
  endif
   IF(dump NE 'y') THEN BEGIN
     TV, image_xz, sf*nya+sf+1,sf*nya+sf+1
   ENDIF ELSE BEGIN
     WRITE_PNG,xz_name, image_xz, red,green,blue
   ENDELSE

  endif

  image_xy = REFORM(a(*,*,iz-1),nxa,nya)
  image_xy = REBIN(image_xy,sf*nxa,sf*nya,/sample)
  image_xy = BYTSCL(image_xy,min=amin,max=amax,top=254)

 if(blocks EQ 'y') then begin
; Linie pionowe
  for pc=0,pdims(0)-1 do begin
    image_xy(0+pc*nx*sf,*)       = white    ; skrajna lewa
    image_xy((pc+1)*nx*sf-1,*)   = white    ; skrajna prawa
    if(nb NE 0) then begin
      image_xy(nbs-1+pc*nx*sf,*)   = white    ; lewe  brzegi bloku
      image_xy((pc+1)*nx*sf-nbs,*) = white    ; prawe brzegi bloku
    endif
  end 

; Linie poziome
  for pc=0,pdims(1)-1 do begin
    image_xy(*,0+pc*ny*sf)       = white    ; skrajna lewa
    image_xy(*,(pc+1)*ny*sf-1)   = white    ; skrajna prawa
    if(nb NE 0) then begin
      image_xy(*,nbs-1+pc*ny*sf)   = white    ; lewe  brzegi bloku
      image_xy(*,(pc+1)*ny*sf-nbs) = white    ; prawe brzegi bloku
    endif
  end 
;  image_xy = rotate(image_xy,-45)
  endif

   IF(dump NE 'y') THEN BEGIN
     TV, image_xy,1,1
   ENDIF ELSE BEGIN
     WRITE_PNG,xy_name, image_xy, red,green,blue
   ENDELSE

  nxa=20
  nya=20
;op = RANDOMU(S,25,50,200)
;FOR I=0, 10 DO op = SMOOTH(op, 3)
;op = BYTSCL(a,TOP=25)
;SHADE_VOLUME, a, 10.0 ,V,P,/LOW
;SCALE3, XRANGE=[0,1], YRANGE=[0,2], ZRANGE=[0,8]
;img = POLYSHADE(V,P,/T3D)
;WINDOW,10
;WSET, 10
;img = PROJECT_VOL(a,100,100,100,DEPTH_Q=0.7,OPAQUE=qp, TRANS=(!P.T))
;print, max(img),min(img)
;TVSCL, img
;VOLUME,a

endif

endfor

;set_plot,'ps'
;device,filename=var+'.ps'

;plot,findgen(nx)*100/(nx-1),REFORM(a(*,1,0),nx), xrange=[-10,110], yrange=[amin,amax], PSYM=4, SYMSIZE=0.5, title=var
;device,/close
;set_plot,'x'
end

pro shock_rj,f=frm,run=run
WINDOW,1, XSIZE=900,YSIZE=900
!P.MULTI=[0,3,3]
plotmhd,zoom=1.0,log=0,ff=frm,lf=frm,var='dens',fix=run
plotmhd,zoom=1.0,log=0,ff=frm,lf=frm,var='eint',fix=run
plotmhd,zoom=1.0,log=0,ff=frm,lf=frm,var='ener',fix=run
plotmhd,zoom=1.0,log=0,ff=frm,lf=frm,var='velx',fix=run
plotmhd,zoom=1.0,log=0,ff=frm,lf=frm,var='vely',fix=run
plotmhd,zoom=1.0,log=0,ff=frm,lf=frm,var='velz',fix=run
plotmhd,zoom=1.0,log=0,ff=frm,lf=frm,var='magy',fix=run
plotmhd,zoom=1.0,log=0,ff=frm,lf=frm,var='magz',fix=run
plotmhd,zoom=1.0,log=0,ff=frm,lf=frm,var='magx',fix=run
end


