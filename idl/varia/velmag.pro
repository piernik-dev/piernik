PRO VELMAG, nvh,nvv,ff,crop 


  COMMON dims, nv,nx,ny,nz,nxd,nyd,nzd,nb
  COMMON mesh, x,y,z
  COMMON state,iter,t
  COMMON params, Eexpl,d0,e0,cx,cy,cz
  COMMON vars, u,b,wa
  ncolors = !D.N_COLORS -1
  ;device,true_color=24,decomposed=0,retain=2
;  device,decomposed=0,retain=2

 write = 'write'  ; set this to 'write' if You want to save images to file
 sf = 1

  xbl = 0
  xbr = 0
  ybl = 0
  ybr = 0
  zbl = 0
  zbr = 0

 plane='XY'
; nvh=15    ; liczba wektorkow w kier. horyzontalnym
; nvv=15    ; liczba wektorkow w kier. wertykalnym 

  dir = '/raid16/pawlaszek/TIposter/TIanomalousHQ/'
  dir = '/raid16/pawlaszek/TIposter/TIohmicHQ/'
  dir = '/mnt/xarth/codes/piernik-trunk/runs/TIanomalousHQ/'
 prefix = 'tearing_tst'


  var1='ccur'
  vecx='magx'
  vecy='magy'
  vecz='magz'

; var1='dens'
; vecx='velx'
; vecy='vely'
; vecz='velz'
  first=ff
  last =ff
  freq =1

  ylog=0
;  if(var EQ 'dens' OR var EQ 'ener') then ylog=1

  if(var1 EQ 'ccur') then begin
    amin= 1.e-2
    amax= 50.0
  endif
  if(var1 EQ 'dens') then begin
    amin= 0.6
    amax= 1.5
  endif

  frame   = string(0, format = '(i2.2)')+'_' $
           +string(0, format = '(i2.2)')+'_' $
           +string(0, format = '(i2.2)')+'_' $
           +string(first, format = '(i4.4)')
  filepref= dir+prefix
  filename = filepref+'_'+frame+'.hdf'


  LOAD_DIMS_HDF, filename, pdims=pdims, pcoords=pcoords, dims=dims, $
                           nxd,nyd,nzd, nxb,nyb,nzb, nb, $
                           xmin, xmax, ymin, ymax, zmin, zmax 
  if(nxb eq dims(0)*pdims(0) and nyb eq dims(1)*pdims(1) and nzb eq dims(2)*pdims(2)) then begin
    nx=nxb
    ny=nyb
    nz=nzb
    nb = 0
    dx = (xmax-xmin)/nxd
    dy = (ymax-ymin)/nyd
    dz = (zmax-zmin)/nzd
    nxds = sf*nxd
    nyds = sf*nyd
    nzds = sf*nzd
    nbs = sf*nb
  endif else begin
    nx = dims(0)*pdims(0)
    ny = dims(1)*pdims(1)
    nz = dims(2)*pdims(2)
    dx = (xmax-xmin)/nxd
    dy = (ymax-ymin)/nyd
    dz = (zmax-zmin)/nzd
    xmin = xmin-nb*dx
    xmax = xmax+nb*dx
    ymin = ymin-nb*dy
    ymax = ymax+nb*dy
    zmin = zmin-nb*dz
    zmax = zmax+nb*dz
    nxds = sf*nx
    nyds = sf*ny
    nzds = sf*nz
    nbs  = 0
  endelse


 case plane of

  'YZ': 	begin
 
        window,1,xsize=sf*ny,ysize=sf*nz, title='YZ axis.'

		end ;-- case 'YZ'
		
  'XY':		begin

;        window,1,xsize=sf*(nx+fix(float(nx)/float(nvh))), $
;                 ysize=sf*(ny+fix(float(ny)/float(nvv))), title='XY axis.'

         window,1,xsize=sf*nx,ysize=sf*ny, title='XY axis.'
         window,1,xsize=700,ysize=700, title='XY axis.'
		end ;-- case 'XY'

   'XZ':	begin
    
        window,1,xsize=sf*nx,ysize=sf*nz, title='XZ axis.'

                end

   endcase


  for step=first,last,freq do begin
;   for ziut=0,20,4 do begin

  framenum   = string(step, format = '(i4.4)')

  frame   = string(0, format = '(i2.2)')+'_' $
           +string(0, format = '(i2.2)')+'_' $
           +string(0, format = '(i2.2)')+'_' $
           +string(step, format = '(i4.4)')




;  filename= dir+'/'+prefix+'_'+frame+'.hdf'  
  filepref= dir+prefix
  filename = filepref+'_'+frame+'.hdf'

  nproc = pdims(2)
  a1   = LOAD_DATA_HDF(dir, prefix, step, var1,  xcoord = x, ycoord = y, zcoord = z, time = t, xbl=xbl, xbr=xbr, ybl=ybl, ybr=ybr, zbl=zbl, zbr=zbr)
  ax   = LOAD_DATA_HDF(dir, prefix, step, vecx,  xcoord = x, ycoord = y, zcoord = z, time = t, xbl=xbl, xbr=xbr, ybl=ybl, ybr=ybr, zbl=zbl, zbr=zbr)
  ay   = LOAD_DATA_HDF(dir, prefix, step, vecy,  xcoord = x, ycoord = y, zcoord = z, time = t, xbl=xbl, xbr=xbr, ybl=ybl, ybr=ybr, zbl=zbl, zbr=zbr)
  az   = LOAD_DATA_HDF(dir, prefix, step, vecz,  xcoord = x, ycoord = y, zcoord = z, time = t, xbl=xbl, xbr=xbr, ybl=ybl, ybr=ybr, zbl=zbl, zbr=zbr)

  timestr = '    t='+strtrim(string(t,format='(f7.2)'),0)

;  if(var1 EQ 'dens') then a1 = alog10(a1)
  print, min(a1),max(a1)
  if(amin EQ 0.0)  then amin = min(a1)
  if(amax EQ 0.0)  then amax = max(a1)

  var=var1
  print,''
  print, 'MIN(',var,') =',amin, '    Max(',var,') =',amax
  print,''

;  if(nxb eq dims(0)*pdims(0) and nyb eq dims(1)*pdims(1) and nzb eq dims(2)*pdims(2)) then begin
;    nx=nxb
;    ny=nyb
;    nz=nzb
;    nb = 0
;    dx = (xmax-xmin)/nxd
;    dy = (ymax-ymin)/nyd
;    dz = (zmax-zmin)/nzd
;    nxds = sf*nxd
;    nyds = sf*nyd
;    nzds = sf*nzd
;    nbs = sf*nb
;  endif else begin
;    nx = dims(0)*pdims(0)
;    ny = dims(1)*pdims(1)
;    nz = dims(2)*pdims(2)
;    dx = (xmax-xmin)/nxd
;    dy = (ymax-ymin)/nyd
;    dz = (zmax-zmin)/nzd
;    xmin = xmin-nb*dx
;    xmax = xmax+nb*dx
;    ymin = ymin-nb*dy
;    ymax = ymax+nb*dy
;    zmin = zmin-nb*dz
;    zmax = zmax+nb*dz
;    nxds = sf*nx
;    nyds = sf*ny
;    nzds = sf*nz
;    nbs  = 0
;  endelse

  cx = nx/2;-25+ziut
  cy = ny/2.
  cz = nz/2

  ix=cx
  iy=cy
  iz=cz

  print, ix,iy,iz

;nxds = sf*nxd
;nyds = sf*nyd
;nzds = sf*nzd
; nbs = sf*nb

lxb = nbs
rxb = nbs+nxds
lyb = nbs
ryb = nbs+nyds
lzb = nbs
rzb = nbs+nzds

white = 255
black = 0

 case plane of

  'YZ': 	begin
 ;-- YZ axis

;  window,1,xsize=sf*ny,ysize=sf*nz, title='YZ axis.'+timestr
  image_yz = REFORM(a1(ix,*,*),ny,nz)
  image_yz = REBIN(image_yz,sf*ny,sf*nz,/sample)
  image_yz = BYTSCL(image_yz,min=amin,max=amax,top=254)
;  image_yz(lyb,*) = white
;  image_yz(ryb-1,*) = white
;  image_yz(*,lzb) = white
;  image_yz(*,rzb-1) = white
;  TV, image_yz,0,0

  bb = sqrt(ax*ax + ay*ay + az*az)
  ay = ay/bb
  az = az/bb
  
  aay=reform(ay(ix,*,*), ny,nz) 
  aaz=reform(az(ix,*,*), ny,nz) 

  vh = congrid(aay,nvh,nvv,/interp)
  vv = congrid(aaz,nvh,nvv,/interp)
  xh = congrid(y,nvh,/interp)
  xv = congrid(z,nvv,/interp) 

  xh = -rotate(xh,2)
  xv = -rotate(xv,2)
  vh = -rotate(vh,2)
  vv = -rotate(vv,2)

  vectors, vh,vv,xh,xv,/noerase, $
           xrange=[ymin,ymax], yrange = [zmin,zmax],$
           position=[1.0,1.0],vthick=1,color=white
		
		end ;-- case 'YZ'
		
  'XY':		begin
;-- XY axis

;  window,1,xsize=sf*nx,ysize=sf*ny, title='XY axis.'+timestr
  image_xy = REFORM(a1(*,*,iz),nx,ny)
;  image_xy = TOTAL(a1,2)
  image_xy = REBIN(image_xy,sf*nx,sf*ny,/sample)
  image_xy = BYTSCL(image_xy,min=amin,max=amax,top=254)
;  image_xy(lxb,*) = white
;  image_xy(rxb-1,*) = white		
;  image_xy(*,lyb) = white
;  image_xy(*,ryb-1) = white
;  image_xy = rotate(image_xy,-45) 
  crop = fix(512./(nvv))-1 
  print, "crop with ",crop 
  TV, REBIN(image_xy,512,512),fix(512./(nvv))-1,fix(512./(nvv))-1

;   aa2=reform(a2(*,*,iz),nx,ny)
;   aay=reform(ay(*,*,iz),nx,ny)

;  bb = sqrt(ax*ax + ay*ay); + az*az)
;  ax = ax/bb
;  ay = ay/bb
  
  aax=reform(ax(*,*,iz), nx,ny) 
  aay=reform(ay(*,*,iz), nx,ny) 

  dx = x(2) - x(1)
  dy = y(2) - y(1)
  vh = congrid(aax,nvh,nvv,/interp)
  vv = congrid(aay,nvh,nvv,/interp)
  xh = congrid(x,nvh,/interp)
  xv = congrid(y,nvv,/interp)

  xh = rotate(xh,1)
  xv = rotate(xv,1)
  vh = -rotate(vh,2)
  vv = -rotate(vv,2)
  vh(nvh-1,nvv-1) = 0.0
  vv(nvh-1,nvv-1) = 0.0

  if(var EQ 'dens') then begin
      cl = black
  endif else begin
      cl = white
  endelse
  
  vectors, vh,vv,xh,xv,/noerase, /noclip, $
           xrange=[xmin,xmax], yrange = [ymin,ymax],$
;          position=[1.0,1.0],vthick=2,color=white
           position=[0.7314,0.7314],vthick=2,color=cl

		end ;-- case 'XY'

   'XZ':	begin
    
;   print, sf*nx, sf*nz
;   window,1,xsize=sf*nx,ysize=sf*nz, title='XZ axis.'+timestr
   image_xz = REFORM(a1(*,iy,*),nx,nz)
   image_xz = REBIN(image_xz,sf*nx,sf*nz,/sample)
   image_xz = BYTSCL(image_xz,min=amin,max=amax,top=254) 
;   image_xz = BYTSCL(image_xz)
  image_xz(lxb,*)   = white
  image_xz(rxb-1,*) = white
  image_xz(*,lzb)   = white
  image_xz(*,rzb-1) = white
;   TV, image_xz,0,0

   bb= sqrt(ax*ax+ay*ay+az*az)
   ax = ax/bb
   az = az/bb

   aax = reform(ax(*,iy,*),nx,nz)
   aaz = reform(az(*,iy,*),nx,nz)

  vh = congrid(aax,nvh,nvv,/interp)
  vv = congrid(aaz,nvh,nvv,/interp)
  xh = congrid(x,nvh,/interp)
  xv = congrid(z,nvv,/interp)

  xh = -rotate(xh,2)
  xv = -rotate(xv,2)
  vh = -rotate(vh,2)
  vv = -rotate(vv,2)
  
  vectors, vh,vv,xh,xv ,/noerase, $
           xrange=[0.85*xmin,0.85*xmax], yrange = [0.8*zmin,zmax*0.8],$
           position=[1.0,1.0],vthick=1,color=white,XSTYLE=4,YSTYLE=4;,$
	   XMARGIN=[0,0];, YMARGIN=[0,0]

		end ;-- case 'XZ'
		
   else:  print, 'Error in vector writing script!'
  endcase	    ;-- plane
  
;  if (write EQ 'write') then begin
;    file1 =var1+'_'+framenum
;    WRITE_PNG, dir+'anim/'+file1+'.png',TVRD(TRUE=1)
;    WRITE_PNG, './'+file1+'.png',TVRD(TRUE=1)
;  endif 

 end
 ;return, crop
 print, timestr
 END
