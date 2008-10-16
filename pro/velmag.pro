PRO VELMAG


  COMMON dims, nv,nx,ny,nz,nxd,nyd,nzd,nb
  COMMON mesh, x,y,z
  COMMON state,iter,t
  COMMON params, Eexpl,d0,e0,cx,cy,cz
  COMMON vars, u,b,wa
  ncolors = !D.N_COLORS -1
  device,true_color=24,decomposed=0,retain=2
;  device,decomposed=0,retain=2

 write = 'write'  ; set this to 'write' if You want to save images to file
 sf = 4

 plane='YZ'
 nvh=30    ; liczba wektorkow w kier. horyzontalnym
 nvv=30    ; liczba wektorkow w kier. wertykalnym 

winnum=2

    dir='/home/wolt/somerun/quadnew/'
    prefix = 'quadrupole_tst'

  var1='maga'
  vecx='magx'
  vecy='magy'
  vecz='magz'

  first=0
  last =0
  freq =1
;  step = 150

  ylog=0
;  if(var EQ 'dens' OR var EQ 'ener') then ylog=1
amin=0.0
amax=1.0e-1
  if(var1 EQ 'dens') then begin
    amin= 0.02;-0.01
    amax= 2.5
  endif
  if(var1 EQ 'eint') then begin
    amin= 0.
    amax= 0.
  endif
  if(var1 EQ 'temp') then begin
    amin= 0.
    amax=15000.
  endif
  if(var1 EQ 'velx') then begin
    amin= 0.
    amax= 0.
  endif
  if(var1 EQ 'esrc') then begin
    amin=  -0.0
    amax=  0.0
  endif
  if(var1 EQ 'magx') then begin
    amin= -30. ;-50.
    amax= 30. ;50.
  endif
  if(var1 EQ 'magy') then begin
    amin= -1.2;-30. ;-50.
    amax= 1.2;30. ;50.
  endif
  if(var1 EQ 'magz') then begin
    amin= -3.0;-30. ;-50.
    amax= 5.0;30. ;50.
  endif
  if(var1 EQ 'esrc') then begin
    amin=  -0.0
    amax=  0.0
  endif
  if(var1 EQ 'curx') then begin
    amin= -0.0
    amax=  0.0
  endif
  if(var1 EQ 'cury') then begin
    amin= -0.0
    amax=  0.0
  endif
  if(var1 EQ 'curz') then begin
    amin= -0.0
    amax=  0.0
  endif
  if(var1 EQ 'ccur') then begin
    amin=  0.0;2.6e-6
    amax=  110.0;15;01; 2.33667
  endif



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

  LOAD_DIMS_HDF, filename, pdims=pdims, pcoords=pcoords, dims=dims, $
                           nxd,nyd,nzd, nxb,nyb,nzb, nb, $
                           xmin, xmax, ymin, ymax, zmin, zmax 
  nproc = pdims(2)


  a1   = LOAD_DATA_HDF(dir, prefix, step, var1,  xcoord = x, ycoord = y, zcoord = z, time = t)
  ax   = LOAD_DATA_HDF(dir, prefix, step, vecx,  xcoord = x, ycoord = y, zcoord = z, time = t)
  ay   = LOAD_DATA_HDF(dir, prefix, step, vecy,  xcoord = x, ycoord = y, zcoord = z, time = t)
  az   = LOAD_DATA_HDF(dir, prefix, step, vecz,  xcoord = x, ycoord = y, zcoord = z, time = t)

  timestr = '    t='+strtrim(string(t,format='(f7.2)'),0)

  if(amin EQ 0.0)  then amin = min(a1)
  if(amax EQ 0.0)  then amax = max(a1)


  var=var1
  print,''
  print, 'MIN(',var,') =',amin, '    Max(',var,') =',amax
  print,''

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

 case plane of

  'YZ': 	begin
 ;-- YZ axis

  window,winnum,xsize=sf*ny,ysize=sf*nz, title='YZ axis.'+timestr
  image_yz = REFORM(a1(ix,*,*),ny,nz)
  image_yz = REBIN(image_yz,sf*ny,sf*nz,/sample)
  image_yz = BYTSCL(image_yz,min=amin,max=amax,top=254)
  image_yz(lyb,*) = white
  image_yz(ryb-1,*) = white
  image_yz(*,lzb) = white
  image_yz(*,rzb-1) = white
  TV, image_yz,0,0

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
           position=[0.997,0.997],vthick=1,color=white
		
		end ;-- case 'YZ'
		
  'XY':		begin
;-- XY axis

  window,1,xsize=sf*nx,ysize=sf*ny, title='XY axis.'+timestr
  image_xy = REFORM(a1(*,*,iz),nx,ny)
  image_xy = REBIN(image_xy,sf*nx,sf*ny,/sample)
  image_xy = BYTSCL(image_xy,min=amin,max=amax,top=254)
  image_xy(lxb,*) = white
  image_xy(rxb-1,*) = white		
  image_xy(*,lyb) = white
  image_xy(*,ryb-1) = white
;  image_xy = rotate(image_xy,-45)
  TV, image_xy,0,0

;   aa2=reform(a2(*,*,iz),nx,ny)
;   aay=reform(ay(*,*,iz),nx,ny)

;  bb = sqrt(ax*ax + ay*ay); + az*az)
;  ax = ax/bb
;  ay = ay/bb
  
  aax=reform(ax(*,*,iz), nx,ny) 
  aay=reform(ay(*,*,iz), nx,ny) 

  vh = congrid(aax,nvh,nvv,/interp)
  vv = congrid(aay,nvh,nvv,/interp)
  xh = congrid(x,nvh,/interp)
  xv = congrid(y,nvv,/interp)

  xh = -rotate(xh,2)
  xv = -rotate(xv,2)
  vh = -rotate(vh,2)
  vv = -rotate(vv,2)
  
  vectors, vh,vv,xh,xv,/noerase, $
           xrange=[xmin,xmax], yrange = [ymin,ymax],$
           position=[0.997,0.997],vthick=1,color=white
		end ;-- case 'XY'

   'XZ':	begin
    
   print, sf*nx, sf*nz
   window,1,xsize=sf*nx,ysize=sf*nz, title='XZ axis.'+timestr
   image_xz = REFORM(a1(*,iy,*),nx,nz)
   image_xz = REBIN(image_xz,sf*nx,sf*nz,/sample)
   image_xz = BYTSCL(image_xz,min=amin,max=amax,top=254) 
;   image_xz = BYTSCL(image_xz)
  image_xz(lxb,*)   = white
  image_xz(rxb-1,*) = white
  image_xz(*,lzb)   = white
  image_xz(*,rzb-1) = white
   TV, image_xz,0,0

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
  
  if (write EQ 'write') then begin
    file1 =var1+'_'+framenum
    WRITE_PNG, dir+'movie/'+file1+'.png',TVRD(TRUE=1)
  endif 

 end

 print, timestr
  END
