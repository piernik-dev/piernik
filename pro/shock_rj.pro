
PRO PLOTMHD,ff=first,lf=last,var=var,zoom=sf,log=llog,min=amin,max=amax,fix=run

; Written by: M. Hanasz, December 2005 - January 2006

  COMMON dims, nv,nx,ny,nz,nxd,nyd,nzd,nb
  COMMON mesh, x,y,z
  COMMON state,iter,t
  COMMON params, Eexpl,d0,e0,cx,cy,cz
  COMMON vars, u,b,wa

  DEVICE,DECOMPOSED=0,RETAIN=2

  ncolors = !D.N_COLORS -1
  white = 255
  plots  ='x'

   dir = '../run/'
    prefix = 'shock1d_r' + run
  freq = 1
  step = first

  frame   = string(0, format = '(i2.2)')+'_' $
           +string(0, format = '(i2.2)')+'_' $
           +string(0, format = '(i2.2)')+'_' $
           +string(first, format = '(i4.4)')

  filename= dir+prefix+'_'+frame+'.hdf'  
  
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
  iz = nzd/2

  filepref= dir+prefix  

  if(var eq 'phi') then begin
    a = LOAD_DATA_HDF(dir,prefix, step, 'magz', $
            xcoord = x, ycoord = y, zcoord = z, $
            nxa=nxa,nya=nya,nza=nza, $
            time = t)
    b = LOAD_DATA_HDF(dir,prefix, step, 'magy', $
            xcoord = x, ycoord = y, zcoord = z, $
            nxa=nxa,nya=nya,nza=nza, $
            time = t)
    a = atan(a/b)
  endif else begin
    a = LOAD_DATA_HDF(dir,prefix, step, var, $
            xcoord = x, ycoord = y, zcoord = z, $
            nxa=nxa,nya=nya,nza=nza, $
            time = t)
  endelse  
  timestr = '    t='+strtrim(string(t,format='(f7.2)'),0)
  if(llog EQ 1)    then a    = alog10(a)
  if(amin EQ 0.0)  then amin = min(a)
  if(amax EQ 0.0)  then amax = max(a)

  yrange = [amin,amax]
  amin = min(a)
  amax = max(a)

  WSET, 1
  WSHOW, 1

  px  =  a(*,iy-1,iz-1)
  px = reform(px)
  PLOT,  x, px, line=0, title=var+timestr, $
    xstyle=1,yrange=yrange,PSYM=4, SYMSIZE=0.5

end

pro shock_rj,f=frm,run=run

; usage:
;   f   -> number of frame (usually 1)
;   run -> run prefix: '1a','1b','2a' etc.


WINDOW,1, XSIZE=900,YSIZE=900
!P.MULTI=[0,3,3]
plotmhd,zoom=1.0,log=0,ff=frm,lf=frm,var='dens',fix=run,min=0.0,max=0.0
plotmhd,zoom=1.0,log=0,ff=frm,lf=frm,var='eint',fix=run,min=0.0,max=0.0
plotmhd,zoom=1.0,log=0,ff=frm,lf=frm,var='ener',fix=run,min=0.0,max=0.0
plotmhd,zoom=1.0,log=0,ff=frm,lf=frm,var='velx',fix=run,min=0.0,max=0.0
plotmhd,zoom=1.0,log=0,ff=frm,lf=frm,var='vely',fix=run,min=0.0,max=0.0
plotmhd,zoom=1.0,log=0,ff=frm,lf=frm,var='velz',fix=run,min=0.0,max=0.0
plotmhd,zoom=1.0,log=0,ff=frm,lf=frm,var='magy',fix=run,min=0.0,max=0.0
plotmhd,zoom=1.0,log=0,ff=frm,lf=frm,var='magz',fix=run,min=0.0,max=0.0
plotmhd,zoom=1.0,log=0,ff=frm,lf=frm,var='phi', fix=run,min=0.0,max=0.0
end


