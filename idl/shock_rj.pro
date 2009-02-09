
PRO PLOTMHD,ff=first,lf=last,var=var,zoom=sf,log=llog,min=amin,max=amax,fix=run

; Written by: M. Hanasz, December 2005 - January 2006

  COMMON dims, nv,nx,ny,nz,nxd,nyd,nzd,nb
  COMMON mesh, x,y,z
  COMMON state,iter,t
  COMMON params, Eexpl,d0,e0,cx,cy,cz
  COMMON vars, u,b,wa

  ncolors = !D.N_COLORS -1
  white = 255
  plots  ='x'

  dir = '../runs/shock1d_rj/'
  prefix = 'shock1d_r' + run
  freq = 1
  step = first

  frame   = string(0, format = '(i2.2)')+'_' $
           +string(0, format = '(i2.2)')+'_' $
           +string(0, format = '(i2.2)')+'_' $
           +string(first, format = '(i4.4)')

  filename= dir+prefix+'_'+frame+'.hdf'  
  
  LOAD_DIMS_HDF, filename, pdims=pdims, pcoords=pcoords, dims=dims, $
                           nxd=nxd,nyd=nyd,nzd=nzd, nxb=nxb,nyb=nyb,nzb=nzb, nb=nb, $
                           xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, zmin=zmin, zmax=zmax 
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
            time = t,allblocks='y')
    b = LOAD_DATA_HDF(dir,prefix, step, 'magy', $
            xcoord = x, ycoord = y, zcoord = z, $
            nxa=nxa,nya=nya,nza=nza, $
            time = t, allblocks='y')
    a = atan(a/(b+1.e-5))
  endif else begin
    a = LOAD_DATA_HDF(dir,prefix, step, var, $
            xcoord = x, ycoord = y, zcoord = z, $
            nxa=nxa,nya=nya,nza=nza, $
            time = t,allblock='y')
  endelse  
  timestr = '    t='+strtrim(string(t,format='(f7.2)'),0)
  if(llog EQ 1)    then a    = alog10(a)
  if(amin EQ 0.0)  then amin = min(a)
  if(amax EQ 0.0)  then amax = max(a)

  yrange = [amin,amax]
  amin = min(a)
  amax = max(a)

  ss = size(a)
  if(ss[0] EQ 1) then begin
     px = reform(a)
  endif else begin
     px  =  a(*,iy-1)
     px = reform(px)
  endelse
  PLOT,  x, px, line=0, title=var+timestr, $
    xstyle=1,yrange=yrange,PSYM=4, SYMSIZE=0.5

end

pro shock_rj,f=frm,run=run

; usage:
;   f   -> number of frame (usually 1)
;   run -> run prefix: '1a','1b','2a' etc.

set_plot,'PS'
device,/ENCAPSULATED, filename='shock1d_r' + run + '.ps', XSIZE=12, YSIZE=9, /INCHES
!P.MULTI=[0,3,3,0,1]
plotmhd,zoom=1.0,log=0,ff=frm,lf=frm,var='den1',fix=run,min=0.0,max=0.0
plotmhd,zoom=1.0,log=0,ff=frm,lf=frm,var='vlx1',fix=run,min=0.0,max=0.0
plotmhd,zoom=1.0,log=0,ff=frm,lf=frm,var='magy',fix=run,min=0.0,max=0.0
plotmhd,zoom=1.0,log=0,ff=frm,lf=frm,var='ein1',fix=run,min=0.0,max=0.0
plotmhd,zoom=1.0,log=0,ff=frm,lf=frm,var='vly1',fix=run,min=0.0,max=0.0
plotmhd,zoom=1.0,log=0,ff=frm,lf=frm,var='magz',fix=run,min=0.0,max=0.0
plotmhd,zoom=1.0,log=0,ff=frm,lf=frm,var='ene1',fix=run,min=0.0,max=0.0
plotmhd,zoom=1.0,log=0,ff=frm,lf=frm,var='vlz1',fix=run,min=0.0,max=0.0
plotmhd,zoom=1.0,log=0,ff=frm,lf=frm,var='phi', fix=run,min=0.0,max=0.0
device,/close
set_plot,'X'
!P.MULTI = 0
end
