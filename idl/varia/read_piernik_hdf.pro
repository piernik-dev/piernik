function READ_PIERNIK_HDF, nhdf=nhdf,var=var,xsize=xsize, ysize=ysize,zsize=zsize,dir=dir,$
                   prefix=prefix,timestr=timestr, xmin=xmin,ymin=ymin,zmin=zmin, $
                   xmax=xmax,ymax=ymax,zmax=zmax, mws=mws, xbl=xbl, xbr=xbr, ybl=ybl, $
                   ybr=ybr, zbl=zbl,zbr=zbr,allblocks=allblocks

  DEVICE, DECOMPOSED=0
  COMMON dims, nv,nx,ny,nz,nxd,nyd,nzd,nb
  COMMON mesh, x,y,z
  COMMON state,iter,t
  COMMON params, Eexpl,d0,e0,cx,cy,cz
  COMMON vars, u,b,wa

  frame   = string(0, format = '(i2.2)')+'_' $
           +string(0, format = '(i2.2)')+'_' $
           +string(0, format = '(i2.2)')+'_' $
           +string(nhdf, format = '(i4.4)')
 
  filename = dir+prefix+'_'+frame+'.hdf'
  
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
    iz = nzd/2
  endelse

  CALC_NSIZE, nxd=nxd, nyd=nyd, nzd=nzd, xsize=xsize, ysize=ysize, zsize=zsize, maxsize=mws

  filepref= dir+prefix  
  a = LOAD_DATA_HDF(dir,prefix, nhdf, var, $
            xcoord=x, ycoord = y, zcoord = z, $
            nxa=nxa,nya=nya,nza=nza, $
            time = t, xbl=xbl, xbr=xbr, $
            ybl=ybl, ybr=ybr, zbl=zbl, zbr=zbr, allblocks=allblocks)
  timestr = '    t='+strtrim(string(t,format='(f7.2)'),0)
  print, timestr
  print, min(a),max(a)
  return,a

end
