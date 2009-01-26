pro show,nhdf,var

;   nhdf   = 48
;   var    = 'den1'
   dir = '../obj/'
   prefix = 'kepler_tst'

;   xbl = 0
;   xbr = 19
;   ybl = 0
;   ybr = 19
;   zbl = 0
;   zbr = 3
   allblocks = 'y'

   max_window_size = 900

   a = read_piernik_hdf(nhdf=nhdf,dir=dir,prefix=prefix,var=var,xsize=xsize,ysize=ysize,zsize=zsize,     $
             xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, zmin=zmin, zmax=zmax, timestr=timestr,     $
             mws=max_window_size,xbl=xbl,ybl=ybl,zbl=zbl,xbr=xbr,ybr=ybr,zbr=zbr,allblocks=allblocks)

   s = size(a)
   nxd = s(1)
   nyd = s(2)
   nzd = s(3)

   ix = 1; nxd / 2
   iy = 1; nyd / 2
   iz = 1
   if(nzd ne 1) then iz = nzd/2

   help,a
   if(s(0) EQ 2) then begin
      xy = a
      window,1,XSIZE=xsize,YSIZE=ysize
      imdisp,xy,/AXIS,/ERASE,xrange=[xmin,xmax], yrange=[ymin,ymax], title="xy slice of "+var+" at "+timestr
   endif else begin
      xy = REFORM(a(*,*,iz))
      xz = REFORM(a(*,iy,*))
      yz = REFORM(a(ix,*,*))
      window,1,XSIZE=xsize,YSIZE=ysize
      imdisp,xy,/AXIS,/ERASE,xrange=[xmin,xmax], yrange=[ymin,ymax], title="xy slice of "+var+" at "+timestr
      window,2,XSIZE=xsize,YSIZE=zsize
      imdisp,xz,/AXIS,/ERASE,xrange=[xmin,xmax], yrange=[zmin,zmax], title="xz slice of "+var+" at "+timestr
      window,3,XSIZE=ysize,YSIZE=zsize
      imdisp,yz,/AXIS,/ERASE,xrange=[ymin,ymax], yrange=[zmin,zmax], title="yz slice of "+var+" at "+timestr
   endelse
   amax = max(a)
   amin = min(a)

   


;  colorbar, range=[-4,3], /VERTICAL, /RIGHT
end
