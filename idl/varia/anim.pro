pro anim,dir=dir,var=var,min=amin,max=amax

n = 52
WINDOW,1, XSIZE=800, YSIZE=961
print, var+string(0, format = '(i3.3)')+'.png'
for i = 1, n do begin
  im = PLOT_MHD(dir=dir,var=var,ff=i,lf=i)
  xy = alog10(im(*,*,79))
  xz = alog10(im(*,399,*))
  TV, BYTSCL(xy,min=amin,max=amax,top=254)
  TV, BYTSCL(xz,min=amin,max=amax,top=254), 0, 801
  WRITE_PNG, var+string(i, format = '(i3.3)')+'.png', TVRD(TRUE=1)
endfor
end
