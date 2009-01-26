pro CALC_NSIZE, nxd=nxd, nyd=nyd, nzd=nzd, xsize=xsize, ysize=ysize, zsize=zsize, maxsize=maxsize

   nd   = max([nxd,nyd,nzd])

   r_x = float(nxd)/float(nd)
   r_y = float(nyd)/float(nd)
   r_z = float(nzd)/float(nd)

   xsize = fix(maxsize * r_x)
   ysize = fix(maxsize * r_y)
   zsize = fix(maxsize * r_z)

end
