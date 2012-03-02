#!/usr/bin/python

import os

def Maclaurin_test(file):

   missing = []
   try: import tables as h5
   except ImportError: missing.append("PyTables")
   try: import numpy as np
   except ImportError: missing.append("NumPy")
   try: import matplotlib
   except ImportError: missing.append("matplotlib")

   if (len(missing) > 0):
      print "You must install the package(s) ",missing
      return

   matplotlib.use('cairo')      # choose output format
   import pylab as P
   from matplotlib.ticker import NullFormatter
   # This may sometimes help with font issues
   # from matplotlib import rc
   # rc('text',usetex=True)

   if (not os.path.isfile(file)):
      print "Cannot find ",file
      return

   try:
      h5f = h5.openFile(file,"r")
   except:
      print "Cannot open '"+file+"' as HDF5"
      return

   xmin = h5f.root._v_attrs.xmin[0]
   xmax = h5f.root._v_attrs.xmax[0]
   ymin = h5f.root._v_attrs.ymin[0]
   ymax = h5f.root._v_attrs.ymax[0]
   zmin = h5f.root._v_attrs.zmin[0]
   zmax = h5f.root._v_attrs.zmax[0]
   fpiG = h5f.root._v_attrs.fpiG[0]
   a    = h5f.root._v_attrs.a1[0]
   x0   = h5f.root._v_attrs.x0[0]
   y0   = h5f.root._v_attrs.y0[0]
   z0   = h5f.root._v_attrs.z0[0]

   try:
      soln    = h5f.root.gpot[:,:,:]
      phi0_3D = h5f.root.apot[:,:,:]
   except:
      print "Cannot find apot or gpot arrays"
      return
   h5f.close()

   nz,ny,nx = soln.shape
   dx,dy,dz = (xmax-xmin)/nx,(ymax-ymin)/ny,(zmax-zmin)/nz

   xl,xr = xmin+dx*0.5, xmax-dx*0.5
   yl,yr = ymin+dy*0.5, ymax-dy*0.5
   zl,zr = zmin+dz*0.5, zmax-dz*0.5

   x2 = (np.arange(xl,xr,dx) - x0)**2
   y2 = (np.arange(yl,yr,dy) - y0)**2
   z2 = (np.arange(zl,zr,dz) - z0)**2

   n = nx*ny*nz

   r     = np.zeros(n)
   phi   = np.zeros(n)
   phi0  = np.zeros(n)

   # *********************************************************************
   ind = 0
   ind_i = range(0,nx-1)
   ind_j = range(0,ny-1)
   ind_k = range(0,nz-1)
   for i in ind_i:
      for j in ind_j:
         for k in ind_k:
            r[ind]    = np.sqrt(x2[i]+y2[j]+z2[k])
            phi[ind]  = soln[k,j,i]
            phi0[ind] = phi0_3D[k,j,i]
            ind=ind+1

   new_r = np.unique(r)
   mean  = np.zeros_like(new_r)
   std   = np.zeros_like(new_r)
   phi_0 = np.zeros_like(new_r)

   for i in range(0,new_r.shape[0]-1):
      ind      = np.where(r == new_r[i])
      temp     = phi[ind]
      temp2    = phi0[ind]
      phi_0[i] = temp2.mean()
      mean[i]  = temp.mean()
      std[i]   = temp.std()
   # *********************************************************************
   # plotting ---------------------------------
   # definitions for the axes
   left, width = 0.15, 0.80

   rect_hi = [left, 0.3 , width, 0.65]
   rect_lo = [left, 0.05, width, 0.2]

   P.figure(1,figsize=(8,8))

   axhi = P.axes(rect_hi)
   axlo = P.axes(rect_lo)

   GM = fpiG/3.

   axhi.set_title("Maclaurin spheroid $e=0,\; a_1=1,\;\\varrho_0=1$")
   axhi.plot(new_r[1:-1],mean[1:-1]/GM,'g.',new_r[1:-1],phi_0[1:-1]/GM,'b')
   axhi.xaxis.set_major_formatter( NullFormatter() )
   axhi.set_ylabel('Gravitational potential / GM')
   axhi.legend( ('Numerical solution - $\\varphi$','Analytical solution - $\\varphi_0$'), loc = 'lower right')

   axlo.plot(new_r[1:-1],(phi_0[1:-1]-mean[1:-1])/GM,'g',new_r[1:-1],(phi_0[1:-1]-mean[1:-1]-std[1:-1])/GM,'r:',new_r[1:-1],(phi_0[1:-1]-mean[1:-1]+std[1:-1])/GM,'r:')
   axlo.legend( ('avg. difference', '+/- deviation'), loc = 'upper right')
   axlo.set_xlabel('Radius')
   axlo.set_ylabel('($\\varphi_0 - \\varphi$) / GM')

   #P.show()
   P.draw()
   P.savefig('maclaurin.png',facecolor='white')

if __name__ == "__main__":
   import sys
   if (len(sys.argv) <= 1):
      print "Usage : ",sys.argv[0]," maclaurin_test_hdf_file"
   else:
      Maclaurin_test(sys.argv[1])
      if (len(sys.argv) > 2):
         print "Ignored arguments: ",sys.argv[2:len(sys.argv)]
