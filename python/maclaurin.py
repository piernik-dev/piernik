#!/usr/bin/python

import tables as h5
import numpy as np
import matplotlib
matplotlib.use('cairo')      # choose output format
import pylab as P
from matplotlib.ticker import NullFormatter

h5f = h5.openFile("maclaurin_sph_0001.h5","r")

xmin = h5f.root._v_attrs.xmin[0]
xmax = h5f.root._v_attrs.xmax[0]
ymin = h5f.root._v_attrs.ymin[0]
ymax = h5f.root._v_attrs.ymax[0]
zmin = h5f.root._v_attrs.zmin[0]
zmax = h5f.root._v_attrs.zmax[0]

soln = h5f.root.gpot[:,:,:]
h5f.close()

nz,ny,nx = soln.shape
dx,dy,dz = (xmax-xmin)/nx,(ymax-ymin)/ny,(zmax-zmin)/nz

xl,xr = xmin+dx*0.5, xmax-dx*0.5
yl,yr = ymin+dy*0.5, ymax-dy*0.5
zl,zr = zmin+dz*0.5, zmax-dz*0.5

x = np.arange(xl,xr,dx)
y = np.arange(yl,yr,dy)
z = np.arange(zl,zr,dz)

n = nx*ny*nz

r   = np.zeros(n)
phi = np.zeros(n)

ind = 0
for i in range(0,nx-1):
   for j in range(0,ny-1):
      for k in range(0,nz-1):
         r[ind] = np.sqrt(x[i]**2+y[j]**2+z[k]**2)
         phi[ind] = soln[k,j,i]
         ind=ind+1

new_r = np.unique(r)
mean = np.zeros_like(new_r)
phi_0 = np.zeros_like(new_r)

for i in range(0,new_r.shape[0]-1):
   temp =  phi[np.where(r == new_r[i])]
   mean[i] = temp.mean()

fpiG = 8.38579043837416375E-007
dpiG = fpiG * 0.5
a   = 1.
rho = 1.

phi_0[np.where(new_r <= a)] = -dpiG * rho * (a**2 - 0.333333*new_r[np.where(new_r <= a)]**2)
phi_0[np.where(new_r > a)] = -fpiG * rho * a**3 / (3.0*new_r[np.where(new_r > a)])

# plotting ---------------------------------
# definitions for the axes
left, width = 0.1, 0.85

rect_hi = [left, 0.3 , width, 0.65]
rect_lo = [left, 0.05, width, 0.2]

P.figure(1,figsize=(8,8))

axhi = P.axes(rect_hi)
axlo   = P.axes(rect_lo)

axhi.set_title("Maclaurin spheroid $e=0,\; a_1=1,\;\\varrho_0=1$")
axhi.plot(new_r[1:-1],mean[1:-1],'g.',new_r[1:-1],phi_0[1:-1],'b')
axhi.xaxis.set_major_formatter( NullFormatter() )
axhi.set_ylabel('Gravitational potential')
axhi.legend( ('Numerical solution - $\\varphi$','Analitycal solution - $\\varphi_0$'), loc = 'lower right')

axlo.plot(new_r[1:-1],phi_0[1:-1]-mean[1:-1],'g.')
axlo.set_xlabel('Radius')
axlo.set_ylabel('$\\varphi_0 - \\varphi$')

#P.show()
P.draw()
P.savefig('maclaurin.png',facecolor='white')
