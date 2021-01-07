#!/usr/bin/env python
import numpy as np

def fsym(vmin,vmax):
   vmx =        np.max([np.abs(vmin),np.abs(vmax)])
   vmn = -1.0 * np.max([np.abs(vmin),np.abs(vmax)])
   if vmn == vmx:
      vmn = vmn - 0.00001
      vmx = vmx + 0.00001
   return vmn, vmx

def scale_manage(sctype, xy, xz, yz, umin, umax, d2min, d2max):

    if (umin == 0.0 and umax == 0.0):
       vmin, vmax = d2min, d2max
    else:
       vmin, vmax = umin, umax

    if (sctype == '1' or sctype == 'symlin'):
           vmin, vmax = fsym(vmin,vmax)

    elif (sctype == '2' or sctype == 'log'):
       if (vmin > 0.0):
           vmin = np.log10(vmin)
       else:
           vmin = np.log10(min(np.min(xz,initial=np.inf,where=(xz>0.0)), np.min(xy,initial=np.inf,where=(xy>0.0)), np.min(yz,initial=np.inf,where=(yz>0.0))))
       if (vmax > 0.0):
           vmax = np.log10(vmax)
       else:
           vmax = -1.
       xy = np.log10(xy)
       xz = np.log10(xz)
       yz = np.log10(yz)
    elif (sctype == '3' or sctype == 'symlog'):
       if (umin > 0.0 and umax > 0.0):
           symmin = umin
           vmax = np.log10(umax/umin)
           vmin = np.log10(vmax)
       else:
           if (d2min*d2max > 0.0):
               smin, smax = d2min, d2max
           else:
               smin = min(np.min(xz, initial=np.inf,where=(xz>0.0)), np.min(xy, initial=np.inf,where=(xy>0.0)), np.min(yz, initial=np.inf,where=(yz>0.0)))
               smax = max(np.min(xz,initial=-np.inf,where=(xz<0.0)), np.max(xy,initial=-np.inf,where=(xy<0.0)), np.max(yz,initial=-np.inf,where=(yz<0.0)))
           symmin = min(np.abs(smin),np.abs(smax))
           vmax = np.log10(max(np.abs(smin),np.abs(smax))/symmin)
       vmin = -vmax
       xy = np.sign(xy)*np.log10(np.maximum(np.abs(xy)/symmin,1.0))
       xz = np.sign(xz)*np.log10(np.maximum(np.abs(xz)/symmin,1.0))
       yz = np.sign(yz)*np.log10(np.maximum(np.abs(yz)/symmin,1.0))

    return xy, xz,yz, vmin, vmax

def labelx():
   return lambda var: '$'+str(var)[2:-1].replace('**','^')+'$'
