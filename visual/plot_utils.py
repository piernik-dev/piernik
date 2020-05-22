#!/usr/bin/env python
import numpy as np

def fsym(vmin,vmax):
   vmx =        np.max([np.abs(vmin),np.abs(vmax)])
   vmn = -1.0 * np.max([np.abs(vmin),np.abs(vmax)])
   if vmn == vmx:
      vmn = vmn - 0.00001
      vmx = vmx + 0.00001
   return vmn, vmx

def labelx():
   return lambda var: '$'+str(var)[2:-1].replace('**','^')+'$'
