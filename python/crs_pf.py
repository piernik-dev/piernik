#!/usr/bin/python
from pylab import *
import matplotlib.pyplot as plt
import matplotlib as mp
import numpy as np
from numpy import log, log10
from math import pi
import os
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from colored_io import die, prtinfo, prtwarn

p_ratios_lo  = []
f_ratios_lo  = []
alpha_tab_lo = []
n_tab_lo     = []

p_ratios_up  = []
f_ratios_up  = []
alpha_tab_up = []
n_tab_up     = []

size = 0


def read_dat_table(table_name):

   sFile = open("../src/fluids/cosmicrays/"+table_name+".dat","r")
   sFile.readline(200)
   sFile.readline(22)
   size = int(sFile.readline(3))

   dataArray =[]
   data_to_plot = [[0.0 for i in range(size)]for j in range(size)]
   i = 0
   with open("../src/fluids/cosmicrays/"+table_name+".dat","r") as sFile:
      next(sFile)
      next(sFile)
      next(sFile)
      for line in sFile:
         data = []
         for item in line.split(' '):
            if item != '': data.append(float(item))
            data_to_plot[i][:] = data
         i = i+1
   table = data_to_plot
   return table

def initialize_pf_arrays(pf_initialized = False):
   global p_ratios_lo, f_ratios_lo, p_ratios_up, f_ratios_up, alpha_tab_lo, n_tab_lo, alpha_tab_up, n_tab_up, size
   p_ratios_lo = read_dat_table("p_ratios_lo")
   f_ratios_lo = read_dat_table("f_ratios_lo")

   p_ratios_up = read_dat_table("p_ratios_up")
   f_ratios_up = read_dat_table("f_ratios_lo")

   size = int(len(p_ratios_lo))
# As in cresp_NR_method
   a_min_lo    = 0.2
   a_max_lo    = 0.999999
   a_min_up    = 1.000005
   a_max_up    = 200.0

   n_min_lo    = 1.0e-11
   n_max_lo    = 5000.0
   n_min_up    = 1.0e-12
   n_max_up    = 1000.0

   alpha_tab_lo = zeros(size) ;  n_tab_lo = zeros(size)
   alpha_tab_up = zeros(size) ;  n_tab_up = zeros(size)

   for i in range(size):
      alpha_tab_lo[i] = ind_to_flog(i, a_min_lo, a_max_lo, size)
      alpha_tab_up[i] = ind_to_flog(i, a_min_up, a_max_up, size)
      n_tab_lo[i]     = ind_to_flog(i, n_min_lo, n_max_lo, size)
      n_tab_up[i]     = ind_to_flog(i, n_min_up, n_max_up, size)

   pf_initialized = True

#initialized
   return

def ind_to_flog(ind, min_in, max_in, length):
   ind_to_flog = min_in * 10.0 ** (( (log10(max_in/min_in)) / float( length-1) ) * float(ind))
   return ind_to_flog

def inverse_f_to_ind(value, min_in, max_in, length):
   inverse_f_to_ind = (log10(value/min_in)/log10(max_in/min_in)) * (length)
   try:
      inverse_f_to_ind = int(inverse_f_to_ind)
   except:
      inverse_f_to_ind = -1 # NaN or inf or -inf or whatnot
   return inverse_f_to_ind

def bilin_interpol(y11, y12, y21, y22, t, u):
#  y11, y12, y21, y22, t, u ! y** - tabularized values of interpolated function, t, u - coefficients
   bilin_interpol = (1.0 - t)*(1.0 - u)*y11  +  t*(1.0 - u)*y12   +  (1.0 - t)*u*y21   +   t*u*y22
   return bilin_interpol

def bl_in_tu(val_left, val_mid, val_right):
   bl_in_tu = (val_mid - val_left) / (val_right - val_left)
   return bl_in_tu

def get_interpolated_ratios(bnd, alpha_val, n_val , interpolation_error, **kwargs):
   verbose = kwargs.get("verbose","False")

   loc1 = [1,1]
   loc2 = [2,2]

   pf_ratio = [1.0, 1.0]
   if bnd == "up":
      p_p = p_ratios_up
      p_f = f_ratios_up
      p_a = alpha_tab_up
      p_n = n_tab_up
   elif bnd == "lo":
      p_p = p_ratios_lo
      p_f = f_ratios_lo
      p_a = alpha_tab_lo
      p_n = n_tab_lo
   else:
      die("Provided bnd not supported: %s" %bnd)

   loc1[0] = inverse_f_to_ind(n_val,     p_n[0], p_n[size-1], size-1)
   loc1[1] = inverse_f_to_ind(alpha_val, p_a[0], p_a[size-1], size-1)

   if min(loc1[:]) > 0 and max(loc1[:]) <= size-1  :
      loc2[0] = loc1[0] + 1 ; loc2[1] = loc1[1] + 1
   else:
      if (verbose): prtwarn("(ERROR) Obtained indices (%s): [%i, %i] -> [%f, %f]" %(bnd, loc1[0], loc1[1], pf_ratio[0], pf_ratio[1]))
      interpolation_error = True
      return pf_ratio, interpolation_error

   pf_ratio[0] = bilin_interpol(p_p[loc1[0]][loc1[1]], p_p[loc2[0]][loc1[1]], p_p[loc1[0]][loc2[1]], \
      p_p[loc2[0]][loc2[1]], bl_in_tu(p_a[loc1[1]], alpha_val, p_a[loc2[1]]), bl_in_tu(p_n[loc1[0]], n_val, p_n[loc2[0]]))

   pf_ratio[1] = bilin_interpol(p_f[loc1[0]][loc1[1]], p_f[loc2[0]][loc1[1]], p_f[loc1[0]][loc2[1]], \
      p_f[loc2[0]][loc2[1]], bl_in_tu(p_a[loc1[1]], alpha_val, p_a[loc2[1]]), bl_in_tu(p_n[loc1[0]], n_val, p_n[loc2[0]]))

   if min(pf_ratio[:]) <= 0.0:
      if (verbose): prtwarn("(ERROR) Obtained p & f (%s) ratios: [%f, %f]" %(bnd, pf_ratio[0], pf_ratio[1]))
      return pf_ratio, interpolation_error

   interpolation_error = False

   return pf_ratio, interpolation_error

def e_small_2_f(e_small, p_outer):
   clight = 1.0
   e_small_2_f = e_small / (4.0 * pi * (clight ** 2.0)  * p_outer ** 3.0)
   return e_small_2_f
