#!/usr/bin/python
from pylab import zeros, sqrt, size
import matplotlib.pyplot as plt
from numpy import log10, log, pi, asfarray, array, linspace
import h5py
import os
import sys

#------ default values of parameters --------
e_small = 1.0e-5
eps     = 1.0e-15
ncre      = 45
p_min_fix = 0.4e0  #
p_max_fix = 1.65e4 #
cre_eff = 0.01
c = 1.0 # 0.3066067E+06  # PSM -> 0.3066067E+06, SI -> 0.2997925E+09

first_run = True
fixed_width = True

def nq2f(n,q,p_l,p_r):
      if p_r> 0.0 and p_l > 0 :
        pr_by_pl = p_r / p_l
        nq2f = n / (4* pi * p_l**3)
        if abs(q-3.0) > eps:
            nq2f = nq2f *(3.0 - q) / (( pr_by_pl )**( 3.0 - q ) - 1.0)
        else:
            nq2f = nq2f / log(p_r/p_l)
      return nq2f
# 1D Newton-Raphson algorithm (to find q):
#
def nr_get_q(alpha, p_ratio):
  q_big = 30.0
  exit_code = True
  iter_limit = 30
  tol_f      = 1.0e-9
  x = 5.0
  df = 1.0
  for i in range(iter_limit):
        if abs(x) >= q_big:
           x = (x/abs(x)) * q_big
           break
        dx = min(max(x*1e-3,1.0e-5),10e-2)
        df = 0.5*(fun(x+dx, alpha, p_ratio)-fun(x-dx, alpha, p_ratio))/dx
        delta = -fun(x, alpha, p_ratio)/df
        if abs(delta) <= tol_f:
            exit_code = False
            return x
        else:
            x = x + delta
  nr_get_q = x
  return nr_get_q

# function used to find q: ----------------------
def fun(x, alpha, p_ratio):
      if abs(x - 3.0) < 1.0e-3:
         fun = -alpha + (-1.0 + p_ratio)/log(p_ratio)
      elif abs(x - 4.0) < 1.0e-3:
         fun = -alpha + p_ratio*log(p_ratio)/(p_ratio - 1.0)
      else:
         fun = -alpha + ((3.0-x)/(4.0-x))*((p_ratio**(4.0-x)-1.0)/(p_ratio**(3.0-x)-1.0))
      return fun

# plot data ------------------------------------
def plot_data(plot_var, pl, pr, fl, fr, q, time, location, i_lo_cut, i_up_cut):
   global first_run
   f_lo_cut = fl[0] ;      f_up_cut = fr[-1]
   p_lo_cut = pl[0] ;   p_up_cut = pr[-1]

   if plot_var  == 'f' :
      plot_var_l      = fl
      plot_var_r      = fr
      plot_var_lo_cut  = f_lo_cut
      plot_var_up_cut = f_up_cut
   elif plot_var == 'n' :
      plot_var_l      = 4*pi*fl*pl**2
      plot_var_r      = 4*pi*fr*pr**2
      plot_var_lo_cut = 4*pi*f_lo_cut*p_lo_cut**2
      plot_var_up_cut = 4*pi*f_up_cut*p_up_cut**2
   elif plot_var == 'e' :
      plot_var_l      = 4*pi*c**2*fl*pl**3
      plot_var_r      = 4*pi*c**2*fr*pr**3
      plot_var_lo_cut = 4*pi*c**2*f_lo_cut*p_lo_cut**3
      plot_var_up_cut = 4*pi*c**2*f_up_cut*p_up_cut**3

   s = plt.subplot(122)
   plt.cla()
   s.set_xscale('log')
   s.set_yscale('log')
   plt.xlabel('p')
   plt.ylabel(plot_var)
   global plot_p_min, plot_p_max, plot_var_min, plot_var_max
   if first_run :
      plot_p_min    =  p_lo_cut
      plot_p_max    =  p_up_cut
      plot_var_min = 0.1*e_small
      first_run = False
      if (plot_var == "e"):
        plot_var_min = 0.1 * e_small
      elif (plot_var == "f" ):
        plot_var_min = e_small / (4*pi * (c ** 2)  * p_max_fix **3) /10.
      elif (plot_var == "n" ):
        plot_var_min = 0.1 / (4*pi * p_max_fix ** 2)

   if fixed_width == True:
      plt.xlim (0.1 * plot_p_min   ,  10.*plot_p_max )
      plt.xlim

   plt.ylim (plot_var_min , plot_ymax)
   plt.grid()

#plot floor value
   p_range = linspace(s.get_xlim()[0],s.get_xlim()[1])
   e_smalls = zeros(len(p_range))
   e_smalls[:] = e_small
   if (plot_var == "e"):
      plt.plot(p_range, e_smalls, color="green", label="$e_{small}$")
   elif(plot_var == "n"):
      plt.plot(p_range, e_small/(c*p_range), color="green",label="$n_{small}$")

   for i in range(0, size(fr)) :
      plt.plot([pl[i],pr[i]],[plot_var_l[i],plot_var_r[i]],color='r')
      plt.plot([pl[i],pl[i]], [plot_var_min,plot_var_l[i]],color='r')
      plt.plot([pr[i],pr[i]], [plot_var_min,plot_var_r[i]],color='r')
   s.set_facecolor('white')
   plt.title("Spectrum of %s(p) \n Time = %7.3f | location: %7.2f %7.2f %7.2f " % (plot_var, time, location[0],location[1],location[2]) )

   return s

def simple_plot_data(plot_var, p, var_array, time, location, i_lo_cut, i_up_cut):
   global first_run, plot_p_min, plot_p_max, fixed_width
   p_lo_cut = p[0] ;   p_up_cut = p[-1]

   s = plt.subplot(122)
   plt.cla()
   s.set_xscale('log')
   s.set_yscale('log')
   plt.xlabel('p')
   plt.ylabel(plot_var)

   if first_run :
      plot_p_min    =  p_lo_cut
      plot_p_max    =  p_up_cut
      plot_var_min = 0.1*e_small
      first_run = False
      if (plot_var == "e"):
        plot_var_min = 0.1 * e_small
      elif (plot_var == "f" ):
        plot_var_min = e_small / (4*pi * (c ** 2)  * p_max_fix **3) /10.
      elif (plot_var == "n" ):
        plot_var_min = 0.1 / (4*pi * p_max_fix ** 2)

   if fixed_width == True:
      plt.xlim (0.25 * plot_p_min   ,  5.*plot_p_max )
      plt.xlim

   plt.ylim (plot_var_min , plot_ymax)
   plt.grid()

   plt.plot(p,var_array,color='r')
   for i in range(i_lo_cut,i_up_cut):
       plt.plot([p[i],p[i]],[plot_var_min, var_array[i]],color="r")

#plot floor values
   p_range = linspace(s.get_xlim()[0],s.get_xlim()[1])
   e_smalls = zeros(len(p_range))
   e_smalls[:] = e_small
   if (plot_var == "e"):
      plt.plot(p_range, e_smalls, color="green", label="$e_{small}$")
   elif(plot_var == "n"):
      plt.plot(p_range, e_small/(c*p_range), color="green",label="$n_{small}$")

   s.set_facecolor('white')
   plt.title(" %s(p) \n Time = %7.3f | location: %7.2f %7.2f %7.2f " % (plot_var, time, location[0],location[1],location[2]) )

   return s
#-----------------------------------------------------------------

def crs_plot_main(parameter_names, parameter_values, plot_var, ncrs, ecrs, field_max, time, location, use_simple):
    global first_run

    try:
        for i in range(len(parameter_names)):
            exec("%s = %s" %(parameter_names[i], parameter_values[i]))
    except:
        print "Exiting: len(names) not equal len(values)"
        sys.exit()

# TODO -------- do it under *args TODO
    fixed_width = True
# TODO -------- do it under *args TODO

    first_run = True
# -------------------
    global plot_ymax
    plot_ymax = field_max * cre_eff
    edges = []
    p_fix = []
    edges[0:ncre] = range(0,ncre+1, 1)
    p_fix[0:ncre] = zeros(ncre+1)

# WARNING !!! cutoff momenta are not precisely described here !!! WARNING
    log_width   = (log10(p_max_fix/p_min_fix))/(ncre-2.0)
# TODO: do it in the first run / calling script TODO
    if first_run:
        for i in range(0,ncre-1):
            p_fix[i+1]  = p_min_fix * 10.0**( log_width * edges[i])
            p_fix_ratio = 10.0 ** log_width
            p_fix[0]    = ( sqrt(p_fix[1] * p_fix[2]) ) / p_fix_ratio
            p_fix[ncre]    = ( sqrt(p_fix[ncre-2]* p_fix[ncre-1]) ) * p_fix_ratio
            p_fix = asfarray(p_fix)

    p_mid_fix = zeros(ncre)
    p_mid_fix[1:ncre-1] = sqrt(p_fix[1:ncre-1]*p_fix[2:ncre])
    p_mid_fix[0]        = p_mid_fix[1] / p_fix_ratio
    p_mid_fix[ncre-1]   = p_mid_fix[ncre-2] * p_fix_ratio
    p_mid_fix = asfarray(p_mid_fix)
    i_lo = 0
    i_up = ncre
    empty_cell = True

#------------ locate cutoff indices
    for i in range(1,ncre):
        i_lo = i
        if ecrs[i] >= e_small:
            empty_cell = False # not plotting if empty
            break
    for i in range(ncre,1,-1):
        i_up = i
        if ecrs[i-1] >= e_small: break

    sys.stdout.write('Time = %6.2f |  i_lo = %2d, i_up = %2d, %11s.'%(time, i_lo if not empty_cell else 0, i_up if not empty_cell else 0, '(empty cell)' if empty_cell else ' '))
    pln = p_fix[i_lo:i_up]
    prn = p_fix[i_lo+1:i_up+1]
    pln = array(pln)
    prn = array(prn)
    ncrs_act = ncrs[i_lo-2:i_up-2]
    ecrs_act = ecrs[i_lo-2:i_up-2]
    q_nr = [] ; fln = [] ; frn = []
    for i in range(0,i_up - i_lo):
        q_nr.append(nr_get_q(ecrs[i+i_lo]/(ncrs[i+i_lo]*c*pln[i]), prn[i]/pln[i], ))
        fln.append(nq2f(ncrs[i+i_lo], q_nr[-1], pln[i], prn[i]))

    q_nr = array(q_nr)
    fln  = array(fln)
    frn  = array(fln)
    frn  = frn * (prn/pln) ** (-q_nr)
    plot = False
    if empty_cell==False:
         if (use_simple):
            var_array = zeros(ncre)
            if (plot_var == "e"):
               var_array = array(ecrs)
            elif (plot_var == "n"):
               var_array = array(ncrs)
            else: # plotvar == "f"
               var_array = array(fln) # TODO not fully prepared
            plot = simple_plot_data(plot_var, p_mid_fix, var_array, time, location, i_lo, i_up)
         else:
            plot = plot_data(plot_var, pln, prn, fln, frn, q_nr, time, location, i_lo, i_up)

    return plot, empty_cell
