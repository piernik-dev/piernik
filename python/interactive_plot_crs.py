#!/usr/bin/python
# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import numpy as np
from math import pi
from numpy import log, log10
import crs_pf
import os, weakref
import sys
import read_h5
import crs_h5
import re
from colored_io import prtinfo, prtwarn, read_var, die
try:
    import yt
    import h5py
    from yt.units import dimensions
except:
    die("You must make yt & h5py available somehow")

plot_field = "cree05"
plot_var = "e"

simple_plot = False # True
plot_vel = False
plot_mag = True

user_limits = True ### True# False # use values below as limits for clickable plot
plot_user_min = 1.0e-3 #1.0e-3 #1.0e-5
plot_user_max = 0.47#1.78#0.15#0.37#0.47#0.034 #0.86#9.93 # 1.6e-3 #9.93
use_logscale  = True
hdf_save_fpq  = False

user_annot_line = False
user_annot_time = False
user_save_spectrum = True
par_epsilon   = 1.0e-15
#####################################

f_run = True
pf_initialized = False
def _total_cree(field,data):
   list_cree = []
   for element in h5ds.field_list:
      if re.search("cree",str(element[1])):
         list_cree.append(element[1])
   cree_tot = data[str(list_cree[0])]
   for element in list_cree[1:]:
      cree_tot = cree_tot + data[element]
   return cree_tot

def _total_cren(field,data):
   list_cren = []
   for element in h5ds.field_list:
      if re.search("cren",str(element[1])):
         list_cren.append(element[1])
   cren_tot = data[str(list_cren[0])]
   for element in list_cren[1:]:
      cren_tot = cren_tot + data[element]
   return cren_tot

def en_ratio(field,data):
   bin_nr = field.name[1][-2:]
   for element in h5ds.field_list:
      if re.search("cree"+str(bin_nr.zfill(2)),str(element[1])):
         cren_data = data["cren"+str(bin_nr.zfill(2))]
         cren_data[cren_data == 0.0 ] = par_epsilon # necessary to avoid FPEs
         cree_data = data["cree"+str(bin_nr.zfill(2))]
         cree_data[cree_data == 0.0 ] = par_epsilon # necessary to avoid FPEs
         en_ratio=cree_data / cren_data
   return en_ratio

def copy_field(field,data):
   field_name_to_copy= field.name[1][:].split("_")[0]
   copied_field = data[field_name_to_copy]
   return copied_field


#---------- reading parameters

filename=sys.argv[-1]
filename_ext= filename.split('.')[-1]
filename_nam  = filename.split('.')[0].split('/')[-1]
if filename_ext[0:2] != 'h5':
    die("Script requires a (list of) hdf5 file(s) on input")

if f_run :
    if not os.path.exists('results'):
        os.makedirs('results')
        prtinfo("Output directory created: %s" %(os.getcwd()+'/results'))

var_array = []
if f_run == True:
    var_names = []
    var_names = [ "ncre", "p_min_fix", "p_max_fix", "e_small", "cre_eff", "q_big", "hdf_save_fpq"]
    if len(var_names) == 0:
        prtwarn("Empty list of parameter names provided: enter names of parameters to read")
        var_names = read_h5.input_names_array()

    var_array = read_h5.read_par(filename, var_names)
    for i in range(len(var_names)):
        exec( "%s=%s" %(var_names[i], var_array[i]))

    if type(hdf_save_fpq) != type(False): hdf_save_fpq = False # if parameter not in problem.par - makes sure it's type is bool and set to False

    print ""
    prtinfo("*** Values read from problem.par@hdf5 file: *** ")
    for i in range(len(var_names)):
        prtinfo ( " %20s =  %8s ( %15s  ) " %(var_names[i], var_array[i], type(var_array[i])))
    print ""

    crs_pf.initialize_pf_arrays(pf_initialized)
#---------- Open file
    h5File = h5py.File(filename,"r") # sorry, I'm not sure how to access timestep via yt
    h5ds = yt.load(filename)
#---------- bounds on domain size
    grid_dim = h5ds.domain_dimensions
    dim_map  = {'x':0,'y':1,'z':2}
    dom_l = np.array(h5ds.domain_left_edge[0:3])
    dom_r = np.array(h5ds.domain_right_edge[0:3])

#----------- Loading other data
    t = h5ds.current_time[0]
    dt= h5File.attrs["timestep"]
    time = t.in_units('Myr')

#------------ Organizing domain data
    length_unit = 'pc'
    prtinfo("Domain shape of in provided file            (i, j, k): [%i, %i, %i] \033[0m" %(grid_dim[0], grid_dim[1], grid_dim[2]) )
    prtinfo("Domain physical dimensions in provided file (x, y, z): [%9.3f,%9.3f,%9.3f]:[%9.3f,%9.3f,%9.3f] %s \033[0m" %(dom_l[0], dom_l[1], dom_l[2], dom_r[0], dom_r[1], dom_r[2], length_unit))

    avail_dim = [0,1,2]
    avail_dims_by_slice = [[1,2],[0,2],[0,1]]

    if len(grid_dim) == 3 and min(grid_dim) != 1:
        slice_ax = ""
        slice_coord = ""
        while slice_ax not in dim_map.keys():
            slice_ax    = read_var("Choose slice ax (x, y, z)      : ")
        while (slice_coord < dom_l[dim_map[slice_ax]]) or (slice_coord > dom_r[dim_map[slice_ax]]): # or slice_coord < -10000:
            try:
                slice_coord = float(read_var("Choose slice coordinate (%d:%d %s ) (if empty, middle is assumed): \033[0m" % (dom_l[dim_map[slice_ax]],dom_r[dim_map[slice_ax]], length_unit) ))
            except:
                slice_coord = (dom_l[dim_map[slice_ax]] + dom_r[dim_map[slice_ax]]) / 2.
                prtwarn(" (empty / incorrect input): Setting slice coordinate to %s %s.\033[0m" %(slice_coord, length_unit))
    elif min(grid_dim) == 1:
        slice_coord = 0.0
        if   grid_dim[0] == 1:
            slice_ax = 'x'
        elif grid_dim[1] == 1:
            slice_ax = 'y'
        else:
            slice_ax = 'z'
    avail_dim = avail_dims_by_slice[dim_map[slice_ax]]
    prtinfo("Slice ax set to %s, coordinate = %f %s \033[0m" %(slice_ax, slice_coord, length_unit))
    resolution = [grid_dim[avail_dim[0]],grid_dim[avail_dim[1]]]

#--------- Preparing clickable image
    s = plt.figure(figsize=(12,8),dpi=100)
    s1 = plt.subplot(121)
    dsSlice = h5ds.slice(slice_ax, slice_coord)

    if (plot_field == "cree_tot" ):
      if str(dsSlice["cree01"].units) is "dimensionless":      #  DEPRECATED
         try:
            h5ds.add_field(("gdf","cree_tot"), units="",function=_total_cree, display_name="Total cr electron energy density", sampling_type="cell")
         except:
            die("Failed to construct field %s" % plot_field)
      else:
         try:
            h5ds.add_field(("gdf","cree_tot"), units="Msun/(Myr**2*pc)",function=_total_cree, display_name="Total cr electron energy density", dimensions=dimensions.energy/dimensions.volume, sampling_type="cell")
         except:
            die("Failed to construct field %s" % plot_field)

    if (plot_field == "cren_tot" ):
      if str(dsSlice["cren01"].units) is "dimensionless":      #  DEPRECATED
         try:
            h5ds.add_field(("gdf","cren_tot"), units="",function=_total_cren, display_name="Total cr electron number density", sampling_type="cell")
         except:
            die("Failed to construct field %s" % plot_field)

      else: # dimensions defined
         try:
            h5ds.add_field(("gdf","cren_tot"), units="1/(pc**3)",function=_total_cren, display_name="Total cr electron number density", dimensions=dimensions.energy/dimensions.volume, sampling_type="cell")
         except:
            die("Failed to construct field %s" % plot_field)

    if ( plot_field[0:-2] == "en_ratio"):
      if str(dsSlice["cren01"].units) is "dimensionless":      #  DEPRECATED
         try:
            h5ds.add_field(("gdf",plot_field), units="", function=en_ratio, display_name="Ratio e/n in %i-th bin" %int(plot_field[-2:]), sampling_type="cell")
         except:
            die("Failed to construct field %s" % plot_field)
      else:
         try:
            h5ds.add_field(("gdf",plot_field), units="Msun*pc**2/Myr**2", function=en_ratio, display_name="Ratio e/n in %i-th bin" %int(plot_field[-2:]),dimensions=dimensions.energy, sampling_type="cell")
         except:
            die("Failed to construct field %s" % plot_field)

    if (plot_field[-3:] != "tot" and plot_field[0:3] == "cre"):
       if (plot_field[0:4] == "cree"):
          disp_name="CR electron energy density (bin %2i)" %int(plot_field[-2:])
          new_field_dimensions = dimensions.energy/dimensions.volume
          new_field_units      = dsSlice[plot_field].units
       elif (plot_field[0:4] == "cren"):
          disp_name="CR electron number density (bin %2i)" %int(plot_field[-2:])
          new_field_dimensions = 1./dimensions.volume
          new_field_units      = dsSlice[plot_field].units
       prtinfo( "Adding display name: %s" %disp_name)
       plot_field_new = str(plot_field+"_updated")
       h5ds.add_field(("gdf",plot_field_new), units=new_field_units, function=copy_field, display_name=disp_name,dimensions=new_field_dimensions, sampling_type="cell")
       plot_field = plot_field_new

    click_coords = [0, 0]
    image_number = 0

    field_max  = float(h5ds.find_max("cr01")[0]) # WARNING - this makes field_max unitless

    w = dom_r[avail_dim[0]] + abs(dom_l[avail_dim[0]])
    h = dom_r[avail_dim[1]] + abs(dom_l[avail_dim[1]])
    if (plot_field[0:-2] != "en_ratio"):
      frb = np.array(dsSlice.to_frb(w, resolution, height=h)[plot_field])
      plot_max   = h5ds.find_max(plot_field)[0]
      plot_units = str(plot_max).split(' ')[1]
      plot_min   = h5ds.find_min(plot_field)[0]
    else:
      frb = np.array(dsSlice.to_frb(w, resolution, height=h)[plot_field])
      plot_min = h5ds.find_min(plot_field)[0]
      plot_max = h5ds.find_max(plot_field)[0]
      h5ds.find_max("cre"+plot_var+str(plot_field[-2:]))[0]
      plot_units = "Msun*pc**2/Myr**2"

    plt.xlabel("Domain cooridnates "+dim_map.keys()[dim_map.values().index(avail_dim[0])]+" ("+length_unit+")" )
    plt.ylabel("Domain cooridnates "+dim_map.keys()[dim_map.values().index(avail_dim[1])]+" ("+length_unit+")" )
    plt.colormap="plasma"
    if (user_limits == True): # Overwrites previously found values
       plot_min = plot_user_min
       plot_max = plot_user_max

    plot_max   = float(plot_max)
    plot_min   = max(float(plot_min), plot_user_min)

    if (use_logscale):
        plt.imshow(frb,extent=[dom_l[avail_dim[0]], dom_r[avail_dim[0]], dom_l[avail_dim[1]], dom_r[avail_dim[1]] ], origin="lower" ,norm=LogNorm(vmin=plot_min, vmax=plot_max))
    else:
        plt.imshow(frb,extent=[dom_l[avail_dim[0]], dom_r[avail_dim[0]], dom_l[avail_dim[1]], dom_r[avail_dim[1]] ], origin="lower")
    plt.title("Component: "+plot_field+" | t = %9.3f Myr"  %time)

    try:
       cbar = plt.colorbar(shrink=0.9, pad=0.01,label=plot_units)
    except:
       die("An empty field might have been picked.")

    yt_data_plot = yt.SlicePlot(h5ds, slice_ax, plot_field)

    colormap_my  = plt.cm.viridis
    colormap_my.set_bad(color=colormap_my(par_epsilon)) # masks bad values
    yt_data_plot.set_cmap(field=plot_field, cmap = colormap_my)

    if (user_annot_line == True):
       prtinfo("Marking line on yt.plot at (0 0 0) : (500 500 0)")
       yt_data_plot.annotate_line((0., 0., 0.), (500., 500.0, 0), plot_args={'color':'white',"lw":2.0} )#, coord_system='axis')
    else:
       prtinfo("Not marking line on yt.plot (user_annot_line = %s)" %(user_annot_line))


    if (plot_vel): yt_data_plot.annotate_velocity(factor=48)
    if (plot_mag): yt_data_plot.annotate_magnetic_field(factor=48)
    yt_data_plot.set_zlim(plot_field,plot_min,plot_max)
    marker_l   = ["x", "+", "*", "X", ".","^", "v","<",">","1"]
    m_size_l   = [500, 500, 400, 400, 500, 350, 350, 350, 350, 500]
    m_e_width  = 5
    marker_index = 0

    plt.subplots_adjust(left=0.075,right=0.975, hspace=0.12)

    print ("")
    prtinfo("\033[44mClick on the colormap to display spectrum")
#---------
    def read_click_and_plot(event):
        global click_coords, image_number, f_run, marker_index
        exit_code = True
        if (marker_index == len(marker_l) or marker_index == len(m_size_l)): marker_index = 0
        click_coords = [ event.xdata, event.ydata ]
        coords = [slice_coord, slice_coord, slice_coord]
        if slice_ax == "x":
            coords[1] = click_coords[0]
            coords[2] = click_coords[1]
        elif slice_ax == "y":
            coords[0] = click_coords[0]
            coords[2] = click_coords[1]
        else: # slice_ax = "z"
            coords[0] = click_coords[0]
            coords[1] = click_coords[1]

# ------------ preparing data and passing -------------------------
        position = h5ds.r[coords:coords]
        if ( plot_field[0:-2] != "en_ratio"):
           prtinfo ("Value of %s at point [%f, %f, %f] = %f " %(plot_field, coords[0], coords[1], coords[2], position[plot_field]))
        else:
           prtinfo("Value of %s at point [%f, %f, %f] = %f " %(plot_field, coords[0], coords[1], coords[2], position["cree"+str(plot_field[-2:])]/position["cren"+str(plot_field[-2:])]))
           plot_max   = h5ds.find_max("cre"+plot_var+str(plot_field[-2:]))[0] # once again appended - needed as ylimit for the plot
        if (hdf_save_fpq != True):
           ecrs = [] ; ncrs = []
           for ind in range(1,ncre+1):
               ecrs.append(float(position['cree'+str(ind).zfill(2)][0].v))
               ncrs.append(float(position['cren'+str(ind).zfill(2)][0].v))

           fig2,exit_code = crs_h5.crs_plot_main(var_names, var_array, plot_var, ncrs, ecrs, field_max, time, coords, simple_plot, marker=marker_l[marker_index])
        else:
           fcrs = [] ; qcrs = [] ; pcut = [ 0., 0.]
           ecrs = [] ; ncrs = []

           for ind in range(1,ncre+1):
               ecrs.append(float(position['cree'+str(ind).zfill(2)][0].v))
               ncrs.append(float(position['cren'+str(ind).zfill(2)][0].v))

           for ind in range(1, ncre+2):
               fcrs.append(float(position['cref'+str(ind).zfill(2)][0].v))
           for ind in range(1, ncre+1):
               qcrs.append(float(position['creq'+str(ind).zfill(2)][0].v))
           pcut[:] = [ position['crep01'][0].v , position['crep02'][0].v ]

           fig2,exit_code = crs_h5.crs_plot_main_fpq(var_names, var_array, plot_var, fcrs, qcrs, pcut, field_max, time, coords, marker=marker_l[marker_index])

        if (exit_code != True):
            point = s1.plot(event.xdata,event.ydata, marker=marker_l[marker_index], color="red")         # plot point only if cell is not empty
            s.savefig('results/'+filename_nam+'_'+plot_var+'_%04d.png' % image_number, transparent ='True')
            prtinfo("  --->  Saved plot to: %s " %str('results/'+filename_nam+'_'+plot_var+'_%04d.png' %image_number))

            yt_data_plot.annotate_marker(coords, marker=marker_l[marker_index],  plot_args={'color':'red','s':m_size_l[marker_index],"lw":4.5}) # cumulatively annotate all clicked coordinates
            marker_index = marker_index + 1

            image_number=image_number+1
#------------- saving just the spectrum
            if (user_save_spectrum):
               extent = fig2.get_window_extent().transformed(s.dpi_scale_trans.inverted())
               s.savefig('results/'+filename_nam+'_'+plot_var+'_spectrum_%04d.pdf' % image_number, transparent ='True', bbox_inches=extent.expanded(1.4, 1.195),quality=95,dpi=150)
               prtinfo("  --->  Saved plot to: %s.\n\033[44mPress 'q' to quit and save yt.SlicePlot with marked coordinates." %str('results/'+filename_nam+'_'+plot_var+'_spectrum_%04d.pdf' %image_number))
        else:
            prtwarn("Empty cell - not saving.")

        if (f_run): f_run = False
    cid = s.canvas.mpl_connect('button_press_event',read_click_and_plot)

    plt.show()
    yt_data_plot.set_font({'size':30})

    text_coords = [0., 0., 0.]; text_coords[dim_map.get(slice_ax)] = slice_coord; text_coords[avail_dim[0]] = dom_l[avail_dim[0]]; text_coords[avail_dim[1]] = dom_l[avail_dim[1]]
    text_coords = [ item * 0.9 for item in text_coords]

    if (user_annot_time == True):
      yt_data_plot.annotate_text(text_coords , 'T = {:0.2f} Myr'.format(float(t.in_units('Myr'))), text_args={'fontsize':30,'color':'white','alpha':'0.0'},inset_box_args={'boxstyle':'round','pad':0.2,'alpha':0.8})
    else:
      prtinfo("Not marking line on yt.plot (user_annot_time = %s)" %(user_annot_time))

    yt_data_plot.save('results/'+filename_nam+'_'+plot_field+'_sliceplot.pdf')  # save image when "q" pressed
    s.canvas.mpl_disconnect(cid)
