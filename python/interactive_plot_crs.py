#!/usr/bin/python
# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
import matplotlib
from matplotlib.colors import LogNorm
import numpy as np
from math import pi
from numpy import log, log10
import os, weakref
import sys
import read_h5
import crs_h5
try:
    import yt
    import h5py
except:
    print "You must make yt & h5py available somehow"
    exit(-1)

f_run = True

plot_vel = False
plot_mag = True
logscale_colors = False
#---------- reading parameters
filename=sys.argv[-1]
exten = filename.split('.')
exten = exten[-1]
if exten[0:2] != 'h5':
    sys.exit("Script requires (a list of) hdf5 file(s) on input")

var_array = []
if f_run == True:
    var_names = []
    var_names = [ "ncre", "p_min_fix", "p_max_fix", "e_small"] # TODO: add cre_eff
    if len(var_names) == 0:
        print ("Empty list of parameter names provided: enter names of parameters to read")
        var_names = read_h5.input_names_array()
    
    var_array = read_h5.read_par(filename, var_names)
    for i in range(len(var_names)):
        exec( "%s=%s" %(var_names[i], var_array[i]))

    print ""
    print "*** Values read from problem.par@hdf5 file: *** "
    for i in range(len(var_names)):
        print ( " %20s =  %8s ( %15s  ) " %(var_names[i], var_array[i], type(var_array[i])))
    print ""
#---------- Open file
    h5File = h5py.File(filename,'r')
    h5_field_list = []
    for grid in h5File['data']:
        for data in h5File['data'][grid]:
            if data not in h5_field_list: h5_field_list.append(data)

    h5ds = yt.load(filename)
#----------- Loading other data
    t = h5ds.current_time[0]
    dt= h5File.attrs["timestep"]
    time = t.in_units('Myr')    

#---------- bounds on domain size
    grid_dim = h5File['grid_dimensions'][0][:]
    dim_map  = {'x':0,'y':1,'z':2}
    avail_dim = [0,1,2]
    avail_dims_by_slice = [[1,2],[0,2],[0,1]]
    dom_l = np.array(h5ds.domain_left_edge[0:3])
    dom_r = np.array(h5ds.domain_right_edge[0:3])

    print "Domain size of provided file (x, y, z): ", grid_dim[:]
    if len(grid_dim) == 3 and min(grid_dim) != 1:
        slice_ax = ''
        slice_coord = -1
        while slice_ax not in dim_map.keys():
            slice_ax    = raw_input("Choose slice ax (x, y, z)      : ")
        while slice_coord > grid_dim[dim_map[slice_ax]] or slice_coord < 0:
            slice_coord = int(raw_input("Choose slice coordinate (%d:%d): " % (0,grid_dim[dim_map[slice_ax]]) ))
    elif min(grid_dim) == 1:
        slice_coord = 0
        if   grid_dim[0] == 1: 
            slice_ax = 'x'
        elif grid_dim[1] == 1: 
            slice_ax = 'y'
        else:                  
            slice_ax = 'z'
    avail_dim = avail_dims_by_slice[dim_map[slice_ax]]
    print "Slice ax set to", slice_ax, ", ind = ", slice_coord #, "available_dims: ", avail_dim
    min_ran = [0,0,0]
    min_ran[dim_map[slice_ax]] = slice_coord
    max_ran = grid_dim
    max_ran[dim_map[slice_ax]] = slice_coord

#---------
    s = plt.figure(figsize=(12,6),dpi=80)
    s1 = plt.subplot(121)

    plot_field = "cr01"
    dset = h5File['data']['grid_0000000000'][plot_field]

    if (slice_ax == "x"):
      fig1 = plt.imshow(dset[:,:,slice_coord], origin="lower") # TODO provide the right coordinates
      field_max = np.max(dset[:,:,slice_coord])
    elif (slice_ax == "y"):
      fig1 = plt.imshow(dset[:,slice_coord,:], origin="lower") # TODO provide the right coordinates
      field_max = np.max(dset[:,slice_coord,:])
    else:
      fig1 = plt.imshow(dset[slice_coord,:,:], origin="lower") # TODO provide the right coordinates
      field_max = np.max(dset[slice_coord,:,:])

    plt.title("Component name: "+plot_field+" | Time = %f Myr"  %time,y=1.07)
    plt.ylabel("Physical domain ("+dim_map.keys()[dim_map.values().index(avail_dim[1])]+") [kpc]" )
    plt.xlabel("Physical domain ("+dim_map.keys()[dim_map.values().index(avail_dim[0])]+") [kpc]" )
    if logscale_colors == True:
            print "WARNING: logscale_colors = True, if negative values are encountered, the script will crash."
            plt.pcolor(dset[:,:,0], norm=LogNorm(vmin=dset[:,:,0].min(), vmax=dset[:,:,0].max()), cmap='viridis') # TODO provide the right cooridinates
    plt.colorbar(shrink=0.9, pad=0.18)

    valuesx = tuple(np.arange(0,grid_dim[avail_dim[0]], grid_dim[avail_dim[0]]/5.))# + (grid_dim[avail_dim[0]],)
    labelsx = tuple(np.arange(dom_l[avail_dim[0]], dom_r[avail_dim[0]], ((abs(dom_l[avail_dim[0]])+abs(dom_r[avail_dim[0]]))/5.)))# + (dom_r[avail_dim[0]],)

    s1.set_xticks(valuesx)
    s1.set_xticklabels(labelsx)

    valuesy = tuple(np.arange(0,grid_dim[avail_dim[1]], grid_dim[avail_dim[1]]/5.))# + (grid_dim[avail_dim[1]],)
    labelsy = tuple(np.arange(dom_l[avail_dim[1]], dom_r[avail_dim[1]], ((abs(dom_l[avail_dim[1]])+abs(dom_r[avail_dim[1]]))/5.)))# + (dom_r[avail_dim[1]],)
    
    s1.set_yticks(valuesy)
    s1.set_yticklabels(labelsy)

    click_coords = [0, 0]
    image_number = 0
#---------
    def read_click_and_plot(event):
        global click_coords, image_number, first_run
        image_number ++1
        #print('click: x=%d, y=%d, xdata=%d, ydata=%d' % (event.x, event.y, int(round(event.xdata)), int(round(event.ydata))))
        click_coords = [ int(round(event.xdata)), int(round(event.ydata)) ]
        point = s1.plot(event.xdata,event.ydata, marker="+", color="red")
        coords = [slice_coord, slice_coord, slice_coord]
        coords[dim_map[slice_ax]] = slice_coord
        if slice_ax == "x":
            coords[1] = click_coords[0]
            coords[0] = click_coords[1]
        elif slice_ax == "y":
            coords[0] = click_coords[1]
            coords[2] = click_coords[0]
        else: # slice_ax = "z"
            coords[2] = click_coords[0]
            coords[1] = click_coords[1]
        print "Coordinates clicked:", coords[::-1]
# ------------ preparing data and passing -------------------------
        time = h5File.attrs['time'][0]
        dt   = h5File.attrs['timestep'][0]
        ecrs = [] ; ncrs = []
        for ind in range(1,ncre+1):
            ecrs.append(h5File['data']['grid_0000000000']['cree'+str(ind).zfill(2)].value[coords[0],coords[1],coords[2]])
            ncrs.append(h5File['data']['grid_0000000000']['cren'+str(ind).zfill(2)].value[coords[0],coords[1],coords[2]])
        plot_var = "e"
        fig2 = crs_h5.crs_plot_main(var_names, var_array, plot_var, ncrs, ecrs, field_max, time, dt, image_number)
        first_run = False
    cid = s.canvas.mpl_connect('button_press_event',read_click_and_plot)
    
    plt.show()
    s.canvas.mpl_disconnect(cid)
