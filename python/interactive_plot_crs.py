#!/usr/bin/python
# -*- coding: utf-8 -*-
from colored_io import die, prtinfo, prtwarn, read_var
import crs_h5
import read_h5
from crs_pf import initialize_pf_arrays
from math import isnan, pi
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from numpy import array as np_array, log, log10, mean, rot90
import os
from optparse import OptionParser
from re import search
from sys import argv, version
import warnings
try:
    import yt
    from yt.units import dimensions
except:
    die("You must make yt available somehow")
if (version[0:3] != "2.7"):
    prtwarn("Using python version higher than 2.7!")
    raw_input = input
    not_py27 = True
else:
    not_py27 = False
# ------- Parse arguments
parser = OptionParser("Usage: %prog FILE [options] [args] or %prog [options] [args] -F FILENAME")
parser.add_option("-F", "--file", dest="filename", default="None", help=u"File to use", type="str")
parser.add_option("-v", "--var", dest="var_name", default="e", help=u"Variable to plot the spectrum (default: e)")
parser.add_option("-f", "--field", dest="fieldname", default="cree_tot", help=u"DS fieldname to image (default:cree_tot)")
parser.add_option("-z", "--zlim", dest="plot_range", default=("x", "y"), help=u"Plot image with this range", nargs=2, metavar="ZMIN ZMAX")
parser.add_option("-s", "--slice", dest="slice_info", default=("a", "Nan"), help=u"DS slice coords to image (default:cree_tot)", metavar="AX COORDINATE", nargs=2)
parser.add_option("-d", "--def", dest="default_range", default=False, help=u"Use min/max on yt.ds(fieldname) for clickable image", action="store_true")
parser.add_option("-l", "--lin", dest="use_linscale", default=False, help=u"Use linear scale for clickable plot (default: log)", action="store_true")
parser.add_option("-V", "--vel", dest="plot_vel", default=False, help=u"Plot velocity field vectors ", action="store_true")
parser.add_option("-m", "--mag", dest="plot_mag", default=False, help=u"Plot magnetic field vectors ", action="store_true")
parser.add_option("-t", "--time", dest="annotate_time", default=False, help=u"Annotate time on resulting DS plot", action="store_true")
parser.add_option("", "--nosave", dest="not_save_spec", default=False, help=u"Do not save output spectrum ", action="store_true")
parser.add_option("-a", "--average", "--avg", dest="avg_layer", default=False, help=u"Plot mean spectrum at pointed layer at ordinate axis", action="store_true")
parser.add_option("-O", "--overlap", dest="overlap_layer", default=False, help=u"Overlap all spectra at pointed layer at ordinate axis", action="store_true")
parser.add_option("", "--verbose", dest="yt_verbose", default=40, help=u"Append yt verbosity (value from 10 [high] to 50 [low])")
parser.add_option("-q", "--quiet", dest="py_quiet", default=False, help=u"Suppress ALL python warnings (not advised)", action="store_true")
parser.add_option("-c", "--coords", dest="coords_dflt", default=("x", "y", "z"), help=u"Provides coordinates for the first plot (non-clickable)", nargs=3, metavar="xc yc zc")
parser.add_option("-k", "--keep", dest="clean_plot", default=True, help=u"Keep spectrum plot (do not clean spectrum field)", action="store_false")
parser.add_option("", "--fontsize", dest="fontsize", default=18, help=u"Set fontsize for SlicePlot (default is 18)")
parser.add_option("", "--noxlabels", dest="no_xlabels", default=False, help=u"Do not show labels of X axis on resulting SlicePlot", action="store_true")
parser.add_option("", "--noylabels", dest="no_ylabels", default=False, help=u"Do not show labels of Y axis on resulting SlicePlot", action="store_true")
parser.add_option("", "--nocbar", dest="no_cbar", default=False, help=u"Do not show colorbar on the resulting SlicePlot", action="store_true")
parser.add_option("", "--noaxes", dest="no_axes", default=False, help=u"Hide axes on the spectrum plot (useful when combining spectra)", action="store_true")
parser.add_option("", "--rectangle", dest="annotate_rect", default=False, help=u"Annotate and average over width/height rectangle surface", action="store_true")
parser.add_option("", "--width", dest="usr_width", default=0., help=u"Set custom frb width")
parser.add_option("", "--height", dest="usr_height", default=0., help=u"Set custom frb width")
parser.add_option("", "--center", dest="usr_center", default=(0., 0.), help=u"Set custom frb center", nargs=2, metavar="XC YC")

(options, args) = parser.parse_args(argv[1:])  # argv[1] is filename
yt.mylog.setLevel(int(options.yt_verbose))    # Reduces the output to desired level, 50 - least output
if (options.py_quiet == True):
    warnings.simplefilter(action='ignore', category=FutureWarning)
    warnings.simplefilter(action='ignore', category=Warning)
    # warnings.filterwarnings("ignore")
    prtwarn("Python warnings are turned off from here (-q, --quiet switch)")
plot_var = options.var_name
user_draw_timebox = options.annotate_time
user_limits = (options.default_range != True)
save_spectrum = (options.not_save_spec != True)
use_logscale = (options.use_linscale != True)
plot_field = options.fieldname
plot_var = options.var_name
plot_vel = options.plot_vel
plot_mag = options.plot_mag
if (plot_vel == True):
    plot_mag = False
if (plot_mag == True):
    plot_vel = False
user_annot_line = False
user_annot_time = True
plot_layer = options.avg_layer
plot_ovlp = options.overlap_layer
options.fontsize = int(options.fontsize)
display_bin_no = False

user_coords_provided = ((options.coords_dflt[0] != "x") and (options.coords_dflt[1] != "y") and (options.coords_dflt[2] != "z"))
if user_coords_provided:
    try:
        user_coords = [float(options.coords_dflt[0]), float(options.coords_dflt[1]), float(options.coords_dflt[2])]
    except:
        die("Got spectrum coordinates %s, but failed to convert it to float." % str(options.coords_dflt))

#####################################
# ------- Parameters ----------------
par_epsilon = 1.0e-15
f_run = True
pf_initialized = False

# ------- Local functions -----------


def _total_cree(field, data):
    list_cree = []
    for element in h5ds.field_list:
        if search("cree", str(element[1])):
            list_cree.append(element[1])
    cree_tot = data[str(list_cree[0])]
    for element in list_cree[1:]:
        cree_tot = cree_tot + data[element]
    return cree_tot


def _total_cren(field, data):
    list_cren = []
    for element in h5ds.field_list:
        if search("cren", str(element[1])):
            list_cren.append(element[1])
    cren_tot = data[str(list_cren[0])]
    for element in list_cren[1:]:
        cren_tot = cren_tot + data[element]
    return cren_tot


def _total_B(field, data):
    b_tot = 2.85 * (data["mag_field_x"]**2 + data["mag_field_y"]**2 + data["mag_field_y"]**2)**0.5
    return b_tot


def en_ratio(field, data):  # DEPRECATED (?)
    bin_nr = field.name[1][-2:]
    for element in h5ds.field_list:
        if search("cree" + str(bin_nr.zfill(2)), str(element[1])):
            cren_data = data["cren" + str(bin_nr.zfill(2))]
            cren_data[cren_data <= par_epsilon**2] = par_epsilon  # necessary to avoid FPEs
            cree_data = data["cree" + str(bin_nr.zfill(2))]
            en_ratio = cree_data / cren_data
    return en_ratio


def copy_field(field, data):
    field_name_to_copy = field.name[1][:].split("_")[0]
    copied_field = data[field_name_to_copy]
    return copied_field


def add_cren_tot_to(h5_dataset):
    try:
        if (h5ds.all_data()["cren01"].units is "dimensionless"):
            h5ds.add_field(("gdf", "cren_tot"), units="", function=_total_cren, display_name="Total CR electron number density", sampling_type="cell")
        else:
            h5ds.add_field(("gdf", "cren_tot"), units="1/(pc**3)", function=_total_cren, display_name="Total CR electron number density", dimensions=dimensions.energy / dimensions.volume, sampling_type="cell", take_log=True)
    except:
        die("Failed to construct field 'cren_tot'")
    return h5_dataset


def add_cree_tot_to(h5_dataset):
    try:
        if (h5ds.all_data()["cree01"].units is "dimensionless"):
            h5ds.add_field(("gdf", "cree_tot"), units="", function=_total_cree, display_name="Total CR electron energy density", sampling_type="cell")
        else:
            h5ds.add_field(("gdf", "cree_tot"), units="Msun/(Myr**2*pc)", function=_total_cree, display_name="Total CR electron energy density", dimensions=dimensions.energy / dimensions.volume, sampling_type="cell")
    except:
        die("Failed to construct field 'cree_tot'")
    return h5_dataset


def add_tot_fields(h5_dataset):
    h5_dataset = add_cree_tot_to(h5_dataset)
    h5_dataset = add_cren_tot_to(h5_dataset)
    return h5_dataset


# ---------- reading parameters
if (options.filename != "None"):
    filename = options.filename
else:
    filename = argv[1]

filename_trimmed = filename.split("/")[-1]
filename_ext = filename_trimmed.split('.')[-1]
filename_nam = filename_trimmed.split('.')[0].split('/')[-1]
if filename_ext[0:2] != 'h5':
    die("Script requires a (list of) hdf5 file(s) on input")

if f_run:
    if not os.path.exists('results'):
        os.makedirs('results')
        prtinfo("Output directory created: %s" % (os.getcwd() + '/results'))

var_array = []
if f_run == True:
    var_names = []
    var_names = ["ncre", "p_min_fix", "p_max_fix", "e_small", "cre_eff", "q_big"]
    var_def = [20, 10., 1.e5, 1.e-6, 0.01, 30., ]
    if len(var_names) == 0:
        prtwarn("Empty list of parameter names provided: enter names of parameters to read")
        var_names = read_h5.input_names_array()

    var_array = read_h5.read_par(filename, var_names, var_def)
    for i in range(len(var_names)):
        exec("%s=%s" % (var_names[i], var_array[i]))

    prtinfo("\n*** Values read from problem.par@hdf5 file: *** \n")
    for i in range(len(var_names)):
        prtinfo(" %15s =  %10s ( %15s  ) " % (var_names[i], var_array[i], type(var_array[i])))

    initialize_pf_arrays(pf_initialized)
# ---------- Open file
    h5ds = yt.load(filename)
# ---------- bounds on domain size
    grid_dim = h5ds.domain_dimensions
    dim_map = {'x': 0, 'y': 1, 'z': 2}
    dom_l = np_array(h5ds.domain_left_edge[0:3])
    dom_r = np_array(h5ds.domain_right_edge[0:3])

# ----------- Loading other data
    t = h5ds.current_time[0]
    time = t.in_units('Myr')
# ----------- Checking user image limits
    try:
        plot_user_min = float(options.plot_range[0])
        plot_user_max = float(options.plot_range[1])
    except:
        prtwarn("No provided ZLIM or failed to convert it into float -- using default (min/max of ds(field) )")
        user_limits = False
# ------------ Organizing domain data
    length_unit = 'pc'
    prtinfo("Domain shape of in provided file            (i, j, k): [%i, %i, %i] \033[0m" % (grid_dim[0], grid_dim[1], grid_dim[2]))
    prtinfo("Domain physical dimensions in provided file (x, y, z): [%9.3f,%9.3f,%9.3f]:[%9.3f,%9.3f,%9.3f] %s \033[0m" % (dom_l[0], dom_l[1], dom_l[2], dom_r[0], dom_r[1], dom_r[2], length_unit))

    avail_dim = [0, 1, 2]
    avail_dims_by_slice = [[1, 2], [0, 2], [0, 1]]

    if len(grid_dim) == 3 and min(grid_dim) != 1:
        slice_ax = str(options.slice_info[0])
        slice_coord = float(options.slice_info[1])
        while slice_ax not in dim_map.keys():
            slice_ax = read_var("Choose slice ax (x, y, z)      : ")
        while (slice_coord < dom_l[dim_map[slice_ax]]) or (slice_coord > dom_r[dim_map[slice_ax]] or isnan(slice_coord) == True):  # or slice_coord < -10000:
            try:
                slice_coord = float(read_var("Choose slice coordinate (%d:%d %s ) (if empty, middle is assumed): \033[0m" % (dom_l[dim_map[slice_ax]], dom_r[dim_map[slice_ax]], length_unit)))
            except:
                slice_coord = (dom_l[dim_map[slice_ax]] + dom_r[dim_map[slice_ax]]) / 2.
                prtwarn(" (empty / incorrect input): Setting slice coordinate to %s %s.\033[0m" % (slice_coord, length_unit))
    elif min(grid_dim) == 1:
        slice_coord = 0.0
        if grid_dim[0] == 1:
            slice_ax = 'x'
        elif grid_dim[1] == 1:
            slice_ax = 'y'
        else:
            slice_ax = 'z'
    avail_dim = avail_dims_by_slice[dim_map[slice_ax]]
    prtinfo("Slice ax set to %s, coordinate = %f %s \033[0m" % (slice_ax, slice_coord, length_unit))
    resolution = [grid_dim[avail_dim[0]], grid_dim[avail_dim[1]]]

# --------- Preparing clickable image
    s = plt.figure(figsize=(12, 8), dpi=100)
    s1 = plt.subplot(121)
    dsSlice = h5ds.slice(slice_ax, slice_coord)

    click_coords = [0, 0]
    image_number = 0

    if (plot_field == "b_tot"):
        try:
            h5ds.add_field(("gdf", plot_field), dimensions=dimensions.magnetic_field, units="", function=_total_B, display_name="Total magnetic field ($\mu$G)", sampling_type="cell")
        except:
            die("Failed to construct field %s" % plot_field)

    if (plot_field[0:-2] == "en_ratio"):
        try:
            if str(dsSlice["cren01"].units) is "dimensionless":  # DEPRECATED
                h5ds.add_field(("gdf", plot_field), units="", function=en_ratio, display_name="Ratio e/n in %i-th bin" % int(plot_field[-2:]), sampling_type="cell")
            else:
                h5ds.add_field(("gdf", plot_field), units="Msun*pc**2/Myr**2", function=en_ratio, display_name="Ratio e/n in %i-th bin" % int(plot_field[-2:]), dimensions=dimensions.energy, sampling_type="cell", take_log=True)
        except:
            die("Failed to construct field %s" % plot_field)

    dsSlice = add_tot_fields(dsSlice)

# For elegant labels when plot_field is cree?? or cren??
    if (plot_field[-3:] != "tot" and plot_field[0:3] == "cre" and plot_field[3:-2] != "ratio"):
        if (plot_field[0:4] == "cree"):
            disp_name = "energy"
            new_field_dimensions = dimensions.energy / dimensions.volume
        elif (plot_field[0:4] == "cren"):
            disp_name = "number"
            new_field_dimensions = 1. / dimensions.volume
        prtinfo("Adding display name: %s density" % disp_name)
        if (display_bin_no):
            disp_name = "CR electron %s density (bin %2i)" % (disp_name, int(plot_field[-2:]))
        else:
            disp_name = "CR electron %s density" % (disp_name)
        new_field = str(plot_field + "_updated")
        new_field_units = dsSlice[plot_field].units
        h5ds.add_field(("gdf", new_field), units=new_field_units, function=copy_field, display_name=disp_name, dimensions=new_field_dimensions, sampling_type="cell")
        plot_field = new_field

    try:
        field_max = h5ds.find_max("cr01")[0].v  # WARNING - this makes field_max unitless
    except:
        field_max = h5ds.find_max("cr1")[0].v  # WARNING - this makes field_max unitless

# prepare limits for framebuffer
    # if (options.usr_width == 0.):
    frb_w = dom_r[avail_dim[0]] + abs(dom_l[avail_dim[0]])
    if (options.usr_width != 0. and not options.annotate_rect):
        frb_w = float(options.usr_width)

    # if (options.usr_height == 0.):
    frb_h = dom_r[avail_dim[1]] + abs(dom_l[avail_dim[1]])
    if (options.usr_height != 0. and not options.annotate_rect):
        frb_h = float(options.usr_height)

    frb_center = None
    if (options.usr_center == [0., 0.] and not options.annotate_rect):
        frb_center = None
    elif (options.usr_center != [0., 0.] and not options.annotate_rect):
        frb_center = [0, 0]
        frb_center[0] = float(options.usr_center[0])
        frb_center[1] = float(options.usr_center[1])

    slice_center = "c"
    slice_center = [0, 0, 0]
    if (options.usr_center != [0., 0.] and not options.annotate_rect):
        slice_center = [0, 0, 0]
        slice_center[avail_dim[0]] = frb_center[0]
        slice_center[avail_dim[1]] = frb_center[1]
        slice_center[dim_map[slice_ax]] = slice_coord

# construct framebuffer
    if (slice_ax == "y"):
        frb = np_array(dsSlice.to_frb(width=frb_h, resolution=resolution, center=slice_center, height=frb_w)[plot_field])
        frb = rot90(frb)
    else:
        frb = np_array(dsSlice.to_frb(width=frb_w, resolution=resolution, center=slice_center, height=frb_h, periodic=False)[plot_field])
    if (not user_limits):
        plot_max = h5ds.find_max(plot_field)[0]
    if (not user_limits):
        plot_min = h5ds.find_min(plot_field)[0]
    plot_units = str(h5ds.all_data()[plot_field].units)

    if (user_limits == True):  # Overwrites previously found values
        plot_min = plot_user_min
        plot_max = plot_user_max

    if (not_py27):
        plt.xlabel("Domain cooridnates " + list(dim_map.keys())[list(dim_map.values()).index(avail_dim[0])] + " (" + length_unit + ")")
        plt.ylabel("Domain cooridnates " + list(dim_map.keys())[list(dim_map.values()).index(avail_dim[1])] + " (" + length_unit + ")")
    else:
        plt.xlabel("Domain cooridnates " + dim_map.keys()[dim_map.values().index(avail_dim[0])] + " (" + length_unit + ")")
        plt.ylabel("Domain cooridnates " + dim_map.keys()[dim_map.values().index(avail_dim[1])] + " (" + length_unit + ")")

    if (options.annotate_rect):
        yt_data_plot = yt.SlicePlot(h5ds, slice_ax, plot_field, width=(dom_r[avail_dim[0]] + abs(dom_l[avail_dim[0]]), dom_r[avail_dim[1]] + abs(dom_l[avail_dim[1]])), center=slice_center)
    else:
        yt_data_plot = yt.SlicePlot(h5ds, slice_ax, plot_field, width=(frb_w, frb_h), center=slice_center)
    yt_data_plot.set_font({'size': options.fontsize})

    plt.colormap = "plasma"
    encountered_nans = False
    plot_max = float(plot_max)
    plot_min = max(float(plot_min), par_epsilon)
    if (isnan(plot_min) == True or isnan(plot_max) == True):
        encountered_nans = True
        prtwarn("Invalid data encountered (NaN), ZLIM will be adjusted")

    im_orig = "lower"
    if (use_logscale):
        plt.imshow(frb, extent=[dom_l[avail_dim[0]], dom_r[avail_dim[0]], dom_l[avail_dim[1]], dom_r[avail_dim[1]]], origin=im_orig, norm=LogNorm(vmin=plot_min, vmax=plot_max) if (encountered_nans == False) else LogNorm())
    elif (use_linscale):
        plt.imshow(frb, extent=[dom_l[avail_dim[0]], dom_r[avail_dim[0]], dom_l[avail_dim[1]], dom_r[avail_dim[1]]], origin=im_orig, vmin=plot_min if (encountered_nans == False) else None, vmax=plot_max if (encountered_nans == False) else None)

    plt.title("Component: " + plot_field + " | t = %9.3f Myr" % time)

    try:
        cbar = plt.colorbar(shrink=0.9, pad=0.01, label=plot_units)
    except:
        die("An empty field might have been picked.")

    colormap_my = plt.cm.viridis
    colormap_my.set_bad(color=colormap_my(par_epsilon))     # masks bad values
    yt_data_plot.set_cmap(field=plot_field, cmap=colormap_my)

    if (user_annot_line == True):
        prtinfo("Marking line on yt.plot at (0 0 0) : (500 500 0)")
        yt_data_plot.annotate_line((0., 0., 0.), (500., 500.0, 0), plot_args={'color': 'white', "lw": 2.0})

    if (plot_vel):
        yt_data_plot.annotate_velocity(factor=32, scale=3e7)
    if (plot_mag):
        yt_data_plot.annotate_magnetic_field(factor=32, scale=40)
    yt_data_plot.set_zlim(plot_field, plot_min, plot_max)
    marker_l = ["x", "+", "*", "X", ".", "^", "v", "<", ">", "1"]
    m_size_l = [350, 500, 400, 400, 500, 350, 350, 350, 350, 500]
    m_e_width = 5
    marker_index = 0

    plt.subplots_adjust(left=0.075, right=0.975, hspace=0.12)
    plt.tight_layout()

    print("")

    crs_h5.crs_initialize(var_names, var_array)

    mplot = yt_data_plot.plots[plot_field]

    xticklabels = mplot.axes.xaxis.get_ticklabels()
    yticklabels = mplot.axes.yaxis.get_ticklabels()

# ---------
    def read_click_and_plot(event):
        global click_coords, image_number, f_run, marker_index
        exit_code = True
        if (marker_index == len(marker_l) or marker_index == len(m_size_l)):
            marker_index = 0
        click_coords = [event.xdata, event.ydata]
        coords = [slice_coord, slice_coord, slice_coord]
        if slice_ax == "x":
            coords[1] = click_coords[0]
            coords[2] = click_coords[1]
        elif slice_ax == "y":
            coords[0] = click_coords[0]
            coords[2] = click_coords[1]
        else:  # slice_ax = "z"
            coords[0] = click_coords[0]
            coords[1] = click_coords[1]

        mark_plot_save(coords)

    def mark_plot_save(coords):
        global click_coords, image_number, f_run, marker_index
# ------------ preparing data and passing -------------------------
        position = h5ds.r[coords:coords]
        if (plot_field[0:-2] != "en_ratio"):
            prtinfo(">>>>>>>>>>>>>>>>>>> Value of %s at point [%f, %f, %f] = %f " % (plot_field, coords[0], coords[1], coords[2], position[plot_field]))
        else:
            prtinfo("Value of %s at point [%f, %f, %f] = %f " % (plot_field, coords[0], coords[1], coords[2], position["cree" + str(plot_field[-2:])] / position["cren" + str(plot_field[-2:])]))
            plot_max = h5ds.find_max("cre" + plot_var + str(plot_field[-2:]))[0]  # once again appended - needed as ylimit for the plot

        btot = (position["mag_field_x"].v**2 + position["mag_field_y"].v**2 + position["mag_field_z"].v**2)**0.5
        btot_uG = 2.85 * btot  # WARNING magic number @btot - conversion factor
        prtinfo("B_tot = %f = %f (uG)" % (btot, btot_uG))
        if (True):   # TODO DEPRECATED save_fqp
            ecrs = []
            ncrs = []
            if (plot_ovlp != True):  # overwrites position
                if (plot_layer is True):
                    prtinfo("Plotting layer...")
                    position = h5ds.r[[coords[0], dom_l[avail_dim[0]], coords[2]]: [coords[0], dom_r[avail_dim[0]], coords[2]]]

                elif (options.annotate_rect):
                    # define borders of the selected region (plane)
                    usr_w = float(options.usr_width)
                    usr_h = float(options.usr_height)
                    usr_c = [0, 0]
                    usr_c[:] = [float(options.usr_center[0]), float(options.usr_center[1])]
                    lb = [0, 0, 0]
                    lb[avail_dim[0]] = usr_c[0] - usr_w * 0.5
                    lb[avail_dim[1]] = usr_c[1] - usr_h * 0.5
                    lb[dim_map[slice_ax]] = slice_coord - (dom_r[avail_dim[0]] - dom_l[avail_dim[0]]) / int(h5ds.domain_dimensions[avail_dim[0]])
                    rb = [0, 0, 0]
                    rb[avail_dim[0]] = usr_c[0] + usr_w * 0.5
                    rb[avail_dim[1]] = usr_c[1] - usr_h * 0.5
                    rb[dim_map[slice_ax]] = slice_coord
                    lt = [0, 0, 0]
                    lt[avail_dim[0]] = usr_c[0] - usr_w * 0.5
                    lt[avail_dim[1]] = usr_c[1] + usr_h * 0.5
                    lt[dim_map[slice_ax]] = slice_coord
                    rt = [0, 0, 0]
                    rt[avail_dim[0]] = usr_c[0] + usr_w * 0.5
                    rt[avail_dim[1]] = usr_c[1] + usr_h * 0.5
                    rt[dim_map[slice_ax]] = slice_coord + (dom_r[avail_dim[0]] - dom_l[avail_dim[0]]) / int(h5ds.domain_dimensions[avail_dim[0]])

                    # select region as plane spreading between lb and rt corners
                    position = yt.data_objects.selection_data_containers.YTRegion(left_edge=lb, right_edge=rt, center=usr_c, ds=h5ds)  # (center, left_edge, right_edge
                    coords[avail_dim[0]:avail_dim[1]] = usr_c[:]
                    coords[dim_map[slice_ax]] = slice_coord

                for ind in range(1, ncre + 1):
                    ecrs.append(float(mean(position['cree' + str(ind).zfill(2)][0].v)))
                    ncrs.append(float(mean(position['cren' + str(ind).zfill(2)][0].v)))

                fig2, exit_code = crs_h5.crs_plot_main(plot_var, ncrs, ecrs, time, coords, marker=marker_l[marker_index], clean_plot=options.clean_plot, hide_axes=options.no_axes)

            elif (plot_ovlp == True):  # for overlap_layer
                prtinfo("Plotting layer with overlap...")
                dnum = int(h5ds.domain_dimensions[avail_dim[0]])
                dl = (dom_r[avail_dim[0]] - dom_l[avail_dim[0]]) / float(dnum)
                for j in range(dnum):
                    position = position = h5ds.r[[coords[0], dom_l[avail_dim[0]] + dl * j, coords[2]]: [coords[0], dom_l[avail_dim[0]] + dl * j, coords[2]]]
                    for ind in range(1, ncre + 1):
                        ecrs.append(position['cree' + str(ind).zfill(2)][0].v)
                        ncrs.append(position['cren' + str(ind).zfill(2)][0].v)
                    fig2, exit_code_tmp = crs_h5.crs_plot_main(plot_var, ncrs, ecrs, time, coords, marker=marker_l[marker_index], i_plot=image_number, clean_plot=options.clean_plot, hide_axes=options.no_axes)
                    if (exit_code_tmp == False):
                        exit_code = exit_code_tmp  # Just one plot is allright
                    ecrs = []
                    ncrs = []
        else:     # for fqp, DEPRECATED probably
            fcrs = []
            qcrs = []
            pcut = [0., 0.]
            ecrs = []
            ncrs = []

            for ind in range(1, ncre + 1):
                ecrs.append(float(position['cree' + str(ind).zfill(2)][0].v))
                ncrs.append(float(position['cren' + str(ind).zfill(2)][0].v))

            for ind in range(1, ncre + 2):
                fcrs.append(float(position['cref' + str(ind).zfill(2)][0].v))
            for ind in range(1, ncre + 1):
                qcrs.append(float(position['creq' + str(ind).zfill(2)][0].v))
            pcut[:] = [position['crep01'][0].v, position['crep02'][0].v]

            fig2, exit_code = crs_h5.crs_plot_main_fpq(var_names, var_array, plot_var, fcrs, qcrs, pcut, field_max, time, coords, marker=marker_l[marker_index], clean_plot=options.clean_plot)

        if (exit_code != True):
            if ((plot_layer == True) or (plot_ovlp == True)):
                line = s1.plot([dom_l[avail_dim[0]], dom_r[avail_dim[0]]], [coords[2], coords[2]], color="white")    # plot line (layer) if cell not empty WARNING - for now only works with mcrwind
            else:
                point = s1.plot(coords[avail_dim[0]], coords[avail_dim[1]], marker=marker_l[marker_index], color="red")  # plot point if cell not empty
            s.savefig('results/' + filename_nam + '_' + plot_var + '_%04d.png' % image_number, transparent='True')
            #prtinfo("  --->  Saved plot to: %s " %str('results/'+filename_nam+'_'+plot_var+'_%04d.png' %image_number))

            if (plot_layer == True):  # Mark averaged level
                yt_data_plot.annotate_line([coords[0], dom_l[avail_dim[0]], coords[2]], [coords[0], dom_r[avail_dim[0]], coords[2]], plot_args={'color': 'white', "lw": 10.0})
            else:
                if (not options.annotate_rect and marker_index > 0):
                    yt_data_plot.annotate_marker(coords, marker=marker_l[marker_index - 1], plot_args={'color': 'red', 's': m_size_l[marker_index - 1], "lw": 4.5})  # cumulatively annotate all clicked coordinates
            if (options.annotate_rect):
                yt_data_plot.annotate_line(lb, lt, plot_args={'color': 'white', "lw": 2.0})
                yt_data_plot.annotate_line(lt, rt, plot_args={'color': 'white', "lw": 2.0})
                yt_data_plot.annotate_line(rt, rb, plot_args={'color': 'white', "lw": 2.0})
                yt_data_plot.annotate_line(rb, lb, plot_args={'color': 'white', "lw": 2.0})

            marker_index = marker_index + 1
            image_number = image_number + 1
# ------------- saving just the spectrum
            if (save_spectrum):
                extent = fig2.get_window_extent().transformed(s.dpi_scale_trans.inverted())
                s.savefig('results/' + filename_nam + '_' + plot_var + '_spec_%03d.pdf' % image_number, transparent='True', bbox_inches="tight", quality=95, dpi=150)  # bbox not working in py27 FIXME
                prtinfo("  --->  Saved plot to: %s.\n\033[44mPress 'q' to quit and save yt.SlicePlot with marked coordinates." % str('results/' + filename_nam + '_' + plot_var + '_spectrum_%04d.pdf' % image_number))
        else:
            prtwarn("Empty cell - not saving.")

        if (f_run):
            f_run = False

    def plot_with_coords_provided(coords_in):
        mark_plot_save(coords_in)

    if (user_coords_provided):
        prtinfo("Provided coordintates %s (clickable map will not be shown) , processing image." % str(user_coords))
        plot_with_coords_provided(user_coords)
    else:
        prtinfo("\033[44mClick LMB on the colormap to display spectrum ('q' to exit)")
        cid = s.canvas.mpl_connect('button_press_event', read_click_and_plot)
        plt.show()

    text_coords = [0., 0., 0.]
    text_coords[dim_map.get(slice_ax)] = slice_coord
    text_coords[avail_dim[0]] = dom_l[avail_dim[0]]
    text_coords[avail_dim[1]] = dom_l[avail_dim[1]]
    text_coords = [item * 0.9 for item in text_coords]

    if (user_annot_time):
        if (user_draw_timebox == True):
            yt_data_plot.annotate_text(text_coords, 'T = {:0.2f} Myr'.format(float(t.in_units('Myr'))), text_args={'fontsize': options.fontsize, 'color': 'white', 'alpha': '0.0'}, inset_box_args={'boxstyle': 'round', 'pad': 0.2, 'alpha': 0.8})
        else:
            prtinfo("Not marking line on yt.plot (user_draw_timebox = %s)" % (user_draw_timebox))
            yt_data_plot.annotate_title('T = {:0.2f} Myr'.format(float(t.in_units('Myr'))))

    if (options.no_cbar):
        yt_data_plot.hide_colorbar()
    if (options.no_xlabels and options.no_ylabels):
        yt_data_plot.hide_axes()

    yt_data_plot.save('results/' + filename_nam + '_' + plot_field + '_sliceplot.pdf')  # save image (spectrum already saved) when finished.
    if (not user_coords_provided):
        s.canvas.mpl_disconnect(cid)
