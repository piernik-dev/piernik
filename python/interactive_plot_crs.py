#!/usr/bin/python
# -*- coding: utf-8 -*-
from colored_io import die, prtinfo, prtwarn, read_var
from copy import copy
from crs_h5 import crs_initialize, crs_plot_main, crs_plot_main_fpq, crs_plot_ratio
from crs_pf import initialize_pf_arrays
from math import isnan, pi
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from numpy import array as np_array, log, log10, mean, rot90, shape
from os import getcwd, makedirs, path
from optparse import OptionParser
from re import search
from read_h5 import read_par, input_names_array
from sys import argv, version
from warnings import simplefilter
try:
    import yt
    from yt.units import dimensions
except:
    die("You must make yt available somehow")
if (version[0:3] != "2.7"):
    raw_input = input
    not_py27 = True
else:
    not_py27 = False

# ------- Parse arguments
parser = OptionParser(
    "Usage: %prog FILE [options] [args] or %prog [options] [args] -F FILENAME")
parser.add_option("-F", "--file", dest="filename",
                  default="None", help=u"File to use", type="str")
parser.add_option("-v", "--var", dest="var_name", default="n",
                  help=u"Variable to plot the spectrum (default: n)")
parser.add_option("-f", "--field", dest="fieldname", default="",
                  help=u"DS fieldname to image (default:cree_tot)")
parser.add_option("-z", "--zlim", dest="plot_range", default=("x", "y"),
                  help=u"Plot image with this range", nargs=2, metavar="ZMIN ZMAX")
parser.add_option("-s", "--slice", dest="slice_info", default=("a", "Nan"),
                  help=u"DS slice coords to image (default:cree_tot)", metavar="AX COORDINATE", nargs=2)
parser.add_option("-S", "--species", dest="crspecies", default="e-",
                  help=u"DS fieldname of cr species (default:cr_e-)", metavar="CR FIELDNAME", nargs=1)
parser.add_option("-d", "--def", dest="default_range", default=False,
                  help=u"Use min/max on yt.ds(fieldname) for clickable image", action="store_true")
parser.add_option("-l", "--lin", dest="use_linscale", default=False,
                  help=u"Use linear scale for clickable plot (default: log)", action="store_true")
parser.add_option("-V", "--vel", dest="plot_vel", default=False,
                  help=u"Plot velocity field vectors ", action="store_true")
parser.add_option("-m", "--mag", dest="plot_mag", default=False,
                  help=u"Plot magnetic field vectors ", action="store_true")
parser.add_option("-t", "--time", dest="annotate_time", default=False,
                  help=u"Annotate time on resulting DS plot", action="store_true")
parser.add_option("", "--nosave", dest="not_save_spec", default=False,
                  help=u"Do not save output spectrum ", action="store_true")
parser.add_option("-a", "--average", "--avg", dest="avg_layer", default=False,
                  help=u"Plot mean spectrum at pointed layer at ordinate axis", action="store_true")
parser.add_option("-O", "--overlap", dest="overlap_layer", default=False,
                  help=u"Overlap all spectra at pointed layer at ordinate axis", action="store_true")
parser.add_option("", "--verbose", dest="yt_verbose", default=40,
                  help=u"Append yt verbosity (value from 10 [high] to 50 [low])")
parser.add_option("-q", "--quiet", dest="py_quiet", default=False,
                  help=u"Suppress ALL python warnings (not advised)", action="store_true")
parser.add_option("-c", "--coords", dest="coords_dflt", default="",
                  help=u"Provides coordinates for the first plots", nargs=1, metavar="x1,y1,z1:x2,y2,z2:...")
parser.add_option("-k", "--keep", dest="clean_plot", default=True,
                  help=u"Keep spectrum plot (do not clean spectrum field)", action="store_false")
parser.add_option("", "--fontsize", dest="fontsize", default=18,
                  help=u"Set fontsize for SlicePlot (default is 18)")
parser.add_option("", "--noxlabels", dest="no_xlabels", default=False,
                  help=u"Do not show labels of X axis on resulting SlicePlot", action="store_true")
parser.add_option("", "--noylabels", dest="no_ylabels", default=False,
                  help=u"Do not show labels of Y axis on resulting SlicePlot", action="store_true")
parser.add_option("", "--nocbar", dest="no_cbar", default=False,
                  help=u"Do not show colorbar on the resulting SlicePlot", action="store_true")
parser.add_option("", "--noaxes", dest="no_axes", default=False,
                  help=u"Hide axes on the spectrum plot (useful when combining spectra)", action="store_true")
parser.add_option("", "--rectangle", dest="annotate_rect", default=False,
                  help=u"Annotate and average over width/height rectangle surface", action="store_true")
parser.add_option("", "--width", dest="usr_width",
                  default=0., help=u"Set custom frb width")
parser.add_option("", "--height", dest="usr_height",
                  default=0., help=u"Set custom frb width")
parser.add_option("", "--center", dest="usr_center", default=(0., 0.),
                  help=u"Set custom frb center", nargs=2, metavar="XC YC")
(options, args) = parser.parse_args(argv[1:])  # argv[1] is filename
# Reduces the output to desired level, 50 - least output
yt.mylog.setLevel(int(options.yt_verbose))
if (options.py_quiet is True):
    simplefilter(action='ignore', category=FutureWarning)
    simplefilter(action='ignore', category=Warning)
    prtwarn("Python warnings are turned off from here (-q, --quiet switch)")
plot_var = options.var_name
user_draw_timebox = options.annotate_time
user_limits = (options.default_range is not True)
save_spectrum = (options.not_save_spec is not True)
use_logscale = (options.use_linscale is not True)
use_linscale = (options.use_linscale)

spc_label = options.crspecies.replace("cr_", "")
spc_n_lab = "cr_" + spc_label + "n"  # -> "cr_e-n" default
spc_e_lab = "cr_" + spc_label + "e"  # -> "cr_e-e" default
if (options.crspecies[0:3] != "cr_"):
    options.crspecies = "cr_" + options.crspecies
plot_field = options.fieldname if len(
    options.fieldname) > 1 else options.crspecies + "n_tot"  # DEFAULT cr_e-e_tot

# plot_var = options.var_name
plot_vel = options.plot_vel
plot_mag = options.plot_mag
spc_label = options.crspecies
if (plot_vel is True):
    plot_mag = False
if (plot_mag is True):
    plot_vel = False
user_annot_line = False
user_annot_time = True
plot_layer = options.avg_layer
plot_ovlp = options.overlap_layer
options.fontsize = int(options.fontsize)
display_bin_no = False
plot_CRisotope = False
user_coords_provided = not options.coords_dflt == ""
if user_coords_provided:
    try:
        user_coords = []
        aux = options.coords_dflt.split(":")
        for auxi in aux:
            user_coords.append(list(float(item) for item in auxi.split(",")))
    except:
        die("Got spectrum coordinates %s, but failed to convert it to float." %
            str(options.coords_dflt))

#####################################
# ------- Parameters ----------------
par_epsilon = 1.0e-15
f_run = True
pf_initialized = False

proton_field_names = ["cr_p+n01", "cr01", "cr1", "cr_p+"]

# ------- Local functions -----------


def _total_cr_e(field, data):
    list_cr_e = []
    for element in h5ds.field_list:
        # print(element, spc_e_lab, str(element[1]), spc_e_lab in element)
        if search(spc_e_lab.replace("+", "\+"), str(element[1])):
            list_cr_e.append(element[1])
    cr_e_tot = data[str(list_cr_e[0])]
    for element in list_cr_e[1:]:
        cr_e_tot = cr_e_tot + data[element]
    return cr_e_tot


def _total_cr_n(field, data):
    list_cr_n = []
    for element in h5ds.field_list:
        if search(spc_n_lab.replace("+", "\+"), str(element[1])):
            list_cr_n.append(element[1])
    cr_n_tot = data[str(list_cr_n[0])]
    for element in list_cr_n[1:]:
        cr_n_tot = cr_n_tot + data[element]
    return cr_n_tot


def _total_B(field, data):
    b_tot = 2.85 * (data["mag_field_x"]**2 +
                    data["mag_field_y"]**2 + data["mag_field_y"]**2)**0.5
    return b_tot


def en_ratio(field, data):  # DEPRECATED (?)
    bin_nr = field.name[1][-2:]
    for element in h5ds.field_list:
        if search(spc_n_lab + str(bin_nr.zfill(2)), str(element[1])):
            cren_data = data[spc_n_lab + str(bin_nr.zfill(2))]
            # necessary to avoid FPEs
            cren_data[cren_data <= par_epsilon**2] = par_epsilon
            cree_data = data[spc_e_lab + str(bin_nr.zfill(2))]
            en_ratio = cree_data / cren_data
    return en_ratio


def BC_ratio(field, data):  # Boron to Carbon
    bin_nr = field.name[1][-2:]
    BC_ratio = []

    for ind in range(1, ncrb + 1):
        print(ind)
        for element1 in h5ds.field_list:
            if ("cr_B11n" + str(ind).zfill(2) == str(element1[1])):
                for element2 in h5ds.field_list:
                    if ("cr_C12n" + str(ind).zfill(2) == str(element2[1])):
                        Bn_data = data["cr_B11n" + str(ind).zfill(2)]
                        Cn_data = data["cr_C12n" + str(ind).zfill(2)]
                        BC_ratio.append(Bn_data / Cn_data)

    print('BC_ratio : ')
    print(BC_ratio)
    print('shape :')
    print(shape(BC_ratio))
    return BC_ratio


def copy_field(field, data):
    field_name_to_copy = field.name[1][:].split("_")[0]
    copied_field = data[field_name_to_copy]
    return copied_field


def add_cr_n_tot_to(h5_dataset, name):
    try:
        if (h5ds.all_data()[name + "01"].units == "dimensionless"):
            h5ds.add_field(("gdf", name + "_tot"), units="", function=_total_cr_n,
                           display_name="Total CR electron number density", sampling_type="cell")
        else:
            h5ds.add_field(("gdf", name + "_tot"), units="Msun/(Myr**2*pc)", function=_total_cr_n, display_name="Total CR electron number density",
                           dimensions=dimensions.energy / dimensions.volume, sampling_type="cell", take_log=True)
    except:
        die("Failed to construct field '%s_tot'" % name)
    return h5_dataset


def _total_cr_species_n(field, data):
    global field_name_total_n
    list_crn = []
    for element in h5ds.field_list:
        if search(field_name_total_n, str(element[1])):
            list_crn.append(element[1])
    CR_species_n_tot = data[str(list_crn[0])]
    for element in list_crn[1:]:
        CR_species_n_tot = CR_species_n_tot + data[element]
    return CR_species_n_tot


def add_cr_e_tot_to(h5_dataset, name):
    try:
        if (h5ds.all_data()[name + '01'].units == "dimensionless"):
            h5ds.add_field(("gdf", name + "_tot"), units="", function=_total_cr_e,
                           display_name="Total CR species energy density", sampling_type="cell")
        else:
            h5ds.add_field(("gdf", name + "_tot"), units="Msun/(Myr**2*pc)", function=_total_cr_e,
                           display_name="Total CR species energy density", dimensions=dimensions.energy / dimensions.volume, sampling_type="cell")
    except:
        die("Failed to construct field '" + name + "e_tot'")
    return h5_dataset


def add_total_n_to(h5_dataset, name):
    try:
        if (h5ds.all_data()[name].units == "dimensionless"):
            h5ds.add_field(("gdf", name), units="", function=_total_cr_species_n,
                           display_name="Total CR " + name + " number density", sampling_type="cell")
        else:
            h5ds.add_field(("gdf", name), units="Msun/(Myr**2*pc)", function=_total_cr_species_n, display_name="Total CR " + name + " number density",
                           dimensions=dimensions.energy / dimensions.volume, sampling_type="cell", take_log=True)
            # TODO BUG units should be "Msun/(Myr**2*pc)"; fix it after fixing it in PIERNIK!
    except:
        die("Failed to construct field '" + name + "n_tot'")
    return h5_dataset


def add_tot_fields(h5_dataset):
    global plot_CRisotope
    if (plot_field[-4:] == "_tot"):
        print("add_tot_fields, plot_field is:", plot_field)
        h5_dataset = add_cr_e_tot_to(h5_dataset, spc_e_lab)
        h5_dataset = add_cr_n_tot_to(h5_dataset, spc_n_lab)
    else:
        # h5_dataset = add_total_n_to(h5_dataset, plot_field)
        plot_CRisotope = True
    return h5_dataset


# ---------- reading parameters
if (options.filename != "None"):
    filename = options.filename
else:
    filename = argv[1]

filename_trimmed = filename.split("/")[-1]
filename_ext = filename_trimmed.split('.')[-1]
filename_nam = filename_trimmed.split('.')[0].split('/')[-1]
if (filename_ext != 'h5'):
    die("Script requires a (list of) hdf5 file(s) on input")

if f_run:
    output_path = path.realpath(
        "./" + filename.strip(filename_trimmed) + "/") + "/results"
    if not path.exists(output_path):
        makedirs(output_path)
        prtinfo("Output directory created: %s" % output_path)

var_array = []
if f_run is True:
    var_names = []
    var_names = ["ncrb", "p_min_fix", "p_max_fix",
                 "e_small", "cre_eff", "q_big"]
    var_def = [20, 10., 1.e5, 1.e-6, 0.01, 30., ]
    if len(var_names) == 0:
        prtwarn(
            "Empty list of parameter names provided: enter names of parameters to read")
        var_names = input_names_array()

    var_array = read_par(filename, var_names, var_def)
    for i in range(len(var_names)):
        exec("%s=%s" % (var_names[i], var_array[i]))

    prtinfo("\n*** Values read from problem.par@hdf5 file: *** \n")
    for i in range(len(var_names)):
        prtinfo(" %15s =  %10s ( %15s  ) " %
                (var_names[i], var_array[i], type(var_array[i])))

# ---------- Open file
    h5ds = yt.load(filename)

    initialize_pf_arrays(filename, pf_initialized)
# ---------- bounds on domain size
    grid_dim = h5ds.domain_dimensions
    dim_map = {'x': 0, 'y': 1, 'z': 2}
    dom_l = np_array(h5ds.domain_left_edge[0:3])
    dom_r = np_array(h5ds.domain_right_edge[0:3])
    prtinfo("Max level of refinement is %i." % (h5ds.max_level))

# ----------- Loading other data
    t = h5ds.current_time[0]
    time = t.in_units('Myr')
# ----------- Checking user image limits
    try:
        plot_user_min = float(options.plot_range[0])
        plot_user_max = float(options.plot_range[1])
    except:
        prtwarn(
            "No provided ZLIM or failed to convert it into float -- using default (min/max of ds(field) )")
        user_limits = False
# ------------ Organizing domain data
    length_unit = 'pc'
    prtinfo("Domain shape of in provided file            (i, j, k): [%i, %i, %i] \033[0m" % (
        grid_dim[0], grid_dim[1], grid_dim[2]))
    prtinfo("Domain physical dimensions in provided file (x, y, z): [%9.3f,%9.3f,%9.3f]:[%9.3f,%9.3f,%9.3f] %s \033[0m" % (
        dom_l[0], dom_l[1], dom_l[2], dom_r[0], dom_r[1], dom_r[2], length_unit))

    avail_dim = [0, 1, 2]
    avail_dims_by_slice = [[1, 2], [0, 2], [0, 1]]

    if len(grid_dim) == 3 and min(grid_dim) != 1:
        slice_ax = str(options.slice_info[0])
        slice_coord = float(options.slice_info[1])
        while slice_ax not in dim_map.keys():
            slice_ax = read_var("Choose slice ax (x, y, z)      : ")
        # or slice_coord < -10000:
        while (slice_coord < dom_l[dim_map[slice_ax]]) or (slice_coord > dom_r[dim_map[slice_ax]] or isnan(slice_coord) is True):
            try:
                slice_coord = float(read_var("Choose slice coordinate (%d:%d %s ) (if empty, middle is assumed): \033[0m" % (
                    dom_l[dim_map[slice_ax]], dom_r[dim_map[slice_ax]], length_unit)))
            except:
                slice_coord = (dom_l[dim_map[slice_ax]] +
                               dom_r[dim_map[slice_ax]]) / 2.
                prtwarn(
                    " (empty / incorrect input): Setting slice coordinate to %s %s.\033[0m" % (slice_coord, length_unit))
    elif min(grid_dim) == 1:
        slice_coord = 0.0
        if grid_dim[0] == 1:
            slice_ax = 'x'
        elif grid_dim[1] == 1:
            slice_ax = 'y'
        else:
            slice_ax = 'z'
    avail_dim = avail_dims_by_slice[dim_map[slice_ax]]
    prtinfo("Slice ax set to %s, coordinate = %f %s \033[0m" % (
        slice_ax, slice_coord, length_unit))
    resolution = [grid_dim[avail_dim[0]] * 2**h5ds.max_level,
                  grid_dim[avail_dim[1]] * 2**h5ds.max_level]

# --------- Preparing clickable image
    s = plt.figure("Displayed file: %s at %s %s of axis %s" % (
        filename, str(slice_coord), length_unit, slice_ax), figsize=(12, 8), dpi=100)
    s1 = plt.subplot(121)
    dsSlice = h5ds.slice(slice_ax, slice_coord)

    click_coords = [0, 0]
    image_number = 0

    if (plot_field == "b_tot"):
        try:
            h5ds.add_field(("gdf", plot_field), dimensions=dimensions.magnetic_field, units="",
                           function=_total_B, display_name="Total magnetic field ($\mu$G)", sampling_type="cell")
        except:
            die("Failed to construct field %s" % plot_field)

    if (plot_field[0:-2] == "en_ratio"):
        try:
            if str(dsSlice["cren01"].units) == "dimensionless":  # DEPRECATED
                h5ds.add_field(("gdf", plot_field), units="", function=en_ratio,
                               display_name="Ratio e/n in %i-th bin" % int(plot_field[-2:]), sampling_type="cell")
            else:
                h5ds.add_field(("gdf", plot_field), units="Msun*pc**2/Myr**2", function=en_ratio, display_name="Ratio e/n in %i-th bin" %
                               int(plot_field[-2:]), dimensions=dimensions.energy, sampling_type="cell", take_log=True)
        except:
            die("Failed to construct field %s" % plot_field)

    if (plot_field == "cr_B11n_"):
        print('hello ! ')
        # try:
        print(dsSlice["cr_B11n01"].units)
        """
        if str(dsSlice["cr_B11n01"].units) == "":  # DEPRECATED
            print('case 1 ')
            h5ds.add_field(("gdf", plot_field), units="", function=BC_ratio,
                           display_name="Ratio B/C in %i-th bin" % int(plot_field[-2:]), sampling_type="cell")
        else:
        """
        h5ds.add_field(("gdf", plot_field), units="", function=BC_ratio, display_name="Ratio B/C in %i-th bin",
                       dimensions=dimensions.energy, sampling_type="cell", take_log=True)
        # except:
        # die("Failed to construct field %s" % plot_field)

    dsSlice = add_tot_fields(dsSlice)

# For elegant labels when plot_field is cree?? or cren??
    if (plot_field[-3:] != "tot" and plot_field[0:3] == "cre" and plot_field[3:-2] != "ratio"):
        if (plot_field[0:4] == spc_e_lab):
            disp_name = "energy"
            new_field_dimensions = dimensions.energy / dimensions.volume
        elif (plot_field[0:4] == spc_n_lab):
            disp_name = "number"
            new_field_dimensions = 1. / dimensions.volume
        prtinfo("Adding display name: %s density" % disp_name)
        if (display_bin_no):
            disp_name = "CR electron %s density (bin %2i)" % (
                disp_name, int(plot_field[-2:]))
        else:
            disp_name = "CR electron %s density" % (disp_name)
        new_field = str(plot_field + "_updated")
        new_field_units = dsSlice[plot_field].units
        h5ds.add_field(("gdf", new_field), units=new_field_units, function=copy_field,
                       display_name=disp_name, dimensions=new_field_dimensions, sampling_type="cell")
        plot_field = new_field

    # WARNING - this makes field_max unitless
    for proton_field in proton_field_names:
        try:
            # WARNING - this makes field_max unitless
            field_max = h5ds.find_max(proton_field)[0].v
            break
        except:
            print("MAX for proton field ", proton_field, " not found")

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
    frbuffer_plot_field = plot_field
    # if (plot_CRisotope): frbuffer_plot_field = "cr_"+plot_field+"n_tot"
    if (slice_ax == "y"):
        frb = np_array(dsSlice.to_frb(width=frb_h, resolution=resolution,
                       center=slice_center, height=frb_w)[frbuffer_plot_field])
        frb = rot90(frb)
    else:
        frb = np_array(dsSlice.to_frb(width=frb_w, resolution=resolution,
                       center=slice_center, height=frb_h, periodic=False)[frbuffer_plot_field])
    if (not user_limits):
        plot_max = h5ds.find_max(frbuffer_plot_field)[0]
    if (not user_limits):
        plot_min = h5ds.find_min(frbuffer_plot_field)[0]
    plot_units = str(h5ds.all_data()[frbuffer_plot_field].units)

    if (user_limits is True):  # Overwrites previously found values
        plot_min = plot_user_min
        plot_max = plot_user_max

    if (not_py27):
        plt.xlabel("Domain coordinates " + list(dim_map.keys())
                   [list(dim_map.values()).index(avail_dim[0])] + " (" + length_unit + ")")
        plt.ylabel("Domain coordinates " + list(dim_map.keys())
                   [list(dim_map.values()).index(avail_dim[1])] + " (" + length_unit + ")")
    else:
        plt.xlabel("Domain coordinates " + dim_map.keys()
                   [dim_map.values().index(avail_dim[0])] + " (" + length_unit + ")")
        plt.ylabel("Domain coordinates " + dim_map.keys()
                   [dim_map.values().index(avail_dim[1])] + " (" + length_unit + ")")

    if (options.annotate_rect):
        yt_data_plot = yt.SlicePlot(h5ds, slice_ax, plot_field, width=(dom_r[avail_dim[0]] + abs(
            dom_l[avail_dim[0]]), dom_r[avail_dim[1]] + abs(dom_l[avail_dim[1]])), center=slice_center)
    else:
        yt_data_plot = yt.SlicePlot(h5ds, slice_ax, frbuffer_plot_field, width=(
            frb_w, frb_h), center=slice_center)
    yt_data_plot.set_font({'size': options.fontsize})

    encountered_nans = False
    plot_max = float(plot_max)
    plot_min = max(float(plot_min), par_epsilon)
    if (isnan(plot_min) is True or isnan(plot_max) is True):
        encountered_nans = True
        prtwarn("Invalid data encountered (NaN), ZLIM will be adjusted")

    colormap_my = copy(plt.cm.viridis)
    colormap_my.set_bad(colormap_my(par_epsilon))

    im_orig = "lower"
    if (use_logscale):
        plt.imshow(frb, extent=[dom_l[avail_dim[0]], dom_r[avail_dim[0]], dom_l[avail_dim[1]], dom_r[avail_dim[1]]], origin=im_orig,
                   cmap=colormap_my, norm=LogNorm(vmin=plot_min, vmax=plot_max) if (encountered_nans is False) else LogNorm())
    elif (use_linscale):
        plt.imshow(frb, extent=[dom_l[avail_dim[0]], dom_r[avail_dim[0]],
                   dom_l[avail_dim[1]], dom_r[avail_dim[1]]], origin=im_orig, cmap=colormap_my)

    plt.title("Component: " + plot_field + " | t = %9.3f Myr" % time)
    try:
        cbar = plt.colorbar(shrink=0.9, pad=0.01, label=plot_units)
    except:
        die("An empty field might have been picked.")

    if (user_annot_line is True):
        prtinfo("Marking line on yt.plot at (0 0 0) : (500 500 0)")
        yt_data_plot.annotate_line((0., 0., 0.), (500., 500.0, 0), plot_args={
                                   'color': 'white', "lw": 2.0})

    if (plot_vel):
        yt_data_plot.annotate_velocity(factor=32, scale=3e7)
    if (plot_mag):
        yt_data_plot.annotate_magnetic_field(factor=32, scale=40)

    yt_data_plot.set_cmap(field=frbuffer_plot_field, cmap=colormap_my)
    yt_data_plot.set_zlim(frbuffer_plot_field, plot_min, plot_max)

    marker_l = ["x", "+", "*", "X", ".", "^", "v", "<", ">", "1"]
    m_size_l = [350, 500, 400, 400, 500, 350, 350, 350, 350, 500]
    m_e_width = 5
    marker_index = 0

    plt.subplots_adjust(left=0.075, right=0.975, hspace=0.12)
    plt.tight_layout()

    print("")

    crs_initialize(var_names, var_array, plot_field)

    mplot = yt_data_plot.plots[frbuffer_plot_field]

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
        fieldname = spc_label
        if (fieldname[-3] == "e" or fieldname[-3] == "n"):
            # If just one bin is plotted on clickable field, strip the bin number + quantity from fieldname
            fieldname = plot_field[0:-3]
        plot_field_click = frbuffer_plot_field
        """
        if (plot_field[0:-2] != "en_ratio"):
            prtinfo(">>>>>>>>>>>>>>>>>>> Value of %s at point [%f, %f, %f] = %f " % (
                plot_field_click, coords[0], coords[1], coords[2], position[plot_field_click]))
        """
        print('Here it comes : ')
        print(plot_field)

        if (plot_field == "cr_B11n_tot"):  # Plot B to C ratio rather than spectra

            BC_ratio = []

            for ind in range(1, ncrb + 1):

                for element1 in h5ds.field_list:
                    if ("cr_B11n" + str(ind).zfill(2) == str(element1[1])):
                        for element2 in h5ds.field_list:
                            if ("cr_C12n" + str(ind).zfill(2) == str(element2[1])):
                                Bn_data = position["cr_B11n" +
                                                   str(ind).zfill(2)]
                                Cn_data = position["cr_C12n" +
                                                   str(ind).zfill(2)]
                                BC_ratio.append(Bn_data / Cn_data)

            BC_ratio = np_array(BC_ratio)
            print('BC_ratio : ')
            print(BC_ratio)
            print('shape :')
            print(shape(BC_ratio))

        else:

            if (plot_field[-3:] != 'tot'):
                prtinfo("Value of %s at point [%f, %f, %f] = %f " % (plot_field_click, coords[0], coords[1],
                        coords[2], position["cree" + str(plot_field_click[-2:])] / position["cren" + str(plot_field_click[-2:])]))
                # once again appended - needed as ylimit for the plot
                plot_max = h5ds.find_max(
                    "cre" + plot_var + str(plot_field_click[-2:]))[0]

        btot = (position["mag_field_x"].v**2 + position["mag_field_y"].v **
                2 + position["mag_field_z"].v**2)**0.5
        btot_uG = 2.85 * btot  # WARNING magic number @btot - conversion factor
        prtinfo("B_tot = %f = %f (uG)" % (btot, btot_uG))
        if (True):   # TODO DEPRECATED save_fqp
            ecrs = []
            ncrs = []
            if (plot_ovlp is not True):  # overwrites position
                if (plot_layer is True):
                    prtinfo("Plotting layer...")
                    position = h5ds.r[[coords[0], dom_l[avail_dim[0]], coords[2]]: [
                        coords[0], dom_r[avail_dim[0]], coords[2]]]

                elif (options.annotate_rect):
                    # define borders of the selected region (plane)
                    usr_w = float(options.usr_width)
                    usr_h = float(options.usr_height)
                    usr_c = [0, 0]
                    usr_c[:] = [float(options.usr_center[0]),
                                float(options.usr_center[1])]
                    lb = [0, 0, 0]
                    lb[avail_dim[0]] = usr_c[0] - usr_w * 0.5
                    lb[avail_dim[1]] = usr_c[1] - usr_h * 0.5
                    lb[dim_map[slice_ax]] = slice_coord - \
                        (dom_r[avail_dim[0]] - dom_l[avail_dim[0]]) / \
                        int(h5ds.domain_dimensions[avail_dim[0]])
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
                    rt[dim_map[slice_ax]] = slice_coord + \
                        (dom_r[avail_dim[0]] - dom_l[avail_dim[0]]) / \
                        int(h5ds.domain_dimensions[avail_dim[0]])

                    # select region as plane spreading between lb and rt corners
                    position = yt.data_objects.selection_data_containers.YTRegion(
                        left_edge=lb, right_edge=rt, center=usr_c, ds=h5ds)  # (center, left_edge, right_edge
                    coords[avail_dim[0]:avail_dim[1]] = usr_c[:]
                    coords[dim_map[slice_ax]] = slice_coord

                for ind in range(1, ncrb + 1):
                    ecrs.append(
                        float(mean(position[spc_e_lab + str(ind).zfill(2)][0].v)))
                    ncrs.append(
                        float(mean(position[spc_n_lab + str(ind).zfill(2)][0].v)))

                fig2, exit_code = crs_plot_main(
                    plot_var, ncrs, ecrs, time, coords, marker=marker_l[marker_index], clean_plot=options.clean_plot, hide_axes=options.no_axes)

            elif (plot_ovlp is True):  # for overlap_layer
                prtinfo("Plotting layer with overlap...")
                dnum = int(h5ds.domain_dimensions[avail_dim[0]])
                dl = (dom_r[avail_dim[0]] - dom_l[avail_dim[0]]) / float(dnum)
                for j in range(dnum):
                    position = position = h5ds.r[[coords[0], dom_l[avail_dim[0]] + dl * j, coords[2]]: [
                        coords[0], dom_l[avail_dim[0]] + dl * j, coords[2]]]
                    for ind in range(1, ncrb + 1):
                        ecrs.append(
                            position['cr_' + fieldname + 'e' + str(ind).zfill(2)][0].v)
                        ncrs.append(
                            position['cr_' + fieldname + 'n' + str(ind).zfill(2)][0].v)
                    fig2, exit_code_tmp = crs_plot_main(
                        plot_var, ncrs, ecrs, time, coords, marker=marker_l[marker_index], i_plot=image_number, clean_plot=options.clean_plot, hide_axes=options.no_axes)
                    if (exit_code_tmp is False):
                        exit_code = exit_code_tmp  # Just one plot is allright
                    ecrs = []
                    ncrs = []
        else:     # for fqp, DEPRECATED probably
            fcrs = []
            qcrs = []
            pcut = [0., 0.]
            ecrs = []
            ncrs = []

            for ind in range(1, ncrb + 1):
                ecrs.append(float(position['cree' + str(ind).zfill(2)][0].v))
                ncrs.append(float(position['cren' + str(ind).zfill(2)][0].v))

            for ind in range(1, ncrb + 2):
                fcrs.append(float(position['cref' + str(ind).zfill(2)][0].v))
            for ind in range(1, ncrb + 1):
                qcrs.append(float(position['creq' + str(ind).zfill(2)][0].v))
            pcut[:] = [position['crep01'][0].v, position['crep02'][0].v]

            fig2, exit_code = crs_plot_main_fpq(var_names, var_array, plot_var, fcrs, qcrs, pcut,
                                                field_max, time, coords, marker=marker_l[marker_index], clean_plot=options.clean_plot)

        if (exit_code is not True):
            if ((plot_layer is True) or (plot_ovlp is True)):
                # plot line (layer) if cell not empty WARNING - for now only works with mcrwind
                line = s1.plot([dom_l[avail_dim[0]], dom_r[avail_dim[0]]], [
                               coords[2], coords[2]], color="white")
            else:
                # plot point if cell not empty
                point = s1.plot(coords[avail_dim[0]], coords[avail_dim[1]],
                                marker=marker_l[marker_index], color="red")

            if (plot_layer is True):  # Mark averaged level
                yt_data_plot.annotate_line([coords[0], dom_l[avail_dim[0]], coords[2]], [
                                           coords[0], dom_r[avail_dim[0]], coords[2]], plot_args={'color': 'white', "lw": 10.0})
            else:
                if (not options.annotate_rect and marker_index > 0):
                    yt_data_plot.annotate_marker(coords, marker=marker_l[marker_index - 1], plot_args={
                                                 'color': 'red', 's': m_size_l[marker_index - 1], "lw": 4.5})  # cumulatively annotate all clicked coordinates
            if (options.annotate_rect):
                yt_data_plot.annotate_line(
                    lb, lt, plot_args={'color': 'white', "lw": 2.0})
                yt_data_plot.annotate_line(
                    lt, rt, plot_args={'color': 'white', "lw": 2.0})
                yt_data_plot.annotate_line(
                    rt, rb, plot_args={'color': 'white', "lw": 2.0})
                yt_data_plot.annotate_line(
                    rb, lb, plot_args={'color': 'white', "lw": 2.0})

            marker_index = marker_index + 1
            image_number = image_number + 1
# ------------- saving just the spectrum
            if (save_spectrum):
                extent = fig2.get_window_extent().transformed(s.dpi_scale_trans.inverted())
                spectrum_file_out = str(output_path + '/' + filename_nam + '_' + 'slice_' +
                                        slice_ax + '_' + plot_var + '_' + spc_label + '_spec_%03d.pdf' % image_number)
                # bbox not working in py27 FIXME
                s.savefig(spectrum_file_out, transparent='True',
                          bbox_inches="tight", dpi=150)
                prtinfo(
                    "  --->  Saved plot to: %s. \n\033[44mPress 'q' to quit and save yt.SlicePlot with marked coordinates." % spectrum_file_out)
        else:
            prtwarn("Empty cell - not saving.")

        if (f_run):
            f_run = False

    def plot_with_coords_provided(coords_in):
        mark_plot_save(coords_in)

    if (user_coords_provided):
        for cindex in range(len(user_coords)):
            prtinfo("Provided coordintates %s (clickable map will not be shown), processing image." % str(
                user_coords[cindex][:]))
            plot_with_coords_provided(user_coords[cindex][:])
    else:
        prtinfo(
            "\033[44mClick LMB on the colormap to display spectrum ('q' to exit)")
        cid = s.canvas.mpl_connect('button_press_event', read_click_and_plot)
        plt.show()

    text_coords = [0., 0., 0.]
    text_coords[dim_map.get(slice_ax)] = slice_coord
    text_coords[avail_dim[0]] = dom_l[avail_dim[0]]
    text_coords[avail_dim[1]] = dom_l[avail_dim[1]]
    text_coords = [item * 0.9 for item in text_coords]

    if (user_annot_time):
        if (user_draw_timebox is True):
            yt_data_plot.annotate_text(text_coords, 'T = {:0.2f} Myr'.format(float(t.in_units('Myr'))), text_args={
                                       'fontsize': options.fontsize, 'color': 'white', 'alpha': '0.0'}, inset_box_args={'boxstyle': 'round', 'pad': 0.2, 'alpha': 0.8})
        else:
            prtinfo("Not marking line on yt.plot (user_draw_timebox = %s)" %
                    (user_draw_timebox))
            yt_data_plot.annotate_title(
                'T = {:0.2f} Myr'.format(float(t.in_units('Myr'))))

    if (options.no_cbar):
        yt_data_plot.hide_colorbar()
    if (options.no_xlabels and options.no_ylabels):
        yt_data_plot.hide_axes()

    # save image (spectrum already saved) when finished.
    yt_plot_file_out = str(output_path + "/" + filename_nam + '_' +
                           plot_field + '_sliceplot_' + slice_ax + '.pdf')
    yt_data_plot.save(yt_plot_file_out)

    if (not user_coords_provided):
        s.canvas.mpl_disconnect(cid)
