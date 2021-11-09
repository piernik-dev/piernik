#!/usr/bin/python
from pylab import zeros, sqrt, size
import matplotlib.pyplot as plt
from numpy import log10, log, pi, asfarray, array, linspace, sign, around
import h5py
import os
import sys
import crs_pf
import matplotlib.lines as mlines
import matplotlib.markers as mmark
from colored_io import prtinfo, prtwarn, read_var, die

# ------ default values of parameters --------
# TODO - treat these as default, but try to read from problem par
e_small = 1.e-6
eps = 1.0e-15
ncre = 45
p_min_fix = 1.0
p_max_fix = 1.0e6
cre_eff = 0.01
q_big = 30.
f_init = 0.00235  # 0.019# 0.025 # 1.0
q_init = 4.1

q_eps = 0.001
arr_dim_q = 1000
helper_arr_dim = int(arr_dim_q / 20)

c = 1.0  # PSM -> 0.3066067E+06, SI -> 0.2997925E+09

first_run = True
got_q_tabs = False
q_explicit = False
interpolate_cutoffs = True
highlighted = False
plotted_init_slope = False
verbosity_1 = True
verbosity_2 = False

global par_visible_gridx, par_visible_gridy, par_vis_all_borders, par_visible_title, par_simple_title, par_alpha, par_plot_legend, \
    par_plot_e_small, par_plot_color, par_plot_width, par_fixed_dims, i_plot, xkcd_colors, use_color_list, tightened, par_legend_loc
xkcd_colors = ['xkcd:blue', 'xkcd:darkblue',
               'xkcd:red', 'xkcd:crimson', 'xkcd:green']
xkcd_colorsh = ['xkcd:azure', 'xkcd:blue',
                'xkcd:coral', 'xkcd:orangered', 'xkcd:yellowgreen']
plot_colors = ["crimson", "xkcd:azure", "green",
               "xkcd:purple", "xkcd:orange", "indigo", "chartreuse"]

linestyles = ["solid", (0, (2, 0.75)), (0, (5, 1)), (0, (3, 1, 1, 1))]

i_plot = 0
colors = xkcd_colors
par_plot_color = colors[1]
par_plot_linestyle = linestyles[0]
use_color_list = True
par_visible_gridx = False
par_visible_gridy = True
par_vis_all_borders = False
par_visible_title = False
par_simple_title = True
par_alpha = 0.85
par_plot_width = 2.0
par_fixed_dims = True
par_plot_legend = True
par_plot_e3 = False
par_plotted_src = False
par_plot_init_slope = True
# (0.150, 0.925)  # (0.220+index_t*0.24,0.85)#(0.350,0.965- 0.0525*index_t*1.5)# 0.867 (0.490,0.925 )
par_legend_loc = (-2, -2)
default_legend_loc = 1
par_test_name = ""
highlight_bins = []  # [1,3,5,7,9,11,13,15] #[10,16]
par_highlight_bins_rect = False
tightened = False
hide_axes = False

fontsize_axlabels = 18
fontsize_legend = 10


def set_plot_color(plot_color, index, color_list):
    global use_color_list
    if use_color_list:
        plot_color = color_list[index % len(color_list)]
    else:
        plot_color = par_plot_color
    return plot_color


def nq2f(n, q, p_l, p_r):
    if p_r > 0.0 and p_l > 0:
        pr_by_pl = p_r / p_l
        nq2f = n / (4 * pi * p_l**3)
        if abs(q - 3.0) > eps:
            nq2f = nq2f * (3.0 - q) / ((pr_by_pl)**(3.0 - q) - 1.0)
        else:
            nq2f = nq2f / log(p_r / p_l)
    return nq2f
# 1D Newton-Raphson algorithm (to find q):
#


def nr_get_q(q_start, alpha, p_ratio, exit_code):
    iter_limit = 30
    tol_f = 1.0e-9
    x = q_start
    df = 1.0
    for i in range(iter_limit):
        if abs(x) >= q_big:
            x = (x / abs(x)) * q_big
            exit_code = True
            break
        dx = min(x * 1e-3, 10e-2)
        dx = sign(dx) * max(abs(dx), 1.0e-10)
        df = 0.5 * (fun(x + dx, alpha, p_ratio) -
                    fun(x - dx, alpha, p_ratio)) / dx
        delta = -fun(x, alpha, p_ratio) / df
        if abs(delta) <= tol_f:
            exit_code = False
            return x, exit_code
        else:
            x = x + delta
    nr_get_q = x

    return nr_get_q, exit_code

# function used to find q: ----------------------


def fun(x, alpha, p_ratio):
    if abs(x - 3.0) < q_eps:
        fun = -alpha + (-1.0 + p_ratio) / log(p_ratio)
    elif abs(x - 4.0) < q_eps:
        fun = -alpha + p_ratio * log(p_ratio) / (p_ratio - 1.0)
    else:
        fun = -alpha + ((3.0 - x) / (4.0 - x)) * \
            ((p_ratio**(4.0 - x) - 1.0) / (p_ratio**(3.0 - x) - 1.0))
    return fun


def prepare_q_tabs():
    global alpha_tab_q, q_grid, q_big, q_space
    q_grid = zeros(arr_dim_q)  # for later interpolation
    q_space = zeros(helper_arr_dim)  # for start values

    q_grid[:] = q_big
    q_grid[int(arr_dim_q / 2):] = -q_big

    def ln_eval_array_val(i, arr_min, arr_max, min_i, max_i):
        b = (log(float(max_i)) - log(float(min_i))) / (arr_max - arr_min)
        ln_eval_array_val = (
            arr_min - log(float(min_i)) / b) + log(float(i)) / b
        return ln_eval_array_val

    for i in range(1, int(0.5 * helper_arr_dim)):
        q_space[i - 1] = ln_eval_array_val(i, q_big,
                                           float(0.05), 1, int(0.5 * helper_arr_dim - 1))

    for i in range(0, int(0.5 * helper_arr_dim) + 1):
        q_space[int(0.5 * helper_arr_dim) + i - 1] = - \
            q_space[int(0.5 * helper_arr_dim) - i]

    a_max_q = (1.0 + 0.1) * p_fix_ratio
    a_min_q = 1.00000005
    alpha_tab_q = zeros(arr_dim_q)
    alpha_tab_q[:] = a_min_q

    j = min(arr_dim_q - int(arr_dim_q / (arr_dim_q / 100.)), arr_dim_q - 1)
    while (q_grid[j] <= (-q_big) and (q_grid[arr_dim_q - 1] <= (-q_big))):
        a_max_q = a_max_q - a_max_q * 0.005
        for i in range(0, arr_dim_q):
            alpha_tab_q[i] = a_min_q * \
                10.0**((log10(a_max_q / a_min_q)) /
                       float(arr_dim_q) * float(i))
        # computing q_grid takes so little time, that saving the grid is not necessary.
        fill_q_grid()
    return


def fill_q_grid():
    global q_grid, alpha_tab_q, p_fix_ratio, arr_dim_q, helper_arr_dim, q_space
    previous_solution = q_grid[int(len(q_grid) / 2)]
    exit_code = True
    x = previous_solution
    for i in range(1, arr_dim_q, 1):
        x, exit_code = nr_get_q(
            previous_solution, alpha_tab_q[i], p_fix_ratio, exit_code)
        if exit_code is True:
            for j in range(1, helper_arr_dim, 1):
                x = q_space[j]
                x, exit_code = nr_get_q(
                    x, alpha_tab_q[i], p_fix_ratio, exit_code)
                if exit_code is False:
                    q_grid[i] = x
                    prev_solution = x
        else:  # exit_code == false
            q_grid[i] = x
            previous_solution = x
    return


def interpolate_q(alpha):
    global arr_dim_q, alpha_tab_q, q_grid
    index = int((log10(alpha / alpha_tab_q[0]) / log10(
        alpha_tab_q[-1] / alpha_tab_q[0])) * (arr_dim_q - 1))  # + 1
    if (index < 0 or index > arr_dim_q - 1):
        index = max(0, min(arr_dim_q - 1, index))
        q_out = q_grid[index]
    else:
        index2 = index + 1
        q_out = q_grid[index] + (alpha - alpha_tab_q[index]) * (
            q_grid[index] - q_grid[index2]) / (alpha_tab_q[index] - alpha_tab_q[index2])

    return q_out
# plot data ------------------------------------


def plot_data(plot_var, pl, pr, fl, fr, q, time, location, i_lo_cut, i_up_cut):
    global first_run, e_small, i_plot, par_plot_color, par_plot_linestyle, s, clean_plot
    global plot_p_min, plot_p_max, plot_var_min, plot_var_max, use_color_list, i_plot, handle_list, tightened, highlighted, plotted_init_slope

    f_lo_cut = fl[0]
    f_up_cut = fr[-1]
    p_lo_cut = pl[0]
    p_up_cut = pr[-1]

    if plot_var == 'f':
        plot_var_l = fl
        plot_var_r = fr
        plot_var_lo_cut = f_lo_cut
        plot_var_up_cut = f_up_cut
    elif plot_var == 'n':
        plot_var_l = 4 * pi * fl * pl**2
        plot_var_r = 4 * pi * fr * pr**2
        plot_var_lo_cut = 4 * pi * f_lo_cut * p_lo_cut**2
        plot_var_up_cut = 4 * pi * f_up_cut * p_up_cut**2
    elif plot_var == 'e':
        plot_var_l = 4 * pi * c**2 * fl * pl**3
        plot_var_r = 4 * pi * c**2 * fr * pr**3
        plot_var_lo_cut = 4 * pi * c**2 * f_lo_cut * p_lo_cut**3
        plot_var_up_cut = 4 * pi * c**2 * f_up_cut * p_up_cut**3
    if (first_run):
        s = plt.subplot(122)

    if clean_plot:
        s.cla()

    s.set_xscale('log')
    s.set_yscale('log')

    plt.xlabel('$p/m_e c$', labelpad=0.2, fontsize=fontsize_axlabels)
    plt.ylabel('d$' + plot_var + ' / $d$p$',
               fontsize=fontsize_axlabels, labelpad=-0.)
    plt.tick_params(axis='both', which='major', labelsize=fontsize_axlabels)

    if first_run:
        plot_p_min = p_lo_cut
        plot_p_max = p_up_cut

        handle_list = []

        if (plot_var == "e"):
            plot_var_min = 0.1 * e_small
        elif (plot_var == "f"):
            plot_var_min = 0.1 * e_small / (4 * pi * (c ** 2) * p_max_fix ** 3)
        elif (plot_var == "n"):
            plot_var_min = 0.1 * e_small / (c * p_max_fix)

    if par_fixed_dims:  # overwrite
        if (plot_var != "e"):
            plt.ylim(9.5e-13, 1.e-3)
            plt.xlim(p_fix[1], p_fix[-2] * 0.5)
        else:
            plt.ylim(e_small, 1.e-2)
            plt.xlim(p_fix[1], p_fix[-2] * 0.5)

    if (par_plot_e3):
        plt.ylim(10. * plot_var_min, 10. *
                 max(plot_var_r) * max(pr)**3)  # override

    if (par_vis_all_borders):
        plt.grid()
    else:
        s.spines['top'].set_visible(False)
        s.spines['right'].set_visible(False)
        s.spines['bottom'].set_linewidth(1.5)
        s.spines['left'].set_linewidth(1.5)

    if (par_visible_gridy):
        plt.grid(True, 'major', 'y', ls='--', lw=.5, c='k', alpha=.3)
    if (par_visible_gridx):
        plt.grid(True, 'major', 'x', ls='--', lw=.5, c='k', alpha=.3)

# plot floor value
    p_range = linspace(s.get_xlim()[0], s.get_xlim()[1])
    e_smalls = zeros(len(p_range))
    e_smalls[:] = e_small
    if (plot_var == "e"):
        plt.plot(p_range, e_smalls, color="xkcd:azure", label="$e_{small}$")
    elif(plot_var == "n"):
        plt.plot(p_range, e_small / (c * p_range),
                 color="xkcd:azure", label="$n_{small}$")

    par_plot_color = set_plot_color(par_plot_color, i_plot, colors)
    # par_plot_linestyle = set_plot_color(par_plot_linestyle, i_plot, linestyles)      ### WARNING temporary trick

    spectrum_label = ("d$%s$(p)/d$p$ %s, \n[%3.1f, %3.1f, %3.1f] kpc " % (
        plot_var, par_test_name, location[0] / 1000., location[1] / 1000., location[2] / 1000.))
    spectrum_label = ("d$%s$(p)/d$p$ [%3.1f, %3.1f, %3.1f] kpc " % (
        plot_var, location[0] / 1000., location[1] / 1000., location[2] / 1000.))
    # spectrum_label  = ("d$%s$/d$p$, %s (  )" % (plot_var, par_test_name) ) #
    # spectrum_label  = (" %s (z=%3.1fkpc)" % ( par_test_name , location[2]/1000.) )

    for i in range(0, size(fr)):
        if (par_plot_e3):  # multiply times gamma**3
            plt.plot([pl[i], pr[i]], [(pl[i]**3) * plot_var_l[i], (pr[i]**3) *
                     plot_var_r[i]], lw=par_plot_width, color=par_plot_color, alpha=par_alpha)
            plt.plot([pl[i], pl[i]], [plot_var_min, (pl[i]**3) * plot_var_l[i]],
                     lw=par_plot_width, color=par_plot_color, alpha=par_alpha)
            plt.plot([pr[i], pr[i]], [plot_var_min, (pr[i]**3) * plot_var_r[i]],
                     lw=par_plot_width, color=par_plot_color, alpha=par_alpha)
        else:
            plt.plot([pl[i], pr[i]], [plot_var_l[i], plot_var_r[i]], lw=2 * par_plot_width,
                     solid_capstyle='round', color=par_plot_color, alpha=par_alpha, linestyle=par_plot_linestyle)
            plt.plot([pl[i], pl[i]], [plot_var_r[i - 1], plot_var_l[i]], lw=2 * par_plot_width,
                     solid_capstyle='round', color=par_plot_color, alpha=par_alpha, linestyle=par_plot_linestyle)
            plt.plot([pl[i], pl[i]], [plot_var_min, plot_var_l[i]], lw=par_plot_width,
                     solid_capstyle='round', color="xkcd:gray", alpha=par_alpha * 0.2)
            plt.plot([pr[i], pr[i]], [plot_var_min, plot_var_r[i]], lw=par_plot_width,
                     solid_capstyle='round', color="xkcd:gray", alpha=par_alpha * 0.2)
    if (not par_plot_e3):
        plt.plot([pr[size(fr) - 1], pr[size(fr) - 1]], [plot_var_r[size(fr) - 1], plot_var_min], lw=2 *
                 par_plot_width, solid_capstyle='round', color=par_plot_color, alpha=par_alpha)  # rightmost edge
    spectrum = mlines.Line2D([], [], color=par_plot_color, solid_capstyle='round',
                             lw=par_plot_width, alpha=par_alpha, linestyle=par_plot_linestyle, label=spectrum_label)

    if (not highlighted):
        if (len(highlight_bins) > 0):
            par_plot_color = set_plot_color(
                par_plot_color, i_plot, xkcd_colorsh)
            for ind in highlight_bins:
                i = ind
                i1 = i + 1
                plt.fill([p_fix[i], p_fix[i1], p_fix[i1], p_fix[i]], [
                         e_small, e_small, 10., 10.], color="mediumseagreen", alpha=0.20)
            if (not (clean_plot is True)):
                highlighted = True

    if ((par_plot_init_slope is True) and (plotted_init_slope is False)):
        if (plot_var == 'n'):
            init_spec = plt.plot(p_range, (1.0 + 2.e-1) * f_init * 4 * pi * p_range**(-(q_init - 2)), color='gray',
                                 linestyle=":", alpha=0.75, label=r"d$n(p,t)$/d$p$, $E<1/bt$", lw=3)     # initial spectrum
        if (plot_var == 'e'):
            init_spec = plt.plot(p_range, (1.0 + 2.e-1) * f_init * 4 * pi * p_range**(-(q_init - 3)), color='gray',
                                 linestyle=":", alpha=0.45, label=r"d$e(p,t)$/d$p$, $E<1/bt$", lw=3)     # initial spectrum
        if (not (clean_plot is True)):
            # if cleaning plot is on, init slope must be replotted each iteration
            plotted_init_slope = True

    if (par_visible_title):
        if (par_simple_title):
            plt.title("Spectrum of %s(p), Time = %7.3f" % (plot_var, time))
        else:
            plt.title("Spectrum of %s(p) \n Time = %7.3f | location: %7.2f %7.2f %7.2f " % (
                plot_var, time, location[0], location[1], location[2]))
    if (tightened is not True):
        plt.tight_layout()
        tightened = True

    if (par_plot_legend):
        handle_list.append(spectrum)
        plt.legend(handles=handle_list, loc=default_legend_loc if par_legend_loc == (-2, -2)
                   else par_legend_loc, edgecolor="gray", facecolor="white", framealpha=0.65, fontsize=fontsize_legend)

    if (clean_plot):
        handle_list = []

    if (first_run is True):
        first_run = False

    if (hide_axes is True):
        # allows one to hide all axes for the plot, useful for combining mulitple plots.
        s.axis('off')

    return s

# -----------------------------------------------------------------


def detect_active_bins_new(n_in, e_in):
    global num_active_bins

    num_active_bins = 0

    ne_gt_zero = []
    f_gt_zero = []
    q_gt_zero = []
    e_ampl_l = []
    e_ampl_r = []
    active_bins_new = []
    pln = []
    prn = []
    i_lo_tmp = 0
    i_up_tmp = ncre

    for i in range(0, ncre):
        # returns nonzero bin numbers reduced by one compared to CRESP fortran
        if (n_in[i] > 0.0 and e_in[i] > 0.0):
            ne_gt_zero.append(i)
            num_active_bins = num_active_bins + 1
    if num_active_bins == 0:
        return active_bins_new, ncre, 0

    i_lo_tmp = max(ne_gt_zero[0], 0)
    i_up_tmp = min(ne_gt_zero[-1], ncre)
    pln = p_fix[0:ncre - 1]
    prn = p_fix[1:ncre]
    num_active_bins = 0

    for i in range(0, i_up_tmp - i_lo_tmp + 1):
        q_tmp = 3.5
        exit_code = False
        if (q_explicit is True):
            q_tmp, exit_code = nr_get_q(
                q_tmp, e_in[i + i_lo_tmp] / (n_in[i + i_lo_tmp] * c * pln[i + i_lo_tmp]), prn[i + i_lo_tmp] / pln[i + i_lo_tmp], exit_code)
        else:
            # this instruction is duplicated, TODO return it via detect_active_bins_new()
            q_tmp = interpolate_q(
                e_in[i + i_lo_tmp] / (n_in[i + i_lo_tmp] * c * pln[i + i_lo_tmp]))

        q_gt_zero.append(q_tmp)
        f_gt_zero.append(
            nq2f(n_in[i + i_lo_tmp], q_gt_zero[-1], pln[i + i_lo_tmp], prn[i + i_lo_tmp]))
        e_ampl_l.append(4 * pi * c**2 * f_gt_zero[-1] * pln[i + i_lo_tmp]**3)
        e_ampl_r.append(
            4 * pi * c**2 * f_gt_zero[-1] * ((prn[i + i_lo_tmp] / pln[i + i_lo_tmp])**(-q_tmp)) * prn[i + i_lo_tmp] ** 3)
        if ((e_ampl_l[-1] > e_small or e_ampl_r[-1] > e_small) and e_in[i + i_lo_tmp] > e_small):
            active_bins_new.append(ne_gt_zero[i])
            num_active_bins = num_active_bins + 1

    if num_active_bins == 0:
        return active_bins_new, i_lo_tmp, i_up_tmp

    i_lo_tmp = max(active_bins_new[0], 0)
    i_up_tmp = min(active_bins_new[-1], ncre)

    active_bins_new = [i for i in range(i_lo_tmp, i_up_tmp + 1)]
    num_active_bins = len(active_bins_new)

    prtinfo("Active_bins: " + str(active_bins_new))
    return active_bins_new, i_lo_tmp, i_up_tmp

# ------------------------------------------


def crs_initialize(parameter_names, parameter_values):

    try:
        for i in range(len(parameter_names)):
            exec("%s = %s" %
                 (parameter_names[i], parameter_values[i]), globals())
    except:
        die("Exiting: len(names) not equal len(values)")

    global p_fix_ratio, p_fix

    edges = []
    p_fix = []
    edges[0:ncre] = range(0, ncre + 1, 1)
    p_fix[0:ncre] = zeros(ncre + 1)
    log_width = (log10(p_max_fix / p_min_fix)) / (ncre - 2.0)
    for i in range(0, ncre - 1):
        p_fix[i + 1] = p_min_fix * 10.0**(log_width * edges[i])
        p_fix_ratio = 10.0 ** log_width
        p_fix[0] = (sqrt(p_fix[1] * p_fix[2])) / p_fix_ratio
        p_fix[ncre] = (sqrt(p_fix[ncre - 2] * p_fix[ncre - 1])) * p_fix_ratio
        p_fix = asfarray(p_fix)

    p_mid_fix = zeros(ncre)
    p_mid_fix[1:ncre - 1] = sqrt(p_fix[1:ncre - 1] * p_fix[2:ncre])
    p_mid_fix[0] = p_mid_fix[1] / p_fix_ratio
    p_mid_fix[ncre - 1] = p_mid_fix[ncre - 2] * p_fix_ratio
    p_mid_fix = asfarray(p_mid_fix)

    p_fix = tuple(p_fix)
    p_mid_fix = tuple(p_mid_fix)
    global clean_plot
    clean_plot = True


def crs_plot_main(plot_var, ncrs, ecrs, time, location, **kwargs):
    global first_run, got_q_tabs, e_small, p_min_fix, p_max_fix, ncre, cre_eff, i_plot, marker, clean_plot, hide_axes

    marker = kwargs.get("marker", "x")
    clean_plot = kwargs.get("clean_plot", "True")
    hide_axes = kwargs.get("hide_axes", False)

    i_lo = 0
    i_up = ncre
    active_bins = []
    empty_cell = True

    if (not got_q_tabs and not q_explicit):
        prepare_q_tabs()
        got_q_tabs = True

    active_bins, i_lo, i_up = detect_active_bins_new(ncrs, ecrs)
    if (num_active_bins > 1):
        empty_cell = False

    # i_lo = max(i_lo,1) # temporarily do not display the leftmost bin # FIXME

    prtinfo("\033[44mTime = %6.2f |  i_lo = %2d, i_up = %2d %s" % (time, i_lo if not empty_cell else 0,
            i_up if not empty_cell else 0, '(empty cell / failed to construct spectrum)' if empty_cell else ' '))

    if (verbosity_1):   # Display number density and energy density before exiting
        ncrs1e3 = []
        for item in ncrs:
            ncrs1e3.append(float('%1.3e' % item))
        prtinfo("n = " + str(ncrs1e3))
        ecrs1e3 = []
        for item in ecrs:
            ecrs1e3.append(float('%1.3e' % item))
        prtinfo("e = " + str(ecrs1e3))
        enpc1e3 = []
        for i in range(len(ecrs)):
            enpc1e3.append(float('%1.3e' % (ecrs[i] / (ncrs[i] * p_fix[i]))))
        prtinfo("e/(npc) = " + str(enpc1e3))

    if (empty_cell):
        return plt.subplot(122), empty_cell

    exit_code = False
    if interpolate_cutoffs:
        exit_code_lo = True
        pf_ratio_lo = [0., 0.]
        pf_ratio_lo, exit_code_lo = crs_pf.get_interpolated_ratios(
            "lo", ecrs[i_lo] / (ncrs[i_lo] * c * p_fix[i_lo + 1]), ncrs[i_lo], exit_code_lo, verbose=verbosity_2)

        exit_code_up = True
        pf_ratio_up = [0., 0.]
        if (i_up == ncre):
            i_up = i_up - 1
        pf_ratio_up, exit_code_up = crs_pf.get_interpolated_ratios("up", ecrs[i_up] / (
            ncrs[i_up] * c * p_fix[i_up]), ncrs[i_up], exit_code_up, verbose=verbosity_2)

    pln = p_fix[0:ncre - 1]
    prn = p_fix[1:ncre]
    pln = array(pln)
    prn = array(prn)

    if interpolate_cutoffs:
        if exit_code_lo is True:
            if (verbosity_1):
                prtwarn("Failed to extract boundary (lo) p and f from e, n: pf_ratio_lo = %.6f. Assuming p_fix value." %
                        pf_ratio_lo[0])  # p_fix assumed
        else:
            pln[0] = p_fix[i_lo + 1] / pf_ratio_lo[0]
        fl_lo = crs_pf.e_small_2_f(e_small, pln[i_lo])
        fr_lo = fl_lo * pf_ratio_lo[1]

        if exit_code_up is True:
            if (verbosity_1):
                prtwarn("Failed to extract boundary (up) p and f from e, n: pf_ratio_up = %.6f. Assuming p_fix value." %
                        pf_ratio_up[0])  # p_fix assumed
        else:
            prn[i_up] = p_fix[i_up] * pf_ratio_up[0]
        fr_up = crs_pf.e_small_2_f(e_small, prn[i_up])
        fl_up = fr_up / pf_ratio_up[1]

    if (not q_explicit):
        if (verbosity_2):
            prtinfo("Spectral indices q will be interpolated")
    else:
        if (verbosity_2):
            prtinfo("Spectral indices q will be obtained explicitly")

    q_nr = []
    fln = []
    frn = []
    for i in range(0, i_up - i_lo + 1):
        if (q_explicit is True):
            q_tmp = 3.5
            exit_code = False
            # this instruction is duplicated, TODO return it via detect_active_bins_new()
            q_tmp, exit_code = nr_get_q(
                q_tmp, ecrs[i + i_lo] / (ncrs[i + i_lo] * c * pln[i + i_lo]), prn[i + i_lo] / pln[i + i_lo], exit_code)
        else:
            # this instruction is duplicated, TODO return it via detect_active_bins_new()
            q_tmp = interpolate_q(
                ecrs[i + i_lo] / (ncrs[i + i_lo] * c * pln[i + i_lo]))
        q_nr.append(q_tmp)
        fln.append(nq2f(ncrs[i + i_lo], q_nr[-1],
                   pln[i + i_lo], prn[i + i_lo]))

    q_nr = array(q_nr)
    fln = array(fln)
    frn = array(fln)
    frn = frn * (prn[i_lo:i_up + 1] / pln[i_lo:i_up + 1]) ** (-q_nr)
    plot = False

    # retrieve slopes and f values for cutoffs
    fln[0] = fl_lo
    frn[0] = fr_lo
    fln[-1] = fl_up
    frn[-1] = fr_up
    q_nr[0] = -log10(frn[0] / fln[0]) / log10(prn[0] / pln[0])
    q_nr[0] = sign(q_nr[0]) * min(abs(q_nr[0]), q_big)
    q_nr[-1] = -log10(frn[-1] / fln[-1]) / log10(prn[-1] / pln[-1])
    q_nr[-1] = sign(q_nr[-1]) * min(abs(q_nr[-1]), q_big)

    if (verbosity_1):
        prtinfo("q = " + str(around(q_nr, 3)))
        fln1e3 = []
        for item in fln:
            fln1e3.append(float('%1.3e' % item))
        prtinfo("f = " + str(fln1e3))

    if (verbosity_1):
        prtinfo("Cutoff indices obtained (lo, up): %i, %i || momenta (lo, up): %f, %f " % (
            i_lo, i_up, pln[i_lo], prn[i_up]))

    if (verbosity_2):
        dummyCRSfile = open("crs.dat", "a")
        p_all = list(p_fix)    # fixes list binding problem - appending p_fix
        p_all[i_lo] = pln[0]
        p_all[i_up] = prn[-1]
        q_all = zeros(ncre)
        q_all[i_lo:i_up] = q_nr
        f_all = zeros(ncre + 1)
        f_all[i_lo:i_up] = fln
        string = "%15.3e %4.2f %2d %2d %2d" % (time, 0.0, ncre, i_lo, i_up)
        string = string + " " + str(p_all).strip("[").strip("]") + " " + str(
            f_all).strip("[").strip("]") + " " + str(q_all).strip("[").strip("]")
        dummyCRSfile.write("%s " % string.replace(
            "\n", "").replace("nan", "0.0") + "\n")
        dummyCRSfile.close()

    if empty_cell is not True:
        plot = plot_data(plot_var, pln[i_lo:i_up + 1], prn[i_lo:i_up + 1], fln, frn,
                         q_nr, time, location, i_lo, i_up)
        i_plot = i_plot + 1

    return plot, empty_cell


def crs_plot_main_fpq(parameter_names, parameter_values, plot_var, fcrs, qcrs, pcrs, field_max, time, location, **kwargs):
    global first_run, got_q_tabs, e_small, p_min_fix, p_max_fix, ncre, cre_eff, i_plot, marker, par_plotted_src

    marker = kwargs.get("marker", "x")
    i_plot = kwargs.get("i_plot", 0)

    try:
        for i in range(len(parameter_names)):
            exec("%s = %s" %
                 (parameter_names[i], parameter_values[i]), globals())
    except:
        die(" len(names) not equal len(values) at input.")

# TODO -------- do it under *args TODO
    fixed_width = True
# -------------------
    global plot_ymax, p_fix_ratio, p_fix
    plot_ymax = field_max * cre_eff

    if first_run:
        dummyCRSfile = open("crs.dat", "w+")
        dummyCRSfile.close()
        p_fix = []
        edges = []
        edges[0:ncre] = range(0, ncre + 1, 1)
        p_fix[0:ncre] = zeros(ncre + 1)
        log_width = (log10(p_max_fix / p_min_fix)) / (ncre - 2.0)
        for i in range(0, ncre - 1):  # organize p_fix
            p_fix[i + 1] = p_min_fix * 10.0**(log_width * edges[i])
            p_fix_ratio = 10.0 ** log_width
            p_fix[0] = (sqrt(p_fix[1] * p_fix[2])) / p_fix_ratio
            p_fix[ncre] = (
                sqrt(p_fix[ncre - 2] * p_fix[ncre - 1])) * p_fix_ratio
            p_fix = asfarray(p_fix)

    i_lo = 0
    i_up = ncre
    empty_cell = True

    for i in range(ncre):
        if fcrs[i] > 0.0:
            i_lo = i - 2
            empty_cell = False
            break
    for i in range(ncre, 1, -1):
        if fcrs[i] > 0.0:
            i_up = max(0, i + 1)
            break

    i_cor = 0
    if (fcrs[i_lo] == 0.0):
        i_cor = 1
    fln = array(fcrs[i_lo + i_cor:i_up])
    q = array(qcrs[i_lo + i_cor:i_up])
    pln = p_fix[i_lo + i_cor:i_up]
    prn = p_fix[i_lo + 1 + i_cor:i_up + 1]

    pln[0] = pcrs[0]
    prn[-1] = pcrs[-1]
    plot = False  # dummy variable until plot is made

    prtinfo("\033[44mTime = %6.2f |  i_lo = %2d, i_up = %2d %s" % (time, i_lo + i_cor if not empty_cell else 0,
            i_up if not empty_cell else 0, '(empty cell / failed to construct spectrum)' if empty_cell else ' '))
    if (i_lo == ncre or i_up == 0):
        empty_cell = True
        return plot, empty_cell

    frn = fln * (prn / pln) ** (-q)
    prtinfo("Cutoff indices obtained (lo, up): %i, %i || momenta (lo, up): %f, %f " % (
        i_lo + i_cor, i_up, pcrs[0], pcrs[1]))
    if (verbosity_1):
        fcrs1e3 = []
        for item in fcrs:
            # : print ("n = %1.3e" %item)
            fcrs1e3.append(float('%1.3e' % item))
        prtinfo("f = " + str(fcrs1e3))
        prtinfo("q = " + str(around(qcrs, 3)))

    if (empty_cell is False):

        plot = plot_data(plot_var, pln, prn, fln, frn,
                         q, time, location, i_lo, i_up)
        i_plot = i_plot + 1

    return plot, empty_cell
