#!/usr/bin/env python
# -*- coding: utf-8 -*-

import getopt
import h5py
import os
import sys
import plot_compose as pc
import pvf_settings as ps

exten = ps.f_exten
plotdir = ps.f_plotdir
cmap = ps.plot2d_colormap
sctype = ps.plot2d_sctype
psize = ps.particles_size
pcolor = 'default'
gcolor = ''
linstyl = ''
axcuts = 'default'
cu, center = False, [0.0, 0.0, 0.0]
zoom = False, [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]
zmin, zmax = 0.0, 0.0
draw_part = False
draw_data = False
draw_grid = False
draw_uni, draw_amr = False, False
plotlevels, gridlist = '', ''
dnames = ''
uaxes = ''
nbins = 1
player = True, '0', '0', '0'

print('PIERNIK VISUALIZATION FACILITY')


def print_usage():
    print('Usage: ./pvf.py <HDF5 files> [options]')
    print('')
    print('Usage for grid structure:  ./pvf.py <HDF5 files> -g COLORS [options]')
    print('Usage for grid datafields: ./pvf.py <HDF5 files> -d VARS [options]')
    print('Usage for particles:       ./pvf.py <HDF5 files> -p [options]')
    print('')
    print('Options:')
    print(' -h, \t\t\t--help \t\t\t\tprint this help')
    print(' -a CUT1[,CUT2], \t--axes CUT1[,CUT2] \t\tselect plot cuts from the following: 1x, 1y, 1z, xy, xz, yz, 1D, 2D [default: all 2D cuts, otherwise all 1D]')
    print('\t\t\t--amr\t\t\t\tcollect all refinement levels of grid to plot [default: True while AMR refinement level structure exists]')
    print(' -b BINS, \t\t--bins BINS \t\t\tmake a 2D histogram plot using BINS number instead of scattering particles [default: 1, which leads to scattering]')
    print(' -c CX,CY,CZ, \t\t--center CX,CY,CZ \t\tplot cuts across given point coordinates CX, CY, CZ [default: computed domain center]')
    print(' -d VAR[,VAR2], \t--dataset VAR[,VAR2] \t\tspecify one or more datafield(s) to plot [default: print available datafields; all or _all_ to plot all available datafields]')
    print(' -D COLORMAP, \t\t--colormap COLORMAP \t\tuse COLORMAP palette [default: %s]' % ps.plot2d_colormap)
    print(' -e EXTENSION, \t\t--extension EXTENSION \t\tsave plot in file using filename extension EXTENSION [default: %s]' % ps.f_exten[1:])
    print(' -g COLOR, \t\t--gridcolor COLOR \t\tshow grids in color COLOR; possible list of colors for different grid refinement levels [default: none]')
    print('\t\t\t--grid-list GRID1[,GRID2] \tplot only selected numbered grid blocks [default: all existing blocks]')
    print(' -l LEVEL1[,LEVEL2], \t--level LEVEL1[,LEVEL2] \tplot only requested grid levels [default: 0 for --uniform, all for --amr]')
    print('\t\t\t--linestyle STYLELIST \t\tline styles list for different refinement levels in 1D plots [default: %s]' % ps.plot1d_linestyle)
    print(' -o OUTPUT, \t\t--output OUTPUT \t\tdump plot files into OUTPUT directory [default: %s]' % ps.f_plotdir)
    print(' -p,\t\t\t--particles\t\t\tscatter particles onto slices [default: switched-off]')
    print(' -P,\t\t\t--particle-color\t\tuse color for particles scattering or colormap for particles histogram plot [default: %s or %s]' % (ps.particles_color, ps.hist2d_colormap))
    print(' -r W1[,W2,W3],\t\t--particle-slice W1[,W2,W3]\tread particles from layers +/-W1 around center; uses different widths for different projections if W1,W2,W3 requested [default: all particles]')
    print(' -R W1[,W2,W3],\t\t--particle-space W1[,W2,W3]\tread particles from square +/-W1 around center or cuboid if W1,W2,W3 requested [default: no limits]')
    print(' -s,\t\t\t--particle-sizes\t\tmarker sizes for scattering particles onto slices [default: switched-off]')
    print(' -t SCALETYPE, \t\t--scale SCALETYPE \t\tdump use SCALETYPE scale type for displaying data (possible values: 0 | linear, 1 | symlin, 2 | log, 3 | symlog) [default: %s]' % ps.plot2d_sctype)
    print(' -u UNIT, \t\t--units UNIT \t\t\tscale plot axes with UNIT [default: dataset units]')
    print('\t\t\t--uniform\t\t\treconstruct uniform grid to plot [default: True while no AMR refinement level structure exists]')
    print(' -z ZMIN,ZMAX, \t\t--zlim ZMIN,ZMAX \t\tlimit colorscale to ZMIN and ZMAX [default: computed data maxima symmetrized]')
    print('\t\t\t--zoom XL,XR,YL,YR,ZL,ZR \tset plot axes ranges [default: domain edges]')


def cli_params(argv):
    try:
        opts, args = getopt.getopt(argv, "a:b:c:d:D:e:g:hl:o:pP:r:R:s:t:u:z:", ["help", "amr", "axes=", "bins=", "center=", "colormap=", "dataset=", "extension=", "gridcolor=", "grid-list=", "level=", "linestyle=", "output=", "particles", "particle-color=", "particle-space=", "particle-sizes=", "particle-slice=", "scale=", "uniform", "units=", "zlim=", "zoom="])
    except getopt.GetoptError:
        print("Unrecognized options: %s \n" % argv)
        print_usage()
        sys.exit(2)
    for opt, arg in opts:
        if opt in ("-h", "--help"):
            print_usage()
            sys.exit()

        elif opt in ("-a", "--axes"):
            global axcuts
            axcuts = arg.split(',')

        elif opt in ("-b", "--bins"):
            global nbins
            nbins = int(arg)

        elif opt in ("-c", "--center"):
            global center, cu
            cx, cy, cz = arg.split(',')
            cu, center = True, [float(cx), float(cy), float(cz)]

        elif opt in ("-d", "--dataset"):
            global dnames
            global draw_data
            dnames = str(arg)
            draw_data = True

        elif opt in ("-D", "--colormap"):
            global cmap
            cmap = str(arg)

        elif opt in ("-e", "--extension"):
            global exten
            exten = '.' + str(arg)
            print(exten)

        elif opt in ("-g", "--gridcolor"):
            global gcolor, draw_grid
            gcolor = str(arg)
            draw_grid = True

        elif opt in ("-l", "--level"):
            global plotlevels
            plotlevels = [int(i) for i in arg.split(',')]

        elif opt in ("--linestyle"):
            global linstyl
            linstyl = arg.split(',')

        elif opt in ("-o", "--output"):
            global plotdir
            plotdir = str(arg)
            print('PLOTDIR: ', plotdir)

        elif opt in ("-p", "--particles"):
            global draw_part
            draw_part = True

        elif opt in ("-P", "--particle-color"):
            global pcolor
            pcolor = str(arg)

        elif opt in ("-r", "--particle-slice"):
            global player
            aux = arg.split(',')
            if len(aux) >= 3:
                player = True, aux[0], aux[1], aux[2]
            else:
                player = True, aux[0], aux[0], aux[0]

        elif opt in ("-R", "--particle-space"):
            aux = arg.split(',')
            if len(aux) >= 3:
                player = False, aux[0], aux[1], aux[2]
            else:
                player = False, aux[0], aux[0], aux[0]

        elif opt in ("-s", "--particle-sizes"):
            global psize
            psize = float(arg)

        elif opt in ("-t", "--scale"):
            global sctype
            sctype = str(arg)

        elif opt in ("-u", "--units"):
            global uaxes
            uaxes = str(arg)

        elif opt in ("-z", "--zlim"):
            global zmin, zmax
            zmin, zmax = arg.split(',')
            zmin = float(zmin)
            zmax = float(zmax)
            print("zmin, zmax = ", zmin, zmax)

        elif opt in ("--amr",):
            global draw_amr
            draw_amr = True

        elif opt in ("--grid-list",):
            global gridlist
            gridlist = [int(i) for i in arg.split(',')]

        elif opt in ("--uniform",):
            global draw_uni
            draw_uni = True

        elif opt in ("--zoom",):
            global zoom
            aux = arg.split(',')
            zoom = True, [float(aux[0]), float(aux[2]), float(aux[4])], [float(aux[1]), float(aux[3]), float(aux[5])]
            print("ZOOM: xmin, xmax = ", zoom[1][0], zoom[2][0], 'ymin, ymax = ', zoom[1][1], zoom[2][1], 'zmin, zmax = ', zoom[1][2], zoom[2][2])


if (len(sys.argv) < 2):
    print_usage()
    exit()

files_list = []
optilist = []
for word in sys.argv[1:]:
    if word.split('.')[-1] == 'h5':
        files_list.append(word)
    else:
        optilist.append(word)

if files_list == []:
    print('No h5 files selected. See ./pvf.py -h for help.')

cli_params(optilist)

if pcolor == 'default':
    if nbins > 1:
        pcolor = ps.hist2d_colormap
    else:
        pcolor = ps.particles_color

p1x, p1y, p1z, p2xy, p2xz, p2yz = False, False, False, False, False, False
if 'all' in axcuts:
    p1x, p1y, p1z, p2xy, p2xz, p2yz = True, True, True, True, True, True
if '2D' in axcuts:
    p2xy, p2xz, p2yz = True, True, True
if 'xy' in axcuts:
    p2xy = True
if 'xz' in axcuts:
    p2xz = True
if 'yz' in axcuts:
    p2yz = True
if '1D' in axcuts:
    p1x, p1y, p1z = True, True, True
if '1x' in axcuts:
    p1x = True
if '1y' in axcuts:
    p1y = True
if '1z' in axcuts:
    p1z = True
axc = [p1x, p1y, p1z], [p2yz, p2xz, p2xy]

options = axc, zmin, zmax, cmap, pcolor, player, psize, sctype, cu, center, draw_grid, draw_data, draw_uni, draw_amr, draw_part, nbins, uaxes, zoom, plotlevels, gridlist, gcolor, linstyl
if not os.path.exists(plotdir):
    os.makedirs(plotdir)

for pthfilen in files_list:
    print('')
    file_exists = os.path.exists(pthfilen)
    if not file_exists:
        print('The file %s does not exist!' % pthfilen)
        continue
    h5f = h5py.File(pthfilen, 'r')
    particles_in_file = 'particle_types' in list(h5f)
    if not (draw_data or draw_part or draw_grid) or (draw_data and dnames == '') or (not draw_data and not draw_grid and draw_part and not particles_in_file):
        partincl = ''
        if particles_in_file:
            partincl = 'and particles'
        else:
            if draw_part:
                print('Particles not available in the file!')
        print('Available datafields in the file %s: \n' % pthfilen, list(h5f['field_types'].keys()), partincl)
        h5f.close()
        continue
    filen = pthfilen.split('/')[-1]

    print("Reading file: %s" % pthfilen)
    prd, prp, prg = '', '', ''
    if draw_data:
        if dnames == "_all_" or dnames == "all":
            varlist = h5f['field_types'].keys()
        else:
            varlist = dnames.split(',')
        prd = 'datasets: %s' % varlist
        if draw_part:
            if particles_in_file:
                prp = 'particles and '
            else:
                print('Particles not available in the file!')
    elif particles_in_file:
        varlist = [ps.particles_output]
        prp = 'particles only'
    elif draw_grid:
        varlist = [ps.grid_output]
        prg = 'grid only'
    else:
        varlist = []
    if varlist != []:
        print('Going to read ' + prp + prd + prg)

    for var in varlist:
        if (not draw_data or var in list(h5f['field_types'].keys())):
            # output = plotdir+'/'+filen.split('/')[-1].replace('.h5',"_%s.png" % var)
            fnl = filen.split('/')[-1]
            output = [plotdir + '/' + '_'.join(fnl.split('_')[:-1]) + '_' + var + '_', fnl.split('_')[-1].replace('.h5', exten)]
            pc.plotcompose(pthfilen, var, output, options)
        else:
            print(var, ' is not available in the file ', pthfilen)

    h5f.close()
