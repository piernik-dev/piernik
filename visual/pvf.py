#!/usr/bin/env python
# -*- coding: utf-8 -*-

import getopt
import h5py
import os
import sys
import plot_compose as pc

cmap = 'viridis'
pcolor = 'default'
gcolor = ''
plotdir = 'frames'
sctype = 'linear'
cu, center = False, [0.0, 0.0, 0.0]
zoom = False, [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]
zmin, zmax = 0.0, 0.0
draw_part = False
draw_data = False
draw_uni, draw_amr = False, False
plotlevels, gridlist = '', ''
dnames = ''
uaxes = ''
nbins = 1
psize = 0
exten = '.png'

print('PIERNIK VISUALIZATION FACILITY')


def print_usage():
    print('Usage: ./pvf.py <file> [options]')
    print('')
    print('Options:')
    print(' -h, \t\t\t--help \t\t\t\tprint this help')
    print('\t\t\t--amr\t\t\t\tcollect all refinement levels of grid to plot [default: True while AMR refinement level structure exists]')
    print(' -a UNIT, \t\t--axes UNIT \t\t\tscale plot axes with UNIT [default: dataset units]')
    print(' -b BINS, \t\t--bins BINS \t\t\tmake a 2D histogram plot using BINS number instead of scattering particles [default: 1, which leads to scattering]')
    print(' -c CX,CY,CZ, \t\t--center CX,CY,CZ \t\tplot cuts across given point coordinates CX, CY, CZ [default: computed domain center]')
    print(' -d VAR[,VAR2], \t--dataset VAR[,VAR2] \t\tspecify one or more datafield(s) to plot [default: print available datafields; all or _all_ to plot all available datafields]')
    print(' -D COLORMAP, \t\t--colormap COLORMAP \t\tuse COLORMAP palette [default: viridis]')
    print(' -e EXTENSION, \t\t--extension EXTENSION \t\tsave plot in file using filename extension EXTENSION [default: png]')
    print(' -g COLOR, \t\t--gridcolor COLOR \t\tshow grids in color COLOR [default: none]')
    print('\t\t\t--grid-list GRID1[,GRID2] \tplot only selected numbered grid blocks [default: all existing blocks]')
    print(' -l LEVEL1[,LEVEL2], \t--level LEVEL1[,LEVEL2] \tplot only requested grid levels [default: 0 for --uniform, all for --amr]')
    print(' -o OUTPUT, \t\t--output OUTPUT \t\tdump plot files into OUTPUT directory [default: frames]')
    print(' -p,\t\t\t--particles\t\t\tscatter particles onto slices [default: switched-off]')
    print(' -P,\t\t\t--particle-color\t\tuse color for particles scattering or colormap for particles histogram plot [default: #1f77b4 (blue) or viridis]')
    print(' -s,\t\t\t--particle-sizes\t\tmarker sizes for scattering particles onto slices [default: switched-off]')
    print(' -t SCALETYPE, \t\t--scale SCALETYPE \t\tdump use SCALETYPE scale type for displaying data (possible values: 0 | linear, 1 | symlin, 2 | log, 3 | symlog) [default: linear]')
    print('\t\t\t--uniform\t\t\treconstruct uniform grid to plot [default: True while no AMR refinement level structure exists]')
    print(' -z ZMIN,ZMAX, \t\t--zlim ZMIN,ZMAX \t\tlimit colorscale to ZMIN and ZMAX [default: computed data maxima symmetrized]')
    print('\t\t\t--zoom XL,XR,YL,YR,ZL,ZR \tset plot axes ranges [default: domain edges]')


def cli_params(argv):
    try:
        opts, args = getopt.getopt(argv, "a:b:c:d:D:e:g:hl:o:pP:s:t:z:", ["help", "amr", "axes=", "bins=", "center=", "colormap=", "dataset=", "extension=", "gridcolor=", "grid-list=", "level=", "output=", "particles", "particle-color=", "particle-sizes=", "scale=", "uniform", "zlim=", "zoom="])
    except getopt.GetoptError:
        print("Unrecognized options: %s \n" % argv)
        print_usage()
        sys.exit(2)
    for opt, arg in opts:
        if opt in ("-h", "--help"):
            print_usage()
            sys.exit()

        elif opt in ("-a", "--axes"):
            global uaxes
            uaxes = str(arg)

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
            global gcolor
            gcolor = str(arg)

        elif opt in ("-l", "--level"):
            global plotlevels
            plotlevels = [int(i) for i in arg.split(',')]

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

        elif opt in ("-s", "--particle-sizes"):
            global psize
            psize = float(arg)

        elif opt in ("-t", "--scale"):
            global sctype
            sctype = str(arg)

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

iw = 1
for word in sys.argv[1:]:
    if word[0] == '-':
        break
    iw += 1

cli_params(sys.argv[iw:])
if pcolor == 'default':
    if nbins > 1:
        pcolor = 'viridis'
    else:
        pcolor = '#1f77b4'

options = zmin, zmax, cmap, pcolor, psize, sctype, cu, center, draw_data, draw_uni, draw_amr, draw_part, nbins, uaxes, zoom, plotlevels, gridlist, gcolor
if not os.path.exists(plotdir):
    os.makedirs(plotdir)


files_list = sys.argv[1:iw]
for pthfilen in files_list:
    print('')
    file_exists = os.path.exists(pthfilen)
    if not file_exists:
        print('The file %s does not exist!' % pthfilen)
        continue
    if pthfilen.split('.')[-1] == 'h5':
        h5f = h5py.File(pthfilen, 'r')
        particles_in_file = 'particle_types' in list(h5f)
    else:
        continue
    if not (draw_data or draw_part) or (draw_data and dnames == '') or (not draw_data and draw_part and not particles_in_file):
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
    prd, prp = '', '',
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
    else:
        if particles_in_file:
            varlist = ['part']
            prp = 'particles only'
        else:
            varlist = []
    if varlist != []:
        print('Going to read ' + prp + prd)

    for var in varlist:
        if (not draw_data or var in list(h5f['field_types'].keys())):
            # output = plotdir+'/'+filen.split('/')[-1].replace('.h5',"_%s.png" % var)
            fnl = filen.split('/')[-1]
            output = plotdir + '/' + '_'.join(fnl.split('_')[:-1]) + '_' + var + '_' + fnl.split('_')[-1].replace('.h5', exten)
            pc.plotcompose(pthfilen, var, output, options)
        else:
            print(var, ' is not available in the file ', pthfilen)

    h5f.close()
