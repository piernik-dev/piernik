#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import getopt
import h5py
import os
import sys
import plot_compose as pc
import pvf_settings as ps
import plot_utils as pu

exten = ps.f_exten
plotdir = ps.f_plotdir
cmap = ps.plot2d_colormap
sctype = ps.plot2d_sctype
pstype = ps.hist2d_sctype
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
cmpr, cmprb, cmprf, cmprd, cmprl, cmprn, cmprt = False, False, '', '', '', '', ps.plot2d_comparetype
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
    print(' -h, \t\t\t--help \t\t\t\t\tprint this help')
    print(' -a CUT1[,CUT2], \t--axes CUT1[,CUT2] \t\t\tselect plot cuts from the following: 1x, 1y, 1z, xy, xz, yz, 1D, 2D [default: all 2D cuts, otherwise all 1D]')
    print('\t\t\t--amr\t\t\t\t\tcollect all refinement levels of grid to plot [default: True while AMR refinement level structure exists]')
    print(' -b BINS, \t\t--bins BINS \t\t\t\tmake a 2D histogram plot using BINS number instead of scattering particles [default: 1, which leads to scattering]')
    print(' -c CX,CY,CZ, \t\t--center CX,CY,CZ \t\t\tplot cuts across given point coordinates CX, CY, CZ (give "max" or "min" key to find CX,CY,CZ automatically) [default: computed domain center]')
    print(' -C\t\t\t--compare-adjusted-grids \t\tadjust ranges of grids to compare [default: switched-off]')
    print(' -d VAR[,VAR2], \t--dataset VAR[,VAR2] \t\t\tspecify one or more datafield(s) to plot [default: print available datafields; all or _all_ to plot all available datafields]')
    print(' -D COLORMAP, \t\t--colormap COLORMAP \t\t\tuse COLORMAP palette [default: %s]' % ps.plot2d_colormap)
    print('\t\t\t--compare-datafield VAR \t\tcompare chosen datafields to another VAR')
    print('\t\t\t--compare-type TYPE \t\t\toperation executed as a comparison: 1 - subtraction, 2 - division, 3 - relative error [default: %s]' % ps.plot2d_comparetype)
    print(' -e EXTENSION, \t\t--extension EXTENSION \t\t\tsave plot in file using filename extension EXTENSION [default: %s]' % ps.f_exten[1:])
    print(' -F FILE\t\t--compare-file FILE \t\t\tcompare chosen datafields to another FILE')
    print(' -g COLOR, \t\t--gridcolor COLOR \t\t\tshow grids in color COLOR; possible list of colors for different grid refinement levels [default: none]')
    print('\t\t\t--grid-list GRID1[,GRID2] \t\tplot only selected numbered grid blocks [default: all existing blocks]')
    print(' -l LEVEL1[,LEVEL2], \t--level LEVEL1[,LEVEL2] \t\tplot only requested grid levels [default: all]')
    print(' -L LEVEL1[,LEVEL2]\t--compare-level LEVEL1[,LEVEL2] \tspecify different grid levels to compare accross files [default: the same levels]')
    print('\t\t\t--linestyle STYLELIST \t\t\tline styles list for different refinement levels in 1D plots [default: %s]' % ps.plot1d_linestyle)
    print(' -o OUTPUT, \t\t--output OUTPUT \t\t\tdump plot files into OUTPUT directory [default: %s]' % ps.f_plotdir)
    print(' -p,\t\t\t--particles\t\t\t\tscatter particles onto slices [default: switched-off]')
    print(' -P,\t\t\t--particle-color\t\t\tuse color for particles scattering or colormap for particles histogram plot [default: %s or %s]' % (ps.particles_color, ps.hist2d_colormap))
    print(' -r W1[,W2,W3],\t\t--particle-slice W1[,W2,W3]\t\tread particles from layers +/-W1 around center; uses different widths for different projections if W1,W2,W3 requested [default: all particles]')
    print(' -R W1[,W2,W3],\t\t--particle-space W1[,W2,W3]\t\tread particles from square +/-W1 around center or cuboid if W1,W2,W3 requested [default: no limits]')
    print(' -s,\t\t\t--particle-sizes\t\t\tmarker sizes for scattering particles onto slices [default: switched-off]')
    print(' -T LOG[,MIN,MAX]\t--particle-h2d-scale LOG[,MIN,MAX]\tscaling particle 2D histogram [default: %s %s %s]' % ps.hist2d_sctype)
    print(' -t SCALETYPE, \t\t--scale SCALETYPE \t\t\tdump use SCALETYPE scale type for displaying data (possible values: 0 | linear, 1 | symlin, 2 | log, 3 | symlog) [default: %s]' % ps.plot2d_sctype)
    print(' -u UNIT, \t\t--units UNIT \t\t\t\tscale plot axes with UNIT [default: dataset units]')
    print('\t\t\t--uniform\t\t\t\treconstruct uniform grid to plot [default: True while no AMR refinement level structure exists]')
    print(' -z ZMIN,ZMAX, \t\t--zlim ZMIN,ZMAX \t\t\tlimit colorscale to ZMIN and ZMAX [default: computed data maxima symmetrized]')
    print('\t\t\t--zoom XL,XR,YL,YR,ZL,ZR | LEVEL \tset plot axes ranges or take ranges from LEVEL [default: domain edges | 0]')


def cli_params(argv):
    try:
        opts, args = getopt.getopt(argv, "a:b:c:Cd:D:e:F:g:hl:L:o:pP:r:R:s:t:T:u:z:", ["help", "amr", "axes=", "bins=", "center=", "colormap=", "compare-adjusted-grids", "compare-datafield=", "compare-file=", "compare-level=", "compare-type=", "dataset=", "extension=", "gridcolor=", "grid-list=", "level=", "linestyle=", "output=", "particles", "particle-color=", "particle-h2d-scale=", "particle-space=", "particle-sizes=", "particle-slice=", "scale=", "uniform", "units=", "zlim=", "zoom="])
    except getopt.GetoptError:
        print("Unrecognized options: %s \n" % argv)
        print_usage()
        sys.exit(2)
    for opt, arg in opts:
        if pu.recognize_opt(opt, ("-h", "--help")):
            print_usage()
            sys.exit()

        elif pu.recognize_opt(opt, ("-a", "--axes")):
            global axcuts
            axcuts = arg.split(',')

        elif pu.recognize_opt(opt, ("-b", "--bins")):
            global nbins
            nbins = int(arg)

        elif pu.recognize_opt(opt, ("-c", "--center")):
            global center, cu
            if arg == 'max':
                cu, center = True, [False, True, ]
            elif arg == 'min':
                cu, center = True, [True, False, ]
            else:
                cx, cy, cz = arg.split(',')
                cu, center = True, [float(cx), float(cy), float(cz)]

        elif pu.recognize_opt(opt, ("-d", "--dataset")):
            global dnames
            global draw_data
            dnames = str(arg)
            draw_data = True

        elif pu.recognize_opt(opt, ("-D", "--colormap")):
            global cmap
            cmap = str(arg)

        elif pu.recognize_opt(opt, ("--compare-datafield",)):
            global cmpr, cmprd, cmprn
            cmpr = True
            cmprd = str(arg)
            if cmprn == '':
                cmprn = '_vs'
            cmprn = cmprn + '_' + cmprd

        elif pu.recognize_opt(opt, ("-F", "--compare-file")):
            global cmprf
            cmpr = True
            cmprf = str(arg)
            if not os.path.exists(cmprf):
                print('The compare file %s does not exist!' % cmprf)
                sys.exit()
            if cmprn == '':
                cmprn = '_vs'
            cmprn = cmprn + '_' + cmprf.split('/')[-1]
            cmprn = ''.join(cmprn.split('.')[:-1])

        elif pu.recognize_opt(opt, ("-C", "--compare-adjusted-grids")):
            global cmprb
            cmprb = True

        elif pu.recognize_opt(opt, ("--compare-type",)):
            global cmprt
            cmprt = int(arg)
            if cmprt != 1 and cmprt != 2 and cmprt != 3:
                print('Warning: Unrecognized type of comparison: %s. Taking 1 (subtraction).' % arg)

        elif pu.recognize_opt(opt, ("-e", "--extension")):
            global exten
            exten = '.' + str(arg)
            print(exten)

        elif pu.recognize_opt(opt, ("-g", "--gridcolor")):
            global gcolor, draw_grid
            gcolor = str(arg)
            draw_grid = True

        elif pu.recognize_opt(opt, ("-l", "--level")):
            global plotlevels
            plotlevels = [int(i) for i in arg.split(',')]

        elif pu.recognize_opt(opt, ("-L", "--compare-level")):
            global cmprl
            cmpr = True
            cmprl = [int(i) for i in arg.split(',')]

        elif pu.recognize_opt(opt, ("--linestyle",)):
            global linstyl
            linstyl = arg.split(',')

        elif pu.recognize_opt(opt, ("-o", "--output")):
            global plotdir
            plotdir = str(arg)
            print('PLOTDIR: ', plotdir)

        elif pu.recognize_opt(opt, ("-p", "--particles")):
            global draw_part
            draw_part = True

        elif pu.recognize_opt(opt, ("-P", "--particle-color")):
            global pcolor
            pcolor = str(arg)

        elif pu.recognize_opt(opt, ("-r", "--particle-slice")):
            global player
            aux = arg.split(',')
            if len(aux) >= 3:
                player = True, aux[0], aux[1], aux[2]
            else:
                player = True, aux[0], aux[0], aux[0]

        elif pu.recognize_opt(opt, ("-R", "--particle-space")):
            aux = arg.split(',')
            if len(aux) >= 3:
                player = False, aux[0], aux[1], aux[2]
            else:
                player = False, aux[0], aux[0], aux[0]

        elif pu.recognize_opt(opt, ("-s", "--particle-sizes")):
            global psize
            psize = float(arg)

        elif pu.recognize_opt(opt, ("-T", "--particle-h2d-scale")):
            global pstype
            aux = arg.split(',')
            if len(aux) > 2:
                pstype = aux[0], aux[1], aux[2]
            elif len(aux) > 1:
                pstype = aux[0], aux[1], pstype[2]
            elif len(aux) > 0:
                pstype = aux[0], pstype[1], pstype[2]
            l1 = aux[0].lower() in ("yes", "true", "t", "1")
            l2, l3 = pstype[1:3]
            pstype = l1, l2, l3

        elif pu.recognize_opt(opt, ("-t", "--scale")):
            global sctype
            sctype = str(arg)

        elif pu.recognize_opt(opt, ("-u", "--units")):
            global uaxes
            uaxes = str(arg)

        elif pu.recognize_opt(opt, ("-z", "--zlim")):
            global zmin, zmax
            zmin, zmax = arg.split(',')
            zmin = float(zmin)
            zmax = float(zmax)
            print("zmin, zmax = ", zmin, zmax)

        elif pu.recognize_opt(opt, ("--amr",)):
            global draw_amr
            draw_amr = True

        elif pu.recognize_opt(opt, ("--grid-list",)):
            global gridlist
            gridlist = [int(i) for i in arg.split(',')]

        elif pu.recognize_opt(opt, ("--uniform",)):
            global draw_uni
            draw_uni = True

        elif pu.recognize_opt(opt, ("--zoom",)):
            global zoom
            aux = arg.split(',')
            if len(aux) > 1:
                zoom = True, [float(aux[0]), float(aux[2]), float(aux[4])], [float(aux[1]), float(aux[3]), float(aux[5])]
                print("ZOOM: xmin, xmax = ", zoom[1][0], zoom[2][0], 'ymin, ymax = ', zoom[1][1], zoom[2][1], 'zmin, zmax = ', zoom[1][2], zoom[2][2])
            else:
                zoom = True, int(aux[0])


if (len(sys.argv) < 2):
    print_usage()
    exit()

files_list = []
optilist = []
opt_cmpf = False
for word in sys.argv[1:]:
    if word.split('.')[-1] == 'h5' and not opt_cmpf:
        files_list.append(word)
    else:
        optilist.append(word)
    opt_cmpf = False
    if pu.recognize_opt(word, ("-F", "--compare-file")):
        opt_cmpf = True

if files_list == []:
    print('No h5 files selected. See ./pvf.py -h for help.')

cli_params(optilist)

if pcolor == 'default':
    if nbins > 1:
        pcolor = ps.hist2d_colormap
    else:
        pcolor = ps.particles_color

if pstype[1] == 'auto':
    pstype = pstype[0], None, pstype[2]
if pstype[2] == 'auto':
    pstype = pstype[0], pstype[1], None

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

compare = cmpr, cmprb, cmprf, cmprd, cmprl, cmprt, False

options = axc, zmin, zmax, cmap, pcolor, player, psize, sctype, pstype, cu, center, compare, draw_grid, draw_data, draw_uni, draw_amr, draw_part, nbins, uaxes, zoom, plotlevels, gridlist, gcolor, linstyl
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
        if pu.recognize_opt(dnames, ("_all_", "all")):
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
            output = [plotdir + '/' + '_'.join(fnl.split('_')[:-1]) + '_' + var + cmprn + '_', fnl.split('_')[-1].replace('.h5', exten)]
            pc.plotcompose(pthfilen, var, output, options)
        else:
            print(var, ' is not available in the file ', pthfilen)

    h5f.close()
