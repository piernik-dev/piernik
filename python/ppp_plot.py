#!/usr/bin/env python3

import sys
import argparse


t_bias = 10  # I'd prefer to reduce the timers initially to 0, but sometimes slaves start before master


class PPP_Node:
    """A single event and links"""

    parent = None
    proc = None
    label = ""

    def __init__(self, label, time):
        self.label = label
        self.start = []
        self.stop = []
        self._set_time(time)
        self.children = {}  # this flattens the tree a bit ut we can recontruct it from timings anyway, if needed
        self.parent = None

    def _set_time(self, time):
        if time > 0:
            if (len(self.start) != len(self.stop)):
                sys.stderr.write("Warning: orphaned start for '" + self.path() + "' " + str(time) + " vs " + str(self.start[-1]) + "\n")
                self.stop.append(None)
            self.start.append(time)
        elif time < 0:
            if (len(self.start) != len(self.stop) + 1):
                sys.stderr.write("Warning: orphaned stop for '" + self.path() + "' " + str(time) + " vs " + str(self.stop[-1]) + "\n")
                self.start.append(None)
            self.stop.append(-time)
        else:
            sys.stderr.write("Warning: time == 0 for '" + self.path() + "'\n")

    def _add(self, label, time):
        if time > 0:  # create/update child
            if label not in self.children:
                self.children[label] = PPP_Node(label, time)
                self.children[label].parent = self
            else:
                self.children[label]._set_time(time)
            return self.children[label]
        else:
            if label == self.label:
                self._set_time(time)
            else:
                sys.stderr.write("Warning: unfinished for '" + self.path() + "' " + str(time) + "\n")
            return self.parent  # most likely won't recover properly

    def path(self):
        return (str(self.proc) if self.parent is None else self.parent.path()) + "/" + self.label

    def shortpath(self):
        return ("" if self.parent is None else self.parent.shortpath() + "/") + self.label

    def print(self, indent=1):  # need to filter through parent range
        try:
            for _ in range(min(len(self.start), len(self.stop))):
                print("  " * indent + "'" + self.label + "' %.6f %.6f" % (self.start[_], self.stop[_] - self.start[_]))
        except TypeError:
            print("  " * indent + "'" + self.label + "' TypeError: ", self.start, self.stop)
        for i in self.children:
            self.children[i].print(indent=indent + 1)

    def get_all_ev(self):
        evl = [[self.shortpath(), self.start, self.stop]]
        for i in self.children:
            evl += self.children[i].get_all_ev()
        return evl


class PPP_Tree:
    """An event tree for one Piernik process"""

    def __init__(self, name):
        self.thread = name
        self._last = None
        self.root = {}

    def _add(self, label, time):
        if self._last is None:
            self.root[label] = PPP_Node(label, time)
            self.root[label].proc = self.thread
            self._last = self.root[label]
        else:
            self._last = self._last._add(label, time)

    def print(self):
        print("Process: " + str(self.thread))
        for i in self.root:
            self.root[i].print()


class PPP:
    """A collection of event trees from one or many Piernik processes"""

    def __init__(self, name):
        self.name = name
        self.trees = {}

    def print(self):
        print(self.name)
        for i in sorted(self.trees):
            self.trees[i].print()

    def _decode_text(self, fname):
        self.descr = ""
        self.nthr = 0
        self.bigbang = None
        try:
            file = open(fname, 'r')
            for line in file:
                l = line.split()
                if self.bigbang is None:
                    self.bigbang = float(l[1]) - t_bias
                self._add(int(l[0]), line[line.index(l[2]):].strip(), float(l[1]) + self.bigbang * (-1. if float(l[1]) > 0. else 1.))  # proc, label, time
                if int(l[0]) >= self.nthr:
                    self.nthr = int(l[0]) + 1
            self.descr = "'%s' (%d thread%s)" % (fname, self.nthr, "s" if self.nthr > 1 else "")
        except IOError:
            sys.stderr.write("Error: cannot open '" + fname + "'\n")
            exit(4)

    def _add(self, proc, label, time):
        if proc not in self.trees:
            self.trees[proc] = PPP_Tree(proc)
        self.trees[proc]._add(label, time)


class PPPset:
    """A collection of event trees from one or many Piernik runs"""

    def __init__(self, name):
        self.name = name
        self.evt = []

    def decode(self, fnamelist):
        """Let's focus on ASCII decoding for a while. Don't bother with HDF5 unltil we implement this type of PPP dump"""

        for fname in fnamelist:
            self.evt.append(PPP(fname))
            self.evt[-1]._decode_text(fname)

    def print(self, otype):
        if otype == "tree":
            for _ in range(len(self.evt)):
                self.evt[_].print()
        elif otype == "gnu":
            self.print_gnuplot()
        elif otype == "summary":
            print("ARGS ", args)
            from pprint import pprint
            for _ in range(len(self.evt)):
                pprint(self.evt[_].__dict__)

    def print_gnuplot(self):
        self.descr = ""
        ev = {}
        peff = 0
        for f in range(len(self.evt)):
            for p in self.evt[f].trees:
                evlist = []
                for r in self.evt[f].trees[p].root:
                    evlist += self.evt[f].trees[p].root[r].get_all_ev()
                for e in evlist:
                    e_base = e[0].split('/')[-1].split()[0]
                    if e_base not in ev:
                        ev[e_base] = {}
                    if p + peff not in ev[e_base]:
                        ev[e_base][p + peff] = []
                    ev[e_base][p + peff].append(e)
            peff = peff + self.evt[f].nthr + (1 if len(self.evt) > 1 else 0)  # spacing only when we have multiple files
            self.descr = self.evt[f].descr + ("\\n" if len(self.descr) > 0 else "") + self.descr

        print("#!/usr/bin/gnuplot\n$PPPdata << EOD\n")
        i_s = 1
        i_e = 2
        for e in ev:
            print("# __" + e + "__ '" + ev[e][next(iter(ev[e].keys()))][0][0] + "'")
            for p in ev[e]:
                for t in range(len(ev[e][p])):
                    depth = ev[e][p][t][0].count("/")
                    for _ in range(min(len(ev[e][p][t][i_s]), len(ev[e][p][t][i_e]))):
                        try:
                            if (ev[e][p][t][i_e][_] - ev[e][p][t][i_s][_]) > 0.:
                                print("0. 0. %.6f %.6f %.6f %.6f" % (ev[e][p][t][i_s][_] - t_bias, ev[e][p][t][i_e][_] - t_bias,
                                                                     depth + float(p) / peff, depth + float(p + 1) / peff), depth, p)
                        except TypeError:
                            sys.stderr.write("Warning: inclomplete event '", e, "' @", p, " #", str(t), ev[e][p][t])
                    print("")
            print("")
        print("EOD\n\n# Suggested gnuplot commands:\nset key outside horizontal")
        print("set xlabel 'time (walltime seconds)'\nset ylabel 'timer depth + proc/nproc'")
        print('set title "%s"' % self.descr.replace('_', "\\\\_"))
        pline = "plot $PPPdata "
        nocomma = True
        fs = "solid 0.1"
        i = 0
        gp_colors = 8
        for e in ev:
            ind = e.split("/")[-1]
            if nocomma:
                nocomma = False
            else:
                pline += ', "" '
            pline += ' index "__' + ind + '__" with boxxy title "' + ind.replace('_', "\\\\\\_") + '" fs  ' + fs
            i += 1
            if i >= gp_colors:  # gnuplot seems to have only gp_colors colors for solid filling
                npat = int(i / gp_colors)
                if npat >= 3:  # pattern #3 does not show borders
                    npat += 1
                fs = "pat %d" % npat
        print(pline)
        print("pause mouse close")


parser = argparse.ArgumentParser(description="Piernik Precise Profiling Presenter")
parser.add_argument("filename", nargs='+', help="PPP ascii file to process")
parser.add_argument("-o", "--output", nargs=1, help="processed output file")

# parser.add_argument("-e", "--exclude", help="do not show TIMER(s)")  # multiple excudes
# parser.add_argument("--root", help="show only ROOT and its children")
# parser.add_argument("-m", "--maxoutput", nargs=1, default=50000, help="limit output to MAXOUTPUT enries")

pgroup = parser.add_mutually_exclusive_group()
pgroup.add_argument("-g", "--gnuplot", action="store_const", dest="otype", const="gnu", default="gnu", help="gnuplot output (default)")
pgroup.add_argument("-t", "--tree", action="store_const", dest="otype", const="tree", help="tree output")

# pgroup.add_argument("-r", "--reduced_tree", action="store_const", dest="otype", const="rtree", help="reduced tree output")
pgroup.add_argument("-s", "--summary", action="store_const", dest="otype", const="summary", help="short summary")

args = parser.parse_args()

evt = PPPset("Collection")
evt.decode(args.filename)
evt.print(args.otype)


# collect ev_tree[] into ev_summary
# ev_summary.print()
# ev_summary.make_image()  # dump an .SVG
# ev_summary.browse()  # create browsable image with zoom

# for j  in *.ppprofile.ascii ; do for i in `awk '{print $2}' $j |  sort |  uniq` ; do in=`grep "$i *  -" $j | wc -l` ; out=`grep "$i *  [1-9]" $j | wc -l` ; [ $in != $out ] && echo $j $i $in $out ; done |  column -t ; done

# awk '{print $2}' *.ppprofile.ascii |  sort |  uniq -c |  sort -n
