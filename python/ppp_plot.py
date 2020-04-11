#!/usr/bin/env python3

import sys
import magic
import argparse


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
        self.children = {}
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

    def print(self, bigbang, indent=1):
        try:
            for _ in range(min(len(self.start), len(self.stop))):
                print("  " * indent + "'" + self.label + "' %.6f %.6f" % (self.start[_] - bigbang, self.stop[_] - self.start[_]))
        except TypeError:
            print("  " * indent + "'" + self.label + "' TypeError: ", self.start, self.stop)
        for i in self.children:
            self.children[i].print(bigbang, indent=indent + 1)

    def get_all_ev(self):
        evl = [[self.shortpath(), self.start, self.stop]]
        for i in self.children:
            evl += self.children[i].get_all_ev()
        return evl


class PPP_Tree:
    """An event tree for one Piernik process"""

    def __init__(self, name):
        self.label = name
        self._last = None
        self.root = {}

    def _add(self, label, time):
        if self._last is None:
            self.root[label] = PPP_Node(label, time)
            self.root[label].proc = self.label
            self._last = self.root[label]
        else:
            self._last = self._last._add(label, time)

    def get_bigbang(self):
        bigbang = None
        for i in self.root:
            bb = min(self.root[i].start)
            bigbang = bb if bigbang is None else min(bb, bigbang)
        return bigbang

    def print(self, bigbang):
        print("Process: " + str(self.label))
        for i in self.root:
            self.root[i].print(bigbang)


class PPP:
    """A collection of event trees from one or many Piernik processes"""

    def __init__(self, name):
        self.name = name
        self.trees = {}

    def get_bigbang_from_master(self):
        try:
            bigbang = self.trees[0].get_bigbang()
        except KeyError:
            bigbang = 0.
        return bigbang

    def print(self):
        print(self.name)
        bigbang = self.get_bigbang_from_master()
        for i in sorted(self.trees):
            self.trees[i].print(bigbang)

    def print_gnuplot(self):
        ev = {}
        np = 0  # number of processes
        bigbang = self.get_bigbang_from_master()
        for p in self.trees:
            if p > np:
                np = p
            evlist = []
            for r in self.trees[p].root:
                evlist += self.trees[p].root[r].get_all_ev()
            for e in evlist:
                e_base = e[0].split('/')[-1].split()[0]
                if e_base not in ev:
                    ev[e_base] = {}
                if p not in ev[e_base]:
                    ev[e_base][p] = []
                ev[e_base][p].append(e)

        print("#!/usr/bin/gnuplot\n$PPPdata << EOD\n")
        i_s = 1
        i_e = 2
        for e in ev:
            print("# __" + e + "__ '" + ev[e][0][0][0] + "'")
            for p in ev[e]:
                for t in range(len(ev[e][p])):
                    depth = ev[e][p][t][0].count("/")
                    for _ in range(min(len(ev[e][p][t][i_s]), len(ev[e][p][t][i_e]))):
                        try:
                            if (ev[e][p][t][i_e][_] - ev[e][p][t][i_s][_]) > 0.:
                                print("0. 0. %.6f %.6f %.6f %.6f" % (ev[e][p][t][i_s][_] - bigbang, ev[e][p][t][i_e][_] - bigbang,
                                                                     depth + float(p) / (np + 1), depth + float(p + 1) / (np + 1)), depth, p)
                        except TypeError:
                            sys.stderr.write("Warning: inclomplete event '", e, "' @", p, " #", str(t), ev[e][p][t])
                    print("")
            print("")
        print("EOD\n\n# Suggested gnuplot commands:\nset key outside horizontal")
        print("set xlabel 'time (walltime seconds)'\nset ylabel 'timer depth + proc/nproc'\nset title 'nproc = %d'" % (np + 1))
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

    def decode(self, fname):
        """Try to extract the data from a given file"""
        try:
            ftype = magic.detect_from_filename(fname).mime_type
        except:
            sys.stderr.write("Error: cannot determine file type or file cannot be read\n")
            exit(2)
        if (ftype == "text/plain"):
            self._decode_text(fname)
        elif (ftype == "application/x-hdf"):
            self._decode_hdf5(fname)
        else:
            sys.stderr.write("Error: don't know what to do with " + ftype + "\n")
            exit(3)

    def _decode_text(self, fname):
        self.name += " of text events"
        try:
            file = open(fname, 'r')
            for line in file:
                l = line.split()
                self._add(int(l[0]), line[line.index(l[1]):line.index(l[-1])].strip(), float(l[-1]))  # proc, label, time
        except IOError:
            sys.stderr.write("Error: cannot open '" + fname + "'\n")
            exit(4)

    def _decode_hdf5(self, fname):
        self.name += " of HDF5 events"
        sys.stderr.write("Error: don't know how to decode HDF5 yet\n")
        exit(5)

    def _add(self, proc, label, time):
        if proc not in self.trees:
            self.trees[proc] = PPP_Tree(proc)
        self.trees[proc]._add(label, time)


parser = argparse.ArgumentParser(description="Piernik Precise Profiling Presenter")
parser.add_argument("filename", nargs=1, help="PPP ascii file to process")
parser.add_argument("-o", "--output", nargs=1, help="processed output file")

# parser.add_argument("-e", "--exclude", help="do not show TIMER(s)")  # multiple excudes
# parser.add_argument("--root", help="show only ROOT and its children")
# parser.add_argument("-m", "--maxoutput", nargs=1, default=50000, help="limit output to MAXOUTPUT enries")

pgroup = parser.add_mutually_exclusive_group()
pgroup.add_argument("-g", "--gnuplot", action="store_const", dest="otype", const="gnu", default="gnu", help="gnuplot output (default)")
pgroup.add_argument("-t", "--tree", action="store_const", dest="otype", const="tree", help="tree output")

# pgroup.add_argument("-r", "--reduced_tree", action="store_const", dest="otype", const="rtree", help="reduced tree output")
# pgroup.add_argument("-s", "--summary", action="store_const", dest="otype", const="summary", help="short summary")

args = parser.parse_args()

evt = PPP("Collection")
evt.decode(args.filename[0])
if args.otype == "tree":
    evt.print()
elif args.otype == "gnu":
    evt.print_gnuplot()


# collect ev_tree[] into ev_summary
# ev_summary.print()
# ev_summary.make_image()  # dump an .SVG
# ev_summary.browse()  # create browsable image with zoom

# for j  in *.ppprofile.ascii ; do for i in `awk '{print $2}' $j |  sort |  uniq` ; do in=`grep "$i *  -" $j | wc -l` ; out=`grep "$i *  [1-9]" $j | wc -l` ; [ $in != $out ] && echo $j $i $in $out ; done |  column -t ; done

# awk '{print $2}' *.ppprofile.ascii |  sort |  uniq -c |  sort -n
