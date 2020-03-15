#!/usr/bin/env python3

import sys
import magic


class PPP_Node:
    """A single event and links"""

    start = None
    stop = None

    def __init__(self, label, time):
        self.label = label
        self._set_time(time)
        self.children = {}
        self.parent = None

    def _set_time(self, time):
        if time > 0:
            self._set_start(time)
        elif time < 0:
            self._set_stop(-time)
        else:
            sys.stderr.write("Warning: time == 0 for '" + self.path() + "'\n")

    def _set_start(self, time):
        if self.start is None:
            self.start = time
        else:
            sys.stderr.write("Warning: self.start already set for '" + self.path() + "' " + str(time) + " vs " + str(self.start) + "\n")

    def _set_stop(self, time):  # a bit of spaghetti here (see _set_start)
        if self.stop is None:
            self.stop = time
        else:
            sys.stderr.write("Warning: self.stop already set for '" + self.path() + "' " + str(time) + " vs " + str(self.stop) + "\n")

    def _add(self, label, time):
        if label == self.label:
            self._set_time(time)
            return self.parent if self.stop is not None else self
        else:
            if time > 0:  # create/update child
                if label not in self.children:
                    self.children[label] = PPP_Node(label, time)
                    self.children[label].parent = self
                    return self.children[label]
                else:
                    return self.children[label]._add(label, time)
            else:
                if self.stop is None:
                    sys.stderr.write("Warning: unfinished for '" + self.path() + "' " + str(time) + "\n")
                return self.parent._add(label, time)

    def path(self):
        return (str(self.proc) if self.parent is None else self.parent.path()) + "/" + self.label

    def shortpath(self):
        return ("" if self.parent is None else self.parent.shortpath() + "/") + self.label

    def print(self, bigbang, indent=1):
        print("  " * indent + "'" + self.label + "' %.6f %.6f" % (self.start - bigbang, self.stop - self.start))
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
            bb = self.root[i].start
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

    def print(self):
        print(self.name)
        for i in sorted(self.trees):
            bigbang = self.trees[0].get_bigbang()
            self.trees[i].print(bigbang)

    def print_gnuplot(self):
        ev = {}
        bigbang = self.trees[0].get_bigbang()
        for p in self.trees:
            evlist = []
            for r in self.trees[p].root:
                evlist += self.trees[p].root[r].get_all_ev()
            for e in evlist:
                if e[0] not in ev:
                    ev[e[0]] = {}
                if p not in ev[e[0]]:
                    ev[e[0]][p] = e[1:]
                else:
                    sys.stderr.write("Warning: repeated event '" + e[0] + "' @ " + str(p) + "\n")
        i = 1
        for e in ev:
            print("#%d " % i + "'" + e + "'")
            for p in ev[e]:
                for t in range(len(ev[e][p])):
                    print(e.count("/"), p, "%.6f" % (ev[e][p][t] - bigbang), i)
                print("")
            i += 1
            print("")

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


if (len(sys.argv) != 2):
    sys.stderr.write("Error: need exactly one file to process\n")
    exit(1)


evt = PPP("Collection")
evt.decode(sys.argv[1])
# evt.print()
evt.print_gnuplot()

# collect ev_tree[] into ev_summary
# ev_summary.print()
# ev_summary.make_image()  # dump an .SVG
# ev_summary.browse()  # create browsable image with zoom
