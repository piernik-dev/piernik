#!/usr/bin/env python3

import sys
import argparse


t_bias = 1.e-6  # For global bigbang this may need to be increased to some large values. For per-thread bigbang it can be some small positive value.
included_threads = []  # The list of threads to analyze (same for each file); all enabled when list is empty.
t_min, t_max = -sys.float_info.max, sys.float_info.max  # "Infinite" time interval by default


def niceprint_set(s):
    outstr = ""
    prev = None
    first = None
    for p in sorted(s):
        if first is None:
            first = p
            outstr += str(p)
        else:
            if p - prev > 1:
                if prev != first:
                    outstr += "-" + str(prev)
                first = p
                outstr += "," + str(p)
        prev = p
    if prev != first:
        outstr += "-" + str(prev)
    if len(outstr) < 1:
        outstr = "None"
    return outstr


def plural(n):
    return "" if n == 1 else "s"


class PPP_Node:
    """A single node on tree of events with all event enties and links"""

    parent = None
    proc = None
    label = ""

    def __init__(self, label, time=None):
        self.label = label
        self.start = []
        self.stop = []
        if time is not None:
            self._set_time(time)
        self.children = {}  # this flattens the tree a bit but we can recontruct it from timings anyway, if needed
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
                sys.stderr.write("Warning: unfinished for '" + self.path() +
                                 "' '" + str(label) + "' " + str(time) + "\n")
            return self.parent  # most likely won't recover properly

    def path(self):
        return str(self.proc) if self.parent is None else self.parent.path() + "/" + self.label

    def shortpath(self):
        return "" if self.parent is None else self.parent.shortpath() + "/" + self.label

    def print(self, indent=0):  # need to filter through parent range
        try:
            for _ in range(min(len(self.start), len(self.stop))):
                print("  " * indent + "'" + self.label + "' %.7f %.7f" % (self.start[_] - t_bias, self.stop[_] - self.start[_]))
        except TypeError:
            print("  " * indent + "'" + self.label + "' TypeError: ", self.start, self.stop)
        for i in self.children:
            self.children[i].print(indent=indent + 1)

    def get_all_ev(self):
        evl = [[self.shortpath(), self.start, self.stop, self.own_time, self.child_time]] if len(self.start) > 0 else []
        for i in self.children:
            evl += self.children[i].get_all_ev()
        return evl

    def calculate_time(self):
        self.own_time = 0.
        for i in range(len(self.start)):
            self.own_time += self.stop[i] - self.start[i]
        self.child_time = 0.
        for _ in self.children:
            self.child_time += self.children[_].calculate_time()
        return self.own_time

    def apply_exclusions(self):
        if args.exclude is not None:
            for _ in list(self.children.keys()):
                if self.children[_].label in args.exclude:
                    del self.children[_]

        if args.maxdepth is not None:
            for _ in list(self.children.keys()):
                if self.shortpath().count("/") >= args.maxdepth:
                    del self.children[_]

        for _ in self.children:
            self.children[_].apply_exclusions()

    def time_filter(self):
        for _ in self.children:
            for i in reversed(range(len(self.children[_].start))):
                if (self.children[_].start[i] > t_max) or (self.children[_].stop[i] < t_min):
                    self.children[_].start.pop(i)
                    self.children[_].stop.pop(i)

        for _ in list(self.children.keys()):
            if len(self.children[_].start) < 1:  # clean up empty entries that may occur when using -T option
                del self.children[_]

        for _ in self.children:
            self.children[_].time_filter()

    def root_matches(self):
        if self.label in args.root:
            return [self]
        else:
            match = []
            for _ in self.children:
                match += self.children[_].root_matches()
            return match

    def shift_time(self, dt):
        for i in range(len(self.start)):
            self.start[i] -= dt
            self.stop[i] -= dt
        for _ in self.children:
            self.children[_].shift_time(dt)


class PPP_Tree:
    """An event tree for one Piernik process"""

    def __init__(self, name):
        self.thread = name
        self.root = PPP_Node("/")  # root node is special and does not have time
        self.root.proc = self.thread
        self._last = self.root

    def _add(self, label, time):
        self._last = self._last._add(label, time)

    def print(self):
        print("Process: " + str(self.thread))
        self.root.print()

    def rebase_root(self):
        self.newroot = PPP_Node("/")
        self.newroot.proc = self.thread
        for r in self.root.root_matches():
            nn = r.parent.path() if r != self.root else "/"
            newname = nn
            i = 1
            while newname in self.newroot.children:  # find some unique label
                newname = nn + " %d" % i
                i += 1
            self.newroot.children[newname] = r
        self.root = self.newroot


class PPP:
    """A collection of event trees from one or many Piernik processes from a single run"""

    def __init__(self, name):
        self.name = name
        self.trees = {}  # separate event tree for each process

    def print(self):
        print(self.name)
        for i in sorted(self.trees):
            self.trees[i].print()

    def _decode_text(self):
        self.nthr = 0
        self.bigbang = None
        file = open(self.name, 'r')
        for line in file:
            ln = line.split()
            if int(ln[0]) >= self.nthr:
                self.nthr = int(ln[0]) + 1
                self.bigbang = None  # force per-thread bigbang
            if self.bigbang is None:
                self.bigbang = float(ln[1]) - t_bias
            self._add(int(ln[0]), line[line.index(ln[2]):].strip(), float(ln[1]) + self.bigbang * (-1. if float(ln[1]) > 0. else 1.))  # proc, label, time
        self.sel_proc = []
        for p in self.trees:
            self.sel_proc.append(p)
        outof_str = " (%d out of %d thread%s)" % (len(set(self.sel_proc)), self.nthr, plural(self.nthr)) if self.nthr != len(set(self.sel_proc)) else ""
        self.descr = "'%s', thread%s %s" % (self.name, plural(set(self.sel_proc)), niceprint_set(self.sel_proc)) + outof_str

    def _add(self, proc, label, time):
        if (len(included_threads) == 0) or (proc in included_threads):
            if proc not in self.trees:
                self.trees[proc] = PPP_Tree(proc)
            self.trees[proc]._add(label, time)

    def calculate_time(self):
        self.total_time = 0.
        for _ in self.trees:
            self.trees[_].root.calculate_time()
            self.total_time += self.trees[_].root.child_time  # root does not have time

    def apply_exclusions(self):
        if args.exclude is not None or args.maxdepth is not None:
            for p in self.trees:
                self.trees[p].root.apply_exclusions()

    def rebase_root(self):
        if args.root is not None:
            earliest = None
            for p in self.trees:
                self.trees[p].rebase_root()
                for i in self.trees[p].root.children:
                    if earliest is None:
                        earliest = self.trees[p].root.children[i].start[0]
                    else:
                        earliest = min(earliest, self.trees[p].root.children[i].start[0])
            if earliest is not None:
                for p in self.trees:
                    self.trees[p].root.shift_time(earliest - t_bias)

    def apply_time_window(self):
        for p in self.trees:
            self.trees[p].root.time_filter()


class PPPset:
    """A collection of event trees from one or many Piernik runs"""

    def __init__(self, fnamelist):
        self.run = {}
        for fname in fnamelist:
            self.run[fname] = PPP(fname)
        for _ in self.run:
            self.run[_]._decode_text()
            self.run[_].apply_exclusions()
            self.run[_].rebase_root()
            self.run[_].apply_time_window()  # has to be called after shift_time was applied
            self.run[_].calculate_time()

    def print(self, otype):
        self.out = ""
        if otype == "tree":
            for _ in self.run:
                print("")
                self.run[_].print()
        elif otype == "gnu":
            self.print_gnuplot()
            if args.output is None:
                print(self.out)
            else:
                ofile = open(args.output[0], "w")
                ofile.write(self.out)
                ofile.close()
        elif otype == "summary":
            # ToDo: add "called from"
            print("ARGS: ", args)
            for f in self.run:
                print("\n## File: " + self.run[f].descr)
                ed = {}
                for p in self.run[f].trees:
                    for e in self.run[f].trees[p].root.get_all_ev():
                        e_base = e[0].split('/')[-1]
                        if e_base not in ed:
                            ed[e_base] = [0, 0., 0.]  # calls, cumul. own time, cumul. child. time
                        ed[e_base][0] += len(e[1])
                        ed[e_base][1] += e[3]
                        ed[e_base][2] += e[4]
                ml = len(max(ed, key=len)) if len(ed) > 0 else 0
                print("# %-*s %20s %8s %15s \033[97m%20s\033[0m %10s" % (ml - 2, "label", "avg. CPU time (s)", "%time", "avg. occurred", "avg. non-child", "% of total"))
                print("# %-*s %20s %8s %15s \033[97m%20s\033[0m %10s" % (ml - 2, "", "(per thread)", "in child", "(per thread)", "time (s)", "time"))
                skip_cnt, skip_val = 0, 0.
                for e in sorted(ed.items(), key=lambda x: x[1][1] - x[1][2], reverse=True):
                    if (e[1][1] - e[1][2]) / self.run[f].total_time > args.cutsmall[0] / 100.:
                        print("%-*s %20.7f %8s %15d%s \033[97m%20.7f\033[0m %10.2f" % (ml, e[0], e[1][1] / len(set(self.run[f].sel_proc)),
                                                                                       ("" if e[1][2] == 0. else "%8.2f" % ((100 * e[1][2] / e[1][1]) if e[1][1] > 0. else 0.)),
                                                                                       e[1][0] / len(set(self.run[f].sel_proc)),
                                                                                       " " if e[1][0] % len(set(self.run[f].sel_proc)) == 0 else "+",
                                                                                       (e[1][1] - e[1][2]) / len(set(self.run[f].sel_proc)),
                                                                                       100. * (e[1][1] - e[1][2]) / self.run[f].total_time))
                    else:
                        skip_cnt += 1
                        skip_val += e[1][1] - e[1][2]
                if skip_cnt > 0:
                    print("# (skipped %d timers that contributed %.7f seconds of non-child time = %.2f%% of total time)" % (skip_cnt, skip_val / len(set(self.run[f].sel_proc)), 100. * skip_val / self.run[f].total_time))

    def print_gnuplot(self):
        self.descr = ""
        ev = {}
        peff = 0
        gcnt = 0
        gomit = 0
        for f in self.run:
            pcnt = 0
            for p in self.run[f].trees:
                evlist = self.run[f].trees[p].root.get_all_ev()
                for e in evlist:
                    e_base = e[0].split('/')[-1].split()[0]
                    if e_base not in ev:
                        ev[e_base] = {}
                    if pcnt + peff not in ev[e_base]:
                        ev[e_base][pcnt + peff] = []
                    if e[3] > args.cutsmall[0] / 100. * self.run[f].total_time / len(set(self.run[f].sel_proc)):
                        ev[e_base][pcnt + peff].append(e)
                pcnt += 1
            peff = peff + len(set(self.run[f].sel_proc)) + (1 if len(self.run) > 1 else 0)  # spacing only when we have multiple files
            self.descr = self.run[f].descr + ("\\n" if len(self.descr) > 0 else "") + self.descr

        ndel = []
        for e in ev:
            n = 0
            for p in ev[e]:
                n += len(ev[e][p])
            if n == 0:
                ndel.append(e)
        for i in ndel:
            del (ev[i])
        if len(ndel) > 0:
            self.descr += "\\n(omitted %d timer%s below %g%% threshold)" % (len(ndel), "" if len(ndel) == 1 else "s", args.cutsmall[0])

        self.out += "#!/usr/bin/gnuplot\n$PPPdata << EOD\n\n"
        i_s = 1
        i_e = 2
        ev_cnt = {}
        for e in ev:
            ev_cnt[e] = 0
            for p in ev[e]:
                for t in range(len(ev[e][p])):
                    ev_cnt[e] += len(ev[e][p][t][i_s])
        for e in sorted(ev_cnt, key=lambda x: ev_cnt[x]):  # if we need to skip something due to MAXOUPUT, let it be the most abundant thing
            ftname = None
            for p in ev[e]:
                if ftname is None and len(ev[e][p]) > 0:
                    ftname = ev[e][p][0][0]
            self.out += "# __" + e + "__ '" + ftname + "'\n"
            for p in ev[e]:
                for t in range(len(ev[e][p])):
                    depth = ev[e][p][t][0].count("/") - 1
                    for _ in range(min(len(ev[e][p][t][i_s]), len(ev[e][p][t][i_e]))):
                        try:
                            if (ev[e][p][t][i_e][_] - ev[e][p][t][i_s][_]) > 0.:
                                if gcnt < args.maxoutput[0]:
                                    self.out += "%.7f %.7f %.7f %.7f %.7f %.7f %d %d\n" % (ev[e][p][t][i_s][_] - t_bias, depth,
                                                                                           ev[e][p][t][i_s][_] - t_bias, ev[e][p][t][i_e][_] - t_bias,
                                                                                           depth + float(p) / peff, depth + float(p + 1) / peff, depth, p)
                                else:
                                    gomit += 1
                                gcnt += 1
                        except TypeError:
                            sys.stderr.write("Warning: inclomplete event '", e, "' @", p, " #", str(t), ev[e][p][t])
                    self.out += "\n"
            self.out += "\n"
        if gomit > 0:
            self.descr += "\\n(omitted %d entries above maxoutput limit of %d)" % (gomit, args.maxoutput[0])
        self.out += "EOD\n\n# Suggested gnuplot commands:\nset key outside horizontal\n"
        self.out += "set xlabel 'time (walltime seconds)'\nset ylabel 'timer depth + proc/nproc'\n"
        self.out += 'set title "%s"\n' % self.descr.replace('_', "\\\\_")
        if len(ev) > 0:
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
            self.out += pline + "\n"
            if len(ndel) > 0:
                self.out += 'print "Timers below threshold: ' + " ".join(ndel) + '"\n'
            self.out += "pause mouse close\n"
        else:
            self.out += 'show title; print "No data to plot"'


parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                 description="Piernik Precise Profiling Presenter", epilog="""
examples:

plot profile of 64th step from file.ascii with gnuplot (hopefully in the interactive mode):
    ppp_plot.py file.ascii -r 'step 64'| gnuplot
    ppp_plot.py file.ascii -r 'step 64' -o file.gnu; gnuplot file.gnu

when gnuplot fails to set up desired teminal by default, try to set $GNUTERM (qt or x11 are recommended):
    ppp_plot.py file.ascii | GNUTERM=qt gnuplot

print list of top-level timers (steps) present in file.ascii:
    ppp_plot.py file.ascii -t -d 1 -p 0
(the step names are followed by their time offset and length)

find most time-consuming timers:
    ppp_plot.py file.ascii -s
(this also helps to identify most often called timers to be excluded in case of performance problems with gnuplot)

plot profile without the identified too-often called timer (e.g. "Loechner_mark") that makes interactive gnuplot to choke:
    ppp_plot.py file.ascii -e Loechner_mark | gnuplot
same as above but don't filter out timers that are contributing less than 0.1%:
    ppp_plot.py file.ascii -e Loechner_mark -% 0 | gnuplot

For massively parallel runs first try to plot only some processes (to avoid slowdowns due to large number of elements shown):
    ppp_plot.py file.ascii -p 0-3,128,255
then it is easier to apply filters like -r or -e and you may also narrow down the output to specified time interval:
    ppp_plot.py file.ascii -T 1.3 2.5 -p 0-12,255 -r "step 3" -e MPI_Waitall:restrict_1v
""")

parser.add_argument("filename", nargs='+', help="PPP ascii file(s) to process")
parser.add_argument("-o", "--output", nargs=1, help="processed output file (gnuplot only)")
parser.add_argument("-%", "--cutsmall", nargs=1, default=[.1], type=float, help="skip contributions below CUTSMALL%% (default = 0.1%%)")
parser.add_argument("-e", "--exclude", nargs='+', help="do not show EXCLUDEd timer(s)")  # multiple excudes
parser.add_argument("-r", "--root", nargs='+', help="show only ROOT and their children")
parser.add_argument("-d", "--maxdepth", type=int, help="limit output to MAXDEPTH")
parser.add_argument("-m", "--maxoutput", nargs=1, default=[50000], type=int, help="limit output to MAXOUTPUT enries (gnuplot only, default = 50000)")
parser.add_argument("-p", "--processes", nargs=1, help="list of threads to display in each file, syntax example: 1,3,7-10, default = all")
parser.add_argument("-T", "--timerange", nargs=2, help="show only events overlapping with specified time interval")
# parser.add_argument("-c", "--check", help="do a formal check only")

pgroup = parser.add_mutually_exclusive_group()
pgroup.add_argument("-g", "--gnuplot", action="store_const", dest="otype", const="gnu", default="gnu", help="gnuplot output (default)")
pgroup.add_argument("-t", "--tree", action="store_const", dest="otype", const="tree", help="tree output")
pgroup.add_argument("-s", "--summary", action="store_const", dest="otype", const="summary", help="short summary")

args = parser.parse_args()

# Process --processes
if args.processes is not None:
    for p in args.processes[0].split(','):
        pp = p.split('-')
        if len(pp) == 1:
            included_threads.append(int(pp[0]))
        elif len(pp) == 2:
            if int(pp[1]) >= int(pp[0]):  # ValueError will occur here if you try to use negative process numbers
                included_threads.extend(range(int(pp[0]), int(pp[1]) + 1))
            else:
                print("Reversed process order: ", pp[0], " > ", pp[1])
                raise ValueError
        else:
            print("I don't know how to interpret '" + p + "'")
            raise ValueError
included_threads = set(included_threads)

if args.timerange is not None:
    t_min = float(args.timerange[0])
    t_max = float(args.timerange[1])
    if t_min > t_max:
        print("t_min > t_max makes an empty interval")
        raise ValueError

all_events = PPPset(args.filename)
all_events.print(args.otype)

# An attempt to find unbalanced timers:
# for j  in *.ppprofile.ascii ; do for i in `awk '{print $1}' $j |  sort |  uniq` ; do in=`grep "$i *  -" $j | wc -l` ; out=`grep "$i *  [1-9]" $j | wc -l` ; [ $in != $out ] && echo $j $i $in $out ; done |  column -t ; done
# No output is good
# An output in the form:
#     filename thread_number closed_timers opened_timers
# means that something bad happened and the data is truncated or corrupted in some other way
