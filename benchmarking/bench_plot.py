#!/usr/bin/env python

sedov_weak, sedov_strong, maclaurin_weak, maclaurin_strong, crtest_weak, crtest_strong = range(6)
make_prep, make_11, make_1n, make_2n, make_4n, make_8n = range(6)


def extr_make_t(columns):
    return float(columns[len(columns)-4]), float(columns[len(columns)-1].replace('%', ''))


def read_timings(file):

    import re  # overkill, I know

    data = {}
    data["filename"] = file
    with open(file, "r") as f:
        make_real = [0 for x in range(make_8n + 1)]
        make_load = [0 for x in range(make_8n + 1)]
        timings = {}

        b_type = -1
        for line in f:
            columns = line.split()
            if re.match("Preparing objects", line):
                make_real[make_prep], make_load[make_prep] = extr_make_t(columns)
            if re.match("Single-thread make object", line):
                make_real[make_11], make_load[make_11] = extr_make_t(columns)
            if re.match("Multi-thread make object", line):
                make_real[make_1n], make_load[make_1n] = extr_make_t(columns)
            if re.match("Multi-thread make two objects", line):
                make_real[make_2n], make_load[make_2n] = extr_make_t(columns)
            if re.match("Multi-thread make four objects", line):
                make_real[make_4n], make_load[make_4n] = extr_make_t(columns)
            if re.match("Multi-thread make eight objects", line):
                make_real[make_8n], make_load[make_8n] = extr_make_t(columns)
            if re.match("(.*)sedov, weak", line):
                b_type = sedov_weak
            if re.match("(.*)sedov, strong", line):
                b_type = sedov_strong
            if re.match("(.*)maclaurin, weak", line):
                b_type = maclaurin_weak
            if re.match("(.*)maclaurin, strong", line):
                b_type = maclaurin_strong
            if re.match("(.*)crtest, weak", line):
                b_type = crtest_weak
            if re.match("(.*)crtest, strong", line):
                b_type = crtest_strong
            if (b_type in (crtest_weak, crtest_strong)):
                d_col = 3
            elif (b_type in (sedov_weak, sedov_strong)):
                d_col = 6
            elif (b_type in (maclaurin_weak, maclaurin_strong)):
                d_col = 5
            if (len(columns) > 0):
                try:
                    nthr = int(columns[0])
                    if (nthr > 2**20): # crude protection against eating too much memory due to bad data lines
                        print "Ignoring bogus thread number: ", columns
                    elif (nthr > 0):
                        if (nthr not in timings):
                            timings[nthr] = [None for x in range(crtest_strong + 1)]
                        if (len(columns) >= d_col+1):
                            timings[nthr][b_type] = float(columns[d_col])
                        else:
                            timings[nthr][b_type] = None
                except ValueError:
                    continue

    data["make_real"] = make_real
    data["make_load"] = make_load
    data["timings"] = timings
    return data


def mkplot(data):
    import matplotlib.pyplot as plt
    import math as m

    plt.figure(figsize=(20, 15))

    m_labels = ["setup", "serial make", "parallel make", "parallel make\n2 objects", "parallel make\n4 objects", "parallel make\n8 objects"]
    t_labels = ["sedov, weak scaling", "sedov, strong scaling", "maclaurin, weak scaling", "maclaurin, strong scaling", "crtest, weak scaling", "crtest, strong scaling"]

    exp = 0.25
    sub = 1
    lines = []
    plt.subplot(4, 2, sub)
    for d in data:
        l, = plt.plot(d["make_real"])
        lines.append(l)
    plt.ylabel("time [s]")
    plt.xticks(range(len(d["make_real"])), m_labels)
    plt.annotate("compilation time", xy=(0.5, 0.1), xycoords="axes fraction", horizontalalignment='center')
    plt.ylim(ymin=0.)
    plt.xlim(-exp, len(m_labels)-1+exp)

    sub = 2
    plt.subplot(4, 2, sub)
    for d in data:
        plt.plot(d["make_load"])
    plt.ylabel("CPU load [%]")
    plt.xticks(range(len(d["make_load"])), m_labels)
    plt.annotate("compilation CPU usage", xy=(0.5, 0.1), xycoords="axes fraction", horizontalalignment='center')
    plt.ylim(ymin=0.)
    plt.xlim(-exp, len(m_labels)-1+exp)

    ntm = 0
    for d in data:
        for k in d["timings"].keys():
            ntm = max(ntm, k)

    for test in (sedov_weak, sedov_strong, maclaurin_weak, maclaurin_strong, crtest_weak, crtest_strong):
        sub += 1
        plt.subplot(4, 2, sub)
        for d in data:
            n = d["timings"].keys()
            y = []
            for x in n:
                y.append(d["timings"][x][test])
            if (test in (sedov_strong, maclaurin_strong, crtest_strong)):
                for i in range(len(y)):
                    if (y[i]):
                        y[i] *= n[i]
            plt.plot(n, y)
        plt.xlabel("N_threads", verticalalignment='center')
        if (test in (sedov_strong, maclaurin_strong, crtest_strong)):
            plt.ylabel("time * N_threads [s]")
        else:
            plt.ylabel("time [s]")
        plt.annotate(t_labels[test], xy=(0.5, 0.1), xycoords="axes fraction", horizontalalignment='center')
        plt.ylim(ymin=0.)
        plt.xlim(1-exp, ntm+exp)

        if (ntm >= 10):
            xf, xi = m.modf(m.log10(ntm))
            xf = pow(10, xf)
            if (xf >= 5.):
                xf = 1
                xi += 1
            elif (xf >= 2.):
                xf = 5
            else:
                xf = 2
            xtstep = int(xf * m.pow(10, xi-1))
            x_ticks = range(0, ntm+xtstep, xtstep)
        else:
            x_ticks = range(1, ntm+1)
        plt.xticks(x_ticks)

    names = []
    for d in data:
        names.append(d["filename"])

    plt.subplots_adjust(top=0.95, bottom=0.05+0.025*int((len(data)-1)/2+1), left=0.05, right=0.95, wspace=0.1)
    plt.figlegend((lines), names, loc="lower center", ncol=2, frameon=False)
    plt.annotate("Piernik benchmarks", xy=(0.5, 0.97), xycoords="figure fraction", horizontalalignment='center', size=20)

#    plt.show()

def mkrplot(rdata):
    import matplotlib.pyplot as plt
    import math as m
    import numpy as np

    plt.figure(figsize=(20, 15))

    m_labels = ["setup", "serial make", "parallel make", "parallel make\n2 objects", "parallel make\n4 objects", "parallel make\n8 objects"]
    t_labels = ["sedov, weak scaling", "sedov, strong scaling", "maclaurin, weak scaling", "maclaurin, strong scaling", "crtest, weak scaling", "crtest, strong scaling"]

    alph = 0.2
    exp = 0.25
    sub = 1
    lines = []
    ld = {}
    plt.subplot(4, 2, sub)
    for d in rdata:
        l, = plt.plot(rdata[d]["avg"]["make_real"])
        plt.fill_between(range(len(rdata[d]["avg"]["make_real"])), rdata[d]["min"]["make_real"], rdata[d]["max"]["make_real"], alpha=alph, color=l.get_color())
        lines.append(l)
        ld[d] = l
    plt.ylabel("time [s]")
    plt.xticks(range(len(rdata[d]["avg"]["make_real"])), m_labels)
    plt.annotate("compilation time", xy=(0.5, 0.1), xycoords="axes fraction", horizontalalignment='center')
    plt.ylim(ymin=0.)
    plt.xlim(-exp, len(m_labels)-1+exp)

    sub = 2
    plt.subplot(4, 2, sub)
    for d in rdata:
        plt.plot(rdata[d]["avg"]["make_load"])
        plt.fill_between(range(len(rdata[d]["avg"]["make_load"])), rdata[d]["min"]["make_load"], rdata[d]["max"]["make_load"], alpha=alph, color=ld[d].get_color())
    plt.ylabel("CPU load [%]")
    plt.xticks(range(len(rdata[d]["avg"]["make_load"])), m_labels)
    plt.annotate("compilation CPU usage", xy=(0.5, 0.1), xycoords="axes fraction", horizontalalignment='center')
    plt.ylim(ymin=0.)
    plt.xlim(-exp, len(m_labels)-1+exp)

    ntm = 0
    for d in rdata:
        for k in rdata[d]["avg"]["timings"].keys():
            ntm = max(ntm, k)

    for test in (sedov_weak, sedov_strong, maclaurin_weak, maclaurin_strong, crtest_weak, crtest_strong):
        sub += 1
        plt.subplot(4, 2, sub)
        for d in rdata:
            n = rdata[d]["avg"]["timings"].keys()
            y = []
            ymin = []
            ymax = []
            for x in n:
                y.append(rdata[d]["avg"]["timings"][x][test])
                ymin.append(rdata[d]["min"]["timings"][x][test])
                ymax.append(rdata[d]["max"]["timings"][x][test])
            if (test in (sedov_strong, maclaurin_strong, crtest_strong)):
                for i in range(len(y)):
                    if (y[i]):
                        y[i] *= n[i]
                        ymin[i] *= n[i]
                        ymax[i] *= n[i]
            ywhere = np.empty_like(y, dtype=bool)
            for i in range(len(y)):
                ywhere[i] = ymin[i] and ymax[i]
                if (not ywhere[i]):
                    ymin[i] = 0.
                    ymax[i] = 0.
            plt.plot(n, y)
            plt.fill_between(n, ymin, ymax, alpha=alph, color=ld[d].get_color(), where=ywhere)
        plt.xlabel("N_threads", verticalalignment='center')
        if (test in (sedov_strong, maclaurin_strong, crtest_strong)):
            plt.ylabel("time * N_threads [s]")
        else:
            plt.ylabel("time [s]")
        plt.annotate(t_labels[test], xy=(0.5, 0.1), xycoords="axes fraction", horizontalalignment='center')
        plt.ylim(ymin=0.)
        plt.xlim(1-exp, ntm+exp)

        if (ntm >= 10):
            xf, xi = m.modf(m.log10(ntm))
            xf = pow(10, xf)
            if (xf >= 5.):
                xf = 1
                xi += 1
            elif (xf >= 2.):
                xf = 5
            else:
                xf = 2
            xtstep = int(xf * m.pow(10, xi-1))
            x_ticks = range(0, ntm+xtstep, xtstep)
        else:
            x_ticks = range(1, ntm+1)
        plt.xticks(x_ticks)

    names = []
    for d in rdata:
        names.append(d)

    plt.subplots_adjust(top=0.95, bottom=0.05+0.025*int((len(rdata)-1)/2+1), left=0.05, right=0.95, wspace=0.1)
    plt.figlegend((lines), names, loc="lower center", ncol=2, frameon=False)
    plt.annotate("Piernik benchmarks", xy=(0.5, 0.97), xycoords="figure fraction", horizontalalignment='center', size=20)

    plt.show()

def reduce(data):
    import os.path
    import pprint
    import numpy as np
    from copy import deepcopy

    rd = {}
    dirnames = set()
    for d in data:
        d["dirname"]=os.path.dirname(d["filename"])
        dirnames.add(d["dirname"])
        if (d["dirname"] not in rd):
            rd[d["dirname"]] = {}
            rd[d["dirname"]]["n"] = 1
            rd[d["dirname"]]["avg"] = {}
            for i in ("make_real", "make_load", "timings"):
                rd[d["dirname"]]["avg"][i] = d[i]
            for i in ("min", "max"):
                rd[d["dirname"]][i] = deepcopy(rd[d["dirname"]]["avg"])
        else:
            rd[d["dirname"]]["n"] += 1
            for i in ("make_real", "make_load"):
                rd[d["dirname"]]["min"][i] = np.minimum(rd[d["dirname"]]["min"][i], d[i])
                rd[d["dirname"]]["max"][i] = np.maximum(rd[d["dirname"]]["max"][i], d[i])
                rd[d["dirname"]]["avg"][i] = np.add(rd[d["dirname"]]["avg"][i], d[i])
            for p in d["timings"]:
                rd[d["dirname"]]["min"]["timings"][p] = np.minimum(rd[d["dirname"]]["min"]["timings"][p], d["timings"][p])
                rd[d["dirname"]]["max"]["timings"][p] = np.maximum(rd[d["dirname"]]["max"]["timings"][p], d["timings"][p])
                try:
                    print d["filename"],p,rd[d["dirname"]]["avg"]["timings"][p], d["timings"][p]
                    rd[d["dirname"]]["avg"]["timings"][p] = np.add(rd[d["dirname"]]["avg"]["timings"][p], d["timings"][p])
                except TypeError:
                    for i in range(len(rd[d["dirname"]]["avg"]["timings"][p])):
                        if (rd[d["dirname"]]["avg"]["timings"][p][i] == None or d["timings"][p][i] == None):
                            rd[d["dirname"]]["avg"]["timings"][p][i] = None
                        else:
                            rd[d["dirname"]]["avg"]["timings"][p][i] = rd[d["dirname"]]["avg"]["timings"][p][i] + d["timings"][p][i]

    for d in rd:
        if (rd[d]["n"] > 1):
            for i in ("make_real", "make_load"):
                rd[d]["avg"][i] /= rd[d]["n"]
            for p in rd[d]["avg"]["timings"]:
                try:
                    rd[d]["avg"]["timings"][p] /= rd[d]["n"]
                except TypeError:
                    print rd[d]["avg"]["timings"][p]
                    for i in range(len(rd[d]["avg"]["timings"][p])):
                        if (rd[d]["avg"]["timings"][p][i] != None):
                            rd[d]["avg"]["timings"][p][i] /= rd[d]["n"]

    print dirnames
    pp = pprint.PrettyPrinter(indent=6)

    pp.pprint(data)
    pp.pprint(rd)

    return rd


from sys import argv

if (len(argv) < 2):
    print "Usage: ", argv[0], " benchmark_file [benchmark_file ...]"
    exit(-1)

data = []
for f in argv[1:]:
    data.append(read_timings(f))

mkplot(data)

rdata = []
rdata = reduce(data)

mkrplot(rdata)
