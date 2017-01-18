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
        if ("min" in rdata[d]):
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
        if ("min" in rdata[d]):
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
                if ("min" in rdata[d]):
                    ymin.append(rdata[d]["min"]["timings"][x][test])
                    ymax.append(rdata[d]["max"]["timings"][x][test])
            if (test in (sedov_strong, maclaurin_strong, crtest_strong)):
                for i in range(len(y)):
                    if (y[i]):
                        y[i] *= n[i]
                        if ("min" in rdata[d]):
                            ymin[i] *= n[i]
                            ymax[i] *= n[i]
            ywhere = np.empty_like(y, dtype=bool)
            if ("min" in rdata[d]):
                for i in range(len(y)):
                    ywhere[i] = ymin[i] and ymax[i]
                    if (not ywhere[i]):
                        ymin[i] = 0.
                        ymax[i] = 0.
            plt.plot(n, y)
            if ("min" in rdata[d]):
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


def singlesample(data):
    import os.path
    rd = {}
    for d in data:
        d["dname"]=d["filename"]
        rd[d["dname"]] = {}
        rd[d["dname"]]["avg"] = {}
        for i in ("make_real", "make_load", "timings"):
            rd[d["dname"]]["avg"][i] = d[i]

    return rd


def reduce(data):
    import os.path
    import numpy as np
    from copy import deepcopy

    rd = {}
    for d in data:
        d["dname"]=os.path.dirname(d["filename"])
        if (d["dname"] not in rd):
            rd[d["dname"]] = {}
            rd[d["dname"]]["nt"] = 1
            rd[d["dname"]]["nm"] = 0
            if (np.product(d["make_real"]) * np.product(d["make_load"]) != 0):
                rd[d["dname"]]["nm"] = 1
            rd[d["dname"]]["avg"] = {}
            for i in ("make_real", "make_load", "timings"):
                rd[d["dname"]]["avg"][i] = d[i]
            for i in ("min", "max"):
                rd[d["dname"]][i] = deepcopy(rd[d["dname"]]["avg"])
        else:
            if (np.product(d["make_real"]) * np.product(d["make_load"]) != 0):
                rd[d["dname"]]["nm"] += 1
                for i in ("make_real", "make_load"):
                    if (rd[d["dname"]]["nm"] > 1):
                        rd[d["dname"]]["min"][i] = np.minimum(rd[d["dname"]]["min"][i], d[i])
                        rd[d["dname"]]["max"][i] = np.maximum(rd[d["dname"]]["max"][i], d[i])
                        rd[d["dname"]]["avg"][i] = np.add(rd[d["dname"]]["avg"][i], d[i])
                if (rd[d["dname"]]["nm"] == 1):
                    for i in ("make_real", "make_load"):
                        rd[d["dname"]]["avg"][i] = d[i]
                    for i in ("min", "max"):
                        for j in ("make_real", "make_load"):
                            rd[d["dname"]][i][j] = deepcopy(rd[d["dname"]]["avg"][j])
            rd[d["dname"]]["nt"] += 1
            for p in d["timings"]:
                rd[d["dname"]]["min"]["timings"][p] = np.minimum(rd[d["dname"]]["min"]["timings"][p], d["timings"][p])
                rd[d["dname"]]["max"]["timings"][p] = np.maximum(rd[d["dname"]]["max"]["timings"][p], d["timings"][p])
                try:
                    rd[d["dname"]]["avg"]["timings"][p] = np.add(rd[d["dname"]]["avg"]["timings"][p], d["timings"][p])
                except TypeError:
                    for i in range(len(rd[d["dname"]]["avg"]["timings"][p])):
                        if (rd[d["dname"]]["avg"]["timings"][p][i] == None or d["timings"][p][i] == None):
                            rd[d["dname"]]["avg"]["timings"][p][i] = None
                        else:
                            rd[d["dname"]]["avg"]["timings"][p][i] = rd[d["dname"]]["avg"]["timings"][p][i] + d["timings"][p][i]

    for d in rd:
        if (rd[d]["nm"] > 1):
            for i in ("make_real", "make_load"):
                rd[d]["avg"][i] /= rd[d]["nm"]
        if (rd[d]["nt"] > 1):
            for p in rd[d]["avg"]["timings"]:
                try:
                    rd[d]["avg"]["timings"][p] /= rd[d]["nt"]
                except TypeError:
                    for i in range(len(rd[d]["avg"]["timings"][p])):
                        if (rd[d]["avg"]["timings"][p][i] != None):
                            rd[d]["avg"]["timings"][p][i] /= rd[d]["nt"]

    return rd


from sys import argv

if (len(argv) < 2):
    print "Usage: ", argv[0], " benchmark_file [benchmark_file ...]"
    exit(-1)

data = []
for f in argv[1:]:
    data.append(read_timings(f))

rdata = []
rdata = reduce(data)
mkrplot(rdata)

sdata = []
sdata = singlesample(data)
mkrplot(sdata)
