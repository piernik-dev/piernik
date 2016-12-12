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
                            timings[nthr] = [0 for x in range(crtest_strong + 1)]
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
        plt.xticks(range(1,ntm+1))

    names = []
    for d in data:
        names.append(d["filename"])

    plt.subplots_adjust(top=0.95, bottom=0.05+0.025*int((len(data)-1)/2+1), left=0.05, right=0.95, wspace=0.1)
    plt.figlegend((lines), names, loc="lower center", ncol=2, frameon=False)
    plt.annotate("Piernik benchmarks", xy=(0.5, 0.97), xycoords="figure fraction", horizontalalignment='center', size=20)

    plt.show()


from sys import argv

if (len(argv) < 2):
    print "Usage: ", argv[0], " benchmark_file [benchmark_file ...]"
    exit(-1)

data = []
for f in argv[1:]:
    data.append(read_timings(f))

mkplot(data)
