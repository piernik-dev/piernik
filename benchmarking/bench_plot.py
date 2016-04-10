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
            if (len(columns) > 0):
                try:
                    nthr = int(columns[0])
                    if (nthr > 0):
                        if (nthr not in timings):
                            timings[nthr] = [0 for x in range(crtest_strong + 1)]
                        timings[nthr][b_type] = float(columns[len(columns)-1])
                except ValueError:
                    continue

    data["make_real"] = make_real
    data["make_load"] = make_load
    data["timings"] = timings
    return data


from sys import argv

if (len(argv) < 2):
    print "Usage: ", argv[0], " benchmark_file [benchmark_file ...]"
    exit(-1)

data = []
for f in argv[1:]:
    data.append(read_timings(f))

print data
