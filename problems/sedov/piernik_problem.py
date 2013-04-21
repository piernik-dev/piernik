from yt.mods import \
    load, SlicePlot, parallel_objects

FIELDS = ['denn', 'enen']

def visualize(files):
    output = []
    for fn in parallel_objects(files, njobs=-1):
        pf = load(fn)
        for field in FIELDS:
            slc = SlicePlot(pf, 'z', field)
            output.append(slc.save(fn.replace('.h5', '_%s.png' % field))[0])
    return output

if __name__ == "__main__":
    import sys
    print visualize(sys.argv[1:])
