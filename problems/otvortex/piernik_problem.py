from yt.mods import \
    load, SlicePlot

FIELDS = ["prei"]

def visualize(files):
    output = []
    for fn in files:
        pf = load(fn)
        for field in FIELDS:
            slc = SlicePlot(pf, 'z', field)
            if field == "prei":
                slc.set_cmap(field, 'gist_stern')
            output.append(slc.save(fn.replace('.h5', '_%s.png' % field))[0])
    return output

if __name__ == "__main__":
    import sys
    print visualize(sys.argv[1:])
