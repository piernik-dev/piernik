from yt.mods import load, SlicePlot, parallel_objects


def visualize(files):
    output = []
    for fn in parallel_objects(files, njobs=-1):
        pf = load(fn)
        slc = SlicePlot(pf, 2, 'denn')
        output += slc.save()
        slc = SlicePlot(pf, 2, 'enen')
        output += slc.save()
    return output
