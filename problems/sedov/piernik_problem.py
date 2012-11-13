import matplotlib
matplotlib.use('cairo')
import matplotlib.pyplot as plt
from yt.mods import load, SlicePlot, parallel_objects

fields = ['denn', 'enen']

def visualize(files):
    output = []
    for fn in parallel_objects(files, njobs=-1):
        pf = load(fn)
        s = pf.h.slice(2, 0.0, fields=fields)
        img = s.to_frb(2.0, (512, 512))
        for field in fields:
            fig = plt.figure(0, figsize=(8, 6))
            ax = plt.subplot(111)
            ax.imshow(img[field])
            plt.draw()
            fn_out = fn.replace('.h5', '_%s.png' % field)
            print fn_out
            plt.savefig(fn_out)
        output.append(fn_out)
    return output

if __name__ == "__main__":
    import sys
    visualize(sys.argv[1:])
