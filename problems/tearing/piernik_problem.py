import matplotlib
matplotlib.use('cairo')
import matplotlib.pyplot as plt
from yt.mods import load, SlicePlot, parallel_objects

field = 'CurrentZ'

def visualize(files):
    output = []
    for fn in parallel_objects(files, njobs=-1):
        pf = load(fn)
        s = pf.h.slice(2, 0.0, fields=[field])
        img = s.to_frb((1, 'unitary'), (512, 512))
        vlim = abs(img[field]).max()
        fig = plt.figure(0, figsize=(8, 6))
        ax = plt.subplot(111)
        ax.set_title(r"$j_z = (\nabla\times\mathbf{B})_z$")
        img = ax.imshow(img[field], cmap='bwr', vmin=-vlim, vmax=vlim)
        plt.colorbar(img)
        plt.draw()
        fn_out = fn.replace('.h5', '_%s.png' % field)
        plt.savefig(fn_out)
        output.append(fn_out)
    return output

if __name__ == "__main__":
    import sys
    visualize(sys.argv[1:])
