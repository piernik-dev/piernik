import matplotlib
matplotlib.use('cairo')
import matplotlib.pyplot as plt
from yt.mods import load, SlicePlot


def visualize(files):
    output = []
    for fn in files:
        fig = plt.figure(0, figsize=(8, 6))
        ax = plt.subplot(111)
        pf = load(fn)
        s = pf.h.slice(2, 0.0, fields=["prei"])
        img = s.to_frb((1.0, "unitary"), (512, 512))
        ax.imshow(img['prei'], cmap='gist_stern')
        plt.draw()
        fn_out = fn.replace('.h5', '.png')
        plt.savefig(fn_out)
        output.append(fn_out)
    return output

if __name__ == "__main__":
    import sys
    print visualize(sys.argv[1:])
