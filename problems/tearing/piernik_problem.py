import numpy as np
from yt.mods import \
    load, SlicePlot, parallel_objects, add_field, ValidateSpatial

FIELDS= ['curz']

def _CurrentZ(field, data):
    new_field = np.zeros(data["magy"].shape, dtype='float64')
    new_field[1:-1,1:-1,:] = (data["magy"][1:-1,1:-1,:] -
                              data["magy"][0:-2,1:-1,:]) \
                              / data["dx"].flat[0]
    new_field[1:-1,1:-1,:] -= (data["magx"][1:-1,1:-1,:] -
                               data["magx"][1:-1,0:-2,:]) \
                               / data["dy"].flat[0]
    return new_field
def _convertCurrent(data):
    return 1.0 / data.convert("cm")
add_field("curz", function=_CurrentZ, convert_function=_convertCurrent,
          validators=[ValidateSpatial(1, ["magx", "magy", "magz"])],
          units=r"\rm{g}/(\rm{cm}\,\rm{s}^{3}\,\rm{gauss})")

def visualize(files):
    output = []
    for fn in parallel_objects(files, njobs=-1):
        pf = load(fn)
        for field in FIELDS:
            slc = SlicePlot(pf, 'z', field)
            if field == 'curz':
                slc.set_cmap(field, 'bwr')
                maxabs = abs(slc._frb[field]).max()
                slc.set_log(field, False)
                slc.set_zlim(field, -maxabs, maxabs)
            output.append(slc.save(fn.replace('.h5', '_%s.png' % field))[0])
    return output

if __name__ == "__main__":
    import sys
    print visualize(sys.argv[1:])
