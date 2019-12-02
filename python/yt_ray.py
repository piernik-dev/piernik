#!/usr/bin/env python

'''If you are in runs/problem directory, you got a bunch of .h5 files
and you want to get values along just one line, you may do:

plot "< ../../python/yt_ray.py file.h5 0 0 0 0 0 0"

it will print you some fields and domain extent, like this:

yt : [WARNING] 'field_units' was overridden by 'dataset_units/density'
yt : [WARNING] 'field_units' was overridden by 'dataset_units/pressure'
yt : [WARNING] 'field_units' was overridden by 'dataset_units/specific_energy'
yt : [WARNING] 'field_units' was overridden by 'dataset_units/velocity_x'
yt : [WARNING] 'field_units' was overridden by 'dataset_units/velocity_y'
yt : [INFO   ] Parameters: current_time              = [ 1.]
yt : [INFO   ] Parameters: domain_dimensions         = [1024 1024    1]
yt : [INFO   ] Parameters: domain_left_edge          = [-1. -1.  0.]
yt : [INFO   ] Parameters: domain_right_edge         = [ 1.  1.  1.]

In columns 1-3, you will get x, y and z coordinates. further columns will
contain "density", "pressure" and so on, in the order as above.
Now plot something:

plot "< ../../python/yt_ray.py file.h5 -1 -.5 0 1 -.5 0" u 1:4

The profile will be density at y=-0.5 coordinate of that 2D data'''

try:
    import yt
except:
    print("you must make yt available somehow")
    exit(-1)


def print_ray(fname, startpoint, endpoint):

    ds = yt.load(fname)
    line = ds.ray(startpoint, endpoint)
    labels = ["x", "y", "z"]
    for f in ds.field_list:
        labels.append(f[1])
    print("Columns: ", labels, file=sys.stderr)

    print("#", end='')
    for l in labels:
        print("{:<20} ".format(l), end='')
    print()

    dlen = len(line[labels[0]].v)
    for i in range(dlen):
        for l in labels:
            print("{:<20.12g} ".format(line[l].v[i]), end='')
        print()

if __name__ == "__main__":
    import sys
    if (len(sys.argv) < 8):
        print("Error: too few arguments.\nUsage: " +
              sys.argv[0] + " hdf_file x_start y_start z_start x_end y_end z_end", file=sys.stderr)
    else:
        ndim = 3
        print_ray(sys.argv[1],
                  (sys.argv[2:2 + ndim]), (sys.argv[2 + ndim:2 + 2 * ndim]))
        if (len(sys.argv) > 8):
            print("Ignored arguments: ", sys.argv[8:len(sys.argv)], file=sys.stderr)
