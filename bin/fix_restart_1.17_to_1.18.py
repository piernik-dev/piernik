#/usr/bin/env python

import sys
import h5py as h5


for f in sys.argv[1:]:
    h5f = h5.File(f)
    try:
       print("Piernik Version %f"%(h5f.attrs['piernik']))
    except:
       h5f.close()
       print("%s is not piernik restart!!!"%(f))
       sys.exit()

    remove_attrs = ['next_t_tsl', 'next_t_log', 'step_res', 'step_hdf']
    add_attrs = ['last_res_time', 'last_plt_time', 'last_tsl_time', 'last_log_time']

    for attr in remove_attrs:
        try:
            del h5f.attrs[attr]
        except:
            pass

    keys = h5f.attrs.keys()
    for attr in add_attrs:
        if not attr in keys:
           h5f.attrs[attr]=h5f.attrs['time']
    if not 'nimg' in keys:
       h5f.attrs['nimg'] = [0]
    h5f.attrs['piernik'] = [1.18]
    h5f.close()

