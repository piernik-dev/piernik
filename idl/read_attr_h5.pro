function read_attr_h5, file

   a = { time:          0.0,  $
         timestep:      0.0,  $
         nstep:         0l,   $
         nxd:           0l,   $
         nyd:           0l,   $
         nzd:           0l,   $
         nb:            0l,   $
         crdsys:        'xyz',$
         xmin:          0.0,  $
         xmax:          0.0,  $
         ymin:          0.0,  $
         ymax:          0.0,  $
         zmin:          0.0,  $
         zmax:          0.0,  $
         pxsize:        0l,   $
         pysize:        0l,   $
         pzsize:        0l,   $
         domain:        '',   $
         problem_name:  '',   $
         run_id:        ''      }

   file_id = H5F_OPEN(file)

   a.time = 0.0d
   attr = H5A_OPEN_NAME(file_id, "time")
   a.time = H5A_READ(attr)
   H5A_CLOSE,attr

   a.timestep = 0.0d
   attr = H5A_OPEN_NAME(file_id, "timestep")
   a.timestep = H5A_READ(attr)
   H5A_CLOSE,attr

   a.nstep = 0
   attr = H5A_OPEN_NAME(file_id, "nstep")
   a.nstep = H5A_READ(attr)
   H5A_CLOSE,attr

   a.nxd = 0
   attr = H5A_OPEN_NAME(file_id, "nxd")
   a.nxd = H5A_READ(attr)
   H5A_CLOSE,attr

   a.nyd = 0
   attr = H5A_OPEN_NAME(file_id, "nyd")
   a.nyd = H5A_READ(attr)
   H5A_CLOSE,attr

   a.nzd = 0
   attr = H5A_OPEN_NAME(file_id, "nzd")
   a.nzd = H5A_READ(attr)
   H5A_CLOSE,attr

   a.nb = 0
   attr = H5A_OPEN_NAME(file_id, "nb")
   a.nb = H5A_READ(attr)
   H5A_CLOSE,attr

   a.xmin = 0
   attr = H5A_OPEN_NAME(file_id, "xmin")
   a.xmin = H5A_READ(attr)
   H5A_CLOSE,attr

   a.xmax = 0
   attr = H5A_OPEN_NAME(file_id, "xmax")
   a.xmax = H5A_READ(attr)
   H5A_CLOSE,attr

   a.ymin = 0
   attr = H5A_OPEN_NAME(file_id, "ymin")
   a.ymin = H5A_READ(attr)
   H5A_CLOSE,attr

   a.ymax = 0
   attr = H5A_OPEN_NAME(file_id, "ymax")
   a.ymax = H5A_READ(attr)
   H5A_CLOSE,attr

   a.zmin = 0
   attr = H5A_OPEN_NAME(file_id, "zmin")
   a.zmin = H5A_READ(attr)
   H5A_CLOSE,attr

   a.zmax = 0
   attr = H5A_OPEN_NAME(file_id, "zmax")
   a.zmax = H5A_READ(attr)
   H5A_CLOSE,attr

   domain = ''
   attr = H5A_OPEN_NAME(file_id, "domain")
   a.domain = H5A_READ(attr)
   H5A_CLOSE,attr

   problem_name = ''
   attr = H5A_OPEN_NAME(file_id, "problem_name")
   a.problem_name = H5A_READ(attr)
   H5A_CLOSE,attr

   run_id = ''
   attr = H5A_OPEN_NAME(file_id, "run_id")
   a.run_id = H5A_READ(attr)
   H5A_CLOSE,attr

   H5F_CLOSE, file_id

   RETURN, a
end
