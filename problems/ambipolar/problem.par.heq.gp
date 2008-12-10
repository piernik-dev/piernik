 $DOMAIN_SIZES
    nxd = 5
    nyd = 5
    nzd = 240
    nb  = 4
 /  

 $MPI_BLOCKS
    pxsize = 1
    pysize = 1
    pzsize = 1    
 /  

 $START_CONTROL
 /

 $RESTART_CONTROL
    nrestart = 0 
    new_id   = ''
 /

 $END_CONTROL
    tend   = 500.0
    nend   = 10000 
 /

 $OUTPUT_CONTROL
    dt_hdf  = 1.0
    dt_res  = 500.0
    dt_log  = 1.0e-6
    dt_tsl  = 1.0e-6
    vars(1) = 'dens'
    vars(2) = 'velx'
    vars(3) = 'vely'
    vars(4) = 'velz'
    vars(5) = 'ener'
    vars(6) = 'magx'
    vars(7) = 'magy'
    vars(8) = 'magz'
    domain  = 'full_domain'
   
 / 

 $BOUNDARIES
    bnd_xl = 'per'  
    bnd_xr = 'per'  
    bnd_yl = 'per'  
    bnd_yr = 'per'  
    bnd_zl = 'outd'  
    bnd_zr = 'outd'
 /

 $DOMAIN_LIMITS
    xmin   =-250.0
    xmax   = 250.0
    ymin   =-1000.0
    ymax   = 1000.0
    zmin   =-500.0
    zmax   = 500.0
 /

 $EQUATION_OF_STATE
    c_si   =  7.
    alpha  = 0.0 
 /  

 $NUMERICAL_SETUP
    cfl    = 0.8
    smalld = 1.e-5
    smallei= 1.e-4
    freezing_speed = 'local'
    integration_order  = 2
  /
 
 $GRAVITY
    r_gc   = 8.5
    nsub   = 100
    tune_zeq     = 1.0  !1.0025
    tune_zeq_bnd = 1.0  !0.9
 /
 
 $THERMAL
 /

 $PROBLEM_CONTROL
    problem_name ='parker_instability'
    run_id =  'gpi'
    d0     =  1.
    nbx0    =  0.
    nby0    = 1.0
    nbz0    =  0.
    a_vp    = 0.0
 /
