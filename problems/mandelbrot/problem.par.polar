 $BASE_DOMAIN
    n_d = 384, 64, 1
    xmin   = -15
    xmax   =  0.5
    ymin   = -3.141592653589793238
    ymax   =  3.141592653589793238
    bnd_xl = "outd"
    bnd_xr = "outd"
    bnd_yl = "outd"
    bnd_yr = "outd"
 /

 $MPI_BLOCKS
    AMR_bsize = 3*16
 /

 $UNITS
    units_set = "cgs"
 /

 $RESTART_CONTROL
    restart  = 'none'
    new_id   = ''
    nrestart = 0
 /

 $END_CONTROL
    tend   = 0.0
    nend   = 0
 /

 $OUTPUT_CONTROL
    problem_name = 'mandelbrot'
    run_id =  'ts0'
    use_v2_io = .true.
    dt_hdf  = 0.5
    vars(1:) = 'mand', 'dist', 'ang', 'level'
 /

 $FLUID_NEUTRAL
 /

 $NUMERICAL_SETUP
 /

 $PROBLEM_CONTROL
    maxiter = 100000
    !smooth_map = .false.
    deref_thr = 0.
    log_polar = .true.
    x_polar = -0.75232384419713483
!-0.91845036919971235

    y_polar = -0.03980694105081919
!0.32967301364153556
 /

 $INTERACTIONS
 /

 $AMR
    level_max = 3
    n_updAMR = 1
 /
