 $BASE_DOMAIN
    n_d = 384, 64, 1
    xmin   = -15
    xmax   =  0.5
    ymin   = -1.5
    ymax   =  4.7832
    bnd_xl = "outd"
    bnd_xr = "outd"
    bnd_yl = "outd"
    bnd_yr = "outd"
 /

 $MPI_BLOCKS
 /

 $UNITS
    units_set = "cgs"
 /

 $RESTART_CONTROL
    restart  = 'none'
    res_id   = ''
 /

 $END_CONTROL
    tend   = 0.0
    nend   = 0
 /

 $OUTPUT_CONTROL
    problem_name = 'mandelbrot_polar'
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
    log_polar = .true.
    x_polar = -0.75232384419713483
    y_polar = -0.03980694105081919
    c_polar = 0.05
 /

 $INTERACTIONS
 /

 $AMR
    bsize = 3*16
    level_max = 3
    n_updAMR = 1
 /
