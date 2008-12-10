 $DOMAIN_SIZES
    nxd =5
    nyd =5
    nzd =120
 /  

 $START_CONTROL
 /

 $RESTART_CONTROL
    nrestart = 0 
    new_id   = ''
 /

 $END_CONTROL
    tend   = 10.0
    nend   = 100000
 /

 $OUTPUT_CONTROL
    dt_hdf  = 1.0
    dt_res  = 0.0
    dt_log  = 5.0
    vars(1) = 'dens'
    vars(2) = 'velx'
    vars(3) = 'vely'
    vars(4) = 'velz'
    vars(5) = 'ener'
    vars(6) = 'magx'
    vars(7) = 'magy'
    vars(8) = 'magz'
    vars(9) = 'esrc'
    vars(10) = 'eint'
    vars(11) = 'temp'
    domain  = 'phys_domain'
   
 / 

 $BOUNDARIES
    bnd_xl = 'per'  
    bnd_xr = 'per'  
    bnd_yl = 'per'  
    bnd_yr = 'per'  
    bnd_zl = 'out'  
    bnd_zr = 'out'
 /

 $DOMAIN_LIMITS
    xmin   =-250.0
    xmax   = 250.0
    ymin   =-1000.0
    ymax   = 1000.0
    zmin   =-600.0
    zmax   = 600.0
 /

 $EQUATION_OF_STATE
    alpha  = 1.0 
 /  

 $NUMERICAL_SETUP
    cfl    = 0.9
    smalld = 1.e-4
    smallei= 1.e-3
    tune_zeq     = 1.05
    tune_zeq_bnd = 0.80
 /
 
 $GRAVITY
    r_gc   = 8.5
    nsub   = 10
    h_gravity_profile =300.
    n_gravity_profile =3
 /

 $THERMAL
    cool_model      = 'c1'
    heat_model      = 'sup1'
    coolheat_active = 'yes'    
    cfl_coolheat  = 0.005
    esrc_upper_lim =  1.e-1
    esrc_lower_lim = -3.e-1
    h_coolheat_profile =300.
    n_coolheat_profile =3
    G_uv1         = 1.0e-24
    G_sup1        = 1.0e-25
    L_C 	  = -1.23848d-55 
    C_heatcond      = 1.e-2
    K_heatcond      = 0.0
    cfl_heatcond    = 0.4
 /
 
 
 
 $PROBLEM_CONTROL
    problem_name = 'parker_thermal'
    run_id =  'tst'
    T0     =  6000.0
    nbx0    =  0.
    nby0    = 1.0
    nbz0    =  0.
    h_pert  = 600.0
    a_pert  = 0.0
 /
