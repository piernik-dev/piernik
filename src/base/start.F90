#include "mhd.def"

module start

! Written by: M. Hanasz, January 2006

  use mpi_setup

  implicit none

  integer nxd, nyd, nzd, nb

  real t,dt
  real dt_mhd, dt_coolheat, dt_visc
  integer nstep
  integer nstep_start

  real tend
  integer nend

  real :: omega, qshear          
  
  character restart*16, new_id*3
  integer   nrestart, resdel

  character domain*16,mag_center*16
  integer, parameter      :: nvarsmx =16
  character, dimension(nvarsmx):: vars*4
  real dt_hdf, dt_res, dt_tsl, dt_log
  integer min_disk_space_MB, sleep_minutes, sleep_seconds
  character*32  user_message_file, system_message_file

  real xmin, xmax, ymin, ymax, zmin, zmax

  real c_si, gamma, alpha 

  real cfl, smalld, smallei, nu_bulk, cfl_visc
  real tune_zeq, tune_zeq_bnd 
  character*16 flux_limiter, freezing_speed, dimensions, magnetic
  integer integration_order, istep
  real rorder
  logical bulk_viscosity
  
  character*16 grav_model 
  real g_z
  real dg_dz
  real r_gc
  real ptmass, ptm_x, ptm_y, ptm_z, r_smooth 
  integer nsub
  real    h_grav,  r_grav
  integer n_gravh, n_gravr,n_gravr2   

  character cool_model*16, heat_model*16, coolheat_active*3,substepping*8
  real G_uv1, G_sup1, cfl_coolheat, esrc_lower_lim, esrc_upper_lim
  real h_coolheat_profile
  integer n_coolheat_profile    
  real L_C, K_heatcond, C_heatcond, dt_heatcond, cfl_heatcond
  logical gravaccel, coolheat, heatcond

  logical magfield, resist
  real    cfl_resist, eta_0, eta_1, j_crit, deint_max
  integer eta_scale

  real  cr_active, gamma_cr, cr_eff, beta_cr, K_cr_paral, K_cr_perp, &
        cfl_cr, amp_cr, smallecr
  real  dt_cr

  real h_sn, r_sn, f_sn_kpc2, t_dw, t_arm, col_dens

! Secondary parameters

  real csi2, csim2, amp_ecr_sn, ethu, f_sn

  real init_mass, mass_loss, mass_loss_tot
  
#ifdef SPLIT
#ifdef ORIG
  real, dimension(2)  :: cn
#endif
#ifdef SSP
  real, dimension(2,3)  :: cn
#endif
#endif

!-------------------------------------------------------------------------------
contains


  subroutine read_params
  
    implicit none

    integer iv
    
  namelist /DOMAIN_SIZES/ nxd, nyd, nzd, nb
    nxd    = 10
    nyd    = 10
    nzd    = 10
    nb     =  4
      
  namelist /START_CONTROL/ nstep, t, dt
    nstep  = 0
    t      = 0.0
    dt     = 0.0

  namelist /RESTART_CONTROL/ restart, new_id, nrestart, resdel
    restart = 'last'   ! 'last': autom. wybor ostatniego 
                       ! niezaleznie od wartosci "nrestart"
		       ! cokolwiek innego: decyduje "nrestart" 
    new_id  = ''
    nrestart=  3
    resdel  = 0

  namelist /END_CONTROL/ tend, nend
    tend   = 1.0
    nend   = 10 


  namelist /OUTPUT_CONTROL/ dt_hdf, dt_res, dt_tsl, dt_log, & 
                            domain, vars, mag_center, &
                            min_disk_space_MB, sleep_minutes, sleep_seconds, &
                            user_message_file, system_message_file
                            
    dt_hdf = 0.0
    dt_res = 0.0
    dt_tsl = 0.0
    dt_log = 0.0
    domain = 'phys_domain'
    vars(:)   = '    '
    mag_center= 'no'
    min_disk_space_MB = 100
    sleep_minutes   = 0
    sleep_seconds   = 0
    user_message_file   = './msg'
    system_message_file = '/etc/ups/user/msg'
      
  namelist /DOMAIN_LIMITS/ xmin, xmax, ymin, ymax, zmin, zmax  
    xmin   = 0.0
    xmax   = 1.0
    ymin   = 0.0
    ymax   = 1.0
    zmin   = 0.0
    zmax   = 1.0
 
  namelist /EQUATION_OF_STATE/ c_si, gamma, alpha
    c_si   = 1.0
    gamma  = 5./3.    ! ignored if ISO
#ifdef ISO
!    gamma  = 1.
#endif    
    alpha  = 1.0   

  namelist /NUMERICAL_SETUP/  cfl, smalld, smallei, &
                              flux_limiter, freezing_speed, &
                              integration_order, &
                              dimensions, magnetic, nu_bulk, cfl_visc

    cfl     = 0.7
    smalld  = 1.e-10
    smallei = 1.e-10
    flux_limiter = 'vanleer'
    freezing_speed = 'local'
    integration_order  = 2
    dimensions = '3d'
    magnetic  = 'yes'
    nu_bulk   = 0.0
    cfl_visc  = 0.4

#ifdef GRAV
  namelist /GRAVITY/ grav_model,  &
                     g_z,   &
                     dg_dz, &
                     r_gc,  &
                     ptmass,ptm_x,ptm_y,ptm_z,r_smooth, &
                     nsub, tune_zeq, tune_zeq_bnd,      &
                     h_grav, r_grav, n_gravr, n_gravr2, n_gravh

    grav_model  = 'null'                     
    g_z     = 0.0 
    dg_dz   = 0.0
    r_gc    = 8500
    ptmass  = 0.0
    ptm_x   = 0.0
    ptm_y   = 0.0
    ptm_z   = 0.0
    r_smooth= 0.0
    nsub    = 10
    tune_zeq     = 1.0
    tune_zeq_bnd = 1.0
    h_grav = 1.e6
    r_grav = 1.e6
    n_gravr = 0
    n_gravr2= 0
    n_gravh = 0
#endif

#ifdef COOL_HEAT
  namelist /THERMAL/ cool_model, heat_model, coolheat_active,  &
                     esrc_lower_lim, esrc_upper_lim, &
                     h_coolheat_profile, n_coolheat_profile, &
                     G_uv1, G_sup1, cfl_coolheat, L_C, &
                     K_heatcond, C_heatcond, cfl_heatcond 

    cool_model   = 'null'
    heat_model   = 'null'
    coolheat_active = 'no'
    esrc_upper_lim =  1.
    esrc_lower_lim = -1.
    h_coolheat_profile = 1.e6
    n_coolheat_profile = 0
    G_uv1        = 0.0
    G_sup1       = 0.0
    cfl_coolheat = 0.005
    L_C          = 0.0
    K_heatcond    = 0.0
    C_heatcond    = 1.5e-5
    cfl_heatcond = 0.45
#endif

#ifdef RESIST
  namelist /RESISTIVITY/ cfl_resist, eta_0, eta_1, eta_scale, j_crit, deint_max
    cfl_resist  =  0.4
    eta_0       =  0.0
    eta_1       =  0.0
    eta_scale   =  4
    j_crit      =  1.0e6
    deint_max   = 0.01
#endif

#ifdef COSM_RAYS
  namelist /COSMIC_RAYS/ cr_active, gamma_cr, cr_eff, beta_cr, &
                         K_cr_paral, K_cr_perp, amp_cr, cfl_cr, smallecr
    cr_active  = 1.0
    gamma_cr   = 4./3.
    beta_cr    = 0.0
    cr_eff     = 0.1       !  canonical conversion rate of SN en.-> CR
                           !  we fix E_SN=10**51 erg
    K_cr_paral = 0.0
    K_cr_perp  = 0.0
    amp_cr     = 0.0
    cfl_cr     = 0.45
    smallecr   = 0.0
#endif

#ifdef GALAXY
  namelist /GALACTIC_PARAMS/ h_sn, r_sn, f_sn_kpc2, t_dw, t_arm,col_dens 
! Default galactic parameters
    h_sn      = 266.0		!  uncertain vertical scaleheight of SN
                                !  from Ferriere 1998
    r_sn       =  10.0       	!  "typical" SNR II radius
    f_sn_kpc2  =  20.0       	!  solar galactic radius SN II freq./kpc**2
    t_dw       = 100.0        	! period of density waves in Myr
    t_arm      = 100.0        	! period of SFR in arms
    col_dens    = 0.0            ! gas column density
#endif GALAXY

#ifdef SHEAR
  namelist /SHEARING/ omega, qshear
    omega  = 0.0
    qshear = 0.0
#endif SHEAR

         
    if(proc .eq. 0) then
          
      open(1,file='problem.par')
        read(unit=1,nml=DOMAIN_SIZES)  
        read(unit=1,nml=START_CONTROL)  
        read(unit=1,nml=RESTART_CONTROL)  
        read(unit=1,nml=END_CONTROL)  
        read(unit=1,nml=OUTPUT_CONTROL)  
        read(unit=1,nml=DOMAIN_LIMITS)
        read(unit=1,nml=EQUATION_OF_STATE)
        read(unit=1,nml=NUMERICAL_SETUP)
#ifdef GRAV	
        read(unit=1,nml=GRAVITY)
#endif
#ifdef COOL_HEAT
        read(unit=1,nml=THERMAL)
#endif
#ifdef RESIST
        read(unit=1,nml=RESISTIVITY)
#endif
#ifdef COSM_RAYS
        read(unit=1,nml=COSMIC_RAYS)
#endif
#ifdef GALAXY
        read(unit=1,nml=GALACTIC_PARAMS)
#endif
#ifdef SHEAR
        read(unit=1,nml=SHEARING)
#endif SHEAR

      close(1)

      open(3, file='tmp.log', position='append')
        write(unit=3,nml=DOMAIN_SIZES)  
        write(unit=3,nml=START_CONTROL)  
        write(unit=3,nml=RESTART_CONTROL)  
        write(unit=3,nml=END_CONTROL)  
        write(unit=3,nml=OUTPUT_CONTROL) 
        write(unit=3,nml=DOMAIN_LIMITS)
        write(unit=3,nml=EQUATION_OF_STATE)
        write(unit=3,nml=NUMERICAL_SETUP)
#ifdef GRAV	
        write(unit=3,nml=GRAVITY)
#endif 	
#ifdef COOL_HEAT
        write(unit=3,nml=THERMAL)
#endif
#ifdef RESIST
        write(unit=3,nml=RESISTIVITY)
#endif
#ifdef COSM_RAYS
        write(unit=3,nml=COSMIC_RAYS)
#endif
#ifdef GALAXY
        write(unit=3,nml=GALACTIC_PARAMS)
#endif 
#ifdef SHEAR
        write(unit=3,nml=SHEARING)
#endif SHEAR
      close(3)
      
    endif 
    
    
    if(proc .eq. 0) then   

!  namelist /DOMAIN_SIZES/ nxd, nyd, nzd, nb
  
      ibuff(1)  = nxd 
      ibuff(2)  = nyd 
      ibuff(3)  = nzd 
      ibuff(4)  = nb 

!  namelist /START_CONTROL/ nstep, t, dt

      ibuff(10) = nstep

      rbuff(10) = t      
      rbuff(11) = dt   

!  namelist /RESTART_CONTROL/ restart, new_id, nrestart, resdel

      cbuff(20) = restart 
      cbuff(21) = new_id

      ibuff(20) = nrestart
      ibuff(21) = resdel

!  namelist /END_CONTROL/ nend, tend 

      ibuff(30) = nend

      rbuff(30) = tend  

!  namelist /OUTPUT_CONTROL/ dt_hdf, dt_res, dt_tsl, domain, vars, mag_center, &
!                            min_disk_space_MB, sleep_minutes, &
!                            user_message_file, system_message_file

      ibuff(40) = min_disk_space_MB
      ibuff(41) = sleep_minutes
      ibuff(42) = sleep_seconds
                                  
      rbuff(40) = dt_hdf 
      rbuff(41) = dt_res 
      rbuff(42) = dt_tsl 
      rbuff(43) = dt_log 

      cbuff(40) = domain 

      do iv = 1, nvarsmx 
        cbuff(40+iv) = vars(iv) 
      enddo
      
      cbuff(60) = mag_center
      cbuff(61) = user_message_file      
      cbuff(62) = system_message_file
            
!  namelist /DOMAIN_LIMITS/ xmin, xmax, ymin, ymax, zmin, zmax  

      rbuff(50) = xmin 
      rbuff(51) = xmax 
      rbuff(52) = ymin 
      rbuff(53) = ymax 
      rbuff(54) = zmin 
      rbuff(55) = zmax 
 
!  namelist /EQUATION_OF_STATE/ c_si, gamma, alpha
      rbuff(70) = c_si  
      rbuff(71) = gamma 
      rbuff(72) = alpha   

!  namelist /NUMERICAL_SETUP/  cfl, smalld, smallei,
!                              flux_limiter, freezing_speed,
!                              integration_order, 
!                              dimensions, magnetic, nu_bulk, cfl_visc

      rbuff(80) = cfl   
      rbuff(83) = smalld
      rbuff(84) = smallei
      rbuff(85) = nu_bulk
      rbuff(86) = cfl_visc

      cbuff(80) = flux_limiter
      cbuff(81) = freezing_speed
      cbuff(82) = dimensions
      cbuff(83) = magnetic
 
      ibuff(80) = integration_order

#ifdef GRAV
!  namelist /GRAVITY/ grav_model, g_z, dg_dz, r_gc, 
!                     ptmass,ptm_x,ptm_y,ptm_z,r_smooth, nsub, 
!                     tune_zeq, tune_zeq_bnd,
!                     h_grav, r_grav, n_gravr, n_gravr2, n_gravh

      ibuff(90) = nsub
      ibuff(91) = n_gravr
      ibuff(92) = n_gravr2
      ibuff(93) = n_gravh

      cbuff(90) = grav_model                 

      rbuff(90)  = g_z    
      rbuff(91)  = dg_dz  
      rbuff(92)  = r_gc   
      rbuff(93)  = ptmass 
      rbuff(94)  = ptm_x  
      rbuff(95)  = ptm_y   
      rbuff(96)  = ptm_z   
      rbuff(97)  = tune_zeq    
      rbuff(98)  = tune_zeq_bnd 
      rbuff(99)  = r_smooth
      rbuff(100) = h_grav
      rbuff(101) = r_grav
#endif

#ifdef COOL_HEAT   
!  namelist /THERMAL/ cool_model, heat_model, coolheat_active, &
!                     esrc_lower_lim, esrc_upper_lim, &
!                     h_coolheat_profile, n_coolheat_profile, &
!                     G_uv1, G_sup1, cfl_coolheat, &
!                     K_heatcond, C_heatcond, cfl_heatcond 
 
       ibuff(110) = n_coolheat_profile  

       cbuff(110) = cool_model
       cbuff(111) = heat_model
       cbuff(112) = coolheat_active

       rbuff(110) = cfl_coolheat
       rbuff(111) = esrc_lower_lim 
       rbuff(112) = esrc_upper_lim 
       rbuff(113) = h_coolheat_profile  

       rbuff(114) = G_uv1
       rbuff(115) = G_sup1
       rbuff(116) = L_C
       rbuff(117) = K_heatcond  
       rbuff(118) = C_heatcond  
       rbuff(119) = cfl_heatcond  
#endif

#ifdef RESIST
!   namelist /RESISTIVITY/ cfl_resist, eta_0, eta_1, j_crit

       ibuff(120) = eta_scale

       rbuff(120) = cfl_resist
       rbuff(121) = eta_0
       rbuff(122) = eta_1
       rbuff(123) = j_crit
       rbuff(124) = deint_max
#endif
             
#ifdef COSM_RAYS
!  namelist /COSMIC_RAYS/ cr_active, gamma_cr, cr_eff, beta_cr, K_cr_paral, K_cr_perp,&
!  			 amp_cr, cfl_cr
       rbuff(130) = cr_active 
       rbuff(131) = gamma_cr 
       rbuff(132) = cr_eff       
       rbuff(133) = beta_cr 
       rbuff(134) = K_cr_paral 
       rbuff(135) = K_cr_perp  
       rbuff(136) = amp_cr    
       rbuff(137) = cfl_cr   
       rbuff(138) = smallecr 
#endif

#ifdef GALAXY
!  namelist /GALACTIC_PARAMS/ h_sn, r_sn, f_sn_kpc2, t_dw, t_arm, col_dens

       rbuff(140) = h_sn  		
       rbuff(141) = r_sn          	
       rbuff(142) = f_sn_kpc2       	
       rbuff(143) = t_dw             	
       rbuff(144) = t_arm   
       rbuff(145) = col_dens              
#endif GALAXY

#ifdef SHEAR
!  namelist /SHEARING/ omega, qshear
       rbuff(160) = omega
       rbuff(161) = qshear
#endif



! Boroadcasting parameters

      call MPI_BCAST(cbuff, 32*buffer_dim, MPI_CHARACTER,        0, comm, ierr)
      call MPI_BCAST(ibuff,    buffer_dim, MPI_INTEGER,          0, comm, ierr)
      call MPI_BCAST(rbuff,    buffer_dim, MPI_DOUBLE_PRECISION, 0, comm, ierr)

    else
    
      call MPI_BCAST(cbuff, 32*buffer_dim, MPI_CHARACTER,        0, comm, ierr)
      call MPI_BCAST(ibuff,    buffer_dim, MPI_INTEGER,          0, comm, ierr)
      call MPI_BCAST(rbuff,    buffer_dim, MPI_DOUBLE_PRECISION, 0, comm, ierr)
      
!  namelist /DOMAIN_SIZES/ nxd, nyd, nzd, nb

      nxd                 = ibuff(1)  
      nyd                 = ibuff(2)  
      nzd                 = ibuff(3)  
      nb                  = ibuff(4)  

!  namelist /START_CONTROL/ nstep, t, dt

      nstep               = ibuff(10) 

      t                   = rbuff(10)      
      dt                  = rbuff(11)   

!  namelist /RESTART_CONTROL/ restart, new_id, nrestart, resdel

      restart             = trim(cbuff(20))  
      new_id              = trim(cbuff(21))

      nrestart            = ibuff(20) 
      resdel              = ibuff(21)

!  namelist /END_CONTROL/ nend, tend 

      nend                = ibuff(30)

      tend                = rbuff(30) 

!  namelist /OUTPUT_CONTROL/ dt_hdf, dt_res, dt_tsl, domain, vars, mag_center, &
!                            min_disk_space_MB, sleep_minutes, &
!                            user_message_file, system_message_file

      min_disk_space_MB   = ibuff(40)
      sleep_minutes       = ibuff(41) 
      sleep_seconds       = ibuff(42) 
                                  
                            
      dt_hdf              = rbuff(40)
      dt_res              = rbuff(41) 
      dt_tsl              = rbuff(42)  
      dt_log              = rbuff(43) 

      domain              = trim(cbuff(40)) 
      do iv=1, nvarsmx
        vars(iv)          = trim(cbuff(40+iv))
      enddo

      mag_center          = trim(cbuff(60)) 
      
      user_message_file   = trim(cbuff(61))    
      system_message_file = trim(cbuff(62))

!  namelist /DOMAIN_LIMITS/ xmin, xmax, ymin, ymax, zmin, zmax  

      xmin                = rbuff(50) 
      xmax                = rbuff(51) 
      ymin                = rbuff(52) 
      ymax                = rbuff(53)  
      zmin                = rbuff(54) 
      zmax                = rbuff(55) 


!  namelist /EQUATION_OF_STATE/ c_si, gamma, alpha

      c_si                = rbuff(70) 
      gamma               = rbuff(71)
      alpha               = rbuff(72)  

!  namelist /NUMERICAL_SETUP/  cfl, smalld, smallei,
!                              flux_limiter, freezing_speed,    
!                              integration_order,
!                              dimensions, magnetic, nu_bulk, cfl_visc

      cfl                 = rbuff(80)   
      smalld              = rbuff(83)
      smallei             = rbuff(84)
      nu_bulk             = rbuff(85) 
      cfl_visc            = rbuff(86) 


      flux_limiter        = trim(cbuff(80)) 
      freezing_speed      = trim(cbuff(81))
      dimensions          = cbuff(82)  
      magnetic            = cbuff(83) 
      
      integration_order   = ibuff(80)
 
#ifdef GRAV
!  namelist /GRAVITY/ grav_model,  g_z, dg_dz, r_gc, 
!                     ptmass,ptm_x,ptm_y,ptm_z,r_smooth
!                     tune_zeq, tune_zeq_bnd,
!                     h_gravity_profile, r_gravity_profile, n_gravity_profile
!                     h_grav, r_grav, n_gravr, n_gravh

      nsub                = ibuff(90)
      n_gravr             = ibuff(91) 
      n_gravr2            = ibuff(92)
      n_gravh             = ibuff(93) 

      grav_model          = trim(cbuff(90))                  

      g_z                 = rbuff(90)    
      dg_dz               = rbuff(91) 
      r_gc                = rbuff(92)  
      ptmass              = rbuff(93) 
      ptm_x               = rbuff(94) 
      ptm_y               = rbuff(95)  
      ptm_z               = rbuff(96)    
      tune_zeq            = rbuff(97)   
      tune_zeq_bnd        = rbuff(98) 
      r_smooth            = rbuff(99) 
      h_grav              = rbuff(100) 
      r_grav              = rbuff(101) 
#endif GRAV
    
#ifdef COOL_HEAT
!  namelist /THERMAL/ cool_model, heat_model, coolheat_active,  &
!                     esrc_lower_lim, esrc_upper_lim, &
!                     h_coolheat_profile, n_coolheat_profile, &
!                     G_uv1, G_sup1, cfl_coolheat, &
!                     K_heatcond, C_heatcond, cfl_heatcond  


      n_coolheat_profile  = ibuff(110) 

      cool_model          = cbuff(110)
      heat_model          = cbuff(111)  
      coolheat_active     = cbuff(112)

      cfl_coolheat        = rbuff(110)  
      esrc_lower_lim      = rbuff(111) 
      esrc_upper_lim      = rbuff(112) 
      h_coolheat_profile  = rbuff(113) 
      G_uv1               = rbuff(114)  
      G_sup1              = rbuff(115)  
      L_C                 = rbuff(116) 
      K_heatcond          = rbuff(117) 
      C_heatcond          = rbuff(118) 
      cfl_heatcond        = rbuff(119) 
#endif COOL_HEAT       

#ifdef RESIST
!   namelist /RESISTIVITY/ cfl_resist, eta_0, eta_1, j_crit

       eta_scale          = ibuff(120)

       cfl_resist         = rbuff(120)  
       eta_0              = rbuff(121)  
       eta_1              = rbuff(122)  
       j_crit             = rbuff(123)  
       deint_max          = rbuff(124)
      
#endif RESIST
             
#ifdef COSM_RAYS
!  namelist /COSMIC_RAYS/ cr_active, gamma_cr, cr_eff, beta_cr, K_cr_paral, K_cr_perp,&
!  			 amp_cr, cfl_cr

       cr_active 	  = rbuff(130) 
       gamma_cr           = rbuff(131) 
       cr_eff             = rbuff(132) 
       beta_cr            = rbuff(133)
       K_cr_paral         = rbuff(134) 
       K_cr_perp          = rbuff(135)
       amp_cr             = rbuff(136) 
       cfl_cr             = rbuff(137)  
       smallecr           = rbuff(138)   
#endif COSM_RAYS


#ifdef GALAXY
!  namelist /GALACTIC_PARAMS/ h_sn, r_sn, f_sn_kpc2, t_dw, t_arm 

       h_sn 		  = rbuff(140)  		
       r_sn 		  = rbuff(141)         	
       f_sn_kpc2 	  = rbuff(142)        	
       t_dw 		  = rbuff(143)             	
       t_arm 		  = rbuff(144)          
       col_dens           = rbuff(145)              
#endif GALAXY

#ifdef SHEAR
!  namelist /SHEARING/ omega, qshear
       omega              = rbuff(160)
       qshear             = rbuff(161)
#endif SHEAR
    endif  ! (proc .eq. 0)    

! Secondary parameters

   csi2  = c_si**2
   csim2 = csi2*(1.+alpha)    ! z-equilibrium defined for fixed
                                ! ratio p_mag/p_gas = alpha = 1/beta
#ifdef COSM_RAYS
   csim2 = csim2 +csi2*beta_cr
#endif COSM_RAYS

   ethu = 7.0**2/(5.0/3.0-1.0) * 1.0    ! thermal energy unit=0.76eV/cm**3
                                        ! for c_si= 7km/s, n=1/cm^3 	 
                                        ! gamma=5.3
#ifdef GALAXY
   amp_ecr_sn = 4.96e6*cr_eff/r_sn**3   ! cosmic ray explosion amplitude 
                                        ! in units: 
        				! e_0 = 1/(5/3-1)*rho_0*c_s0**2
        				! rho_0=1.67e-24g/cm**3,
        				! c_s0 = 7km/s
   f_sn = f_sn_kpc2 * (xmax-xmin)/1000.0 * (ymax-ymin)/1000.0 ! SN frequency per horizontal
                                              ! surface area of the comp. box 
#endif GALAXY

#ifdef SPLIT
#ifdef ORIG
  cn(1:2) = (/0.5 , 1./);
  if(integration_order .eq. 1) then
    cn(1) = 1.
  endif
#endif
#ifdef SSP
  cn(1:2,1) = (/ 1.   , 1.   /)
  cn(1:2,2) = (/ 0.75 , 0.25 /)
  cn(1:2,3) = (/ 1./3., 2./3./)

  if(integration_order .eq. 2) then
    cn(1:2,2) = (/ 0.5 , 0.5 /)
  endif
#endif
#endif


!-------------------------                               
#ifdef ORIG
    if(integration_order .gt. 2) then
      stop 'For "ORIG" scheme integration_order must be 1 or 2'
    endif
#endif
#ifdef SSP
    if(integration_order .gt. 3) then
      stop 'For "SSP" scheme integration_order can`t be greater than 3'
    endif
#endif

    select case (dimensions)
      case('3d','2dxy')
        !do nothing
      case default
        stop '"dimensions" must be one of the following: "3d","2dxy"'
    end select 


#ifdef GRAV
      gravaccel = .true.
#else GRAV
      gravaccel = .false.
#endif GRAV   

#ifdef COOL_HEAT
    if((cool_model .ne. 'null') .or. (heat_model .ne. 'null')) then
      coolheat = .true.
#else COOL_HEAT
      coolheat = .false.
#endif COOL_HEAT
    
#ifdef HEAT_COND
      heatcond = .true.
#else HEAT_COND
      heatcond = .false.
#endif HEAT_COND

    if(magnetic .eq. 'yes') then
      magfield = .true.
    else
      magfield = .false.
    endif
   
#ifdef RESIST
      resist = .true.
      if(eta_scale .lt. 0) then
        write(*,*) 'eta_scale must be greater or equal 0'
        stop
      endif
#else RESIST
      resist = .false.
#endif
   
#ifdef VISC
    if(nu_bulk .ne. 0.0) then
       bulk_viscosity = .true.
#else VISC
       bulk_viscosity = .false.
#endif VISC
       
#ifdef MASS_COMPENS
    mass_loss = 0.0
    mass_loss_tot = 0.0
#endif MASS_COMPENS       

  end subroutine read_params
  
end module start








