! $Id$
#include "piernik.def"

!#ifdef IONIZED
!#define NUMBION IONIZED
!#else /* IONIZED */
!#define NUMBION 0
!#endif /* IONIZED */
!#ifdef NEUTRAL
!#define NUMBNEUT NEUTRAL
!#else /* NEUTRAL */
!#define NUMBNEUT 0
!#endif /* NEUTRAL */
!#ifdef DUST
!#define NUMBDUST DUST
!#else /* DUST */
!#define NUMBDUST 0
!#endif /* DUST */
!#define NUMBFLUID NUMBION+NUMBNEUT+NUMBDUST

module start

! Written by: M. Hanasz, January 2006

  use mpi_setup


  implicit none


  real t,dt
  real collfaq, cfl_colls
  integer nstep
  integer nstep_start

  real tend
  integer nend, maxxyz

  real :: omega, qshear

  real cfl, smalld, smallei, nu_bulk, cfl_visc
#ifdef VZ_LIMITS
  real   :: floor_vz, ceil_vz
#endif /* VZ_LIMITS */

  real tune_zeq, tune_zeq_bnd
  character*16 flux_limiter, freezing_speed, dimensions, magnetic
  integer integration_order, istep
  real rorder

  real, dimension(1) :: gamma   !!! do poprawy

  character*3 gpt_hdf
  real :: g_z, g_y
  real dg_dz
  real r_gc
  real ptmass, ptm_x, ptm_y, ptm_z, r_smooth
  integer nsub
  real    h_grav,  r_grav
  integer n_gravh, n_gravr, n_gravr2

  real G_uv1, G_sup1, cfl_coolheat, esrc_lower_lim, esrc_upper_lim
  logical gravaccel

  logical magfield, resist
  real    cfl_resist, eta_0, eta_1, j_crit, deint_max
  integer eta_scale

!  real  cr_active, gamma_cr, cr_eff, beta_cr, K_cr_paral, K_cr_perp, &
!        cfl_cr, amp_cr, smallecr
  real  dt_cr, dt_colls, dt_supp

  real t_dw, t_arm, col_dens

  real h_sn, r_sn, f_sn_kpc2, amp_dip_sn, snenerg, snemass, sn1time, sn2time, r0sn
  integer howmulti
  character*3 add_mass, add_ener, add_encr, add_magn
! Secondary parameters

  real csi2, csim2, amp_ecr_sn, ethu, f_sn

  real init_mass, mass_loss, mass_loss_tot

  real, dimension(3,2)  :: cn

!-------------------------------------------------------------------------------
contains


  subroutine read_params

    implicit none


  character par_file*(100), tmp_log_file*(100)

  namelist /NUMERICAL_SETUP/  cfl, smalld, smallei, &
#ifdef VZ_LIMITS
                              floor_vz, ceil_vz, &
#endif /* VZ_LIMITS */
#ifdef COLLISIONS
                              cfl_colls, &
#endif /* COLLISIONS */
                              integration_order, &
                              dimensions, magnetic 

#ifdef SHEAR
  namelist /SHEARING/ omega, qshear
#endif /* SHEAR */
#ifdef GRAV
  namelist /GRAVITY/ gpt_hdf,  &
                     g_z,   &
                     g_y,   &
                     dg_dz, &
                     r_gc,  &
                     ptmass,ptm_x,ptm_y,ptm_z,r_smooth, &
                     nsub, tune_zeq, tune_zeq_bnd,      &
                     h_grav, r_grav, n_gravr, n_gravr2, n_gravh
#endif /* GRAV */
#ifdef RESISTIVE
  namelist /RESISTIVITY/ cfl_resist, eta_0, eta_1, eta_scale, j_crit, deint_max
#endif /* RESISTIVE */
!#ifdef COSM_RAYS
!  namelist /COSMIC_RAYS/ cr_active, gamma_cr, cr_eff, beta_cr, &
!                         K_cr_paral, K_cr_perp, amp_cr, cfl_cr, smallecr
!#endif /* COSM_RAYS */
#ifdef GALAXY
  namelist /GALACTIC_PARAMS/ t_dw, t_arm,col_dens
#endif /* GALAXY */
#ifdef SHEAR
  namelist /SHEARING/ omega, qshear
#endif /* SHEAR */
#ifdef SN_SRC
  namelist /SN_PARAMS/ h_sn, r_sn, f_sn_kpc2, amp_dip_sn, howmulti
#endif /* SN_SRC */
#ifdef SNE_DISTR
  namelist /SN_DISTR/ snenerg, snemass, sn1time, sn2time, r0sn, add_mass, add_ener, add_encr, add_magn
#endif /* SNE_DISTR */

    par_file = trim(cwd)//'/problem.par'
    tmp_log_file = trim(cwd)//'/tmp.log'

    nstep  = 0
    t      = 0.0
    dt     = 0.0

    tend   = 1.0
    nend   = 10

    cfl     = 0.7
    smalld  = 1.e-10
    smallei = 1.e-10
    flux_limiter = 'vanleer'
    freezing_speed = 'local'
    integration_order  = 2
    dimensions = '3d'
    magnetic  = 'yes'
#ifdef VZ_LIMITS
    floor_vz  = -1.e99
    ceil_vz   =  1.e99
#endif /* VZ_LIMITS */
#if defined COLLISIONS || defined KEPLER_SUPPRESSION
    cfl_colls = 0.01
#endif /* COLLISIONS || KEPLER_SUPPRESION */

#ifdef GRAV
    gpt_hdf = 'no'
    g_z     = 0.0
    g_y     = 0.0
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
#endif /* GRAV */

#ifdef RESISTIVE
    cfl_resist  =  0.4
    eta_0       =  0.0
    eta_1       =  0.0
    eta_scale   =  4
    j_crit      =  1.0e6
    deint_max   = 0.01
#endif /* RESISTIVE */

#ifdef COSM_RAYS
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
#endif /* COSM_RAYS */

#ifdef GALAXY
! Default galactic parameters
    t_dw       = 100.0           ! period of density waves in Myr
    t_arm      = 100.0           ! period of SFR in arms
    col_dens    = 0.0            ! gas column density
#endif /* GALAXY */

#ifdef SHEAR
    omega  = 0.0
    qshear = 0.0
#endif /* SHEAR */
#ifdef SN_SRC
    h_sn       = 266.0          !  vertical scaleheight of SN from Ferriere 1998
    r_sn       =  10.0          !  "typical" SNR II radius
    f_sn_kpc2  =  20.0          !  solar galactic radius SN II freq./kpc**2
    amp_dip_sn =   1.0e6
    howmulti   = 2              ! 1 for dipols, 2 for quadrupoles
#endif /* SN_SRC */
#ifdef SNE_DISTR
    snenerg    = 1.e51          !  typical energy of supernova explosion [erg]
    snemass    =  10.0          !  typical preSN stellar matter mass injection [Msun]
    sn1time    = 445.0          !  mean time between typ I supernovae explosions [year]
    sn2time    =  52.0          !  mean time between typ II supernovae explosions [year]
    r0sn       =  50.0          !  radius of an area, where mass/ener/encr is added [actually used unit of length]
    add_mass   = 'yes'          !  permission for inserting snemass inside randomly selected areas
    add_ener   = 'yes'          !  permission for inserting snenerg inside randomly selected areas
    add_encr   = 'yes'          !  permission for inserting CR energy inside randomly selected areas
    add_magn   = 'yes'          !  permission for inserting dipolar magnetic field centered at randomly selected areas
#endif /* SNE_DISTR */

    if(proc .eq. 0) then

      open(1,file=par_file)
        read(unit=1,nml=END_CONTROL)
        read(unit=1,nml=NUMERICAL_SETUP)
#ifdef GRAV
        read(unit=1,nml=GRAVITY)
#endif /* GRAV */
#ifdef RESISTIVE
        read(unit=1,nml=RESISTIVITY)
#endif /* RESISTIVE */
!#ifdef COSM_RAYS
!        read(unit=1,nml=COSMIC_RAYS)
!#endif /* COSM_RAYS */
#ifdef GALAXY
        read(unit=1,nml=GALACTIC_PARAMS)
#endif /* GALAXY */
#ifdef SHEAR
        read(unit=1,nml=SHEARING)
#endif /* SHEAR */
#ifdef SN_SRC
        read(unit=1,nml=SN_PARAMS)
#endif /* SN_SRC */
#ifdef SNE_DISTR
        read(unit=1,nml=SN_DISTR)
#endif /* SNE_DISTR */

      close(1)

      open(3, file=tmp_log_file, position='append')
        write(unit=3,nml=END_CONTROL)
        write(unit=3,nml=NUMERICAL_SETUP)
#ifdef GRAV
        write(unit=3,nml=GRAVITY)
#endif /* GRAV */
#ifdef RESISTIVE
        write(unit=3,nml=RESISTIVITY)
#endif /* RESISTIVE */
!#ifdef COSM_RAYS
!        write(unit=3,nml=COSMIC_RAYS)
!#endif /* COSM_RAYS */
#ifdef GALAXY
        write(unit=3,nml=GALACTIC_PARAMS)
#endif /* GALAXY */
#ifdef SHEAR
        write(unit=3,nml=SHEARING)
#endif /* SHEAR */
#ifdef SN_SRC
        write(unit=3,nml=SN_PARAMS)
#endif /* SN_SRC */
#ifdef SNE_DISTR
        write(unit=3,nml=SN_DISTR)
#endif /* SNE_DISTR */
      close(3)

    endif


    if(proc .eq. 0) then

!  namelist /END_CONTROL/ nend, tend

      ibuff(30) = nend

      rbuff(30) = tend

!      do iv = 1, NUMBFLUID
!        rbuff(72+iv) = gamma(iv)
!      enddo
!
!  namelist /NUMERICAL_SETUP/  cfl, smalld, smallei,
!                              flux_limiter, freezing_speed,
!                              integration_order,
!                              dimensions, magnetic, nu_bulk, cfl_visc
!                              floor_vz, ceil_vz, cfl_colls
      rbuff(80) = cfl
      rbuff(83) = smalld
      rbuff(84) = smallei
#ifdef VZ_LIMITS
      rbuff(87) = floor_vz
      rbuff(88) = ceil_vz
#endif /* VZ_LIMITS */
#ifdef COLLISIONS
      rbuff(89) = cfl_colls
#endif /* COLLISIONS */

      cbuff(80) = flux_limiter
      cbuff(81) = freezing_speed
      cbuff(82) = dimensions
      cbuff(83) = magnetic

      ibuff(80) = integration_order

#ifdef GRAV
!  namelist /GRAVITY/ gpt_hdf, g_z, g_y, dg_dz, r_gc,
!                     ptmass,ptm_x,ptm_y,ptm_z,r_smooth, nsub,
!                     tune_zeq, tune_zeq_bnd,
!                     h_grav, r_grav, n_gravr, n_gravr2, n_gravh

      ibuff(90) = nsub
      ibuff(91) = n_gravr
      ibuff(92) = n_gravr2
      ibuff(93) = n_gravh

      cbuff(90) = gpt_hdf

      rbuff(90)  = g_z
      rbuff(185) = g_y
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
#endif /* GRAV */

#ifdef RESISTIVE
!   namelist /RESISTIVITY/ cfl_resist, eta_0, eta_1, j_crit

       ibuff(120) = eta_scale

       rbuff(120) = cfl_resist
       rbuff(121) = eta_0
       rbuff(122) = eta_1
       rbuff(123) = j_crit
       rbuff(124) = deint_max
#endif /* RESISTIVE */

!#ifdef COSM_RAYS
!  namelist /COSMIC_RAYS/ cr_active, gamma_cr, cr_eff, beta_cr, K_cr_paral, K_cr_perp,&
!                         amp_cr, cfl_cr
!       rbuff(130) = cr_active
!       rbuff(131) = gamma_cr
!       rbuff(132) = cr_eff
!       rbuff(133) = beta_cr
!       rbuff(134) = K_cr_paral
!       rbuff(135) = K_cr_perp
!       rbuff(136) = amp_cr
!       rbuff(137) = cfl_cr
!       rbuff(138) = smallecr
!#endif /* COSM_RAYS */

#ifdef GALAXY
!  namelist /GALACTIC_PARAMS/ h_sn, r_sn, f_sn_kpc2, t_dw, t_arm, col_dens

       rbuff(140) = h_sn
       rbuff(141) = r_sn
       rbuff(142) = f_sn_kpc2
       rbuff(143) = t_dw
       rbuff(144) = t_arm
       rbuff(145) = col_dens
#endif /* GALAXY */

#ifdef SHEAR
!  namelist /SHEARING/ omega, qshear
       rbuff(160) = omega
       rbuff(161) = qshear
#endif /* SHEAR */

#ifdef SN_SRC
!  namelist /SN_PARAMS/ h_sn, r_sn, f_sn_kpc2, amp_dip_sn, howmulti
       rbuff(170) = h_sn
       rbuff(171) = r_sn
       rbuff(172) = f_sn_kpc2
       rbuff(173) = amp_dip_sn
       ibuff(170) = howmulti
#endif /* SN_SRC */
#ifdef SNE_DISTR
!  namelist /SN_DISTR/ snenerg, snemass, sn1time, sn2time, r0sn, add_mass, add_ener, add_encr, add_magn
       rbuff(180) = snenerg
       rbuff(181) = snemass
       rbuff(182) = sn1time
       rbuff(183) = sn2time
       rbuff(184) = r0sn
       cbuff(180) = add_mass
       cbuff(181) = add_ener
       cbuff(182) = add_encr
       cbuff(183) = add_magn
#endif /* SNE_DISTR */

! Boroadcasting parameters

      call MPI_BCAST(cbuff, 32*buffer_dim, MPI_CHARACTER,        0, comm, ierr)
      call MPI_BCAST(ibuff,    buffer_dim, MPI_INTEGER,          0, comm, ierr)
      call MPI_BCAST(rbuff,    buffer_dim, MPI_DOUBLE_PRECISION, 0, comm, ierr)

    else

      call MPI_BCAST(cbuff, 32*buffer_dim, MPI_CHARACTER,        0, comm, ierr)
      call MPI_BCAST(ibuff,    buffer_dim, MPI_INTEGER,          0, comm, ierr)
      call MPI_BCAST(rbuff,    buffer_dim, MPI_DOUBLE_PRECISION, 0, comm, ierr)

!  namelist /END_CONTROL/ nend, tend

      nend                = ibuff(30)

      tend                = rbuff(30)

!      do iv=1,NUMBFLUID
!        gamma(iv)         = rbuff(72+iv)
!      enddo

!  namelist /NUMERICAL_SETUP/  cfl, smalld, smallei,
!                              flux_limiter, freezing_speed,
!                              integration_order,
!                              dimensions, magnetic, nu_bulk, cfl_visc
!                              floor_vz, ceil_vz , cfl_colls

      cfl                 = rbuff(80)
      smalld              = rbuff(83)
      smallei             = rbuff(84)
#ifdef VZ_LIMITS
      floor_vz            = rbuff(87)
      ceil_vz             = rbuff(88)
#endif /* VZ_LIMITS */
#ifdef COLLISIONS
      cfl_colls	          = rbuff(89)
#endif /* COLLISIONS */


      flux_limiter        = trim(cbuff(80))
      freezing_speed      = trim(cbuff(81))
      dimensions          = cbuff(82)(1:16)
      magnetic            = cbuff(83)(1:16)

      integration_order   = ibuff(80)

#ifdef GRAV
!  namelist /GRAVITY/ gpt_hdf,  g_z, g_y, dg_dz, r_gc,
!                     ptmass,ptm_x,ptm_y,ptm_z,r_smooth
!                     tune_zeq, tune_zeq_bnd,
!                     h_gravity_profile, r_gravity_profile, n_gravity_profile
!                     h_grav, r_grav, n_gravr, n_gravh

      nsub                = ibuff(90)
      n_gravr             = ibuff(91)
      n_gravr2            = ibuff(92)
      n_gravh             = ibuff(93)

      gpt_hdf             = trim(cbuff(90))

      g_z                 = rbuff(90)
      g_y                 = rbuff(185)
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
#endif /* GRAV */

#ifdef RESISTIVE
!   namelist /RESISTIVITY/ cfl_resist, eta_0, eta_1, j_crit

       eta_scale          = ibuff(120)

       cfl_resist         = rbuff(120)
       eta_0              = rbuff(121)
       eta_1              = rbuff(122)
       j_crit             = rbuff(123)
       deint_max          = rbuff(124)
#endif /* RESISTIVE */

!#ifdef COSM_RAYS
!  namelist /COSMIC_RAYS/ cr_active, gamma_cr, cr_eff, beta_cr, K_cr_paral, K_cr_perp,&
!                         amp_cr, cfl_cr
!
!       cr_active          = rbuff(130)
!       gamma_cr           = rbuff(131)
!       cr_eff             = rbuff(132)
!       beta_cr            = rbuff(133)
!       K_cr_paral         = rbuff(134)
!       K_cr_perp          = rbuff(135)
!       amp_cr             = rbuff(136)
!       cfl_cr             = rbuff(137)
!       smallecr           = rbuff(138)
!#endif /* COSM_RAYS */


#ifdef GALAXY
!  namelist /GALACTIC_PARAMS/ h_sn, r_sn, f_sn_kpc2, t_dw, t_arm

       h_sn 		  = rbuff(140)
       r_sn 		  = rbuff(141)
       f_sn_kpc2 	  = rbuff(142)
       t_dw 		  = rbuff(143)
       t_arm 		  = rbuff(144)
       col_dens           = rbuff(145)
#endif /* GALAXY */

#ifdef SHEAR
!  namelist /SHEARING/ omega, qshear
       omega              = rbuff(160)
       qshear             = rbuff(161)
#endif /* SHEAR */
    endif  ! (proc .eq. 0)

#ifdef SN_SRC
!  namelist /SN_PARAMS/ h_sn, r_sn, f_sn_kpc2, amp_dip_sn, howmulti
       h_sn               = rbuff(170)
       r_sn               = rbuff(171)
       f_sn_kpc2          = rbuff(172)
       amp_dip_sn         = rbuff(173)
       howmulti           = ibuff(170)
#endif /* SN_SRC */
#ifdef SNE_DISTR
!  namelist /SN_DISTR/ snenerg, snemass, sn1time, sn2time, r0sn, add_mass, add_ener, add_encr, add_magn
       snenerg            = rbuff(180)
       snemass            = rbuff(181)
       sn1time            = rbuff(182)
       sn2time            = rbuff(183)
       r0sn               = rbuff(184)
       add_mass           = cbuff(180)
       add_ener           = cbuff(181)
       add_encr           = cbuff(182)
       add_magn           = cbuff(183)
#endif /* SNE_DISTR */

! Secondary parameters

!   csi2  = c_si**2
!   csim2 = csi2*(1.+alpha)    ! z-equilibrium defined for fixed
                                ! ratio p_mag/p_gas = alpha = 1/beta
#ifdef COSM_RAYS
!   csim2 = csim2 +csi2*beta_cr
#endif /* COSM_RAYS */

   ethu = 7.0**2/(5.0/3.0-1.0) * 1.0    ! thermal energy unit=0.76eV/cm**3
                                        ! for c_si= 7km/s, n=1/cm^3
                                        ! gamma=5/3
#ifdef GALAXY
#ifdef COSM_RAYS
   amp_ecr_sn = 4.96e6*cr_eff/r_sn**3   ! cosmic ray explosion amplitude
                                        ! in units:
        				! e_0 = 1/(5/3-1)*rho_0*c_s0**2
        				! rho_0=1.67e-24g/cm**3,
        				! c_s0 = 7km/s
#endif /* COSM_RAYS */
   f_sn = f_sn_kpc2 * (xmax-xmin)/1000.0 * (ymax-ymin)/1000.0 ! SN frequency per horizontal
                                                              ! surface area of the comp. box
#endif /* GALAXY */

  cn(1:3,1) = (/ 1. , 0.5 , 0.0 /)
  cn(1:3,2) = (/ 1. , 1.  , 0.0 /)
  if(integration_order .eq. 1) then
    cn(2,1) = 1.
  endif


!-------------------------
    if(integration_order .gt. 2) then
      stop 'For "ORIG" scheme integration_order must be 1 or 2'
    endif

    select case (dimensions)
      case('3d','2dxy')
        !do nothing
      case default
        stop '"dimensions" must be one of the following: "3d","2dxy"'
    end select


#ifdef GRAV
      gravaccel = .true.
#else /* GRAV */
      gravaccel = .false.
#endif /* GRAV */

    if(magnetic .eq. 'yes') then
      magfield = .true.
    else
      magfield = .false.
    endif

#ifdef RESISTIVE
      resist = .true.
      if(eta_scale .lt. 0) then
        write(*,*) 'eta_scale must be greater or equal 0'
        stop
      endif
#else /* RESISTIVE */
      resist = .false.
#endif /* RESISTIVE */

  end subroutine read_params

end module start








