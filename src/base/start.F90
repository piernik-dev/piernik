! $Id$
#include "piernik.def"

module start

! Written by: M. Hanasz, January 2006

  use mpi_setup


  implicit none


  real t,dt
  integer nstep
  integer nstep_start

  real tend
  integer nend, maxxyz

  real cfl, smalld, smallei

  character*16 dimensions
  integer integration_order, istep
  real rorder

  real, dimension(1) :: gamma   !!! do poprawy

  real    cfl_resist, eta_0, eta_1, j_crit, deint_max
  integer eta_scale

  real  dt_cr

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

    namelist /END_CONTROL/ nend, tend

    namelist /NUMERICAL_SETUP/  cfl, smalld, smallei, &
                              integration_order, &
                              dimensions

#ifdef RESISTIVE
  namelist /RESISTIVITY/ cfl_resist, eta_0, eta_1, eta_scale, j_crit, deint_max
#endif /* RESISTIVE */

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
    integration_order  = 2
    dimensions = '3d'

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
#ifdef RESISTIVE
        read(unit=1,nml=RESISTIVITY)
#endif /* RESISTIVE */
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
#ifdef RESISTIVE
        write(unit=3,nml=RESISTIVITY)
#endif /* RESISTIVE */
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
!                              dimensions, nu_bulk, cfl_visc
!                              floor_vz, ceil_vz, cfl_colls
      rbuff(80) = cfl
      rbuff(83) = smalld
      rbuff(84) = smallei

      cbuff(82) = dimensions

      ibuff(80) = integration_order

#ifdef RESISTIVE
!   namelist /RESISTIVITY/ cfl_resist, eta_0, eta_1, j_crit

       ibuff(120) = eta_scale

       rbuff(120) = cfl_resist
       rbuff(121) = eta_0
       rbuff(122) = eta_1
       rbuff(123) = j_crit
       rbuff(124) = deint_max
#endif /* RESISTIVE */

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


      dimensions          = cbuff(82)(1:16)

      integration_order   = ibuff(80)

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


#ifdef RESISTIVE
      if(eta_scale .lt. 0) then
        write(*,*) 'eta_scale must be greater or equal 0'
        stop
      endif
#endif /* RESISTIVE */

  end subroutine read_params

end module start








