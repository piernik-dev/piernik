!
! PIERNIK Code Copyright (C) 2006 Michal Hanasz
!
!    This file is part of PIERNIK code.
!
!    PIERNIK is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    PIERNIK is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with PIERNIK.  If not, see <http://www.gnu.org/licenses/>.
!
!    Initial implementation of PIERNIK code was based on TVD split MHD code by
!    Ue-Li Pen
!        see: Pen, Arras & Wong (2003) for algorithm and
!             http://www.cita.utoronto.ca/~pen/MHD
!             for original source code "mhd.f90"
!
!    For full list of developers see $PIERNIK_HOME/license/pdt.txt
!
#include "piernik.h"

!>
!! \brief Initialization of Cosmic Ray Energy SPectrum (CRESP) algorithm
!<

module initcrspectrum
! pulled by COSM_RAY_ELECTRONS
   use constants,       only: cbuff_len

   implicit none

   public ! QA_WARN no secrets are kept here

! contains routines reading namelist in problem.par file dedicated to cosmic ray electron spectrum and initializes types used.
! available via namelist COSMIC_RAY_SPECTRUM
   logical            :: use_cresp                   !< determines whether CRESP update is called by fluidupdate
   integer(kind=4)    :: ncre                        !< number of bins
   real(kind=8)       :: p_min_fix                   !< fixed momentum grid lower cutoff
   real(kind=8)       :: p_max_fix                   !< fixed momentum grid upper cutoff
   real(kind=8)       :: p_lo_init                   !< initial lower cutoff momentum
   real(kind=8)       :: p_up_init                   !< initial upper cutoff momentum
   character(len=cbuff_len) :: initial_condition     !< available types: bump, powl, brpl, symf, syme. Description below.
   real(kind=8)       :: p_br_init                   !< initial low energy break
   real(kind=8)       :: f_init                      !< initial value of distr. func. for isolated case
   real(kind=8)       :: q_init                      !< initial value of power law-like spectrum exponent
   real(kind=8)       :: q_br_init                   !< initial q for low energy break
   real(kind=8)       :: q_big                       !< maximal amplitude of q
   real(kind=8)       :: cfl_cre                     !< CFL parameter  for cr electrons
   real(kind=8)       :: cre_eff                     !< fraction of energy passed to cr-electrons by nucleons (mainly protons)
   real(kind=8)       :: K_cre_paral_1               !< maximal parallell diffusion coefficient value
   real(kind=8)       :: K_cre_perp_1                !< maximal perpendicular diffusion coefficient value
   real(kind=8)       :: K_cre_pow                   !< exponent for power law-like diffusion-energy dependance
   integer(kind=4)    :: expan_order                 !< 1,2,3 order of Taylor expansion for p_update (cresp_crspectrum)
   real(kind=8)       :: e_small                     !< lower energy cutoff for energy-approximated cutoff momenta
   logical            :: approx_cutoffs              !< T,F - turns off/on all approximating terms
   integer(kind=1)    :: e_small_approx_p_lo         !< 0,1 - turns off/on energy (e_small) approximated lower cutoff momentum in isolated case
   integer(kind=1)    :: e_small_approx_p_up         !< 0,1 - turns off/on energy (e_small) approximated upper cutoff momentum in isolated case
   integer(kind=1)    :: e_small_approx_init_cond    !< 0,1 - turns off/on energy (e_small) approximated momenta at initialization
   real(kind=8)       :: smallecrn                   !< floor value for CRESP number density
   real(kind=8)       :: smallecre                   !< floor value for CRESP energy density
   real(kind=8)       :: Gamma_min_fix               ! < min of Lorentzs' Gamma factor, lower range of CRESP fixed grid
   real(kind=8)       :: Gamma_max_fix               ! < max of Lorentzs' Gamma factor, upper range of CRESP fixed grid
   real(kind=8)       :: Gamma_lo_init               ! < min of Lorentzs' Gamma factor, lower range of initial spectrum
   real(kind=8)       :: Gamma_up_init               ! < max of Lorentzs' Gamma factor, upper range of initial spectrum
   real(kind=8)       :: max_p_ratio                 !< maximal ratio of momenta for solution grids resolved at initialization via cresp_NR_method
   integer(kind=2)    :: NR_iter_limit               !< maximal number of iterations for NR algorithm
   logical            :: force_init_NR               !< forces resolving new ratio solution grids at initialization
   logical            :: NR_run_refine_pf            !< enables "refine_grids" subroutines that fill empty spaces on the solution grid
   logical            :: NR_refine_solution_q        !< enables NR_1D refinement for value of interpolated "q" value
   logical            :: NR_refine_pf_lo             !< enables NR_2D refinement for interpolated values of "p" and "f" for lower cutoff. Note - algorithm tries to refine values if interpolation was unsuccessful.
   logical            :: NR_refine_pf_up             !< enables NR_2D refinement for interpolated values of "p" and "f" for upper cutoff. Note - algorithm tries to refine values if interpolation was unsuccessful.

   logical            :: nullify_empty_bins          !< nullifies empty bins when entering CRESP module / exiting empty cell.
   logical            :: allow_source_spectrum_break !< allow extension of spectrum to adjacent bins if momenta found exceed set p_fix
   logical            :: synch_active                !< TEST feature - turns on / off synchrotron cooling @ CRESP
   logical            :: adiab_active                !< TEST feature - turns on / off adiabatic   cooling @ CRESP
   logical            :: cre_gpcr_ess                !< electron essentiality for gpcr computation
   real(kind=8)       :: cre_active                  !< electron contribution to Pcr

! NR parameters
   real(kind=8)       :: tol_f                       !< tolerance for f abs. error in NR algorithm
   real(kind=8)       :: tol_x                       !< tolerance for x abs. error in NR algorithm
   real(kind=8)       :: tol_f_1D                    !< tolerance for f abs. error in NR algorithm (1D)
   real(kind=8)       :: tol_x_1D                    !< tolerance for x abs. error in NR algorithm (1D)
   integer(kind=4)    :: arr_dim, arr_dim_q

   real(kind=8), parameter  :: eps = 1.0e-15          !< epsilon parameter for real number comparisons
!----------------------------------
   real(kind=8), allocatable, dimension(:) :: p_fix, p_mid_fix, n_small_bin
   real(kind=8)                            :: w

   real(kind=8), allocatable, dimension(:) :: mom_cre_fix, mom_mid_cre_fix, Gamma_fix, Gamma_mid_fix, gamma_beta_c_fix
   real(kind=8)                            :: Gamma_fix_ratio
   real(kind=8)                            :: G_w

! Types used in module:
   type bin_old
      integer                                :: i_lo
      integer                                :: i_up
      real(kind=8), allocatable,dimension(:) :: p
      real(kind=8), allocatable,dimension(:) :: f
      real(kind=8), allocatable,dimension(:) :: q
      real(kind=8), allocatable,dimension(:) :: e
      real(kind=8), allocatable,dimension(:) :: n
      real(kind=8)                           :: dt
   end type bin_old

   type cr_spectrum
      real(kind=8), allocatable,dimension(:) :: e
      real(kind=8), allocatable,dimension(:) :: n
   end type cr_spectrum

   type(cr_spectrum) cresp
   type(cr_spectrum) norm_init_spectrum
   type(bin_old) crel
! For passing terms to compute energy sources / sinks

   type spec_mod_trms
      real(kind=8) :: ub
      real(kind=8) :: ud
      real(kind=8) :: ucmb
   end type spec_mod_trms

   real(kind=8)     :: total_init_cree
   real(kind=8)     :: p_fix_ratio
   integer, allocatable, dimension(:) :: cresp_all_edges, cresp_all_bins

! CRESP names
   character(len=*), parameter :: nam_cresp_f = "cref" !< helping array for CRESP number density
   character(len=*), parameter :: nam_cresp_p = "crep" !< helping array for CRESP energy density
   character(len=*), parameter :: nam_cresp_q = "creq" !< helping array for CRESP energy density

   logical          :: hdf_save_fpq                    ! diagnostic, if true - adding 'cref', 'crep', 'creq' to hdf_vars must follow

!====================================================================================================
!
 contains
!
!====================================================================================================
   subroutine init_cresp

      use constants,       only: cbuff_len, I_ZERO, zero, one, ten
      use cresp_variables, only: clight ! use units,   only: clight
      use dataio_pub,      only: printinfo, warn, msg, die, nh
      use diagnostics,     only: my_allocate_with_index
      use mpisetup,        only: rbuff, ibuff, lbuff, cbuff, master, slave, piernik_MPI_Bcast
      use units,           only: me

      implicit none

      logical, save            :: first_run = .true.
      integer                  :: i       ! enumerator
      real(kind=8)             :: p_br_def, q_br_def

      namelist /COSMIC_RAY_SPECTRUM/ cfl_cre, p_lo_init, p_up_init, f_init, q_init, q_big, ncre, initial_condition, &
      &                         p_min_fix, p_max_fix, cre_eff, K_cre_paral_1, K_cre_perp_1, cre_active, p_br_init,  &
      &                         K_cre_pow, expan_order, e_small, cre_gpcr_ess, use_cresp, e_small_approx_init_cond, &
      &                         e_small_approx_p_lo, e_small_approx_p_up, force_init_NR, NR_iter_limit, max_p_ratio,&
      &                         synch_active, adiab_active, arr_dim, arr_dim_q, q_br_init, Gamma_min_fix,           &
      &                         Gamma_max_fix, Gamma_lo_init, Gamma_up_init, nullify_empty_bins, approx_cutoffs,    &
      &                         NR_run_refine_pf, NR_refine_solution_q, NR_refine_pf_lo, NR_refine_pf_up, hdf_save_fpq

! Default values
      use_cresp         = .true.
      ncre              = 0
      p_min_fix         = 1.5e1
      p_max_fix         = 1.65e4
      p_lo_init         = 1.5e1
      p_up_init         = 7.5e2
      p_br_def          = p_lo_init
      initial_condition = "powl"
      f_init            = 1.0
      q_init            = 4.1
      q_br_def          = q_init
      q_big             = 30.0d0
      p_br_init         = p_br_def ! < in case it was not provided "powl" is assumed
      q_br_init         = q_br_def ! < in case it was not provided "powl" is assumed
      cfl_cre           = 0.1
      cre_eff           = 0.01
      K_cre_paral_1     = 0.
      K_cre_perp_1      = 0.
      K_cre_pow         = 0.
      expan_order       = 1
      Gamma_min_fix     = 2.5
      Gamma_max_fix     = 1000.0
      Gamma_lo_init     = 10.0
      Gamma_up_init     = 200.0

      approx_cutoffs    = .true.
      e_small           = 1.0e-5
      e_small_approx_p_lo = 1
      e_small_approx_p_up = 1
      e_small_approx_init_cond = 1
      max_p_ratio       = 2.5
      NR_iter_limit     = 100
      force_init_NR     = .false.
      NR_run_refine_pf  = .false.
      NR_refine_solution_q   = .false.
      NR_refine_pf_lo   = .false.
      NR_refine_pf_up   = .false.
      nullify_empty_bins     = .false.
      smallecrn              = 0.0
      smallecre              = 0.0
      allow_source_spectrum_break  = .false.
      synch_active = .true.
      adiab_active = .true.
      cre_gpcr_ess = .false.
      cre_active   = 0.0

! NR parameters
      tol_f    = 1.0e-11
      tol_x    = 1.0e-11
      tol_f_1D = 1.0e-14
      tol_x_1D = 1.0e-14
      arr_dim  = 200
      arr_dim_q = 500

      hdf_save_fpq = .false.

      if (master) then
         if (.not.nh%initialized) call nh%init()
         open(newunit=nh%lun, file=nh%tmp1, status="unknown")
         write(nh%lun,nml=COSMIC_RAY_SPECTRUM)
         close(nh%lun)
         open(newunit=nh%lun, file=nh%par_file)
         nh%errstr=""
         read(unit=nh%lun, nml=COSMIC_RAY_SPECTRUM, iostat=nh%ierrh, iomsg=nh%errstr)
         close(nh%lun)
         call nh%namelist_errh(nh%ierrh, "COSMIC_RAY_SPECTRUM")
         read(nh%cmdl_nml,nml=COSMIC_RAY_SPECTRUM, iostat=nh%ierrh)
         call nh%namelist_errh(nh%ierrh, "COSMIC_RAY_SPECTRUM", .true.)
         open(newunit=nh%lun, file=nh%tmp2, status="unknown")
         write(nh%lun,nml=COSMIC_RAY_SPECTRUM)
         close(nh%lun)
         call nh%compare_namelist()
      endif

      rbuff(:)   = huge(1.)                         ! mark unused entries to allow automatic determination of nn

      if (master) then
         ibuff(1)  = ncre
         ibuff(2)  = expan_order

         ibuff(3)  = e_small_approx_p_lo
         ibuff(4)  = e_small_approx_p_up
         ibuff(5)  = e_small_approx_init_cond

         ibuff(6)  =  NR_iter_limit
         ibuff(7)  =  arr_dim
         ibuff(8)  =  arr_dim_q

         lbuff(1)  =  use_cresp
         lbuff(2)  =  cre_gpcr_ess
         lbuff(3)  =  allow_source_spectrum_break
         lbuff(4)  =  synch_active
         lbuff(5)  =  adiab_active

         lbuff(6)  =  force_init_NR
         lbuff(7)  =  NR_run_refine_pf
         lbuff(8)  =  NR_refine_solution_q
         lbuff(9)  =  NR_refine_pf_lo
         lbuff(10) =  NR_refine_pf_up
         lbuff(11) =  nullify_empty_bins
         lbuff(12) =  approx_cutoffs

         lbuff(13) = hdf_save_fpq

         rbuff(1)  = cfl_cre
         rbuff(2)  = cre_eff
         rbuff(3)  = smallecrn
         rbuff(4)  = smallecre
         rbuff(5)  = cre_active
         rbuff(6)  = p_lo_init
         rbuff(7)  = p_up_init
         rbuff(8)  = f_init
         rbuff(9)  = q_init
         rbuff(10) = q_big
         rbuff(11) = p_min_fix
         rbuff(12) = p_max_fix
         rbuff(13) = K_cre_paral_1
         rbuff(14) = K_cre_perp_1
         rbuff(15) = K_cre_pow

         rbuff(16) = e_small
         rbuff(17) = max_p_ratio

         rbuff(18) = tol_f
         rbuff(19) = tol_x
         rbuff(20) = tol_f_1D
         rbuff(21) = tol_x_1D

         rbuff(22) = Gamma_min_fix
         rbuff(23) = Gamma_max_fix
         rbuff(24) = Gamma_lo_init
         rbuff(25) = Gamma_up_init

         rbuff(26) = p_br_init
         rbuff(27) = q_br_init

         cbuff(1)  = initial_condition
      endif

      call piernik_MPI_Bcast(ibuff)
      call piernik_MPI_Bcast(rbuff)
      call piernik_MPI_Bcast(lbuff)
      call piernik_MPI_Bcast(cbuff, cbuff_len)

!!\deprecated
      open(10, file='crs.dat',status='replace',position='rewind')     ! diagnostic files
      open(11, file='crs_ne.dat',status='replace',position='rewind')  ! diagnostic files

      if (slave) then
         ncre                        = int(ibuff(1),kind=4)
         expan_order                 = int(ibuff(2),kind=4)

         e_small_approx_p_lo         = int(ibuff(3),kind=1)
         e_small_approx_p_up         = int(ibuff(4),kind=1)
         e_small_approx_init_cond    = int(ibuff(5),kind=1)

         NR_iter_limit               = int(ibuff(6),kind=2)
         arr_dim                     = int(ibuff(7),kind=4)
         arr_dim_q                   = int(ibuff(8),kind=4)

         use_cresp                   = lbuff(1)
         cre_gpcr_ess                = lbuff(2)
         allow_source_spectrum_break = lbuff(3)
         synch_active                = lbuff(4)
         adiab_active                = lbuff(5)

         force_init_NR               = lbuff(6)
         NR_run_refine_pf            = lbuff(7)
         NR_refine_solution_q        = lbuff(8)
         NR_refine_pf_lo             = lbuff(9)
         NR_refine_pf_up             = lbuff(10)
         nullify_empty_bins          = lbuff(11)
         approx_cutoffs              = lbuff(12)

         hdf_save_fpq                = lbuff(13)

         cfl_cre                     = rbuff(1)
         cre_eff                     = rbuff(2)
         smallecrn                   = rbuff(3)
         smallecre                   = rbuff(4)
         cre_active                  = rbuff(5)
         p_lo_init                   = rbuff(6)
         p_up_init                   = rbuff(7)
         f_init                      = rbuff(8)
         q_init                      = rbuff(9)
         q_big                       = rbuff(10)
         p_min_fix                   = rbuff(11)
         p_max_fix                   = rbuff(12)
         K_cre_paral_1               = rbuff(13)
         K_cre_perp_1                = rbuff(14)
         K_cre_pow                   = rbuff(15)

         e_small                     = rbuff(16)
         max_p_ratio                 = rbuff(17)

         tol_f                       = rbuff(18)
         tol_x                       = rbuff(19)
         tol_f_1D                    = rbuff(20)
         tol_x_1D                    = rbuff(21)

         Gamma_min_fix               = rbuff(22)
         Gamma_max_fix               = rbuff(23)
         Gamma_lo_init               = rbuff(24)
         Gamma_up_init               = rbuff(25)

         p_br_init                   = rbuff(26)
         q_br_init                   = rbuff(27)

         initial_condition           = trim(cbuff(1))

      endif
      if (first_run .eqv. .true.) then
         if (ncre .ne. I_ZERO)  then
#ifdef VERBOSE
               write (msg,'(A)')            '[initcrspectrum:init_cresp] Initial CRESP parameters read:'
               call printinfo(msg)
               write (msg, '(A, L1)')       '[initcrspectrum:init_cresp] use_cresp   = ', use_cresp
               call printinfo(msg)
               write (msg, '(A, 1I3)')      '[initcrspectrum:init_cresp] ncre        = ', ncre
               call printinfo(msg)
               write (msg, '(A, 1E15.7)')   '[initcrspectrum:init_cresp] p_min_fix   = ', p_min_fix
               call printinfo(msg)
               write (msg, '(A, 1E15.7)')   '[initcrspectrum:init_cresp] p_max_fix   = ', p_max_fix
               call printinfo(msg)
               write (msg, '(A, 1E15.7)')   '[initcrspectrum:init_cresp] p_lo_init   = ', p_lo_init
               call printinfo(msg)
               write (msg, '(A, 1E15.7)')   '[initcrspectrum:init_cresp] p_up_init   = ', p_up_init
               call printinfo(msg)
               write (msg, '(A, 1F15.7)')   '[initcrspectrum:init_cresp] q_init      = ', q_init
               call printinfo(msg)
               write (msg, '(A, 1F15.7)')   '[initcrspectrum:init_cresp] q_big       = ', q_big
               call printinfo(msg)
               write (msg, '(A, 1F15.7)')   '[initcrspectrum:init_cresp] cfl_cre     = ', cfl_cre
               call printinfo(msg)
               write (msg, '(A, 10E15.7)')  '[initcrspectrum:init_cresp] K_cre_paral1 = ', K_cre_paral_1
               call printinfo(msg)
               write (msg, '(A, 10E15.7)')  '[initcrspectrum:init_cresp] K_cre_perp_1 =', K_cre_perp_1
               call printinfo(msg)
               write (msg, '(A, 10E15.7)')  '[initcrspectrum:init_cresp] K_cre_pow    =', K_cre_pow
               call printinfo(msg)
               write (msg, '(A, 1E15.7)')   '[initcrspectrum:init_cresp] e_small      =', e_small
               call printinfo(msg)
               write (msg, '(A, A5)')       '[initcrspectrum:init_cresp] initial_condition =' , initial_condition
               call printinfo(msg)
               write (msg, '(A, L1)')       '[initcrspectrum:init_cresp] cre_gpcr_ess = ', cre_gpcr_ess
               call printinfo(msg)
               write (msg, '(A, L1)')       '[initcrspectrum:init_cresp] cre_active   = ', cre_active
               call printinfo(msg)
               write (msg, '(A, I1)')       '[initcrspectrum:init_cresp] Approximate cutoff momenta at initialization: e_small_approx_init_cond =', e_small_approx_init_cond
               call printinfo(msg)
               write (msg, '(A, I1)')       '[initcrspectrum:init_cresp] Approximate lower momentum cutoff: e_small_approx_p_lo =', e_small_approx_p_lo
               call printinfo(msg)
               write (msg, '(A, I1)')       '[initcrspectrum:init_cresp] Approximate upper momentum cutoff: e_small_approx_p_up =', e_small_approx_p_up
               call printinfo(msg)
               write (msg, '(A, 1F10.5 )')  '[initcrspectrum:init_cresp] max_p_ratio      =', max_p_ratio
               call printinfo(msg)
               write (msg, '(A, L2 )' )     '[initcrspectrum:init_cresp] force_init_NR    = ', force_init_NR
               call printinfo(msg)
               write (msg, '(A, I4)')       '[initcrspectrum:init_cresp] NR_iter_limit    = ', NR_iter_limit
               call printinfo(msg)
               write (msg, '(A, 1E15.7)')   '[initcrspectrum:init_cresp] epsilon(eps)     = ', eps
               call printinfo(msg)
               write (msg, '(A, 1E15.7)')   '[initcrspectrum:init_cresp] Gamma_min_fix    =', Gamma_min_fix
               call printinfo(msg)
               write (msg, '(A, 1E15.7)')   '[initcrspectrum:init_cresp] Gamma_max_fix    =', Gamma_max_fix
               call printinfo(msg)
               write (msg, '(A, 1E15.7)')   '[initcrspectrum:init_cresp] Gamma_lo_init    =', Gamma_lo_init
               call printinfo(msg)
               write (msg, '(A, 1E15.7)')   '[initcrspectrum:init_cresp] Gamma_up_init    =', Gamma_up_init
               call printinfo(msg)
               write (msg,'(A, L1)')        '[initcrspectrum:init_cresp] nullify_empty_bins =', nullify_empty_bins
               call printinfo(msg)
               write (msg, '(A, 1E15.7)')   '[initcrspectrum:init_cresp] smallecrn        = ', smallecrn
               call printinfo(msg)
               write (msg, '(A, 1E15.7)')   '[initcrspectrum:init_cresp] smallecre        = ', smallecre
               call printinfo(msg)
               write (msg, '(A, L1)')       '[initcrspectrum:init_cresp] allow_source_spectrum_break =', allow_source_spectrum_break
               call printinfo(msg)
               write (msg, '(A, L1)')       '[initcrspectrum:init_cresp] synch_active = ', synch_active
               call printinfo(msg)
               write (msg, '(A, L1)')       '[initcrspectrum:init_cresp] adiab_active = ', adiab_active
               call printinfo(msg)

#endif /* VERBOSE */
               if (ncre .lt. 3) then
                  write (msg,'(A)') "[initcrspectrum:init_cresp] CRESP algorithm currently requires at least 3 bins (ncre) in order to work properly, check your parameters."
                  call die(msg)
               endif

               if (approx_cutoffs) then
                  e_small_approx_p_lo = 1; e_small_approx_p_up = 1
                  write (msg,'(A)') "[initcrspectrum:init_cresp] approx_cutoffs = .true. -- will use e_small to approximate spectrum cutoffs and initial state spectrum."
               else
                  e_small_approx_p_lo      = 0 ; e_small_approx_p_up = 0 ! e_small_approx_init_cond stays default, unless user changes.
                  write (msg,'(A)') "[initcrspectrum:init_cresp] approx_cutoffs = .false. -- will not use e_small approximated cutoffs, but still approximate initial state. To turn it off use e_small_approx_init_cond = 0."
               endif
               call printinfo(msg)

               if ( (e_small_approx_p_lo+e_small_approx_p_up) .gt. 0 .and. e_small_approx_init_cond .lt. 1) then
                  e_small_approx_init_cond = 1  !
                  write (msg,'(A)') "[initcrspectrum:init_cresp] Approximation of boundary momenta is active -> modifying e_small_approx_init_cond to 1."
                  if (master) call warn(msg)
                  call sleep(1)
               endif

               if (hdf_save_fpq) then
                  write(msg, '(A)') "[initcrspectrum:init_cresp] hdf_save_fpq is set. Adding 'cref', 'crep', 'creq' to hdf_vars must follow."
                  call warn(msg)
               endif
! countermeasure - in case unrecognized or invalid parameters are provided

               if ( e_small_approx_p_lo .gt. 0 ) then ; e_small_approx_p_lo = 1 ; else ; e_small_approx_p_lo = 0 ; endif
               if ( e_small_approx_p_up .gt. 0 ) then ; e_small_approx_p_up = 1 ; else ; e_small_approx_p_up = 0 ; endif
               if ( e_small_approx_init_cond .gt. 0 ) then ; e_small_approx_init_cond = 1 ; else ; e_small_approx_init_cond = 0 ; endif

               if (e_small_approx_p_lo+e_small_approx_p_up .eq. 0) then
                  NR_refine_solution_q = .true. !< for testing we leave precise solutions of q (especially for outer momenta)
               endif

               if (e_small_approx_init_cond + e_small_approx_p_lo + e_small_approx_p_up .eq. 0) then
                  e_small = zero                !< no threshold energy for bin activation necessary
               endif
! arrays initialization
               call my_allocate_with_index(p_fix,ncre,0)
               call my_allocate_with_index(p_mid_fix,ncre,1)
               call my_allocate_with_index(cresp_all_edges,ncre,0)
               call my_allocate_with_index(cresp_all_bins, ncre,1)
               call my_allocate_with_index(n_small_bin,ncre,1)

               call my_allocate_with_index(Gamma_fix,ncre,0)
               call my_allocate_with_index(Gamma_mid_fix,ncre,1)
               call my_allocate_with_index(mom_cre_fix,ncre,0)
               call my_allocate_with_index(mom_mid_cre_fix,ncre,1)
               call my_allocate_with_index(gamma_beta_c_fix,ncre,0)

               cresp_all_edges = (/ (i,i=0,ncre) /)
               cresp_all_bins  = (/ (i,i=1,ncre) /)

!!\brief for now algorithm requires at least 3 bins
               p_fix = zero
               w  = (log10(p_max_fix/p_min_fix))/real(ncre-2,kind=8)
               p_fix(1:ncre-1)  =  p_min_fix*ten**(w* real((cresp_all_edges(1:ncre-1)-1),kind=8) )
               p_fix(0)    = zero
               p_fix(ncre) = zero
               p_fix_ratio = ten**w

               p_mid_fix = 0.0
               p_mid_fix(2:ncre-1) = sqrt(p_fix(1:ncre-2)*p_fix(2:ncre-1))
               p_mid_fix(1)    = p_mid_fix(2) / p_fix_ratio
               p_mid_fix(ncre) = p_mid_fix(ncre-1) * p_fix_ratio

!> set Gamma arrays, analogically to p_fix arrays, that will be constructed using Gamma arrays
               Gamma_fix            = one             !< Gamma factor obviously cannot be lower than 1
               G_w                  = (log10(Gamma_max_fix/Gamma_min_fix))/real(ncre-2,kind=8)
               Gamma_fix(1:ncre-1)  = Gamma_min_fix * ten**(G_w * real((cresp_all_edges(1:ncre-1)-1),kind=8))
               Gamma_fix_ratio      = ten**w

               Gamma_mid_fix = one
               Gamma_mid_fix(2:ncre-1) = sqrt( Gamma_fix(1:ncre-2)   * Gamma_fix(2:ncre-1) )
               Gamma_mid_fix(1)        = sqrt( Gamma_mid_fix(1)      * Gamma_mid_fix(2))
               Gamma_mid_fix(ncre)     = sqrt( Gamma_mid_fix(ncre-1) * Gamma_mid_fix(ncre-1) * Gamma_fix_ratio )
! compute physical momenta of particles in given unit set
               mom_cre_fix      = (/ (cresp_get_mom(Gamma_fix(i),me),     i=0,ncre ) /)
               mom_mid_cre_fix  = (/ (cresp_get_mom(Gamma_mid_fix(i),me), i=1,ncre ) /)

               gamma_beta_c_fix = mom_cre_fix / me

               n_small_bin(:) = e_small / (p_mid_fix(:) * clight)

#ifdef VERBOSE
               write (msg,'(A, 50I3)')    '[initcrspectrum:init_cresp] fixed all edges: ', cresp_all_edges
               call printinfo(msg)
               write (msg,'(A, 50E15.7)') '[initcrspectrum:init_cresp] Fixed momentum grid: ', p_fix
               call printinfo(msg)
               write (msg,'(A, 50E15.7)') '[initcrspectrum:init_cresp] Bin p-width (log10): ', w
               call printinfo(msg)
               write (msg,'(A, 50E15.7)') '[initcrspectrum:init_cresp] Fixed momentum grid (bin middle):   ',p_mid_fix(1:ncre)
               call printinfo(msg)
               write (msg,'(A, 50F13.2)') '[initcrspectrum:init_cresp] Fixed Gamma      grid: ', Gamma_fix
               call printinfo(msg)
               write (msg,'(A, 50F10.5)') '[initcrspectrum:init_cresp] Gamma bin width(log10): ', G_w
               call printinfo(msg)
               write (msg,'(A, 50F10.5)') '[initcrspectrum:init_cresp] Fixed mid-Gamma     : ', Gamma_mid_fix(1:ncre)
               call printinfo(msg)
               write (msg,'(A, 50E15.7)') '[initcrspectrum:init_cresp] Fixed phys momentum : ', mom_cre_fix
               call printinfo(msg)
               write (msg,'(A, 50E15.7)') '[initcrspectrum:init_cresp] Fixed phys mid mom  : ', mom_mid_cre_fix
               call printinfo(msg)
#endif /* VERBOSE */

! Input parameters check
               else
                  write (msg,"(A,I4,A)") '[initcrspectrum:init_cresp] ncre   = ', ncre, &
                        '; cr-electrons NOT initnialized. If COSM_RAY_ELECTRONS flag is on, please check your parameters.'
                  call die(msg)
               endif
!>
!!\brief Correctness of "initial_condition" is checked here
!!
!! Description of initial_condition keywords: powl - pure power-law like distribution function,
!! brpl - broken power-law like (with break in the first bin, making it easier for NR algorithm
!! to find solution of lower cutoff momentum), bump - gaussian-like spectrum,
!! \deprecated syme - symmetric energy distribution relative to the middle of the initial spectrum,
!! \deprecated symf - similar, but symmetric in distribution function.
!<
               if (initial_condition .ne. 'powl' .and. initial_condition .ne. 'bump' .and. initial_condition .ne. 'brpl' &
                           .and. initial_condition .ne. 'symf' .and. initial_condition .ne. 'syme'  .and. initial_condition .ne. 'brpg' ) then
                  write(msg,"(A,A,A)") "[initcrspectrum:init_cresp] Provided unrecognized initial_condition (",initial_condition,&
                                                   "). Make sure that value is correctly provided."
                  call die(msg)
               endif
               if ( f_init .lt. eps) then
                  if (initial_condition == 'powl' .or. initial_condition == 'brpl') then
                  write (msg,"(A,A,A)") "[initcrspectrum:init_cresp] Provided power law type spectrum (",initial_condition &
                     ,") with initial amplitude f_init ~ zero. Check your parameters."
                  call die(msg)
               endif
         endif

         if (initial_condition .eq. "brpl" ) then ! FIXME TODO
            if (abs(p_br_init - p_br_def) .le. eps) then
               write (msg,"(A)") "[initcrspectrum:init_cresp] Parameter for 'brpl' spectrum: p_br_init has default value (probably unitialized). Assuming p_lo_init value ('powl' spectrum)."
               if (master) call warn(msg)
            else
!>
!! \brief p_br_init should be equal to one of p_fix values
!<
               i = minloc(abs(p_fix - p_br_init),dim=1)-1
               write (msg,"(A,E14.7,1A)") "[initcrspectrum:init_cresp] p_br_init was set, but should be equal to one of p_fix. Assuming p_br_init =", p_fix(i),"."
               p_br_init = p_fix(i)
               if (master) call warn(msg)
            endif
            if (abs(q_br_init - q_br_def) .le. eps) then
               write (msg,"(A)") "[initcrspectrum:init_cresp] Parameter for 'brpl' spectrum: q_br_init has default value (probably unitialized). Assuming q_init value    ('powl' spectrum)."
               if (master) call warn(msg)
            endif
         endif

         call init_cresp_types
      endif
   end subroutine init_cresp

!----------------------------------------------------------------------------------------------------

   subroutine init_cresp_types

      use constants,   only: zero, I_ZERO
      use diagnostics, only: my_allocate_with_index

      implicit none

      if (.not. allocated(crel%p)) call my_allocate_with_index(crel%p,ncre,0)
      if (.not. allocated(crel%f)) call my_allocate_with_index(crel%f,ncre,0)
      if (.not. allocated(crel%q)) call my_allocate_with_index(crel%q,ncre,1)
      if (.not. allocated(crel%n)) call my_allocate_with_index(crel%n,ncre,1)
      if (.not. allocated(crel%e)) call my_allocate_with_index(crel%e,ncre,1)

      if (.not. allocated(cresp%n)) call my_allocate_with_index(cresp%n,ncre,1)
      if (.not. allocated(cresp%e)) call my_allocate_with_index(cresp%e,ncre,1)
      if (.not. allocated(norm_init_spectrum%n)) call my_allocate_with_index(norm_init_spectrum%n,ncre,1)
      if (.not. allocated(norm_init_spectrum%e)) call my_allocate_with_index(norm_init_spectrum%e,ncre,1)

      crel%p = zero
      crel%q = zero
      crel%f = zero
      crel%e = zero
      crel%n = zero
      crel%i_lo = I_ZERO
      crel%i_up = I_ZERO

      cresp%e = zero
      cresp%n = zero

      norm_init_spectrum%n = zero
      norm_init_spectrum%e = zero

   end subroutine init_cresp_types

!----------------------------------------------------------------------------------------------------

   function cresp_get_mom(gamma, particle_mass)

      use constants, only: zero, one
      use units,     only: clight

      implicit none

      real(kind=8)           :: gamma
      real(kind=8), optional :: particle_mass
      real(kind=8)           :: cresp_get_mom

      cresp_get_mom = zero
      if ( (gamma - one) .gt. eps ) then
         cresp_get_mom = gamma * particle_mass * sqrt(one - one/(gamma**2)) * clight
      endif

   end function cresp_get_mom

!----------------------------------------------------------------------------------------------------

   subroutine cleanup_cresp_work_arrays

      use diagnostics, only: my_deallocate

      implicit none

      if (allocated(cresp%n))   call my_deallocate(cresp%n)
      if (allocated(cresp%e))   call my_deallocate(cresp%e)
      if (allocated(norm_init_spectrum%n))   call my_deallocate(norm_init_spectrum%n)
      if (allocated(norm_init_spectrum%e))   call my_deallocate(norm_init_spectrum%e)

      if (allocated(p_fix)) call my_deallocate(p_fix)
      if (allocated(p_mid_fix)) call my_deallocate(p_mid_fix)
      if (allocated(cresp_all_edges)) call my_deallocate(cresp_all_edges)
      if (allocated(cresp_all_bins )) call my_deallocate(cresp_all_bins)
      if (allocated(crel%p))  call my_deallocate(crel%p)
      if (allocated(crel%f)) call my_deallocate(crel%f)
      if (allocated(crel%q)) call my_deallocate(crel%q)
      if (allocated(crel%e)) call my_deallocate(crel%e)
      if (allocated(crel%n)) call my_deallocate(crel%n)

   end subroutine cleanup_cresp_work_arrays

!----------------------------------------------------------------------------------------------------

end module initcrspectrum
