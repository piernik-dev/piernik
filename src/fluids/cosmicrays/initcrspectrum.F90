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
! pulled by CRESP

   use constants, only: cbuff_len, fnamelen

   implicit none

   private
   public :: use_cresp, use_cresp_evol, p_init, initial_spectrum, p_br_init, f_init, q_init, q_br_init, q_big, cfl_cre, cre_eff, expan_order, e_small, e_small_approx_p, e_small_approx_init_cond,  &
           & smallcren, smallcree, max_p_ratio, NR_iter_limit, force_init_NR, NR_run_refine_pf, NR_refine_solution_q, NR_refine_pf, nullify_empty_bins, synch_active, adiab_active,                 &
           & icomp_active, allow_source_spectrum_break, cre_active, tol_f, tol_x, tol_f_1D, tol_x_1D, arr_dim_a, arr_dim_n, arr_dim_q, eps, eps_det, w, p_fix, p_mid_fix, total_init_cree, p_fix_ratio,           &
           & spec_mod_trms, cresp_all_edges, cresp_all_bins, norm_init_spectrum, cresp, crel, dfpq, f_synchIC, init_cresp, cleanup_cresp_sp, check_if_dump_fpq, cleanup_cresp_work_arrays, q_eps,     &
           & u_b_max, def_dtsynchIC, def_dtadiab, NR_smap_file, NR_allow_old_smaps, cresp_substep, n_substeps_max, allow_unnatural_transfer, K_cresp_paral, K_cresp_perp, p_min_fix, p_max_fix, redshift

! contains routines reading namelist in problem.par file dedicated to cosmic ray electron spectrum and initializes types used.
! available via namelist COSMIC_RAY_SPECTRUM
   logical         :: use_cresp                   !< determines whether CRESP routines are called anywhere
   logical         :: use_cresp_evol              !< determines whether CRESP update is called by fluidupdate
   real            :: p_min_fix                   !< fixed momentum grid lower cutoff
   real            :: p_max_fix                   !< fixed momentum grid upper cutoff
   real            :: p_lo_init                   !< initial lower cutoff momentum
   real            :: p_up_init                   !< initial upper cutoff momentum
   real, dimension(2) :: p_init                   !< vector to store p_lo_init and p_up_init
   character(len=cbuff_len) :: initial_spectrum   !< available types: bump, powl, brpl, symf, syme. Description below.
   real            :: p_br_init_lo, p_br_init_up  !< initial low energy break
   real, dimension(2) :: p_br_init                !< vector to store p_br_init_lo and p_br_init_up
   real            :: f_init                      !< initial value of distr. func. for isolated case
   real            :: q_init                      !< initial value of power law-like spectrum exponent
   real            :: q_br_init                   !< initial q for low energy break
   real            :: q_big                       !< maximal amplitude of q
   real            :: cfl_cre                     !< CFL parameter  for cr electrons
   real            :: cre_eff                     !< fraction of energy passed to cr-electrons by nucleons (mainly protons)
   real, dimension(:), allocatable :: K_cresp_paral !< array containing parallel diffusion coefficients of all CR CRESP components (number density and energy density)
   real, dimension(:), allocatable :: K_cresp_perp  !< array containing perpendicular diffusion coefficients of all CR CRESP components (number density and energy density)
   real            :: K_cre_pow                   !< exponent for power law-like diffusion-energy dependence
   real            :: p_diff                      !< momentum to which diffusion coefficients refer to
   integer(kind=4) :: expan_order                 !< 1,2,3 order of Taylor expansion for p_update (cresp_crspectrum)
   real            :: e_small                     !< lower energy cutoff for energy-approximated cutoff momenta
   logical         :: approx_cutoffs              !< T,F - turns off/on all approximating terms
   integer(kind=4), dimension(2) :: e_small_approx_p !< vector to store e_small_approx_p_lo and e_approx_p_up
   integer(kind=4) :: e_small_approx_p_lo         !< 0,1 - turns off/on energy (e_small) approximated lower cutoff momentum in isolated case
   integer(kind=4) :: e_small_approx_p_up         !< 0,1 - turns off/on energy (e_small) approximated upper cutoff momentum in isolated case
   integer(kind=1) :: e_small_approx_init_cond    !< 0,1 - turns off/on energy (e_small) approximated momenta at initialization
   real            :: smallcren                   !< floor value for CRESP number density
   real            :: smallcree                   !< floor value for CRESP energy density
   real            :: Gamma_min_fix               ! < min of Lorentzs' Gamma factor, lower range of CRESP fixed grid
   real            :: Gamma_max_fix               ! < max of Lorentzs' Gamma factor, upper range of CRESP fixed grid
   real            :: max_p_ratio                 !< maximal ratio of momenta for solution grids resolved at initialization via cresp_NR_method
   integer(kind=2) :: NR_iter_limit               !< maximal number of iterations for NR algorithm
   logical         :: force_init_NR               !< forces resolving new ratio solution grids at initialization
   logical         :: NR_run_refine_pf            !< enables "refine_grids" subroutines that fill empty spaces on the solution grid
   logical         :: NR_refine_solution_q        !< enables NR_1D refinement for value of interpolated "q" value
   logical         :: NR_refine_pf_lo             !< enables NR_2D refinement for interpolated values of "p" and "f" for lower cutoff. Note - algorithm tries to refine values if interpolation was unsuccessful.
   logical         :: NR_refine_pf_up             !< enables NR_2D refinement for interpolated values of "p" and "f" for upper cutoff. Note - algorithm tries to refine values if interpolation was unsuccessful.
   logical         :: NR_allow_old_smaps          !< allows to override h5 smap reading in favor of old ".dat" files ! WARNING : parameter not included in the namelist ! WARNING
   logical, dimension(2) :: NR_refine_pf          !< vector to store NR_refine_pf_lo and NR_refine_pf_up
   character(len=fnamelen):: NR_smap_file         !< provides name for NR solution maps to be read from / saved to

   logical         :: nullify_empty_bins          !< nullifies empty bins when entering CRESP module / exiting empty cell.
   logical         :: allow_source_spectrum_break !< allow extension of spectrum to adjacent bins if momenta found exceed set p_fix
   logical         :: allow_unnatural_transfer    !< allows unnatural transfer of n & e with 'manually_deactivate_bins_via_transfer'
   logical         :: synch_active                !< TEST feature - turns on / off synchrotron cooling @ CRESP
   logical         :: adiab_active                !< TEST feature - turns on / off adiabatic   cooling @ CRESP
   logical         :: icomp_active                !< TEST feature - turns on / off Inv-Compton cooling @ CRESP
   real            :: redshift                    !< redshift for chosen epoch WARNING this remains constant
   real            :: cre_active                  !< electron contribution to Pcr

! substepping parameters
   logical         :: cresp_substep               !< turns on / off usage of substepping for each cell independently
   integer(kind=4) :: n_substeps_max              !< maximal allowed number of substeps

! NR parameters
   real            :: tol_f                       !< tolerance for f abs. error in NR algorithm
   real            :: tol_x                       !< tolerance for x abs. error in NR algorithm
   real            :: tol_f_1D                    !< tolerance for f abs. error in NR algorithm (1D)
   real            :: tol_x_1D                    !< tolerance for x abs. error in NR algorithm (1D)
   integer(kind=4) :: arr_dim_a, arr_dim_n, arr_dim_q
   real            :: q_eps                       !< parameter for q tolerance (alpha_to_q)
   real            :: b_max_db, u_b_max           !< parameter limiting maximal value of B and implying maximal MF energy density u_b

   real, parameter :: eps = 1.0e-15          !< epsilon parameter for real number comparisons
   real, parameter :: eps_det = eps * 1.0e-15
!----------------------------------
   real, allocatable, dimension(:) :: p_fix, p_mid_fix, n_small_bin
   real                            :: w

   real, allocatable, dimension(:) :: mom_cre_fix, mom_mid_cre_fix, Gamma_fix, Gamma_mid_fix, gamma_beta_c_fix
   real                            :: Gamma_fix_ratio
   real                            :: G_w

! Types used in module:
   type bin_old
      integer,           dimension(2) :: i_cut
      real, allocatable, dimension(:) :: p
      real, allocatable, dimension(:) :: f
      real, allocatable, dimension(:) :: q
      real, allocatable, dimension(:) :: e
      real, allocatable, dimension(:) :: n
   end type bin_old

   type cr_spectrum
      real, allocatable, dimension(:) :: e
      real, allocatable, dimension(:) :: n
   end type cr_spectrum

   type(cr_spectrum) :: cresp
   type(cr_spectrum) :: norm_init_spectrum
   type(bin_old)     :: crel
! For passing terms to compute energy sources / sinks

   type spec_mod_trms
      real :: ub
      real :: ud
      real :: umag
      real :: ucmb
   end type spec_mod_trms

   real :: total_init_cree
   real :: p_fix_ratio
   integer(kind=4), allocatable, dimension(:) :: cresp_all_edges, cresp_all_bins

! CRESP names
   integer, parameter :: cnlen = 4
   type dump_fpq_type
      character(len=cnlen) :: f_nam = "cref" !< helping array for CRESP number density
      character(len=cnlen) :: p_nam = "crep" !< helping array for CRESP energy density
      character(len=cnlen) :: q_nam = "creq" !< helping array for CRESP energy density
      logical :: any_dump, dump_f, dump_p, dump_q  ! diagnostic, if true - adding 'cref', 'crep', 'creq' to hdf_vars must follow
   end type dump_fpq_type
   type(dump_fpq_type) :: dfpq

   real :: f_synchIC, def_dtadiab, def_dtsynchIC

!====================================================================================================
!
contains
!
!====================================================================================================
   subroutine init_cresp

      use bcast,           only: piernik_MPI_Bcast
      use constants,       only: cbuff_len, I_ZERO, I_ONE, zero, one, three, ten, half, logten, LO, HI
      use cr_data,         only: cr_table, icr_E
      use cresp_variables, only: clight_cresp
      use dataio_pub,      only: printinfo, warn, msg, die, nh
      use diagnostics,     only: my_allocate_with_index
      use global,          only: disallow_CRnegatives
      use func,            only: emag
      use initcosmicrays,  only: ncrb, ncr2b, ncrn, nspc, ncrtot, K_cr_paral, K_cr_perp, K_crs_paral, K_crs_perp, use_smallecr
      use mpisetup,        only: rbuff, ibuff, lbuff, cbuff, master, slave
      use units,           only: clight, me, sigma_T

      implicit none

      integer(kind=4) :: i
      real            :: p_br_def, q_br_def

      namelist /COSMIC_RAY_SPECTRUM/ cfl_cre, p_lo_init, p_up_init, f_init, q_init, q_big, initial_spectrum, p_min_fix, p_max_fix,  &
      &                         cre_eff, cre_active, K_cre_pow, expan_order, e_small, use_cresp, use_cresp_evol,                    &
      &                         e_small_approx_init_cond, p_br_init_lo, e_small_approx_p_lo, e_small_approx_p_up, force_init_NR,    &
      &                         NR_iter_limit, max_p_ratio, synch_active, adiab_active, arr_dim_a, arr_dim_n, arr_dim_q, q_br_init, &
      &                         Gamma_min_fix, Gamma_max_fix, nullify_empty_bins, approx_cutoffs, NR_run_refine_pf, b_max_db,       &
      &                         NR_refine_solution_q, NR_refine_pf_lo, NR_refine_pf_up, smallcree, smallcren, p_br_init_up, p_diff, &
      &                         q_eps, NR_smap_file, cresp_substep, n_substeps_max, allow_unnatural_transfer, icomp_active, redshift

! Default values
      use_cresp         = .true.
      use_cresp_evol    = .true.
      p_min_fix         = 1.5e1
      p_max_fix         = 1.65e4
      p_lo_init         = 1.5e1
      p_up_init         = 7.5e2
      p_br_def          = p_lo_init
      initial_spectrum  = "powl"
      f_init            = 1.0
      q_init            = 4.1
      q_br_def          = q_init
      q_big             = 30.0d0
      p_br_init_lo      = p_br_def ! < in case it was not provided "powl" can be assumed
      p_br_init_up      = p_br_def ! < in case it was not provided "powl" can be assumed
      q_br_init         = q_br_def ! < in case it was not provided "powl" can be assumed
      p_diff            = 10000.0
      cfl_cre           = 0.1
      cre_eff           = 0.01
      K_cre_pow         = 0.
      expan_order       = 1
      Gamma_min_fix     = 2.5
      Gamma_max_fix     = 1000.0

      approx_cutoffs       = .true.
      e_small              = 1.0e-5
      e_small_approx_p_lo  = 1
      e_small_approx_p_up  = 1
      e_small_approx_init_cond = 1
      max_p_ratio          = 2.5
      NR_iter_limit        = 100
      force_init_NR        = .false.
      NR_run_refine_pf     = .false.
      NR_refine_solution_q = .false.
      NR_refine_pf_lo      = .false.
      NR_refine_pf_up      = .false.
      NR_smap_file         = "CRESP_smaps.h5"
      NR_allow_old_smaps   = .false.
      nullify_empty_bins   = .false.
      smallcren            = 0.0
      smallcree            = 0.0
      allow_source_spectrum_break  = .false.
      allow_unnatural_transfer     = .false.
      synch_active         = .true.
      adiab_active         = .true.
      icomp_active         = .false.
      cre_active           = 0.0
      b_max_db             = 10.  ! default value of B limiter
      redshift             = 0.

! NR parameters
      tol_f     = 1.0e-11
      tol_x     = 1.0e-11
      tol_f_1D  = 1.0e-14
      tol_x_1D  = 1.0e-14
      arr_dim_a = 200
      arr_dim_n = 200
      arr_dim_q = 1000
      q_eps     = 0.001

      cresp_substep           = .false.
      n_substeps_max          = 100

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
         ibuff(1)  = expan_order

         ibuff(2)  = e_small_approx_p_lo
         ibuff(3)  = e_small_approx_p_up
         ibuff(4)  = e_small_approx_init_cond

         ibuff(5)  =  NR_iter_limit
         ibuff(6)  =  arr_dim_a
         ibuff(7)  =  arr_dim_n
         ibuff(8)  =  arr_dim_q

         ibuff(9)  =  n_substeps_max

         lbuff(1)  =  use_cresp
         lbuff(2)  =  use_cresp_evol
         lbuff(3)  =  allow_source_spectrum_break
         lbuff(4)  =  synch_active
         lbuff(5)  =  adiab_active
         lbuff(6)  =  icomp_active
         lbuff(7)  =  force_init_NR
         lbuff(8)  =  NR_run_refine_pf
         lbuff(9)  =  NR_refine_solution_q
         lbuff(10) =  NR_refine_pf_lo
         lbuff(11) =  NR_refine_pf_up
         lbuff(12) =  nullify_empty_bins
         lbuff(13) =  approx_cutoffs
         lbuff(14) =  NR_allow_old_smaps

         lbuff(15) =  cresp_substep
         lbuff(16) =  allow_unnatural_transfer

         rbuff(1)  = cfl_cre
         rbuff(2)  = cre_eff
         rbuff(3)  = smallcren
         rbuff(4)  = smallcree
         rbuff(5)  = cre_active
         rbuff(6)  = p_lo_init
         rbuff(7)  = p_up_init
         rbuff(8)  = f_init
         rbuff(9)  = q_init
         rbuff(10) = q_big
         rbuff(11) = p_min_fix
         rbuff(12) = p_max_fix
         rbuff(13) = K_cre_pow

         rbuff(14) = e_small
         rbuff(15) = max_p_ratio

         rbuff(16) = tol_f
         rbuff(17) = tol_x
         rbuff(18) = tol_f_1D
         rbuff(19) = tol_x_1D

         rbuff(20) = Gamma_min_fix
         rbuff(21) = Gamma_max_fix

         rbuff(22) = p_br_init_lo
         rbuff(23) = p_br_init_up
         rbuff(24) = q_br_init
         rbuff(25) = p_diff
         rbuff(26) = q_eps
         rbuff(27) = b_max_db

         rbuff(28) = redshift

         cbuff(1)  = initial_spectrum
      endif

      call piernik_MPI_Bcast(ibuff)
      call piernik_MPI_Bcast(rbuff)
      call piernik_MPI_Bcast(lbuff)
      call piernik_MPI_Bcast(cbuff, cbuff_len)
      call piernik_MPI_Bcast(NR_smap_file, fnamelen)

      if (slave) then
         expan_order                 = int(ibuff(1), kind=4)

         e_small_approx_p_lo         = int(ibuff(2), kind=1)
         e_small_approx_p_up         = int(ibuff(3), kind=1)
         e_small_approx_init_cond    = int(ibuff(4), kind=1)

         NR_iter_limit               = int(ibuff(5), kind=2)
         arr_dim_a                   = int(ibuff(6), kind=4)
         arr_dim_n                   = int(ibuff(7), kind=4)
         arr_dim_q                   = int(ibuff(8), kind=4)

         n_substeps_max              = int(ibuff(9), kind=4)

         use_cresp                   = lbuff(1)
         use_cresp_evol              = lbuff(2)
         allow_source_spectrum_break = lbuff(3)
         synch_active                = lbuff(4)
         adiab_active                = lbuff(5)
         icomp_active                = lbuff(6)
         force_init_NR               = lbuff(7)
         NR_run_refine_pf            = lbuff(8)
         NR_refine_solution_q        = lbuff(9)
         NR_refine_pf_lo             = lbuff(10)
         NR_refine_pf_up             = lbuff(11)
         nullify_empty_bins          = lbuff(12)
         approx_cutoffs              = lbuff(13)
         NR_allow_old_smaps          = lbuff(14)

         cresp_substep               = lbuff(15)
         allow_unnatural_transfer    = lbuff(16)

         cfl_cre                     = rbuff(1)
         cre_eff                     = rbuff(2)
         smallcren                   = rbuff(3)
         smallcree                   = rbuff(4)
         cre_active                  = rbuff(5)
         p_lo_init                   = rbuff(6)
         p_up_init                   = rbuff(7)
         f_init                      = rbuff(8)
         q_init                      = rbuff(9)
         q_big                       = rbuff(10)
         p_min_fix                   = rbuff(11)
         p_max_fix                   = rbuff(12)
         K_cre_pow                   = rbuff(13)

         e_small                     = rbuff(14)
         max_p_ratio                 = rbuff(15)

         tol_f                       = rbuff(16)
         tol_x                       = rbuff(17)
         tol_f_1D                    = rbuff(18)
         tol_x_1D                    = rbuff(19)

         Gamma_min_fix               = rbuff(20)
         Gamma_max_fix               = rbuff(21)

         p_br_init_lo                = rbuff(22)
         p_br_init_up                = rbuff(23)
         q_br_init                   = rbuff(24)
         p_diff                      = rbuff(25)

         q_eps                       = rbuff(26)
         b_max_db                    = rbuff(27)
         redshift                    = rbuff(28)

         initial_spectrum            = trim(cbuff(1))

      endif

      ! since now use only e_small_approx_p, p_init and p_br_init
      e_small_approx_p = [e_small_approx_p_lo, e_small_approx_p_up]
      p_init           = [p_lo_init, p_up_init]
      p_br_init        = [p_br_init_lo, p_br_init_up]
      NR_refine_pf     = [NR_refine_pf_lo, NR_refine_pf_up]

! Input parameters check
      if (use_cresp .and. (ncrb <= I_ZERO .or. nspc == I_ZERO))  then
         write (msg,"(A,I4,A)") '[initcrspectrum:init_cresp] ncrb   = ', ncrb, '; CR bins NOT initnialized. Switching CRESP module off.'
         if (master) call warn(msg)
         use_cresp      = .false.
         use_cresp_evol = .false.
         ncrb           = 0
      endif

      if (.not. use_cresp) then
         if (master) call warn("[initcrspectrum:init_cresp] Switching 'use_cresp_evol' off: superior 'use_cresp' is switched off.")
         use_cresp_evol = .false.
         return
      endif

      if (ncrb < 3) call die("[initcrspectrum:init_cresp] CRESP algorithm currently requires at least 3 bins (ncrb) in order to work properly, check your parameters.")

      if (approx_cutoffs) then
         e_small_approx_p = 1
         write (msg,'(A)') "[initcrspectrum:init_cresp] approx_cutoffs = .true. -- will use e_small to approximate spectrum cutoffs and initial state spectrum."
      else
         e_small_approx_p = 0 ! e_small_approx_init_cond stays default, unless user changes.
         write (msg,'(A)') "[initcrspectrum:init_cresp] approx_cutoffs = .false. -- will not use e_small approximated cutoffs, but still approximate initial state. To turn it off use e_small_approx_init_cond = 0."
      endif
      if (master) call printinfo(msg)

      if (sum(e_small_approx_p) > 0 .and. e_small_approx_init_cond < 1) then
         e_small_approx_init_cond = 1
         if (master) call warn("[initcrspectrum:init_cresp] Approximation of boundary momenta is active -> modifying e_small_approx_init_cond to 1.")
      endif

! countermeasure - in case unrecognized or invalid parameters are provided

      where (e_small_approx_p > 0) ;           e_small_approx_p = 1 ;    elsewhere ; e_small_approx_p = 0 ;      endwhere
      if (e_small_approx_init_cond > 0) then ; e_small_approx_init_cond = 1 ; else ; e_small_approx_init_cond = 0 ; endif

      if (sum(e_small_approx_p) == 0) NR_refine_solution_q = .true. !< for testing we leave precise solutions of q (especially for outer momenta)

      if (e_small_approx_init_cond + sum(e_small_approx_p) == 0) e_small = zero                !< no threshold energy for bin activation necessary

! arrays initialization
      call my_allocate_with_index(p_fix,            ncrb, I_ZERO)
      call my_allocate_with_index(p_mid_fix,        ncrb, I_ONE )
      call my_allocate_with_index(cresp_all_edges,  ncrb, I_ZERO)
      call my_allocate_with_index(cresp_all_bins,   ncrb, I_ONE )
      call my_allocate_with_index(n_small_bin,      ncrb, I_ONE )

      call my_allocate_with_index(Gamma_fix,        ncrb, I_ZERO)
      call my_allocate_with_index(Gamma_mid_fix,    ncrb, I_ONE )
      call my_allocate_with_index(mom_cre_fix,      ncrb, I_ZERO)
      call my_allocate_with_index(mom_mid_cre_fix,  ncrb, I_ONE )
      call my_allocate_with_index(gamma_beta_c_fix, ncrb, I_ZERO)

      call my_allocate_with_index(K_cresp_paral,    ncr2b, I_ONE)
      call my_allocate_with_index(K_cresp_perp,     ncr2b, I_ONE)

      cresp_all_edges = [(i, i = I_ZERO, ncrb)]
      cresp_all_bins  = [(i, i = I_ONE,  ncrb)]

!!\brief for now algorithm requires at least 3 bins
      p_fix = zero
      w  = log10(p_max_fix/p_min_fix) / real(ncrb-2)
      p_fix(1:ncrb-1) = p_min_fix*ten**(w*real(cresp_all_edges(1:ncrb-1)-1))
      p_fix(0)    = zero
      p_fix(ncrb) = zero
      p_fix_ratio = ten**w

      p_mid_fix = 0.0
      p_mid_fix(2:ncrb-1) = sqrt(p_fix(1:ncrb-2)*p_fix(2:ncrb-1))
      p_mid_fix(1)    = p_mid_fix(2) / p_fix_ratio
      p_mid_fix(ncrb) = p_mid_fix(ncrb-1) * p_fix_ratio

!> set Gamma arrays, analogically to p_fix arrays, that will be constructed using Gamma arrays
      Gamma_fix            = one             !< Gamma factor obviously cannot be lower than 1
      G_w                  = log10(Gamma_max_fix/Gamma_min_fix) / real(ncrb-2)
      Gamma_fix(1:ncrb-1)  = Gamma_min_fix * ten**(G_w * real(cresp_all_edges(1:ncrb-1)-1))
      Gamma_fix_ratio      = ten**w

      Gamma_mid_fix = one
      Gamma_mid_fix(2:ncrb-1) = sqrt( Gamma_fix(1:ncrb-2)   * Gamma_fix(2:ncrb-1) )
      Gamma_mid_fix(1)        = sqrt( Gamma_mid_fix(1)      * Gamma_mid_fix(2))
      Gamma_mid_fix(ncrb)     = sqrt( Gamma_mid_fix(ncrb-1) * Gamma_mid_fix(ncrb-1) * Gamma_fix_ratio )
! compute physical momenta of particles in given unit set
      mom_cre_fix      = [(cresp_get_mom(Gamma_fix(i),me),     i = I_ZERO, ncrb )]
      mom_mid_cre_fix  = [(cresp_get_mom(Gamma_mid_fix(i),me), i = I_ONE,  ncrb )]

      gamma_beta_c_fix = mom_cre_fix / me

      n_small_bin(:) = e_small / (p_mid_fix(:) * clight_cresp)

#ifdef VERBOSE
         write (msg,'(A, 50I3)')    '[initcrspectrum:init_cresp] fixed all edges: ', cresp_all_edges
         call printinfo(msg)
         write (msg,'(A, 50E15.7)') '[initcrspectrum:init_cresp] Fixed momentum grid: ', p_fix
         call printinfo(msg)
         write (msg,'(A, 50E15.7)') '[initcrspectrum:init_cresp] Bin p-width (log10): ', w
         call printinfo(msg)
         write (msg,'(A, 50E15.7)') '[initcrspectrum:init_cresp] Fixed momentum grid (bin middle):   ',p_mid_fix(1:ncrb)
         call printinfo(msg)
         write (msg,'(A, 50F13.2)') '[initcrspectrum:init_cresp] Fixed Gamma      grid: ', Gamma_fix
         call printinfo(msg)
         write (msg,'(A, 50F10.5)') '[initcrspectrum:init_cresp] Gamma bin width(log10): ', G_w
         call printinfo(msg)
         write (msg,'(A, 50F10.5)') '[initcrspectrum:init_cresp] Fixed mid-Gamma     : ', Gamma_mid_fix(1:ncrb)
         call printinfo(msg)
         write (msg,'(A, 50E15.7)') '[initcrspectrum:init_cresp] Fixed phys momentum : ', mom_cre_fix
         call printinfo(msg)
         write (msg,'(A, 50E15.7)') '[initcrspectrum:init_cresp] Fixed phys mid mom  : ', mom_mid_cre_fix
         call printinfo(msg)
#endif /* VERBOSE */

      !> check correctness of "initial_spectrum" is checked here
      if (f_init < eps) then
         if (initial_spectrum == 'powl' .or. initial_spectrum == 'brpl') then
            write (msg,"(A,A,A)") "[initcrspectrum:init_cresp] Provided power law type spectrum (",initial_spectrum,") with initial amplitude f_init ~ zero. Check your parameters."
            call die(msg)
         endif
      endif

      if (initial_spectrum == "brpl" ) then
         if (abs(p_br_init(LO) - p_br_def) <= eps) then
            write (msg,"(A)") "[initcrspectrum:init_cresp] Parameter for 'brpl' spectrum: p_br_init_lo has default value (probably unitialized). Assuming p_lo_init value ('powl' spectrum)."
            if (master) call warn(msg)
         else
            !> p_br_init_lo should be equal to one of p_fix values
            i = int(minloc(abs(p_fix - p_br_init(LO)), dim=1), kind=4) - I_ONE
            write (msg,"(A,E14.7,1A)") "[initcrspectrum:init_cresp] p_br_init_lo was set, but should be equal to one of p_fix. Assuming p_br_init_lo =", p_fix(i),"."
            p_br_init(LO) = p_fix(i)
            if (master) call warn(msg)
         endif
         if (abs(q_br_init - q_br_def) <= eps) then
            write (msg,"(A)") "[initcrspectrum:init_cresp] Parameter for 'brpl' spectrum: q_br_init has default value (probably unitialized). Assuming q_init value ('powl' spectrum)."
            if (master) call warn(msg)
         endif
      endif

      if (initial_spectrum == "plpc") then
         if (abs(p_br_init(LO) - p_br_def) <= eps .or. abs(p_br_init(HI) - p_br_def) <= eps) then
            write (msg,"(A)") "[initcrspectrum:init_cresp] Parameters for 'plpc' spectrum: p_br_init_lo or p_br_init_up has default value (probably unitialized). Check spectrum parameters."
            if (master) call die(msg)
         else
            !> p_br_init_lo should be equal to one of p_fix values
            i = int(minloc(abs(p_fix - p_br_init(LO)), dim=1), kind=4) - I_ONE
            write (msg,"(A,E14.7,1A)") "[initcrspectrum:init_cresp] p_br_init_lo was set, but should be equal to one of p_fix. Assuming p_br_init_lo =", p_fix(i),"."
            p_br_init(LO) = p_fix(i)
            if (master) call warn(msg)
            !> p_br_init_up should also be equal to one of p_fix values
            i = int(minloc(abs(p_fix - p_br_init(HI)), dim=1), kind=4) - I_ONE
            write (msg,"(A,E14.7,1A)") "[initcrspectrum:init_cresp] p_br_init_up was set, but should be equal to one of p_fix. Assuming p_br_init_up =", p_fix(i),"."
            p_br_init(HI) = p_fix(i)
            if (master) call warn(msg)
         endif
      endif

      call init_cresp_types

      K_cresp_paral(1:ncrb) = K_cr_paral(cr_table(icr_E)) * (p_mid_fix(1:ncrb) / p_diff)**K_cre_pow
      K_cresp_perp(1:ncrb)  = K_cr_perp(cr_table(icr_E))  * (p_mid_fix(1:ncrb) / p_diff)**K_cre_pow

#ifdef VERBOSE
      write (msg,"(A,*(E14.5))") "[initcrspectrum:init_cresp] K_cresp_paral = ", K_cresp_paral(1:ncrb) ; if (master) call printinfo(msg)
      write (msg,"(A,*(E14.5))") "[initcrspectrum:init_cresp] K_cresp_perp = ",  K_cresp_perp(1:ncrb)  ; if (master) call printinfo(msg)
#endif /* VERBOSE */

      K_cresp_paral(ncrb+1:ncr2b) = K_cresp_paral(1:ncrb)
      K_cresp_perp (ncrb+1:ncr2b) = K_cresp_perp (1:ncrb)
      K_crs_paral(ncrn+1:ncrtot) = K_cresp_paral(1:ncr2b)
      K_crs_perp (ncrn+1:ncrtot) = K_cresp_perp (1:ncr2b)

      f_synchIC = (4. / 3. ) * sigma_T / (me * clight)
      write (msg, *) "[initcrspectrum:init_cresp] 4/3 * sigma_T / ( me * c ) = ", f_synchIC

      def_dtadiab   = cfl_cre * half * three * logten * w
      def_dtsynchIC = cfl_cre * half * w

      if (master) call printinfo(msg)

      u_b_max = f_synchIC * emag(b_max_db, 0., 0.)   !< initializes factor for comparing u_b with u_b_max

      write (msg, "(A,F10.4,A,ES12.5)") "[initcrspectrum:init_cresp] Maximal B_tot =",b_max_db, "mGs, u_b_max = ", u_b_max
      if (master)  call warn(msg)
      write (msg, "(A,ES12.5,A,ES15.8,A,ES15.8)") "[initcrspectrum:init_cresp] dt_synch(p_max_fix = ",p_max_fix,", u_b_max = ",u_b_max,") = ", def_dtsynchIC / (p_max_fix* u_b_max)
      if (master)  call warn(msg)

      if (cresp_substep) then
         if (master) then
            write(msg,"(A, I4)") "[initcrspectrum:init_cresp] Substep for CRESP for each cell is ON, max. substeps: ", n_substeps_max
            call printinfo(msg)
         endif
      else
         n_substeps_max = 1            !< for sanity assuming 1 substep if cresp_substep = .false.
      endif

      if (.not. disallow_CRnegatives) then
         if (.not. use_smallecr) then
            if (master) call warn("[initcrspectrum:init_cresp] Detecting negative values of n,e in CRESP module & performing CFL violation actions related is DISABLED via disallow_CRnegatives.")
            if (master) call warn("[initcrspectrum:init_cresp] as is 'use_smallecr'; should negative values show in CRESP, they will not be fixed.")
         else
            if (master) call warn("[initcrspectrum:init_cresp] Detecting negative values of n,e in CRESP module & performing CFL violation actions related is DISABLED via disallow_CRnegatives.")
         endif
      endif

      if ((q_init < three) .and. any(e_small_approx_p == I_ONE)) then
         call warn("[initcrspectrum:init_cresp] Initial parameters: q_init < 3.0 and approximation of outer momenta is on, approximation of outer momenta with hard energy spectrum might not work.")
      endif

   end subroutine init_cresp

!----------------------------------------------------------------------------------------------------

   subroutine cleanup_cresp_sp

      use diagnostics, only: my_deallocate

      implicit none

      call my_deallocate(p_fix)
      call my_deallocate(p_mid_fix)
      call my_deallocate(cresp_all_edges)
      call my_deallocate(cresp_all_bins)
      call my_deallocate(n_small_bin)

      call my_deallocate(Gamma_fix)
      call my_deallocate(Gamma_mid_fix)
      call my_deallocate(mom_cre_fix)
      call my_deallocate(mom_mid_cre_fix)
      call my_deallocate(gamma_beta_c_fix)

      call my_deallocate(K_cresp_paral)
      call my_deallocate(K_cresp_perp)

   end subroutine cleanup_cresp_sp

!----------------------------------------------------------------------------------------------------

   subroutine init_cresp_types

      use constants,      only: zero, I_ONE
      use diagnostics,    only: my_allocate_with_index
      use initcosmicrays, only: ncrb

      implicit none

      if (.not. allocated(cresp%n)) call my_allocate_with_index(cresp%n, ncrb, I_ONE)
      if (.not. allocated(cresp%e)) call my_allocate_with_index(cresp%e, ncrb, I_ONE)
      if (.not. allocated(norm_init_spectrum%n)) call my_allocate_with_index(norm_init_spectrum%n, ncrb, I_ONE)
      if (.not. allocated(norm_init_spectrum%e)) call my_allocate_with_index(norm_init_spectrum%e, ncrb, I_ONE)

      cresp%e = zero
      cresp%n = zero

      norm_init_spectrum%n = zero
      norm_init_spectrum%e = zero

   end subroutine init_cresp_types

!----------------------------------------------------------------------------------------------------

   subroutine init_crel

      use constants,      only: zero, I_ZERO, I_ONE
      use diagnostics,    only: my_allocate_with_index
      use initcosmicrays, only: ncrb

      implicit none

      if (.not. allocated(crel%p)) call my_allocate_with_index(crel%p, ncrb, I_ZERO)
      if (.not. allocated(crel%f)) call my_allocate_with_index(crel%f, ncrb, I_ZERO)
      if (.not. allocated(crel%q)) call my_allocate_with_index(crel%q, ncrb, I_ONE )
      if (.not. allocated(crel%n)) call my_allocate_with_index(crel%n, ncrb, I_ONE )
      if (.not. allocated(crel%e)) call my_allocate_with_index(crel%e, ncrb, I_ONE )

      crel%p = zero
      crel%q = zero
      crel%f = zero
      crel%e = zero
      crel%n = zero
      crel%i_cut = I_ZERO

   end subroutine init_crel
!----------------------------------------------------------------------------------------------------

   real function cresp_get_mom(gamma, particle_mass)

      use constants, only: zero, one
      use units,     only: clight

      implicit none

      real, intent(in) :: gamma, particle_mass

      cresp_get_mom = zero
      if (gamma > one) cresp_get_mom = particle_mass * sqrt(gamma**2 - one) * clight

   end function cresp_get_mom

!----------------------------------------------------------------------------------------------------
!>
!! \todo try to merge this to common_hdf5:init_hdf5
!<
   subroutine check_if_dump_fpq(vars)

      use constants, only: dsetnamelen

      implicit none

      character(len=dsetnamelen), dimension(:), intent(in) :: vars  !< quantities to be plotted, see dataio::vars
      integer :: i

      do i = lbound(vars, 1), ubound(vars, 1)
         select case (trim(vars(i)))
            case ('cref') !< CRESP distribution function
               dfpq%dump_f   = .true.
               dfpq%any_dump = .true.
            case ('crep') !< CRESP cutoff momenta
               dfpq%dump_p   = .true.
               dfpq%any_dump = .true.
            case ('creq') !< CRESP spectrum index
               dfpq%dump_q   = .true.
               dfpq%any_dump = .true.
         end select
      enddo

      if (dfpq%any_dump) call init_crel

   end subroutine check_if_dump_fpq

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

      if (dfpq%any_dump) then
         if (allocated(crel%p)) call my_deallocate(crel%p)
         if (allocated(crel%f)) call my_deallocate(crel%f)
         if (allocated(crel%q)) call my_deallocate(crel%q)
         if (allocated(crel%e)) call my_deallocate(crel%e)
         if (allocated(crel%n)) call my_deallocate(crel%n)
      endif

   end subroutine cleanup_cresp_work_arrays

!----------------------------------------------------------------------------------------------------

end module initcrspectrum
