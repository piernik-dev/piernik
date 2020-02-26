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
   use constants, only: cbuff_len

   implicit none

   private
   public :: use_cresp, p_init, initial_spectrum, p_br_init, f_init, q_init, q_br_init, q_big, cfl_cre, cre_eff, expan_order, e_small, e_small_approx_p, e_small_approx_init_cond,  &
           & smallcren, smallcree, max_p_ratio, NR_iter_limit, force_init_NR, NR_run_refine_pf, NR_refine_solution_q, NR_refine_pf, nullify_empty_bins, synch_active, adiab_active, &
           & allow_source_spectrum_break, cre_active, tol_f, tol_x, tol_f_1D, tol_x_1D, arr_dim, arr_dim_q, eps, eps_det, w, p_fix, p_mid_fix, total_init_cree, p_fix_ratio,        &
           & spec_mod_trms, cresp_all_edges, cresp_all_bins, norm_init_spectrum, cresp, crel, dfpq, fsynchr, init_cresp, check_if_dump_fpq, cleanup_cresp_work_arrays, q_eps, u_b_max

! contains routines reading namelist in problem.par file dedicated to cosmic ray electron spectrum and initializes types used.
! available via namelist COSMIC_RAY_SPECTRUM
   logical         :: use_cresp                   !< determines whether CRESP update is called by fluidupdate
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
   real            :: K_cre_paral_1               !< maximal parallell diffusion coefficient value
   real            :: K_cre_perp_1                !< maximal perpendicular diffusion coefficient value
   real            :: K_cre_pow                   !< exponent for power law-like diffusion-energy dependance
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
   logical, dimension(2) :: NR_refine_pf          !< vector to store NR_refine_pf_lo and NR_refine_pf_up

   logical         :: nullify_empty_bins          !< nullifies empty bins when entering CRESP module / exiting empty cell.
   logical         :: allow_source_spectrum_break !< allow extension of spectrum to adjacent bins if momenta found exceed set p_fix
   logical         :: synch_active                !< TEST feature - turns on / off synchrotron cooling @ CRESP
   logical         :: adiab_active                !< TEST feature - turns on / off adiabatic   cooling @ CRESP
   real            :: cre_active                  !< electron contribution to Pcr

! NR parameters
   real            :: tol_f                       !< tolerance for f abs. error in NR algorithm
   real            :: tol_x                       !< tolerance for x abs. error in NR algorithm
   real            :: tol_f_1D                    !< tolerance for f abs. error in NR algorithm (1D)
   real            :: tol_x_1D                    !< tolerance for x abs. error in NR algorithm (1D)
   integer(kind=4) :: arr_dim, arr_dim_q
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
      real :: ucmb
   end type spec_mod_trms

   real :: total_init_cree
   real :: p_fix_ratio
   integer, allocatable, dimension(:) :: cresp_all_edges, cresp_all_bins

! CRESP names
   type dump_fpq_type
      character(len=4) :: f_nam = "cref" !< helping array for CRESP number density
      character(len=4) :: p_nam = "crep" !< helping array for CRESP energy density
      character(len=4) :: q_nam = "creq" !< helping array for CRESP energy density
      logical :: any_dump, dump_f, dump_p, dump_q  ! diagnostic, if true - adding 'cref', 'crep', 'creq' to hdf_vars must follow
   end type dump_fpq_type
   type(dump_fpq_type) :: dfpq

   real    :: fsynchr

!====================================================================================================
!
 contains
!
!====================================================================================================
   subroutine init_cresp

      use constants,       only: cbuff_len, I_ZERO, I_ONE, zero, one, three, ten, LO, HI
      use cresp_variables, only: clight_cresp
      use dataio_pub,      only: printinfo, warn, msg, die, nh
      use diagnostics,     only: my_allocate_with_index
      use func,            only: emag
      use initcosmicrays,  only: ncrn, ncre, K_crs_paral, K_crs_perp, K_cre_paral, K_cre_perp
      use mpisetup,        only: rbuff, ibuff, lbuff, cbuff, master, slave, piernik_MPI_Bcast
      use units,           only: clight, me, sigma_T, mGs

      implicit none

      integer :: i
      real    :: p_br_def, q_br_def

      namelist /COSMIC_RAY_SPECTRUM/ cfl_cre, p_lo_init, p_up_init, f_init, q_init, q_big, initial_spectrum, p_min_fix, p_max_fix, &
      &                         cre_eff, K_cre_paral_1, K_cre_perp_1, cre_active, K_cre_pow, expan_order, e_small, use_cresp,      &
      &                         e_small_approx_init_cond, p_br_init_lo, e_small_approx_p_lo, e_small_approx_p_up, force_init_NR,   &
      &                         NR_iter_limit, max_p_ratio, synch_active, adiab_active, arr_dim, arr_dim_q, q_br_init,             &
      &                         Gamma_min_fix, Gamma_max_fix, nullify_empty_bins, approx_cutoffs, NR_run_refine_pf, b_max_db,      &
      &                         NR_refine_solution_q, NR_refine_pf_lo, NR_refine_pf_up, smallcree, smallcren, p_br_init_up, p_diff, q_eps

! Default values
      use_cresp         = .true.
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
      K_cre_paral_1     = 0.
      K_cre_perp_1      = 0.
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
      nullify_empty_bins   = .false.
      smallcren            = 0.0
      smallcree            = 0.0
      allow_source_spectrum_break  = .false.
      synch_active         = .true.
      adiab_active         = .true.
      cre_active           = 0.0
      b_max_db             = 10.  ! default value of B limiter
! NR parameters
      tol_f    = 1.0e-11
      tol_x    = 1.0e-11
      tol_f_1D = 1.0e-14
      tol_x_1D = 1.0e-14
      arr_dim  = 200
      arr_dim_q = 500
      q_eps     = eps

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
         ibuff(6)  =  arr_dim
         ibuff(7)  =  arr_dim_q

         lbuff(1)  =  use_cresp
         lbuff(2)  =  allow_source_spectrum_break
         lbuff(3)  =  synch_active
         lbuff(4)  =  adiab_active

         lbuff(5)  =  force_init_NR
         lbuff(6)  =  NR_run_refine_pf
         lbuff(7)  =  NR_refine_solution_q
         lbuff(8)  =  NR_refine_pf_lo
         lbuff(9)  =  NR_refine_pf_up
         lbuff(10) =  nullify_empty_bins
         lbuff(11) =  approx_cutoffs

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

         rbuff(24) = p_br_init_lo
         rbuff(25) = p_br_init_up
         rbuff(26) = q_br_init
         rbuff(27) = p_diff
         rbuff(28) = q_eps
         rbuff(29) = b_max_db

         cbuff(1)  = initial_spectrum
      endif

      call piernik_MPI_Bcast(ibuff)
      call piernik_MPI_Bcast(rbuff)
      call piernik_MPI_Bcast(lbuff)
      call piernik_MPI_Bcast(cbuff, cbuff_len)

      if (slave) then
         expan_order                 = int(ibuff(1),kind=4)

         e_small_approx_p_lo         = int(ibuff(2),kind=1)
         e_small_approx_p_up         = int(ibuff(3),kind=1)
         e_small_approx_init_cond    = int(ibuff(4),kind=1)

         NR_iter_limit               = int(ibuff(5),kind=2)
         arr_dim                     = int(ibuff(6),kind=4)
         arr_dim_q                   = int(ibuff(7),kind=4)

         use_cresp                   = lbuff(1)
         allow_source_spectrum_break = lbuff(2)
         synch_active                = lbuff(3)
         adiab_active                = lbuff(4)

         force_init_NR               = lbuff(5)
         NR_run_refine_pf            = lbuff(6)
         NR_refine_solution_q        = lbuff(7)
         NR_refine_pf_lo             = lbuff(8)
         NR_refine_pf_up             = lbuff(9)
         nullify_empty_bins          = lbuff(10)
         approx_cutoffs              = lbuff(11)

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

         p_br_init_lo                = rbuff(24)
         p_br_init_up                = rbuff(25)
         q_br_init                   = rbuff(26)
         p_diff                      = rbuff(27)

         q_eps                       = rbuff(28)
         b_max_db                    = rbuff(29)
         initial_spectrum            = trim(cbuff(1))

      endif

      ! since now use only e_small_approx_p, p_init and p_br_init
      e_small_approx_p = [e_small_approx_p_lo, e_small_approx_p_up]
      p_init           = [p_lo_init, p_up_init]
      p_br_init        = [p_br_init_lo, p_br_init_up]
      NR_refine_pf     = [NR_refine_pf_lo, NR_refine_pf_up]

! Input parameters check
      if (ncre < 3) then
         if (ncre <= I_ZERO)  then
            write (msg,"(A,I4,A)") '[initcrspectrum:init_cresp] ncre   = ', ncre, '; cr-electrons NOT initnialized. If COSM_RAY_ELECTRONS flag is on, please check your parameters.'
            call die(msg)
         endif
         call die("[initcrspectrum:init_cresp] CRESP algorithm currently requires at least 3 bins (ncre) in order to work properly, check your parameters.")
      endif

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
      call my_allocate_with_index(p_fix,           ncre, 0)
      call my_allocate_with_index(p_mid_fix,       ncre, 1)
      call my_allocate_with_index(cresp_all_edges, ncre, 0)
      call my_allocate_with_index(cresp_all_bins,  ncre, 1)
      call my_allocate_with_index(n_small_bin,     ncre, 1)

      call my_allocate_with_index(Gamma_fix,        ncre, 0)
      call my_allocate_with_index(Gamma_mid_fix,    ncre, 1)
      call my_allocate_with_index(mom_cre_fix,      ncre, 0)
      call my_allocate_with_index(mom_mid_cre_fix,  ncre, 1)
      call my_allocate_with_index(gamma_beta_c_fix, ncre, 0)

      cresp_all_edges = [(i, i = 0, ncre)]
      cresp_all_bins  = [(i, i = 1, ncre)]

!!\brief for now algorithm requires at least 3 bins
      p_fix = zero
      w  = log10(p_max_fix/p_min_fix) / real(ncre-2)
      p_fix(1:ncre-1) = p_min_fix*ten**(w*real(cresp_all_edges(1:ncre-1)-1))
      p_fix(0)    = zero
      p_fix(ncre) = zero
      p_fix_ratio = ten**w

      p_mid_fix = 0.0
      p_mid_fix(2:ncre-1) = sqrt(p_fix(1:ncre-2)*p_fix(2:ncre-1))
      p_mid_fix(1)    = p_mid_fix(2) / p_fix_ratio
      p_mid_fix(ncre) = p_mid_fix(ncre-1) * p_fix_ratio

!> set Gamma arrays, analogically to p_fix arrays, that will be constructed using Gamma arrays
      Gamma_fix            = one             !< Gamma factor obviously cannot be lower than 1
      G_w                  = log10(Gamma_max_fix/Gamma_min_fix) / real(ncre-2)
      Gamma_fix(1:ncre-1)  = Gamma_min_fix * ten**(G_w * real(cresp_all_edges(1:ncre-1)-1))
      Gamma_fix_ratio      = ten**w

      Gamma_mid_fix = one
      Gamma_mid_fix(2:ncre-1) = sqrt( Gamma_fix(1:ncre-2)   * Gamma_fix(2:ncre-1) )
      Gamma_mid_fix(1)        = sqrt( Gamma_mid_fix(1)      * Gamma_mid_fix(2))
      Gamma_mid_fix(ncre)     = sqrt( Gamma_mid_fix(ncre-1) * Gamma_mid_fix(ncre-1) * Gamma_fix_ratio )
! compute physical momenta of particles in given unit set
      mom_cre_fix      = [(cresp_get_mom(Gamma_fix(i),me),     i = 0, ncre )]
      mom_mid_cre_fix  = [(cresp_get_mom(Gamma_mid_fix(i),me), i = 1, ncre )]

      gamma_beta_c_fix = mom_cre_fix / me

      n_small_bin(:) = e_small / (p_mid_fix(:) * clight_cresp)

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
            i = minloc(abs(p_fix - p_br_init(LO)),dim=1)-1
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
            i = minloc(abs(p_fix - p_br_init(LO)),dim=1)-1
            write (msg,"(A,E14.7,1A)") "[initcrspectrum:init_cresp] p_br_init_lo was set, but should be equal to one of p_fix. Assuming p_br_init_lo =", p_fix(i),"."
            p_br_init(LO) = p_fix(i)
            if (master) call warn(msg)
            !> p_br_init_up should also be equal to one of p_fix values
            i = minloc(abs(p_fix - p_br_init(HI)),dim=1)-1
            write (msg,"(A,E14.7,1A)") "[initcrspectrum:init_cresp] p_br_init_up was set, but should be equal to one of p_fix. Assuming p_br_init_up =", p_fix(i),"."
            p_br_init(HI) = p_fix(i)
            if (master) call warn(msg)
         endif
      endif

      call init_cresp_types

      K_cre_paral(1:ncre) = K_cre_paral_1 * (p_fix(0:ncre-1) / p_diff)**K_cre_pow
      K_cre_paral(1)      = K_cre_paral_1 * (p_fix(1) / p_fix_ratio / p_diff)**K_cre_pow

      K_cre_perp(1:ncre)  = K_cre_perp_1  * (p_fix(0:ncre-1) / p_diff)**K_cre_pow
      K_cre_perp(1)       = K_cre_perp_1  * (p_fix(1) / p_fix_ratio / p_diff)**K_cre_pow
#ifdef VERBOSE
      write (msg,"(A,*(E14.5))") "[initcrspectrum:init_cresp] K_cre_paral = ", K_cre_paral(1:ncre) ; if (master) call printinfo(msg)
      write (msg,"(A,*(E14.5))") "[initcrspectrum:init_cresp] K_cre_perp = ", K_cre_perp(1:ncre)   ; if (master) call printinfo(msg)
#endif /* VERBOSE */

      K_cre_paral(ncre+1:2*ncre)      = K_cre_paral(1:ncre)
      K_cre_perp (ncre+1:2*ncre)      = K_cre_perp (1:ncre)
      K_crs_paral(ncrn+1:ncrn+2*ncre) = K_cre_paral(1:2*ncre)
      K_crs_perp (ncrn+1:ncrn+2*ncre) = K_cre_perp (1:2*ncre)

      fsynchr =  (4. / 3. ) * sigma_T / (me * clight)
      write (msg, *) "[initcrspectrum:init_cresp] 4/3 * sigma_T / ( me * c ) = ", fsynchr

      if (master) call printinfo(msg)

      u_b_max = fsynchr * emag(b_max_db, 0., 0.)   !< initializes factor for comparising u_b with u_b_max

      write (msg, "(A,F10.4,A,ES12.5)") "[initcrspectrum:init_cresp] Maximal B_tot =",b_max_db, "mGs, u_b_max = ", u_b_max
      if (master)  call warn(msg)
      write (msg, "(A,ES12.5,A,ES15.8,A,ES15.8)") "[initcrspectrum:init_cresp] dt_synch(p_max_fix = ",p_max_fix,", u_b_max = ",u_b_max,") = ", cfl_cre * w / (p_max_fix* u_b_max)
      if (master)  call warn(msg)

      if ((q_init < three) .and. any(e_small_approx_p == I_ONE)) then
         call warn("[initcrspectrum:init_cresp] Initial parameters: q_init < 3.0 and approximation of outer momenta is on, approximation of outer momenta with hard energy spectrum might not work.")
      endif

   end subroutine init_cresp

!----------------------------------------------------------------------------------------------------

   subroutine init_cresp_types

      use constants,      only: zero
      use diagnostics,    only: my_allocate_with_index
      use initcosmicrays, only: ncre

      implicit none

      if (.not. allocated(cresp%n)) call my_allocate_with_index(cresp%n,ncre,1)
      if (.not. allocated(cresp%e)) call my_allocate_with_index(cresp%e,ncre,1)
      if (.not. allocated(norm_init_spectrum%n)) call my_allocate_with_index(norm_init_spectrum%n,ncre,1)
      if (.not. allocated(norm_init_spectrum%e)) call my_allocate_with_index(norm_init_spectrum%e,ncre,1)

      cresp%e = zero
      cresp%n = zero

      norm_init_spectrum%n = zero
      norm_init_spectrum%e = zero

   end subroutine init_cresp_types

!----------------------------------------------------------------------------------------------------

   subroutine init_crel

      use constants,      only: zero, I_ZERO
      use diagnostics,    only: my_allocate_with_index
      use initcosmicrays, only: ncre

      implicit none

      if (.not. allocated(crel%p)) call my_allocate_with_index(crel%p,ncre,0)
      if (.not. allocated(crel%f)) call my_allocate_with_index(crel%f,ncre,0)
      if (.not. allocated(crel%q)) call my_allocate_with_index(crel%q,ncre,1)
      if (.not. allocated(crel%n)) call my_allocate_with_index(crel%n,ncre,1)
      if (.not. allocated(crel%e)) call my_allocate_with_index(crel%e,ncre,1)

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
