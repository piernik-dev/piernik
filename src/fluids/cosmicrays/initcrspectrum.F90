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
           & smallcren, smallcree, max_p_ratio, NR_iter_limit, force_init_NR, NR_run_refine_pf, NR_refine_solution_q, NR_refine_pf, nullify_empty_bins, synch_active, adiab_active, &
           & allow_source_spectrum_break, cre_active, tol_f, tol_x, tol_f_1D, tol_x_1D, arr_dim, arr_dim_q, eps, eps_det, w, p_fix, p_mid_fix, total_init_cree, p_fix_ratio,        &
           & spec_mod_trms, cresp_all_edges, cresp_all_bins, cresp, crel, dfpq, fsynchr, init_cresp, cleanup_cresp_sp, check_if_dump_fpq, cleanup_cresp_work_arrays, q_eps,       &
           & u_b_max, def_dtsynch, def_dtadiab, write_cresp_to_restart, NR_smap_file, NR_allow_old_smaps, cresp_substep, n_substeps_max, allow_unnatural_transfer, K_cresp_paral, &
           & K_cresp_perp, norm_init_spectrum_n, norm_init_spectrum_e, bin_old

! contains routines reading namelist in problem.par file dedicated to cosmic ray electron spectrum and initializes types used.
! available via namelist COSMIC_RAY_SPECTRUM
   logical         :: use_cresp                   !< determines whether CRESP routines are called anywhere
   logical         :: use_cresp_evol              !< determines whether CRESP update is called by fluidupdate
   real            :: p_min_fix                   !< fixed momentum grid lower cutoff
   real            :: p_max_fix                   !< fixed momentum grid upper cutoff
   real, dimension(:), allocatable :: p_lo_init                   !< initial lower cutoff momentum
   real, dimension(:), allocatable :: p_up_init                   !< initial upper cutoff momentum
   real, dimension(:,:), allocatable :: p_init                   !< vector to store p_lo_init and p_up_init
   character(len=cbuff_len) :: initial_spectrum   !< available types: bump, powl, brpl, symf, syme. Description below.
   real, dimension(:), allocatable :: p_br_init_lo, p_br_init_up  !< initial low energy break
   real, dimension(:,:), allocatable :: p_br_init                !< vector to store p_br_init_lo and p_br_init_up
   real, dimension(:), allocatable :: f_init                      !< initial value of distr. func. for isolated case
   real, dimension(:), allocatable :: q_init                      !< initial value of power law-like spectrum exponent
   real, dimension(:), allocatable :: q_br_init                   !< initial q for low energy break
   real            :: q_big                       !< maximal amplitude of q
   real, dimension(:), allocatable :: cfl_cre                     !< CFL parameter  for CR spectrally resolved components! TODO FIXME RENAME ME PLEASE!!!!
   real, dimension(:), allocatable :: cre_eff                     !< fraction of energy passed to CR spectrally resolved components by nucleons (mainly protons) ! TODO wat do with cr_eff now?! TODO FIXME RENAME ME PLEASE!!!!
   real, dimension(:,:), allocatable :: K_cresp_paral !< array containing parallel diffusion coefficients of all CR CRESP components (number density and energy density)
   real, dimension(:,:), allocatable :: K_cresp_perp  !< array containing perpendicular diffusion coefficients of all CR CRESP components (number density and energy density)
   real, dimension(:), allocatable :: K_cre_pow   !< exponent for power law-like diffusion-energy dependence ! TODO FIXME RENAME ME PLEASE!!!!
   real, dimension(:), allocatable :: p_diff                      !< momentum to which diffusion coefficients refer to
   integer(kind=4) :: expan_order                 !< 1,2,3 order of Taylor expansion for p_update (cresp_crspectrum)
   real            :: e_small                     !< lower energy cutoff for energy-approximated cutoff momenta
   logical         :: approx_cutoffs              !< T,F - turns off/on all approximating terms
   integer(kind=4), dimension(2) :: e_small_approx_p !< vector to store e_small_approx_p_lo and e_approx_p_up
   integer(kind=4) :: e_small_approx_p_lo         !< 0,1 - turns off/on energy (e_small) approximated lower cutoff momentum in isolated case
   integer(kind=4) :: e_small_approx_p_up         !< 0,1 - turns off/on energy (e_small) approximated upper cutoff momentum in isolated case
   integer(kind=1) :: e_small_approx_init_cond    !< 0,1 - turns off/on energy (e_small) approximated momenta at initialization
   real            :: smallcren                   !< floor value for CRESP number density ! TODO FIXME RENAME ME PLEASE!!!!   MULTIDIMENSIONALITY OPTIONAL
   real            :: smallcree                   !< floor value for CRESP energy density ! TODO FIXME RENAME ME PLEASE!!!!   MULTIDIMENSIONALITY OPTIONAL
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
   logical, dimension(:), allocatable :: synch_active !< TEST feature - turns on / off synchrotron cooling @ CRESP
   logical, dimension(:), allocatable :: adiab_active !< TEST feature - turns on / off adiabatic   cooling @ CRESP
   real,    dimension(:), allocatable :: cre_active   !< electron contribution to Pcr ! TODO FIXME RENAME ME PLEASE!!!!

! substepping parameters
   logical         :: cresp_substep               !< turns on / off usage of substepping for each cell independently
   integer(kind=4) :: n_substeps_max              !< maximal allowed number of substeps

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
   real, dimension(:,:), allocatable :: norm_init_spectrum_n !TODO FIXME change me back into type(cr_spectrum)
   real, dimension(:,:), allocatable :: norm_init_spectrum_e !TODO FIXME change me back into type(cr_spectrum)
   type(bin_old)     :: crel
   type(bin_old), allocatable, dimension(:) :: bins_primary


! For passing terms to compute energy sources / sinks

   type spec_mod_trms
      real :: ub
      real :: ud
      real :: ucmb
   end type spec_mod_trms

   real, dimension(:), allocatable :: total_init_cree
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

   real, allocatable, dimension(:) :: fsynchr, def_dtadiab, def_dtsynch

!====================================================================================================
!
contains
!
!====================================================================================================
   subroutine init_cresp

      use constants,       only: cbuff_len, I_ZERO, I_ONE, zero, one, three, ten, half, logten, LO, HI
      use cr_data,         only: cr_mass, cr_sigma_N, cr_names, icr_spc, icr_H1, cr_spectral
      use cresp_variables, only: clight_cresp
      use dataio_pub,      only: printinfo, warn, msg, die, nh
      use diagnostics,     only: my_allocate_with_index, my_allocate, ma1d
      use global,          only: disallow_CRnegatives
      use func,            only: emag
      use initcosmicrays,  only: ncrb, ncr2b, ncrn, nspc, K_cr_paral, K_cr_perp, K_crs_paral, K_crs_perp, use_smallecr
      use mpisetup,        only: rbuff, ibuff, lbuff, cbuff, master, slave, piernik_MPI_Bcast
      use units,           only: clight, me, amu

      implicit none

      integer(kind=4) :: i, j
      real, dimension(:), allocatable ::  p_br_def, q_br_def

      namelist /COSMIC_RAY_SPECTRUM/ cfl_cre, p_lo_init, p_up_init, f_init, q_init, q_big, initial_spectrum, p_min_fix, p_max_fix, &
      &                         cre_eff, cre_active, K_cre_pow, expan_order, e_small, use_cresp, use_cresp_evol,                   &
      &                         e_small_approx_init_cond, p_br_init_lo, e_small_approx_p_lo, e_small_approx_p_up, force_init_NR,   &
      &                         NR_iter_limit, max_p_ratio, synch_active, adiab_active, arr_dim, arr_dim_q, q_br_init,             &
      &                         Gamma_min_fix, Gamma_max_fix, nullify_empty_bins, approx_cutoffs, NR_run_refine_pf, b_max_db,      &
      &                         NR_refine_solution_q, NR_refine_pf_lo, NR_refine_pf_up, smallcree, smallcren, p_br_init_up, p_diff,&
      &                         q_eps, NR_smap_file, cresp_substep, n_substeps_max, allow_unnatural_transfer

      call allocate_spectral_CRspecies_arrays(nspc, ncrb)
      ma1d = [nspc]
      call my_allocate(p_br_def, ma1d)
      call my_allocate(q_br_def, ma1d)

      allocate(bins_primary(nspc))
      do i = 1, nspc
         call init_crel(bins_primary(i))
      enddo
! Default values
      use_cresp         = .true.
      use_cresp_evol    = .true.
      p_min_fix         = 1.5e1
      p_max_fix         = 1.65e4
      p_lo_init(:)      = 1.5e1
      p_up_init(:)      = 7.5e2
      p_br_def(:)       = p_lo_init(:)
      initial_spectrum  = "powl"
      f_init(:)         = 1.0
      q_init(:)         = 4.1
      q_br_def(:)       = q_init(:)
      q_big             = 30.0d0
      p_br_init_lo(:)   = p_br_def(:) ! < in case it was not provided "powl" can be assumed
      p_br_init_up(:)   = p_br_def(:) ! < in case it was not provided "powl" can be assumed
      q_br_init(:)      = q_br_def(:) ! < in case it was not provided "powl" can be assumed
      p_diff(:)         = 10000.0
      cfl_cre(:)        = 0.1
      cre_eff(:)        = 0.01
      K_cre_pow(:)      = 0.
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
      synch_active(:)      = .true.
      adiab_active(:)      = .true.
      cre_active(:)        = 0.0

      if(cr_spectral(icr_H1)) cre_active(findloc(icr_spc, icr_H1)) = 1.0

      if(size(synch_active) > 1) synch_active(2:) = .false. ! non relevant for hadronic species by default

      b_max_db             = 10.  ! default value of B limiter
! NR parameters
      tol_f    = 1.0e-11
      tol_x    = 1.0e-11
      tol_f_1D = 1.0e-14
      tol_x_1D = 1.0e-14
      arr_dim  = 200
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
         ibuff(6)  =  arr_dim
         ibuff(7)  =  arr_dim_q

         ibuff(8)  =  n_substeps_max

         lbuff(1)  =  use_cresp
         lbuff(2)  =  use_cresp_evol
         lbuff(3)  =  allow_source_spectrum_break
         lbuff(4:3+nspc)        = synch_active(:)
         lbuff(4+nspc:3+2*nspc) = adiab_active(:)
         lbuff(4+2*nspc)        = force_init_NR
         lbuff(5+2*nspc)        = NR_run_refine_pf
         lbuff(6+2*nspc)        = NR_refine_solution_q
         lbuff(7+2*nspc)        = NR_refine_pf_lo
         lbuff(8+2*nspc)        = NR_refine_pf_up
         lbuff(9+2*nspc)        = nullify_empty_bins
         lbuff(10+2*nspc)       = approx_cutoffs
         lbuff(11+2*nspc)       = NR_allow_old_smaps

         lbuff(12+2*nspc)       = cresp_substep
         lbuff(13+2*nspc)       = allow_unnatural_transfer

         rbuff(1:nspc)        = cfl_cre(1:nspc)
         rbuff(1+nspc:2*nspc) = cre_eff(1:nspc)
         rbuff(2*nspc+1)      = smallcren
         rbuff(2*nspc+2)      = smallcree
         rbuff(2*nspc+3:3*nspc+2)      = cre_active
         rbuff(3*nspc+3:4*nspc+2)  = p_lo_init(1:nspc)
         rbuff(4*nspc+3:5*nspc+2)  = p_up_init(1:nspc)
         rbuff(5*nspc+3:6*nspc+2)  = f_init(1:nspc)
         rbuff(6*nspc+3:7*nspc+2)  = q_init(1:nspc)
         rbuff(7*nspc+3)      = q_big
         rbuff(7*nspc+4)      = p_min_fix
         rbuff(7*nspc+5)      = p_max_fix
         rbuff(7*nspc+6:8*nspc+5) = K_cre_pow(1:nspc)

         rbuff(8*nspc+6)  = e_small
         rbuff(8*nspc+7)  = max_p_ratio

         rbuff(8*nspc+8)  = tol_f
         rbuff(8*nspc+9) = tol_x
         rbuff(8*nspc+10) = tol_f_1D
         rbuff(8*nspc+11) = tol_x_1D

         rbuff(8*nspc+12) = Gamma_min_fix
         rbuff(8*nspc+13) = Gamma_max_fix

         rbuff(8*nspc+14:9*nspc+13)   = p_br_init_lo(1:nspc)
         rbuff(9*nspc+14:10*nspc+13)   = p_br_init_up(1:nspc)
         rbuff(10*nspc+14:11*nspc+13)  = q_br_init(1:nspc)
         rbuff(11*nspc+14:12*nspc+13) = p_diff(1:nspc)
         rbuff(12*nspc+15) = q_eps
         rbuff(12*nspc+16) = b_max_db

         cbuff(1)  = initial_spectrum
      endif

      call piernik_MPI_Bcast(ibuff)
      call piernik_MPI_Bcast(rbuff)
      call piernik_MPI_Bcast(lbuff)
      call piernik_MPI_Bcast(cbuff, cbuff_len)
      call piernik_MPI_Bcast(NR_smap_file, fnamelen)

      if (slave) then
         expan_order                 = int(ibuff(1),kind=4)

         e_small_approx_p_lo         = int(ibuff(2),kind=1)
         e_small_approx_p_up         = int(ibuff(3),kind=1)
         e_small_approx_init_cond    = int(ibuff(4),kind=1)

         NR_iter_limit               = int(ibuff(5),kind=2)
         arr_dim                     = int(ibuff(6),kind=4)
         arr_dim_q                   = int(ibuff(7),kind=4)

         n_substeps_max              = int(ibuff(8),kind=4)

         use_cresp                   = lbuff(1)
         use_cresp_evol              = lbuff(2)
         allow_source_spectrum_break = lbuff(3)
         synch_active                = lbuff(4:3+nspc)
         adiab_active                = lbuff(4+nspc:3+2*nspc)
         force_init_NR               = lbuff(4+2*nspc)
         NR_run_refine_pf            = lbuff(5+2*nspc)
         NR_refine_solution_q        = lbuff(6+2*nspc)
         NR_refine_pf_lo             = lbuff(7+2*nspc)
         NR_refine_pf_up             = lbuff(8+2*nspc)
         nullify_empty_bins          = lbuff(9+2*nspc)
         approx_cutoffs              = lbuff(10+2*nspc)
         NR_allow_old_smaps          = lbuff(11+2*nspc)

         cresp_substep               = lbuff(12+2*nspc)
         allow_unnatural_transfer    = lbuff(13+2*nspc)

         cfl_cre(1:nspc)      = rbuff(1:nspc)  !TODO check if i'm correct :)
         cre_eff(1:nspc)      = rbuff(1+nspc:2*nspc)
         smallcren            = rbuff(2*nspc+1)
         smallcree            = rbuff(2*nspc+2)
         cre_active           = rbuff(2*nspc+3:3*nspc+2)
         p_lo_init(1:nspc)    = rbuff(3*nspc+3:4*nspc+2)
         p_up_init(1:nspc)    = rbuff(4*nspc+3:5*nspc+2)
         f_init(1:nspc)       = rbuff(5*nspc+3:6*nspc+2)
         q_init(1:nspc)       = rbuff(6*nspc+3:7*nspc+2)
         q_big                = rbuff(7*nspc+3)
         p_min_fix            = rbuff(7*nspc+4)
         p_max_fix            = rbuff(7*nspc+5)
         K_cre_pow(1:nspc)    = rbuff(7*nspc+6:8*nspc+5)

         e_small              = rbuff(8*nspc+6)
         max_p_ratio          = rbuff(8*nspc+7)

         tol_f                = rbuff(8*nspc+8)
         tol_x                = rbuff(8*nspc+9)
         tol_f_1D             = rbuff(8*nspc+10)
         tol_x_1D             = rbuff(8*nspc+11)

         Gamma_min_fix        = rbuff(8*nspc+12)
         Gamma_max_fix        = rbuff(8*nspc+13)

         p_br_init_lo(1:nspc) = rbuff(8*nspc+14:9*nspc+13)
         p_br_init_up(1:nspc) = rbuff(9*nspc+14:10*nspc+13)
         q_br_init(1:nspc)    = rbuff(10*nspc+14:11*nspc+13)
         p_diff(1:nspc)       = rbuff(11*nspc+14:12*nspc+13)
         q_eps                = rbuff(12*nspc+15)
         b_max_db             = rbuff(12*nspc+16)

         initial_spectrum            = trim(cbuff(1))

      endif

      ! since now use only e_small_approx_p, p_init and p_br_init
      e_small_approx_p = [e_small_approx_p_lo, e_small_approx_p_up]
      p_init(1,:)      = p_lo_init(:) ;      p_init(2,:)      = p_up_init(:)
      p_br_init(1,:)   = p_br_init_lo(:) ;   p_br_init(2,:)   = p_br_init_up(:)
      NR_refine_pf     = [NR_refine_pf_lo, NR_refine_pf_up]

! Input parameters check
      if (use_cresp .and. ncrb <= I_ZERO)  then
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
      if (minval(f_init(1:nspc)) < eps) then ! TODO change the eps into epsilon()
         if (initial_spectrum == 'powl' .or. initial_spectrum == 'brpl') then
            write (msg,"(A,A,A)") "[initcrspectrum:init_cresp] Provided power law type spectrum (",initial_spectrum,") with initial amplitude f_init ~ zero. Check your parameters."
            call die(msg)
         endif
      endif

      do j = 1, nspc
        if (initial_spectrum == "brpl" ) then
           if (abs(p_br_init(LO, j) - p_br_def(j)) <= eps) then
              write (msg,"(A,I2,A)") "[initcrspectrum:init_cresp] Parameter for 'brpl' spectrum (component ", j,"): p_br_init_lo has default value (probably unitialized). Assuming p_lo_init value ('powl' spectrum)."
              if (master) call warn(msg)
           else
              !> p_br_init_lo should be equal to one of p_fix values
              i = int(minloc(abs(p_fix - p_br_init(LO, j)), dim=1), kind=4) - I_ONE
              write (msg,"(A,I2,A,E14.7,1A)") "[initcrspectrum:init_cresp] p_br_init_lo was set, but should be equal to one of p_fix (component", j,"). Assuming p_br_init_lo =", p_fix(i),"."
              p_br_init(LO, j) = p_fix(i)
              if (master) call warn(msg)
           endif
           if (abs(q_br_init(j) - q_br_def(j)) <= eps) then
              write (msg,"(A,I2,A)") "[initcrspectrum:init_cresp] Parameter for 'brpl' spectrum (component ", j,"): q_br_init has default value (probably unitialized). Assuming q_init value ('powl' spectrum)."
              if (master) call warn(msg)
           endif
        endif

        if (initial_spectrum == "plpc") then
           if (abs(p_br_init(LO, j) - p_br_def(j)) <= eps .or. abs(p_br_init(HI, j) - p_br_def(j)) <= eps) then
              write (msg,"(A,I2,A)") "[initcrspectrum:init_cresp] Parameters for 'plpc' spectrum (component ", j,"): p_br_init_lo or p_br_init_up has default value (probably unitialized). Check spectrum parameters."
              if (master) call die(msg)
           else
              !> p_br_init_lo should be equal to one of p_fix values
              i = int(minloc(abs(p_fix - p_br_init(LO, j)), dim=1), kind=4) - I_ONE
              write (msg,"(A,I2,A,E14.7,1A)") "[initcrspectrum:init_cresp] p_br_init_lo was set, but should be equal to one of p_fix (component ", j,"). Assuming p_br_init_lo =", p_fix(i),"."
              p_br_init(LO, j) = p_fix(i)
              if (master) call warn(msg)
              !> p_br_init_up should also be equal to one of p_fix values
              i = int(minloc(abs(p_fix - p_br_init(HI, j)), dim=1), kind=4) - I_ONE
              write (msg,"(A,I2,A,E14.7,1A)") "[initcrspectrum:init_cresp] p_br_init_up was set, but should be equal to one of p_fix (component ", j,"). Assuming p_br_init_up =", p_fix(i),"."
              p_br_init(HI, j) = p_fix(i)
              if (master) call warn(msg)
           endif
        endif
        fsynchr(j) =  (4. / 3. ) * cr_sigma_N(icr_spc(j)) / (cr_mass(icr_spc(j)) * amu * clight)

        write (msg, *) "[initcrspectrum:init_cresp] CR ", cr_names(icr_spc(j)), ": 4/3 * sigma_N / ( m * c ) = ", fsynchr(j)
        if (master) call printinfo(msg)

        K_cresp_paral(j, 1:ncrb) = K_cr_paral(icr_spc(j)) * (p_mid_fix(1:ncrb) / p_diff(j))**K_cre_pow(j)
        K_cresp_perp(j,  1:ncrb) = K_cr_perp(icr_spc(j))  * (p_mid_fix(1:ncrb) / p_diff(j))**K_cre_pow(j)

        K_cresp_paral(j, ncrb+1:ncr2b) = K_cresp_paral(j, 1:ncrb)
        K_cresp_perp (j, ncrb+1:ncr2b) = K_cresp_perp (j, 1:ncrb)

        K_crs_paral(ncrn + 1 + (j - 1) * ncr2b : ncrn + ncr2b * j) = K_cresp_paral(j, 1:ncr2b)
        K_crs_perp (ncrn + 1 + (j - 1) * ncr2b : ncrn + ncr2b * j) = K_cresp_perp (j, 1:ncr2b)
      enddo

      call init_cresp_types

#ifdef VERBOSE
      write (msg,"(A,*(E14.5))") "[initcrspectrum:init_cresp] K_cresp_paral = ", K_cresp_paral(1:ncrb) ; if (master) call printinfo(msg)
      write (msg,"(A,*(E14.5))") "[initcrspectrum:init_cresp] K_cresp_perp = ",  K_cresp_perp(1:ncrb)  ; if (master) call printinfo(msg)
#endif /* VERBOSE */

      def_dtadiab(:) = cfl_cre(:) * half * three * logten * w
      def_dtsynch(:) = cfl_cre(:) * half * w

      u_b_max = maxval(fsynchr(:)) * emag(b_max_db, 0., 0.)   !< initializes factor for comparing u_b with u_b_max

      write (msg, "(A,F10.4,A,ES12.5)") "[initcrspectrum:init_cresp] Maximal B_tot =",b_max_db, "mGs, u_b_max = ", u_b_max
      if (master)  call warn(msg)
      write (msg, "(A,ES12.5,A,ES15.8,A,ES15.8)") "[initcrspectrum:init_cresp] Minimal dt_synch(p_max_fix = ",p_max_fix,", u_b_max = ",u_b_max,") = ", minval(def_dtsynch) / (p_max_fix* u_b_max)
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

      if ( (minval(q_init(:) - three) < eps) .and. any(e_small_approx_p == I_ONE)) then
         call warn("[initcrspectrum:init_cresp] Initial parameters: q_init < 3.0 and approximation of outer momenta is on, approximation of outer momenta with hard energy spectrum might not work.")
      endif

      ! TODO (OPTIONAL) make a clenaup of redundant arrays after initialization

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

      cresp%e = zero
      cresp%n = zero

      norm_init_spectrum_n(:,:) = zero
      norm_init_spectrum_e(:,:) = zero

   end subroutine init_cresp_types

!----------------------------------------------------------------------------------------------------

   subroutine init_crel(bin_old_type)

      use constants,      only: zero, I_ZERO, I_ONE
      use diagnostics,    only: my_allocate_with_index
      use initcosmicrays, only: ncrb

      implicit none

      type(bin_old) :: bin_old_type

      if (.not. allocated(bin_old_type%p)) call my_allocate_with_index(bin_old_type%p, ncrb, I_ZERO)
      if (.not. allocated(bin_old_type%f)) call my_allocate_with_index(bin_old_type%f, ncrb, I_ZERO)
      if (.not. allocated(bin_old_type%q)) call my_allocate_with_index(bin_old_type%q, ncrb, I_ONE )
      if (.not. allocated(bin_old_type%n)) call my_allocate_with_index(bin_old_type%n, ncrb, I_ONE )
      if (.not. allocated(bin_old_type%e)) call my_allocate_with_index(bin_old_type%e, ncrb, I_ONE )

      bin_old_type%p = zero
      bin_old_type%q = zero
      bin_old_type%f = zero
      bin_old_type%e = zero
      bin_old_type%n = zero
      bin_old_type%i_cut = I_ZERO

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

      if (dfpq%any_dump) call init_crel(crel)

   end subroutine check_if_dump_fpq

!----------------------------------------------------------------------------------------------------

   subroutine write_cresp_to_restart(file_id)

      use hdf5,            only: HID_T, SIZE_T
      use h5lt,            only: h5ltset_attribute_int_f, h5ltset_attribute_double_f
      use initcosmicrays,  only: ncrb

      implicit none

      integer(HID_T), intent(in) :: file_id
      integer(SIZE_T)            :: bufsize
      integer(kind=4)            :: error
      integer(kind=4), dimension(1) :: lnsnbuf_i
      real,    dimension(1)      :: lnsnbuf_r

      bufsize = 1
      lnsnbuf_i(bufsize) = ncrb
      call h5ltset_attribute_int_f(file_id,    "/", "ncrb",      lnsnbuf_i, bufsize, error)

      lnsnbuf_r(bufsize) = p_min_fix
      call h5ltset_attribute_double_f(file_id, "/", "p_min_fix", lnsnbuf_r, bufsize, error)

      lnsnbuf_r(bufsize) = p_max_fix
      call h5ltset_attribute_double_f(file_id, "/", "p_max_fix", lnsnbuf_r, bufsize, error)

      lnsnbuf_r(bufsize) = q_big
      call h5ltset_attribute_double_f(file_id, "/", "q_big",     lnsnbuf_r, bufsize, error)

      lnsnbuf_r(bufsize) = e_small
      call h5ltset_attribute_double_f(file_id, "/", "e_small",   lnsnbuf_r, bufsize, error)

   end subroutine write_cresp_to_restart
!----------------------------------------------------------------------------------------------------

   subroutine cleanup_cresp_work_arrays

      use diagnostics, only: my_deallocate

      implicit none

      if (allocated(cresp%n))   call my_deallocate(cresp%n)
      if (allocated(cresp%e))   call my_deallocate(cresp%e)
      if (allocated(norm_init_spectrum_n))   deallocate(norm_init_spectrum_n)
      if (allocated(norm_init_spectrum_e))   deallocate(norm_init_spectrum_e)

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
   subroutine allocate_spectral_CRspecies_arrays(nsp, nb)  ! Allocate arrays for spectrally resolved CR species

      use diagnostics,  only: my_allocate, my_allocate_with_index, ma1d, ma2d

      implicit none

      integer, intent(in)  :: nsp, nb

      ma1d = [nsp]
      call my_allocate(p_lo_init, ma1d)
      call my_allocate(p_up_init, ma1d)
      call my_allocate(p_br_init_lo, ma1d)
      call my_allocate(p_br_init_up, ma1d)
      call my_allocate(f_init, ma1d)
      call my_allocate(q_init, ma1d)
      call my_allocate(q_br_init, ma1d)
      call my_allocate(cfl_cre, ma1d)
      call my_allocate(cre_eff, ma1d)
      call my_allocate(K_cre_pow, ma1d)
      call my_allocate(p_diff, ma1d)
      call my_allocate(fsynchr, ma1d)
      call my_allocate(def_dtadiab, ma1d)
      call my_allocate(def_dtsynch, ma1d)
      call my_allocate(total_init_cree, ma1d)
      call my_allocate(cre_active, ma1d)
      call my_allocate_with_index(synch_active, nsp, 1)
      call my_allocate_with_index(adiab_active, nsp, 1)

      ma2d = [nsp, nb * 2]
      call my_allocate(K_cresp_paral, ma2d)
      call my_allocate(K_cresp_perp, ma2d)

      ma2d = [2, nsp]
      call my_allocate(p_init, ma2d)
      call my_allocate(p_br_init, ma2d)

      ma2d = [nsp, nb]
      if (.not. allocated(norm_init_spectrum_n)) call my_allocate(norm_init_spectrum_n, ma2d)
      if (.not. allocated(norm_init_spectrum_e)) call my_allocate(norm_init_spectrum_e, ma2d)

   end subroutine allocate_spectral_CRspecies_arrays

end module initcrspectrum
