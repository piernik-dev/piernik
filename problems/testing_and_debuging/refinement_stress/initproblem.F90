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
!! @file
!! @brief Problem-specific routines
!!
!! \details This module contains routines that are specific to the problem being solved.
!! It includes routines for reading problem parameters, setting up initial conditions,
!! and pointers to user-defined functions that are called during the simulation.
!!
!! This is based on advection test problem but with different refinement criteria.
!! Here the refinement follows a Lissajous curve (nonperiodic by default), stressing the
!! refinement-related parts of the code.
!<

module initproblem

   use constants,  only: dsetnamelen, ndims, LO, HI
   use fluidtypes, only: component_fluid

   implicit none

   private
   public :: read_problem_par, problem_initial_conditions, problem_pointers

   ! namelist parameters
   real, dimension(ndims) :: pulse_size  !< size of the density pulse
   real, dimension(ndims) :: pulse_off   !< offset of the pulse
   real                   :: pulse_amp   !< amplitude of the density pulse compared to the ambient level
   real, dimension(ndims) :: pulse_vel   !< uniform velocity components
   real                   :: pulse_pres  !< pulse pressure
   integer(kind=4)        :: norm_step   !< how often to calculate the L2-norm
   real, dimension(ndims) :: ref_a       !< amplitude to calculate refinement point coordinate
   real, dimension(ndims) :: ref_om      !< frequency to calculate refinement point coordinate
   logical                :: usedust     !< If .false. then do not set velocity for dust
   real, dimension(ndims) :: B_const     !<  constant-B component strength
   real                   :: divB0_amp   !< Amplitude of the non-divergent component of the magnetic field
   real                   :: divBc_amp   !< Amplitude of constant-divergence component of the magnetic field (has artifacts on periodic domains due to nondifferentiability)
   real                   :: divBs_amp   !< Amplitude of sine-wave divergence component of the magnetic field (should behave well on periodic domains)
   real, dimension(ndims) :: divBb_amp   !< Amplitude of [x, y, z]-magnetic field components for  blob of divergence (should behave well on periodic domains)
   real                   :: divBb_r0    !< divergence defect spatial scaling factor
   integer(kind=4), dimension(ndims) :: divB_k   !< wave numbers for creating initial field

   namelist /PROBLEM_CONTROL/  pulse_size, pulse_off, pulse_vel, pulse_amp, pulse_pres, norm_step, ref_a, ref_om, usedust, &
        &                      divB0_amp, divBc_amp, divBs_amp, divBb_amp, divB_k, B_const, divBb_r0

   ! other private data
   real, dimension(ndims, LO:HI) :: pulse_edge
   real :: pulse_low_density
   character(len=dsetnamelen), parameter :: inid_n = "inid"
   class(component_fluid), pointer :: fl  !< will point either to neutral or ionized fluid

contains

!-----------------------------------------------------------------------------

   subroutine problem_pointers

      use user_hooks,  only: finalize_problem, problem_customize_solution
#ifdef HDF5
      use dataio_user, only: user_vars_hdf5
#endif /* HDF5 */

      implicit none

      finalize_problem           => calculate_error_norm
      problem_customize_solution => calculate_error_norm_wrapper
#ifdef HDF5
      user_vars_hdf5             => inid_var_hdf5
#endif /* HDF5 */

   end subroutine problem_pointers

!-----------------------------------------------------------------------------

   subroutine read_problem_par

      use bcast,            only: piernik_MPI_Bcast
      use constants,        only: I_ONE, xdim, zdim, dpi
      use dataio_pub,       only: nh      ! QA_WARN required for diff_nml
      use dataio_pub,       only: warn, die
      use domain,           only: dom
      use fluidindex,       only: flind
      use func,             only: operator(.notequals.)
      use global,           only: smalld, smallei
      use mpisetup,         only: rbuff, ibuff, lbuff, master, slave, proc, have_mpi, LAST
      use unified_ref_crit_list, only: urc_list
#ifdef MAGNETIC
      use constants,        only: GEO_XYZ
#endif /* MAGNETIC */

      implicit none

      ! namelist default parameter values
      pulse_size(:) = 1.0                  !< size of the pulse
      pulse_off(:)  = 0.0                  !< center of the pulse
      pulse_vel(:)  = 0.0                  !< pulse velocity
      pulse_amp     = 2.0                  !< pulse relative amplitude
      pulse_pres    = smallei * 67.  !< pulse pressure, mimic the old hardcoded default
      norm_step     = 5
      ref_a(:)      = 0.9                  !< nearly whole domain
      ref_om(:)     = dpi * [ 1., sqrt(2.), sqrt(3.) ]  !< nonperiodic Lissajous
      usedust       = .false.
      divB0_amp     = 0.                   !< should be safe to set non-0
      divBc_amp     = 0.                   !< unphysical, only for testing
      divBs_amp     = 0.                   !< unphysical, only for testing
      divBb_amp     = 0.                   !< unphysical, only for testing
      divB_k(:)     = [ I_ONE, I_ONE, I_ONE ]
      B_const       = 0.
      divBb_r0      = 1./sqrt(8.)  ! this should match the blob size of Tricco, Price & Bate when domain size is 2

      if (master) then

         if (.not.nh%initialized) call nh%init()
         open(newunit=nh%lun, file=nh%tmp1, status="unknown")
         write(nh%lun,nml=PROBLEM_CONTROL)
         close(nh%lun)
         open(newunit=nh%lun, file=nh%par_file)
         nh%errstr=""
         read(unit=nh%lun, nml=PROBLEM_CONTROL, iostat=nh%ierrh, iomsg=nh%errstr)
         close(nh%lun)
         call nh%namelist_errh(nh%ierrh, "PROBLEM_CONTROL")
         read(nh%cmdl_nml,nml=PROBLEM_CONTROL, iostat=nh%ierrh)
         call nh%namelist_errh(nh%ierrh, "PROBLEM_CONTROL", .true.)
         open(newunit=nh%lun, file=nh%tmp2, status="unknown")
         write(nh%lun,nml=PROBLEM_CONTROL)
         close(nh%lun)
         call nh%compare_namelist()

         rbuff(1)   = pulse_amp
         rbuff(5)   = divB0_amp
         rbuff(6)   = divBc_amp
         rbuff(7)   = divBs_amp
         rbuff(8:10)= divBb_amp
         rbuff(11)  = divBb_r0
         rbuff(20+xdim:20+zdim) = pulse_size(:)
         rbuff(23+xdim:23+zdim) = pulse_vel(:)
         rbuff(26+xdim:26+zdim) = pulse_off(:)
         rbuff(29)              = pulse_pres
         rbuff(30+xdim:30+zdim) = B_const(:)
         rbuff(40+     xdim:40+  zdim) = ref_a
         rbuff(40+zdim+xdim:40+2*zdim) = ref_om

         ibuff(1)   = norm_step
         ibuff(10+xdim:10+zdim) = divB_k(:)

         lbuff(1)   = usedust

      endif

      call piernik_MPI_Bcast(ibuff)
      call piernik_MPI_Bcast(lbuff)
      call piernik_MPI_Bcast(rbuff)

      if (slave) then

         pulse_amp  = rbuff(1)
         divB0_amp  = rbuff(5)
         divBc_amp  = rbuff(6)
         divBs_amp  = rbuff(7)
         divBb_amp  = rbuff(8:10)
         divBb_r0   = rbuff(11)
         pulse_size = rbuff(20+xdim:20+zdim)
         pulse_vel  = rbuff(23+xdim:23+zdim)
         pulse_off  = rbuff(26+xdim:26+zdim)
         pulse_pres = rbuff(29)
         B_const    = rbuff(30+xdim:30+zdim)
         ref_a      = rbuff(40+     xdim:40+  zdim)
         ref_om     = rbuff(40+zdim+xdim:40+2*zdim)

         norm_step  = int(ibuff(1), kind=4)
         divB_k(:)  = ibuff(10+xdim:10+zdim)

         usedust    = lbuff(1)

      endif

      if (any(pulse_size <= 0. .and. dom%has_dir)) call die("[initproblem:read_problem_par] Pulse size has to be positive")

      if (pulse_amp <= 0.) then
         if (have_mpi) then
            pulse_amp = 1. + proc/real(LAST)
            pulse_size = 1.
            if (master) call warn("[initproblem:read_problem_par] The analytical solution will not be correctly advected (not implemented yet)")
         else
            pulse_amp = 2.
         endif
      endif

      where (dom%has_dir(:))
         pulse_edge(:, LO) = pulse_off(:) - pulse_size/2.
         pulse_edge(:, HI) = pulse_off(:) + pulse_size/2.
      elsewhere
         pulse_edge(:, LO) = -huge(1.)
         pulse_edge(:, HI) =  huge(1.)
      endwhere

      if (associated(flind%neu)) then
         fl => flind%neu
      else if (associated(flind%ion)) then
         fl => flind%ion
      else
         call die("[initproblem] current implementation requires either ionized or neutral fluid to be defined")
      endif

      !BEWARE: hardcoded magic numbers
      pulse_low_density = smalld * 1e5

      if (norm_step <= 0) norm_step = huge(I_ONE)

      ! Create the initial density arrays (it is called before reading restart file, so there is no need to associate user_reg_var_restart)
      call register_user_var

      ! Register user-defined refinement criteria provided by this module
      call urc_list%add_user_urc(set_refinement, .true.)

      if (any([divB0_amp, divBc_amp, divBs_amp, divBb_amp] .notequals. 0.)) then
#ifdef MAGNETIC
         if (dom%geometry_type /= GEO_XYZ) then
            call warn("[initproblem:read_problem_par] Only cartesian formulas for magnetic field is implemented. Forcing all amplitudes to 0.")
            divB0_amp = 0.
            divBc_amp = 0.
            divBs_amp = 0.
            divBb_amp = 0.
         endif
#else
         if (master) call warn("[initproblem:read_problem_par] Ignoring magnetic field amplitudes")
#endif /* !MAGNETIC */
      endif

   end subroutine read_problem_par

!-----------------------------------------------------------------------------

   subroutine problem_initial_conditions

      use cg_list,          only: cg_list_element
      use cg_leaves,        only: leaves
      use constants,        only: xdim, ydim, zdim, GEO_XYZ, GEO_RPZ
      use dataio_pub,       only: die
      use domain,           only: dom
      use fluidindex,       only: flind
      use func,             only: operator(.notequals.), emag
      use global,           only: smallei, t
      use grid_cont,        only: grid_container
      use named_array_list, only: qna
      use non_inertial,     only: get_omega
#ifdef MAGNETIC
#ifdef IONIZED
      use constants,        only: half
#endif /* IONIZED */
      use constants,        only: ndims, I_ONE, I_TWO, I_THREE, dpi
      use dataio_pub,       only: warn
      use div_B,            only: print_divB_norm
      use global,           only: cc_mag
      use mpisetup,         only: master
#endif /* MAGNETIC */

      implicit none

      type(cg_list_element), pointer :: cgl
      type(grid_container), pointer :: cg

      integer :: i, j, k
      real    :: om
#ifdef MAGNETIC
      real, dimension(ndims) :: kk !< wavenumbers to fit one sine wave inside domain
      real :: sx, sy, sz, cx, cy, cz !< precomputed values of sine and cosine at cell centers
      real :: sfx, sfy, sfz, cfx, cfy, cfz !< precomputed values of sine and cosine at cell left faces
      integer :: right_face
      real :: r02, rr02

      kk = 0.
      where (dom%D_ > 0) kk = divB_k * dpi / dom%L_
      right_face = 1
      if (cc_mag) right_face = 0
      r02 = divBb_r0**2
#endif /* MAGNETIC */

      om = get_omega()

      call analytic_solution(t)

      cgl => leaves%first
      do while (associated(cgl))
         cg => cgl%cg

#ifdef MAGNETIC
         if (dom%geometry_type /= GEO_XYZ) then
            call cg%set_constant_b_field([0., 0., 0.])
            if (master) call warn("[initproblem:problem_initial_conditions] noncartesian coordinates not implemented for magnetic field.")
         else
            do k = cg%ks, cg%ke + right_face * dom%D_z
               sz = sin(kk(zdim) * cg%z(k))
               cz = cos(kk(zdim) * cg%z(k))
               sfz = sin(kk(zdim) * (cg%z(k) - cg%dz/2.))
               cfz = cos(kk(zdim) * (cg%z(k) - cg%dz/2.))
               do j = cg%js, cg%je + right_face * dom%D_y
                  sy = sin(kk(ydim) * cg%y(j))
                  cy = cos(kk(ydim) * cg%y(j))
                  sfy = sin(kk(ydim) * (cg%y(j) - cg%dy/2.))
                  cfy = cos(kk(ydim) * (cg%y(j) - cg%dy/2.))
                  do i = cg%is, cg%ie + right_face * dom%D_x
                     sx = sin(kk(xdim) * cg%x(i))
                     cx = cos(kk(xdim) * cg%x(i))
                     sfx = sin(kk(xdim) * (cg%x(i) - cg%dx/2.))
                     cfx = cos(kk(xdim) * (cg%x(i) - cg%dx/2.))

                     cg%b(:, i, j, k) = B_const(:) + divBc_amp * [cg%x(i), cg%y(j), cg%z(k)] ! slight offset between cell- and face-centered is unimportant here

                     ! div B pulse, as described in Tricco, Price & Bate, https://arxiv.org/abs/1607.02394
                     rr02 = sum(([cg%x(i), cg%y(j), cg%z(k)])**2, mask=dom%has_dir)/r02
                     if (rr02 < 1.) then
                        where (dom%has_dir)  cg%b(:, i, j, k) = cg%b(:, i, j, k) + divBb_amp * (rr02**4 - 2 *rr02**2 + 1)
                     endif

                     select case (dom%eff_dim)
                        case (I_ONE) ! can't do anything fancy, just set up something non-zero
                           cg%b(:, i, j, k) = cg%b(:, i, j, k) + divB0_amp
                           if (cc_mag) then
                              if (dom%D_x == 1) then
                                 cg%b(:, i, j, k) = cg%b(:, i, j, k) + divBs_amp * [ kk(xdim)*cx, 1., 1. ]
                              else if (dom%D_y == 1) then
                                 cg%b(:, i, j, k) = cg%b(:, i, j, k) + divBs_amp * [ 1., kk(ydim)*cy, 1. ]
                              else
                                 cg%b(:, i, j, k) = cg%b(:, i, j, k) + divBs_amp * [ 1., 1., kk(zdim)*cz ]
                              endif
                           else
                              if (dom%D_x == 1) then
                                 cg%b(:, i, j, k) = cg%b(:, i, j, k) + divBs_amp * [ kk(xdim)*cfx, 1., 1. ]
                              else if (dom%D_y == 1) then
                                 cg%b(:, i, j, k) = cg%b(:, i, j, k) + divBs_amp * [ 1., kk(ydim)*cfy, 1. ]
                              else
                                 cg%b(:, i, j, k) = cg%b(:, i, j, k) + divBs_amp * [ 1., 1., kk(zdim)*cfz ]
                              endif
                           endif
                        case (I_TWO)
                           ! [sin(x)*sin(y), cos(x)*cos(y), 0] should produce divB == 0. for XY case (curl([0, 0, -sin(x)*cos(y)]))
                           ! The div(B) is really close to numerical noise around 0 only in the case of exactly the same resolution per sine wave in all directions.
                           ! If the resolutions of sine waves don't match, then numerical estimates of mixed derivatives of the vector potential don't cancel out and only high-order estimates of div(b) are close to 0.
                           if (cc_mag) then
                              if (dom%D_z == 0) then
                                 cg%b(:, i, j, k) = cg%b(:, i, j, k) + &
                                      divB0_amp * [ kk(ydim)*sx*sy, kk(xdim)*cx*cy, 1. ] + &
                                      divBs_amp * [ kk(xdim)*cx*sy, kk(ydim)*sx*cy, 1. ]
                              else if (dom%D_y == 0) then
                                 cg%b(:, i, j, k) = cg%b(:, i, j, k) + &
                                      divB0_amp * [ kk(zdim)*sx*sz, 1., kk(xdim)*cx*cz ] + &
                                      divBs_amp * [ kk(xdim)*cx*sz, 1., kk(zdim)*sx*cz ]
                              else
                                 cg%b(:, i, j, k) = cg%b(:, i, j, k) + &
                                      divB0_amp * [ 1., kk(zdim)*sy*sz, kk(ydim)*cy*cz ] + &
                                      divBs_amp * [ 1., kk(ydim)*cy*sz, kk(zdim)*sy*cz ]
                              endif
                           else
                              if (dom%D_z == 0) then
                                 cg%b(:, i, j, k) = cg%b(:, i, j, k) + &
                                      divB0_amp * [ kk(ydim)*sfx*sy, kk(xdim)*cx*cfy, 1. ] + &
                                      divBs_amp * [ kk(xdim)*cfx*sy, kk(ydim)*sx*cfy, 1. ]
                              else if (dom%D_y == 0) then
                                 cg%b(:, i, j, k) = cg%b(:, i, j, k) + &
                                      divB0_amp * [ kk(zdim)*sfx*sz, 1., kk(xdim)*cx*cfz ] + &
                                      divBs_amp * [ kk(xdim)*cfx*sz, 1., kk(zdim)*sx*cfz ]
                              else
                                 cg%b(:, i, j, k) = cg%b(:, i, j, k) + &
                                      divB0_amp * [ 1., kk(zdim)*sfy*sz, kk(ydim)*cy*cfz ] + &
                                      divBs_amp * [ 1., kk(ydim)*cfy*sz, kk(zdim)*sy*cfz ]
                              endif
                           endif
                        case (I_THREE)
                           ! curl([sin(x)*sin(y)*sin(z), sin(x)*sin(y)*sin(z), sin(x)*sin(y)*sin(z)]) should produce div(B) == 0, but see the notes for 2D case.
                           ! setting up a div(B)-free field in flattened domain requires careful choice of kk(:)
                           if (cc_mag) then
                              cg%b(:, i, j, k) = cg%b(:, i, j, k) + divB0_amp * [ &
                                   kk(ydim)*cx*sy*cz - kk(zdim)*cx*cy*sz, &
                                   kk(zdim)*cx*cy*sz - kk(xdim)*sx*cy*cz, &
                                   kk(xdim)*sx*cy*cz - kk(ydim)*cx*sy*cz ] + &
                                   divBs_amp * [ kk(xdim)*cx*sy*sz, kk(ydim)*sx*cy*sz, kk(zdim)*sx*sy*cz ]
                           else
                              cg%b(:, i, j, k) = cg%b(:, i, j, k) + divB0_amp * [ &
                                   kk(ydim)*cfx*sy*cz - kk(zdim)*cfx*cy*sz, &
                                   kk(zdim)*cx*cfy*sz - kk(xdim)*sx*cfy*cz, &
                                   kk(xdim)*sx*cy*cfz - kk(ydim)*cx*sy*cfz ] + &
                                   divBs_amp * [ kk(xdim)*cfx*sy*sz, kk(ydim)*sx*cfy*sz, kk(zdim)*sx*sy*cfz ]
                           endif
                        case default
                           call die("[initproblem:problem_initial_conditions] unsupported dimensionality")
                     end select
                  enddo
               enddo
            enddo
         endif
#else /* !MAGNETIC */
         call cg%set_constant_b_field([0., 0., 0.])
#endif /* !MAGNETIC */

         cg%u(fl%idn, cg%is:cg%ie, cg%js:cg%je, cg%ks:cg%ke) = cg%q(qna%ind(inid_n))%arr(cg%is:cg%ie, cg%js:cg%je, cg%ks:cg%ke)

         select case (dom%geometry_type)
            case (GEO_XYZ)
               if (om .notequals. 0.) then
                  ! Include rotation for pulse_vel = 0., 0., 0. case
                  do i = cg%is, cg%ie
                     do j = cg%js, cg%je
                        cg%u(fl%imx, i, j, :) = - om*cg%y(j) * cg%u(fl%idn, i, j, :)
                        cg%u(fl%imy, i, j, :) = + om*cg%x(i) * cg%u(fl%idn, i, j, :)
                     enddo
                  enddo
                  cg%u(fl%imz, :, :, :) = pulse_vel(zdim) * cg%u(fl%idn, :, :, :)
               else
                  ! Make uniform, completely boring flow
                  cg%u(fl%imx, :, :, :) = pulse_vel(xdim) * cg%u(fl%idn, :, :, :)
                  cg%u(fl%imy, :, :, :) = pulse_vel(ydim) * cg%u(fl%idn, :, :, :)
                  cg%u(fl%imz, :, :, :) = pulse_vel(zdim) * cg%u(fl%idn, :, :, :)
               endif
            case (GEO_RPZ)
               do k = cg%ks, cg%ke
                  do j = cg%js, cg%je
                     do i = cg%is, cg%ie
                        cg%u(fl%imx, i, j, k) = ( pulse_vel(xdim)*cos(cg%y(j)) + pulse_vel(ydim)*sin(cg%y(j))) * cg%u(fl%idn, i, j, k)
                        cg%u(fl%imy, i, j, k) = (-pulse_vel(xdim)*sin(cg%y(j)) + pulse_vel(ydim)*cos(cg%y(j))) * cg%u(fl%idn, i, j, k)
                     enddo
                  enddo
               enddo
               cg%u(fl%imz, :, :, :) = pulse_vel(zdim) * cg%u(fl%idn, :, :, :)
            case default
               call die("[initproblem:problem_initial_conditions] only cartesian and cylindrical geometries are supported")
         end select

         ! Set up the internal energy
         cg%u(fl%ien,:,:,:) = max(smallei, pulse_pres / fl%gam_1 + 0.5 * sum(cg%u(fl%imx:fl%imz,:,:,:)**2,1) / cg%u(fl%idn,:,:,:))

#if defined MAGNETIC && defined IONIZED
         if (cc_mag) then
            cg%u(fl%ien,:,:,:) = cg%u(fl%ien,:,:,:) + emag(cg%b(xdim,:,:,:), cg%b(ydim,:,:,:), cg%b(zdim,:,:,:))
         else
            cg%u(fl%ien, cg%is:cg%ie, cg%js:cg%je, cg%ks:cg%ke) = cg%u(fl%ien, cg%is:cg%ie, cg%js:cg%je, cg%ks:cg%ke) + &
                 emag(half*(cg%b(xdim, cg%is:cg%ie, cg%js:cg%je, cg%ks:cg%ke) + cg%b(xdim, cg%is+dom%D_x:cg%ie+dom%D_x, cg%js        :cg%je,         cg%ks        :cg%ke        )), &
                 &    half*(cg%b(ydim, cg%is:cg%ie, cg%js:cg%je, cg%ks:cg%ke) + cg%b(ydim, cg%is        :cg%ie,         cg%js+dom%D_y:cg%je+dom%D_y, cg%ks        :cg%ke        )), &
                 &    half*(cg%b(zdim, cg%is:cg%ie, cg%js:cg%je, cg%ks:cg%ke) + cg%b(zdim, cg%is        :cg%ie,         cg%js        :cg%je,         cg%ks+dom%D_z:cg%ke+dom%D_z)))
         endif
#endif  /* MAGNETIC && IONIZED */

         if (associated(flind%dst)) then
            cg%u(flind%dst%idn, :, :, :) = cg%u(fl%idn, :, :, :)
            if (usedust) then
               cg%u(flind%dst%imx, :, :, :) = cg%u(fl%imx, :, :, :)
               cg%u(flind%dst%imy, :, :, :) = cg%u(fl%imy, :, :, :)
               cg%u(flind%dst%imz, :, :, :) = cg%u(fl%imz, :, :, :)
            else
               cg%u(flind%dst%imx:flind%dst%imz, :, :, :) = 0.
            endif
         endif

         cgl => cgl%nxt
      enddo

#ifdef MAGNETIC
      call print_divB_norm
#endif /* MAGNETIC */

   end subroutine problem_initial_conditions

!-----------------------------------------------------------------------------

   subroutine inid_var_hdf5(var, tab, ierrh, cg)

      use global,           only: t
      use grid_cont,        only: grid_container
      use named_array_list, only: qna

      implicit none

      character(len=*),              intent(in)    :: var
      real, dimension(:,:,:),        intent(inout) :: tab
      integer,                       intent(inout) :: ierrh
      type(grid_container), pointer, intent(in)    :: cg

      call analytic_solution(t) ! cannot handle this automagically because here we modify it

      ierrh = 0
      if (qna%exists(var)) then
         tab(:,:,:) = real(cg%q(qna%ind(var))%span(cg%ijkse), kind(tab))
      else
         ierrh = -1
      endif

   end subroutine inid_var_hdf5

!-----------------------------------------------------------------------------

   subroutine register_user_var

      use cg_list_global, only: all_cg
      use constants,      only: AT_NO_B

      implicit none

      call all_cg%reg_var(inid_n, restart_mode = AT_NO_B)

   end subroutine register_user_var

!-----------------------------------------------------------------------------

   subroutine calculate_error_norm_wrapper(forward)

      implicit none

      logical, intent(in) :: forward

      call calculate_error_norm
      return
      if (.false. .and. forward) pulse_size = 0.0 ! suppress compiler warnings on unused arguments

   end subroutine calculate_error_norm_wrapper

!-----------------------------------------------------------------------------

   subroutine calculate_error_norm

      use allreduce,        only: piernik_MPI_Allreduce
      use cg_list,          only: cg_list_element
      use cg_leaves,        only: leaves
      use constants,        only: PIERNIK_FINISHED, pSUM, pMIN, pMAX, idlen
      use dataio_pub,       only: code_progress, halfstep, msg, printinfo, warn
      use fluidindex,       only: flind
      use func,             only: operator(.notequals.)
      use global,           only: t, nstep
      use grid_cont,        only: grid_container
      use mpisetup,         only: master
      use named_array_list, only: qna
#ifdef MAGNETIC
      use div_B,            only: print_divB_norm
#endif /* MAGNETIC */

      implicit none

      enum, bind(C)
         enumerator :: N_D, N_2
      end enum
      enum, bind(C)
         enumerator :: GAS, DST_
      end enum
      real, dimension(N_D:N_2, GAS:DST_) :: norm
      real, dimension(GAS:DST_)          :: neg_err, pos_err
      type(cg_list_element),  pointer   :: cgl
      type(grid_container),   pointer   :: cg
      real, dimension(:,:,:), pointer   :: inid
      integer                           :: i, j
      character(len=idlen)              :: descr

      if (code_progress < PIERNIK_FINISHED .and. (mod(nstep, norm_step) /= 0 .or. halfstep)) return

      norm = 0.
      neg_err = huge(1.0)
      pos_err = -neg_err

      call analytic_solution(t)

      cgl => leaves%first
      do while (associated(cgl))
         cg => cgl%cg

         inid => cg%q(qna%ind(inid_n))%arr
         if (.not. associated(inid))then
            if (master) call warn("[initproblem:calculate_error_norm] Cannot compare results with the initial conditions.")
            return
         endif

         cg%wa(cg%is:cg%ie, cg%js:cg%je, cg%ks:cg%ke) = inid(cg%is:cg%ie, cg%js:cg%je, cg%ks:cg%ke) - cg%u(fl%idn, cg%is:cg%ie, cg%js:cg%je, cg%ks:cg%ke)
         norm(N_D, GAS) = norm(N_D, GAS) + sum(cg%wa(cg%is:cg%ie, cg%js:cg%je, cg%ks:cg%ke)**2, mask=cg%leafmap)
         norm(N_2, GAS) = norm(N_2, GAS) + sum(inid( cg%is:cg%ie, cg%js:cg%je, cg%ks:cg%ke)**2, mask=cg%leafmap)
         neg_err(GAS) = min(neg_err(GAS), minval(cg%wa(cg%is:cg%ie, cg%js:cg%je, cg%ks:cg%ke), mask=cg%leafmap))
         pos_err(GAS) = max(pos_err(GAS), maxval(cg%wa(cg%is:cg%ie, cg%js:cg%je, cg%ks:cg%ke), mask=cg%leafmap))

         if (associated(flind%dst)) then
            cg%wa(cg%is:cg%ie, cg%js:cg%je, cg%ks:cg%ke) = inid(cg%is:cg%ie, cg%js:cg%je, cg%ks:cg%ke) - cg%u(flind%dst%idn, cg%is:cg%ie, cg%js:cg%je, cg%ks:cg%ke)
            norm(N_D, DST_) = norm(N_D, DST_) + sum(cg%wa(cg%is:cg%ie, cg%js:cg%je, cg%ks:cg%ke)**2, mask=cg%leafmap)
            norm(N_2, DST_) = norm(N_2, DST_) + sum(inid( cg%is:cg%ie, cg%js:cg%je, cg%ks:cg%ke)**2, mask=cg%leafmap)
            neg_err(DST_) = min(neg_err(DST_), minval(cg%wa(cg%is:cg%ie, cg%js:cg%je, cg%ks:cg%ke), mask=cg%leafmap))
            pos_err(DST_) = max(pos_err(DST_), maxval(cg%wa(cg%is:cg%ie, cg%js:cg%je, cg%ks:cg%ke), mask=cg%leafmap))
         endif

         cgl => cgl%nxt
      enddo

      do i = GAS, DST_
         do j = N_D, N_2
            call piernik_MPI_Allreduce(norm(j, i), pSUM)
         enddo
         call piernik_MPI_Allreduce(neg_err(i), pMIN)
         call piernik_MPI_Allreduce(pos_err(i), pMAX)
      enddo

      if (master) then
         do i = GAS, DST_
            if (i == DST_ .and. .not. usedust) exit
            select case (i)
               case (GAS)
                  descr = "GAS"
               case (DST_)
                  descr = "DST"
            end select
            if (norm(N_2, i) .notequals. 0.) then
               write(msg,'(3a,f12.6,a,2f15.6)')"[initproblem:calculate_error_norm] L2 error norm (",descr,") = ", sqrt(norm(N_D, i)/norm(N_2, i)), &
                    ", min and max error = ", neg_err(i), pos_err(i)
               call printinfo(msg)
            endif
         enddo
      endif

#ifdef MAGNETIC
      call print_divB_norm
#endif /* MAGNETIC */

   end subroutine calculate_error_norm

   !>
   !! \brief Put analytic solution in the inid arrays
   !!
   !! \details Density is shaped as an uniform box and translated according to initial velocity and given time
   !<

   subroutine analytic_solution(t)

      use cg_list,          only: cg_list_element
      use cg_leaves,        only: leaves
      use constants,        only: xdim, ydim, zdim, ndims, GEO_XYZ, GEO_RPZ
      use dataio_pub,       only: warn, die
      use domain,           only: dom
      use func,             only: operator(.notequals.)
      use grid_cont,        only: grid_container
      use non_inertial,     only: get_omega
      use mpisetup,         only: master
      use named_array_list, only: qna

      implicit none

      real, intent(in)                :: t !< time of the solution

      real                            :: dini
      integer                         :: i, j, k, d
      type(cg_list_element),  pointer :: cgl
      type(grid_container),   pointer :: cg
      real, dimension(:,:,:), pointer :: inid !< analytic solution
      real, dimension(ndims)          :: pos
      real, dimension(xdim:ydim)      :: paux
      real                            :: om, cot, sot

      pos = 0. ! suppres compiler warning
      om = get_omega()
      cot = cos(om * t)
      sot = sin(om * t)

      cgl => leaves%first
      do while (associated(cgl))
         cg => cgl%cg

         inid => cg%q(qna%ind(inid_n))%arr
         if (.not. associated(inid))then
            if (master) call warn("[initproblem:analytic_solution] Cannot store the initial conditions.")
            return
         endif

         do k = cg%ks, cg%ke
            do j = cg%js, cg%je
               do i = cg%is, cg%ie

                  select case (dom%geometry_type)
                     case (GEO_XYZ)
                        pos = [cg%x(i), cg%y(j), cg%z(k)]
                        if (om .notequals. 0.) then
                           paux = [ pos(xdim) * cot - pos(ydim) * sot, pos(xdim) * sot + pos(ydim) * cot ]
                           pos(xdim:ydim) = paux
                        endif
                        pos = pos - t * pulse_vel(:)
                     case (GEO_RPZ)
                        pos = [cg%x(i)*cos(cg%y(j)), cg%x(i)*sin(cg%y(j)), cg%z(k)] - t * pulse_vel(:)
                     case default
                        call die("[initproblem:analytic_solution] only cartesian and cylindrical geometries are supported")
                  end select
                  do d = xdim, zdim
                     if ((dom%geometry_type == GEO_XYZ .or. (dom%geometry_type == GEO_RPZ .and. d == zdim)) .and. dom%periodic(d)) then
                        if (pos(d) < dom%edge(d, LO)) then
                           pos(d) = pos(d) + dom%L_(d) * ceiling((dom%edge(d, LO) - pos(d))/dom%L_(d))
                        else if (pos(d) > dom%edge(d, HI)) then
                           pos(d) = pos(d) + dom%L_(d) * floor  ((dom%edge(d, HI) - pos(d))/dom%L_(d))
                        endif
                     endif
                  enddo

                  dini = 0.
                  if (all(pos(:) > pulse_edge(:, LO) - cg%dl(:)/2.).and. all(pos(:) < pulse_edge(:, HI) + cg%dl(:)/2.)) then
                     dini = pulse_low_density * (pulse_amp - 1.)
                     do d = xdim, zdim
                        if (dom%has_dir(d)) then
                           if (abs(pos(d) - pulse_edge(d, LO)) < cg%dl(d)/2.) dini = dini * (0.5 + (pos(d) - pulse_edge(d, LO))/cg%dl(d))
                           if (abs(pos(d) - pulse_edge(d, HI)) < cg%dl(d)/2.) dini = dini * (0.5 - (pos(d) - pulse_edge(d, HI))/cg%dl(d))
                        endif
                     enddo
                  endif

                  inid(i, j, k) = dini + pulse_low_density

               enddo
            enddo
         enddo

         cgl => cgl%nxt
      enddo

   end subroutine analytic_solution

!>
!! \brief set up refinement that varies with time
!! \details This is an exmaple of a user-defined refinement criterion.
!<

   subroutine set_refinement(this, cg)

      use cg_level_finest,       only: finest
      use constants,             only: xdim, zdim, INVALID
      use domain,                only: dom
      use global,                only: t, dt
      use grid_cont,             only: grid_container
      use refinement,            only: n_updAMR
      use unified_ref_crit_user, only: urc_user

      implicit none

      class(urc_user),               intent(inout) :: this  !< an object invoking the type-bound procedure
      type(grid_container), pointer, intent(inout) :: cg    !< current grid piece

      real, dimension(xdim:zdim) :: cp, cf, c
      real :: cos_ang
      real, parameter :: ref_thr = 0.9
      integer :: i, j, k

      ! OPT: beware executed multiple times
      where (dom%has_dir)
         cp = dom%C_ + 0.5 * dom%L_ * ref_a * sin( ref_om * t )  ! position of the point
         cf = dom%C_ + 0.5 * dom%L_ * ref_a * sin( ref_om * (t + 2 * dt * n_updAMR))  ! future position of the point, assumes no drastic changes of dt when n_updAMR > 1
      elsewhere
         cp = dom%C_
         cf = dom%C_
      endwhere

      do k = cg%ks, cg%ke
         do j = cg%js, cg%je
            do i = cg%is, cg%ie

               ! Mark cells for refinement based on the angle between directions to cp and cf.
               ! The angle is calculated using dot product and normalized vectors.
               ! The marked area will follow a line segment.
               c = [ cg%x(i), cg%y(j), cg%z(k) ]
               cos_ang = dot_product(c-cp, cf-c)/sqrt(dot_product(c-cp, c-cp) * dot_product(cf-c, cf-c))
               if (cos_ang > ref_thr) call cg%flag%set(i, j, k)
               if (this%iplot /= INVALID) cg%q(this%iplot)%arr(i, j, k) = cos_ang

               ! Failover for dt = 0: mark at least a point.
               if (all(abs(c-cp) <= cg%dl) .or. all(abs(c-cf) <= cg%dl)) then
                  call cg%flag%set(i, j, k)
                  if (this%iplot /= INVALID) cg%q(this%iplot)%arr(i, j, k) = 1.
               endif

               ! Mark additional pilot cells when the point is moving fast.
               ! 2*cf - cp is the predicted coordinate at t + 4 * dt * n_updAMR (twice the future position).
               ! It is often needed to predict the future position of the required refinement to avoid the situation when the region of interest drifts to coarser cells.
               if (all(abs(c-(2*cf-cp)) <= cg%dl) .and. cg%l%id < finest%level%l%id) then
                  call cg%flag%set(i, j, k)
                  if (this%iplot /= INVALID) cg%q(this%iplot)%arr(i, j, k) = 0.5
               endif
            enddo
         enddo
      enddo

   end subroutine set_refinement

end module initproblem
