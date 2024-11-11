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
   integer(kind=4)        :: nflip       !< how often to call refine/derefine routine
   real                   :: flipratio   !< percentage of blocks on each level to be refined on flip
   real                   :: ref_thr     !< refinement threshold
   logical                :: usedust     !< If .false. then do not set velocity for dust
   real, dimension(ndims) :: B_const     !<  constant-B component strength
   real                   :: divB0_amp   !< Amplitude of the non-divergent component of the magnetic field
   real                   :: divBc_amp   !< Amplitude of constant-divergence component of the magnetic field (has artifacts on periodic domains due to nondifferentiability)
   real                   :: divBs_amp   !< Amplitude of sine-wave divergence component of the magnetic field (should behave well on periodic domains)
   real, dimension(ndims) :: divBb_amp   !< Amplitude of [x, y, z]-magnetic field components for  blob of divergence (should behave well on periodic domains)
   real                   :: divBb_r0    !< divergence defect spatial scaling factor
   integer(kind=4), dimension(ndims) :: divB_k   !< wave numbers for creating initial field

   namelist /PROBLEM_CONTROL/  pulse_size, pulse_off, pulse_vel, pulse_amp, pulse_pres, norm_step, nflip, flipratio, ref_thr, usedust, &
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
      use constants,        only: I_ONE, xdim, zdim
      use dataio_pub,       only: warn, die, nh
      use domain,           only: dom
      use fluidindex,       only: flind, iarr_all_dn
      use func,             only: operator(.notequals.)
      use global,           only: smalld, smallei
      use mpisetup,         only: rbuff, ibuff, lbuff, master, slave, proc, have_mpi, LAST
      use named_array_list, only: wna
      use refinement,       only: set_n_updAMR, n_updAMR
      use unified_ref_crit_list, only: urc_list
      use user_hooks,       only: problem_refine_derefine
#ifdef MAGNETIC
      use constants,        only: GEO_XYZ
#endif /* MAGNETIC */

      implicit none

      integer(kind=4) :: id

      ! namelist default parameter values
      pulse_size(:) = 1.0                  !< size of the pulse
      pulse_off(:)  = 0.0                  !< center of the pulse
      pulse_vel(:)  = 0.0                  !< pulse velocity
      pulse_amp     = 2.0                  !< pulse relative amplitude
      pulse_pres    = smallei * 67.  !< pulse pressure, mimic the old hardcoded default
      norm_step     = 5
      nflip         = 0
      ref_thr       = 0.1
      flipratio     = 1.
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
         rbuff(2)   = ref_thr
         rbuff(4)   = flipratio
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

         ibuff(1)   = norm_step
         ibuff(2)   = nflip
         ibuff(10+xdim:10+zdim) = divB_k(:)

         lbuff(1)   = usedust

      endif

      call piernik_MPI_Bcast(ibuff)
      call piernik_MPI_Bcast(lbuff)
      call piernik_MPI_Bcast(rbuff)

      if (slave) then

         pulse_amp  = rbuff(1)
         ref_thr    = rbuff(2)
         flipratio  = rbuff(4)
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

         norm_step  = int(ibuff(1), kind=4)
         nflip      = ibuff(2)
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

      if (nflip > 0) then
         problem_refine_derefine => flip_flop
         if (n_updAMR /= nflip .and. master) call warn("[initproblem:read_problem_par] Forcing n_updAMR == nflip")
         call set_n_updAMR(nflip)
      else
         ! Automatic refinement criteria
         do id = lbound(iarr_all_dn, dim=1, kind=4), ubound(iarr_all_dn, dim=1, kind=4)
            call urc_list%add_user_urcv(wna%fi, id, ref_thr*pulse_amp, 0., "grad", .true.)
         enddo
      endif

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

         cg%u(fl%idn, :, :, :) = cg%q(qna%ind(inid_n))%arr(:, :, :)

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
         cg%u(fl%ien,RNG) = max(smallei, pulse_pres / fl%gam_1 + 0.5 * sum(cg%u(fl%imx:fl%imz,RNG)**2,1) / cg%u(fl%idn,RNG))

#if defined MAGNETIC && defined IONIZED
         if (cc_mag) then
            cg%u(fl%ien,RNG) = cg%u(fl%ien,RNG) + emag(cg%b(xdim,RNG), cg%b(ydim,RNG), cg%b(zdim,RNG))
         else
            cg%u(fl%ien, RNG) = cg%u(fl%ien, RNG) + &
                 emag(half*(cg%b(xdim, RNG) + cg%b(xdim, cg%is+dom%D_x:cg%ie+dom%D_x, cg%js        :cg%je,         cg%ks        :cg%ke        )), &
                 &    half*(cg%b(ydim, RNG) + cg%b(ydim, cg%is        :cg%ie,         cg%js+dom%D_y:cg%je+dom%D_y, cg%ks        :cg%ke        )), &
                 &    half*(cg%b(zdim, RNG) + cg%b(zdim, cg%is        :cg%ie,         cg%js        :cg%je,         cg%ks+dom%D_z:cg%ke+dom%D_z)))
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

      use global, only: nstep

      implicit none

      logical, intent(in) :: forward

      if (forward .and. mod(nstep, norm_step) == 0) call calculate_error_norm

   end subroutine calculate_error_norm_wrapper

!-----------------------------------------------------------------------------

   subroutine calculate_error_norm

      use allreduce,        only: piernik_MPI_Allreduce
      use cg_list,          only: cg_list_element
      use cg_leaves,        only: leaves
      use constants,        only: pSUM, pMIN, pMAX, idlen
      use dataio_pub,       only: msg, printinfo, warn
      use fluidindex,       only: flind
      use func,             only: operator(.notequals.)
      use global,           only: t
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

         cg%wa(RNG) = inid(RNG) - cg%u(fl%idn, RNG)
         norm(N_D, GAS) = norm(N_D, GAS) + sum(cg%wa(RNG)**2, mask=cg%leafmap)
         norm(N_2, GAS) = norm(N_2, GAS) + sum(inid( RNG)**2, mask=cg%leafmap)
         neg_err(GAS) = min(neg_err(GAS), minval(cg%wa(RNG), mask=cg%leafmap))
         pos_err(GAS) = max(pos_err(GAS), maxval(cg%wa(RNG), mask=cg%leafmap))

         if (associated(flind%dst)) then
            cg%wa(RNG) = inid(RNG) - cg%u(flind%dst%idn, RNG)
            norm(N_D, DST_) = norm(N_D, DST_) + sum(cg%wa(RNG)**2, mask=cg%leafmap)
            norm(N_2, DST_) = norm(N_2, DST_) + sum(inid( RNG)**2, mask=cg%leafmap)
            neg_err(DST_) = min(neg_err(DST_), minval(cg%wa(RNG), mask=cg%leafmap))
            pos_err(DST_) = max(pos_err(DST_), maxval(cg%wa(RNG), mask=cg%leafmap))
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

         inid(:,:,:) = pulse_low_density  ! workaround for use of uninitialized values in problem_initial_conditions

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

!> \brief Periodically refine and derefine whole domain

   subroutine flip_flop

      use cg_level_base,      only: base
      use cg_level_connected, only: cg_level_connected_t
      use cg_list,            only: cg_list_element
      use constants,          only: I_TWO
      use global,             only: nstep

      implicit none

      type(cg_level_connected_t), pointer :: curl
      type(cg_list_element),      pointer :: cgl

      integer :: i

      curl => base%level
      do while (associated(curl))
         cgl => curl%first
         i = 0
         do while (associated(cgl))
            call cgl%cg%flag%clear
            cgl%cg%flag%derefine = .false.
            if (real(i)/curl%cnt <= flipratio) then
               if (mod(nstep, nflip) == 0) then
                  if (mod(nstep, I_TWO*nflip) /= 0) call cgl%cg%flag%set
                  cgl%cg%flag%derefine = .not. cgl%cg%flag%get()
               endif
            endif
            i = i + 1
            cgl => cgl%nxt
         enddo
         curl => curl%finer
      enddo

   end subroutine flip_flop

end module initproblem
