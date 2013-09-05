! $Id$
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
#include "macros.h"

module initproblem

   use constants, only: dsetnamelen, ndims, LO, HI

   implicit none

   private
   public :: read_problem_par, problem_initial_conditions, problem_pointers

   ! namelist parameters
   real, dimension(ndims) :: pulse_size  !< size of the density pulse
   real, dimension(ndims) :: pulse_off   !< offset of the pulse
   real                   :: pulse_amp   !< amplitude of the density pulse compared to the ambient level
   real, dimension(ndims) :: pulse_vel   !< uniform velocity components
   integer(kind=4)        :: norm_step   !< how often to calculate the L2-norm
   integer(kind=4)        :: nflip       !< how often to call refine/derefine routine
   real                   :: flipratio   !< percentage of blocks on each level to be refined on flip
   real                   :: ref_thr     !< refinement threshold
   real                   :: deref_thr   !< derefinement threshold

   namelist /PROBLEM_CONTROL/  pulse_size, pulse_off, pulse_vel, pulse_amp, norm_step, nflip, flipratio, ref_thr, deref_thr

   ! other private data
   real, dimension(ndims, LO:HI) :: pulse_edge
   real :: pulse_low_density
   character(len=dsetnamelen), parameter :: inid_n = "inid"

contains

!-----------------------------------------------------------------------------

   subroutine problem_pointers

      use dataio_user, only: user_vars_hdf5
      use user_hooks,  only: finalize_problem, problem_customize_solution

      implicit none

      finalize_problem           => calculate_error_norm
      user_vars_hdf5             => inid_var_hdf5
      problem_customize_solution => calculate_error_norm_wrapper

   end subroutine problem_pointers

!-----------------------------------------------------------------------------

   subroutine read_problem_par

      use constants,  only: I_ONE, xdim, zdim
      use dataio_pub, only: nh      ! QA_WARN required for diff_nml
      use dataio_pub, only: warn, die
      use domain,     only: dom
      use global,     only: smalld
      use mpisetup,   only: rbuff, ibuff, master, slave, proc, have_mpi, LAST, piernik_MPI_Bcast
      use refinement, only: set_n_updAMR, n_updAMR
      use user_hooks, only: problem_refine_derefine

      implicit none

      ! namelist default parameter values
      pulse_size(:) = 1.0                  !< size of the pulse
      pulse_off(:)  = 0.0                  !< center of the pulse
      pulse_vel(:)  = 0.0                  !< pulse velocity
      pulse_amp     = 2.0                  !< pulse relative amplitude
      norm_step     = 5
      nflip         = 0
      ref_thr       = 0.1
      deref_thr     = 0.01
      flipratio     = 1.

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
         rbuff(3)   = deref_thr
         rbuff(4)   = flipratio
         rbuff(20+xdim:20+zdim) = pulse_size(:)
         rbuff(23+xdim:23+zdim) = pulse_vel(:)
         rbuff(26+xdim:26+zdim) = pulse_off(:)

         ibuff(1)   = norm_step
         ibuff(2)   = nflip

      endif

      call piernik_MPI_Bcast(ibuff)
      call piernik_MPI_Bcast(rbuff)

      if (slave) then

         pulse_amp  = rbuff(1)
         ref_thr    = rbuff(2)
         deref_thr  = rbuff(3)
         flipratio  = rbuff(4)
         pulse_size = rbuff(20+xdim:20+zdim)
         pulse_vel  = rbuff(23+xdim:23+zdim)
         pulse_off  = rbuff(26+xdim:26+zdim)

         norm_step  = int(ibuff(1), kind=4)
         nflip      = ibuff(2)

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
         problem_refine_derefine => mark_surface
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
      use global,           only: t
      use grid_cont,        only: grid_container
      use named_array_list, only: qna

      implicit none

      type(cg_list_element), pointer :: cgl
      type(grid_container), pointer :: cg

      integer :: i, j, k

      call analytic_solution(t)

      cgl => leaves%first
      do while (associated(cgl))
         cg => cgl%cg

         cg%b(:, :, :, :) = 0.

         cg%u(flind%dst%idn, cg%is:cg%ie, cg%js:cg%je, cg%ks:cg%ke) = cg%q(qna%ind(inid_n))%arr(cg%is:cg%ie, cg%js:cg%je, cg%ks:cg%ke)

         select case (dom%geometry_type)
            case (GEO_XYZ)
               ! Make uniform, completely boring flow
               cg%u(flind%dst%imx, :, :, :) = pulse_vel(xdim) * cg%u(flind%dst%idn, :, :, :)
               cg%u(flind%dst%imy, :, :, :) = pulse_vel(ydim) * cg%u(flind%dst%idn, :, :, :)
               cg%u(flind%dst%imz, :, :, :) = pulse_vel(zdim) * cg%u(flind%dst%idn, :, :, :)
            case (GEO_RPZ)
               do k = cg%ks, cg%ke
                  do j = cg%js, cg%je
                     do i = cg%is, cg%ie
                        cg%u(flind%dst%imx, i, j, k) = ( pulse_vel(xdim)*cos(cg%y(j)) + pulse_vel(ydim)*sin(cg%y(j))) * cg%u(flind%dst%idn, i, j, k)
                        cg%u(flind%dst%imy, i, j, k) = (-pulse_vel(xdim)*sin(cg%y(j)) + pulse_vel(ydim)*cos(cg%y(j))) * cg%u(flind%dst%idn, i, j, k)
                     enddo
                  enddo
               enddo
               cg%u(flind%dst%imz, :, :, :) = pulse_vel(zdim) * cg%u(flind%dst%idn, :, :, :)
            case default
               call die("[initproblem:problem_initial_conditions] only cartesian and cylindrical geometries are supported")
         end select

         where (cg%u(flind%dst%idn, :, :, :) < 2.*pulse_low_density)
            cg%u(flind%dst%imx, :, :, :) = 0.
            cg%u(flind%dst%imy, :, :, :) = 0.
            cg%u(flind%dst%imz, :, :, :) = 0.
         endwhere
         cgl => cgl%nxt
      enddo

   end subroutine problem_initial_conditions

!-----------------------------------------------------------------------------

   subroutine inid_var_hdf5(var, tab, ierrh, cg)

      use global,           only: t
      use grid_cont,        only: grid_container
      use named_array_list, only: qna

      implicit none

      character(len=*), intent(in)                    :: var
      real(kind=4), dimension(:,:,:), intent(inout)   :: tab
      integer, intent(inout)                          :: ierrh
      type(grid_container), pointer, intent(in)       :: cg

      call analytic_solution(t) ! cannot handle this automagically because here we modify it

      ierrh = 0
      if (qna%exists(var)) then
         tab(:,:,:) = real(cg%q(qna%ind(var))%span(cg%ijkse), 4)
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

      use cg_list,          only: cg_list_element
      use cg_leaves,        only: leaves
      use constants,        only: PIERNIK_FINISHED, pSUM, pMIN, pMAX
      use dataio_pub,       only: code_progress, halfstep, msg, printinfo, warn
      use fluidindex,       only: flind
      use global,           only: t, nstep
      use grid_cont,        only: grid_container
      use mpisetup,         only: master, piernik_MPI_Allreduce
      use named_array_list, only: qna

      implicit none

      enum, bind(C)
         enumerator :: N_D, N_2
      end enum
      real, dimension(N_D:N_2)          :: norm
      real                              :: neg_err, pos_err
      type(cg_list_element),  pointer   :: cgl
      type(grid_container),   pointer   :: cg
      real, dimension(:,:,:), pointer   :: inid

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

         cg%wa(cg%is:cg%ie, cg%js:cg%je, cg%ks:cg%ke) = inid(cg%is:cg%ie, cg%js:cg%je, cg%ks:cg%ke) - cg%u(flind%dst%idn, cg%is:cg%ie, cg%js:cg%je, cg%ks:cg%ke)
         norm(N_D) = norm(N_D) + sum(cg%wa(cg%is:cg%ie, cg%js:cg%je, cg%ks:cg%ke)**2, mask=cg%leafmap)
         norm(N_2) = norm(N_2) + sum(inid( cg%is:cg%ie, cg%js:cg%je, cg%ks:cg%ke)**2, mask=cg%leafmap)
         neg_err   = min(neg_err, minval(cg%wa(cg%is:cg%ie, cg%js:cg%je, cg%ks:cg%ke), mask=cg%leafmap))
         pos_err   = max(pos_err, maxval(cg%wa(cg%is:cg%ie, cg%js:cg%je, cg%ks:cg%ke), mask=cg%leafmap))

         cgl => cgl%nxt
      enddo

      call piernik_MPI_Allreduce(norm(:), pSUM)
      call piernik_MPI_Allreduce(neg_err, pMIN)
      call piernik_MPI_Allreduce(pos_err, pMAX)

      if (master) then
         write(msg,'(a,f12.6,a,2f15.6)')"[initproblem:calculate_error_norm] L2 error norm (DST) = ", sqrt(norm(N_D)/norm(N_2)), &
              ", min and max error = ", neg_err, pos_err
         call printinfo(msg)
      endif

   end subroutine calculate_error_norm

   !>
   !! \brief Put analytic solution in the inid arrays
   !!
   !! \details Density is shaped as an uniform box and translated according to initial velocity and given time
   !<

   subroutine analytic_solution(t)

      use cg_list,          only: cg_list_element
      use cg_leaves,        only: leaves
      use constants,        only: xdim, zdim, ndims, GEO_XYZ, GEO_RPZ
      use dataio_pub,       only: warn, die
      use domain,           only: dom
      use grid_cont,        only: grid_container
      use mpisetup,         only: master
      use named_array_list, only: qna

      implicit none

      real, intent(in)                :: t !< time of the solution

      real                            :: dini
      integer                         :: i, j, k, d
      type(cg_list_element),  pointer :: cgl
      type(grid_container),   pointer :: cg
      real, dimension(:,:,:), pointer :: inid
      real, dimension(ndims)          :: pos

      pos = 0. ! suppres compiler warning
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
                        pos = [cg%x(i), cg%y(j), cg%z(k)] - t * pulse_vel(:)
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
      use cg_level_connected, only: cg_level_connected_T
      use cg_list,            only: cg_list_element
      use constants,          only: I_TWO
      use global,             only: nstep

      implicit none

      type(cg_level_connected_T), pointer :: curl
      type(cg_list_element),      pointer :: cgl

      integer :: i

      curl => base%level
      do while (associated(curl))
         cgl => curl%first
         i = 0
         do while (associated(cgl))
            cgl%cg%refine_flags%refine   = .false.
            cgl%cg%refine_flags%derefine = .false.
            if (real(i)/curl%cnt <= flipratio) then
               if (mod(nstep, nflip) == 0) then
                  cgl%cg%refine_flags%refine   = (mod(nstep, I_TWO*nflip) /= 0)
                  cgl%cg%refine_flags%derefine = .not. cgl%cg%refine_flags%refine
               endif
            endif
            i = i + 1
            cgl => cgl%nxt
         enddo
         curl => curl%finer
      enddo

   end subroutine flip_flop

!> \brief Request refinement along the surface of the pulse. Derefine inside and outside the pulse if possible.

   subroutine mark_surface

      use cg_leaves,        only: leaves
      use cg_list,          only: cg_list_element
      use constants,        only: LO, HI, xdim, ydim, zdim
      use fluidindex,       only: iarr_all_dn
      use named_array_list, only: wna, qna

      implicit none

      type(cg_list_element), pointer :: cgl
      real :: dmin, dmax
      integer :: id

      ! make sure that density is communicated
      !> \todo set up a flag that tells whether this is required or the data has been recently exchanged
      call leaves%internal_boundaries_4d(wna%fi)

      ! fill cg%wa with its guardcells with values corresponding to cgl%cg%leafmap
      !> \todo Consider extending cgl%cg%leafmap into guardcells if it will be helpful in other places too
      cgl => leaves%first
      do while (associated(cgl))
         cgl%cg%wa = 0.
         where (cgl%cg%leafmap(:,:,:)) cgl%cg%wa(cgl%cg%is:cgl%cg%ie, cgl%cg%js:cgl%cg%je, cgl%cg%ks:cgl%cg%ke) = 1.
         cgl => cgl%nxt
      enddo
      call leaves%internal_boundaries_3d(qna%wai)

      ! Detect the edge of the density pulse using density values relative to initial values
      !> \deprecated this method may refine the whole domain when the pulse gets diffused enough
      !> \todo replace with some slope filter
      cgl => leaves%first
      do while (associated(cgl))
         dmax = -huge(1.)
         dmin =  huge(1.)
         do id = lbound(iarr_all_dn, dim=1), ubound(iarr_all_dn, dim=1)
            ! Look one cell beyond local boundary
            dmax = max(dmax, maxval(cgl%cg%u(id, cgl%cg%lh1(xdim, LO):cgl%cg%lh1(xdim, HI), &
                 &                               cgl%cg%lh1(ydim, LO):cgl%cg%lh1(ydim, HI), &
                 &                               cgl%cg%lh1(zdim, LO):cgl%cg%lh1(zdim, HI)), mask = (cgl%cg%wa( &
                 &                               cgl%cg%lh1(xdim, LO):cgl%cg%lh1(xdim, HI), &
                 &                               cgl%cg%lh1(ydim, LO):cgl%cg%lh1(ydim, HI), &
                 &                               cgl%cg%lh1(zdim, LO):cgl%cg%lh1(zdim, HI)) == 1)))
            dmin = min(dmin, minval(cgl%cg%u(id, cgl%cg%lh1(xdim, LO):cgl%cg%lh1(xdim, HI), &
                 &                               cgl%cg%lh1(ydim, LO):cgl%cg%lh1(ydim, HI), &
                 &                               cgl%cg%lh1(zdim, LO):cgl%cg%lh1(zdim, HI)), mask = (cgl%cg%wa( &
                 &                               cgl%cg%lh1(xdim, LO):cgl%cg%lh1(xdim, HI), &
                 &                               cgl%cg%lh1(ydim, LO):cgl%cg%lh1(ydim, HI), &
                 &                               cgl%cg%lh1(zdim, LO):cgl%cg%lh1(zdim, HI)) == 1)))
         enddo
         cgl%cg%refine_flags%derefine = (dmax < (1+deref_thr)*pulse_low_density .or.  dmin > pulse_low_density * (pulse_amp - deref_thr))
         cgl%cg%refine_flags%refine   = (dmax > (1+  ref_thr)*pulse_low_density .and. dmin < pulse_low_density * (pulse_amp -   ref_thr))
         cgl => cgl%nxt
      enddo

   end subroutine mark_surface

end module initproblem
