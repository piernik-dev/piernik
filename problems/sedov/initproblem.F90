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

! Initial condition for Sedov-Taylor explosion
! Written by: M. Hanasz, March 2006

   implicit none

   private
   public  :: read_problem_par, problem_initial_conditions, problem_pointers

   integer(kind=4) :: n_sn
   real            :: d0, p0, bx0, by0, bz0, Eexpl, x0, y0, z0, r0, dt_sn, r, t_sn, dtrig
   real :: ref_thr !< refinement threshold
   real :: deref_thr !< derefinement threshold

   namelist /PROBLEM_CONTROL/ d0, p0, bx0, by0, bz0, Eexpl, x0, y0, z0, r0, n_sn, dt_sn, ref_thr, deref_thr, dtrig

contains

!-----------------------------------------------------------------------------

   subroutine problem_pointers

      use dataio_user, only: user_tsl
      use user_hooks,  only: problem_refine_derefine, problem_domain_update, late_initial_conditions

      implicit none

      user_tsl       => sedov_tsl
      problem_refine_derefine => mark_overdensity

      problem_domain_update => sedov_dist_to_edge
      late_initial_conditions => sedov_late_init

   end subroutine problem_pointers

!-----------------------------------------------------------------------------

   subroutine read_problem_par

      use constants,  only: DST
      use dataio_pub, only: nh      ! QA_WARN required for diff_nml
      use dataio_pub, only: msg, printinfo, die
      use domain,     only: dom
      use fluidindex, only: flind
      use mpisetup,   only: ibuff, rbuff, master, slave, piernik_MPI_Bcast

      implicit none

      integer :: p

      t_sn = 0.0

      d0      = 1.0
      dtrig   = -1.5
      p0      = 1.e-3
      bx0     =   0.
      by0     =   0.
      bz0     =   0.
      Eexpl   = 1.e5
      x0      = 0.0
      y0      = 0.0
      z0      = 0.0
      r0      = minval(dom%L_(:)/dom%n_d(:), mask=dom%has_dir(:))/2.
      n_sn    = 1
      dt_sn   = 0.0
      ref_thr      = 3.    !< Refine if density is greater than this value
      deref_thr    = 1.5    !< Derefine if density is smaller than this value

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

         rbuff(1) = d0
         rbuff(2) = p0
         rbuff(3) = bx0
         rbuff(4) = by0
         rbuff(5) = bz0
         rbuff(6) = Eexpl
         rbuff(7) = x0
         rbuff(8) = y0
         rbuff(9) = z0
         rbuff(10)= r0
         rbuff(11)= dt_sn
         rbuff(12)= ref_thr
         rbuff(13)= deref_thr
         rbuff(14)= dtrig

         ibuff(1) = n_sn

      endif

      call piernik_MPI_Bcast(ibuff)
      call piernik_MPI_Bcast(rbuff)

      if (slave) then

         d0           = rbuff(1)
         p0           = rbuff(2)
         bx0          = rbuff(3)
         by0          = rbuff(4)
         bz0          = rbuff(5)
         Eexpl        = rbuff(6)
         x0           = rbuff(7)
         y0           = rbuff(8)
         z0           = rbuff(9)
         r0           = rbuff(10)
         dt_sn        = rbuff(11)
         ref_thr      = rbuff(12)
         deref_thr    = rbuff(13)
         dtrig        = rbuff(14)

         n_sn         = ibuff(1)

      endif

      if (flind%energ < flind%fluids) call die("[initproblem:problem_initial_conditions] Not all fluids are adiabatic!")
      do p = 1, flind%energ
         if (flind%all_fluids(p)%fl%tag == DST) call die("[initproblem:problem_initial_conditions] This setup is not suitable for dust!")
         write(msg, '(a,i2)')"Working with fluid#", flind%all_fluids(p)%fl%tag
         if (master) call printinfo(msg)
      enddo

   end subroutine read_problem_par
!-----------------------------------------------------------------------------
   subroutine problem_initial_conditions

      use cg_leaves,  only: leaves
      use cg_list,    only: cg_list_element
      use constants,  only: ION, xdim, ydim, zdim, LO, HI
      use domain,     only: dom
      use fluidindex, only: flind
      use grid_cont,  only: grid_container

      implicit none

      integer, parameter              :: isub = 4
      integer                         :: i, j, k, p, ii, jj, kk
      type(cg_list_element),  pointer :: cgl
      type(grid_container),   pointer :: cg
      real :: x, y, z

      ! BEWARE:
      !  3 triple loop are completely unnecessary here, but this problem serves
      !  as an educational example

      do p = 1, flind%energ
         associate(fl => flind%all_fluids(p)%fl)

! Uniform equilibrium state

         cgl => leaves%first
         do while (associated(cgl))
            cg => cgl%cg

            do k = cg%lhn(zdim,LO), cg%lhn(zdim,HI)
               do j = cg%lhn(ydim,LO), cg%lhn(ydim,HI)
                  do i = cg%lhn(xdim,LO), cg%lhn(xdim,HI)
                     cg%u(fl%idn,i,j,k) = d0
                     cg%u(fl%imx,i,j,k) = 0.0
                     cg%u(fl%imy,i,j,k) = 0.0
                     cg%u(fl%imz,i,j,k) = 0.0
                     cg%u(fl%ien,i,j,k) = p0/(fl%gam_1)
                     cg%u(fl%ien,i,j,k) = cg%u(fl%ien,i,j,k) + 0.5*(cg%u(fl%imx,i,j,k)**2 +cg%u(fl%imy,i,j,k)**2 + cg%u(fl%imz,i,j,k)**2)/cg%u(fl%idn,i,j,k)
                  enddo
               enddo
            enddo

! Explosion

            do k = cg%lhn(zdim,LO), cg%lhn(zdim,HI)
               do j = cg%lhn(ydim,LO), cg%lhn(ydim,HI)
                  do i = cg%lhn(xdim,LO), cg%lhn(xdim,HI)
                     r = sqrt( (cg%x(i)-x0)**2 + (cg%y(j)-y0)**2 + (cg%z(k)-z0)**2 )
                     if ( r**2 < r0**2+2*sum(cg%dl**2, mask=dom%has_dir)) then
                        do ii = 1, isub
                           x = cg%x(i)-x0
                           if (dom%has_dir(xdim)) x = x + cg%dx*(2*ii -isub - 1)/real(2*isub)
                           do jj = 1, isub
                              y = cg%y(j)-y0
                              if (dom%has_dir(ydim)) y = y + cg%dy*(2*jj -isub - 1)/real(2*isub)
                              do kk = 1, isub
                                 z = cg%z(k)-z0
                                 if (dom%has_dir(zdim)) z = z + cg%dx*(2*kk -isub - 1)/real(2*isub)
                                 if (x*x+y*y+z*z < r0**2) cg%u(fl%ien,i,j,k)   = cg%u(fl%ien,i,j,k) + Eexpl/isub**dom%eff_dim
                              enddo
                           enddo
                        enddo
                     endif
                  enddo
               enddo
            enddo

            if (fl%tag == ION) then
               do k = cg%lhn(zdim,LO), cg%lhn(zdim,HI)
                  do j = cg%lhn(ydim,LO), cg%lhn(ydim,HI)
                     do i = cg%lhn(xdim,LO), cg%lhn(xdim,HI)
                        cg%b(xdim,i,j,k) = bx0
                        cg%b(ydim,i,j,k) = by0
                        cg%b(zdim,i,j,k) = bz0
                        cg%u(fl%ien,i,j,k) = cg%u(fl%ien,i,j,k) + 0.5*(cg%b(xdim,i,j,k)**2 + cg%b(ydim,i,j,k)**2 + cg%b(zdim,i,j,k)**2)
                     enddo
                  enddo
               enddo
            endif

            cgl => cgl%nxt
         enddo

#ifdef TRACER
         do k = cg%lhn(zdim,LO), cg%lhn(zdim,HI)
            do j = cg%lhn(ydim,LO), cg%lhn(ydim,HI)
               do i = cg%lhn(xdim,LO), cg%lhn(xdim,HI)
                  r = sqrt( (cg%x(i)-x0)**2 + (cg%y(j)-y0)**2 + (cg%z(k)-z0)**2 )
                  if ( r**2 < r0**2) then
                     cg%u(flind%trc%pos, i, j, k)   = 1.0
                  else
                     cg%u(flind%trc%pos, i, j, k)   = 0.0
                  endif
               enddo
            enddo
         enddo
#endif /* TRACER */
         end associate
      enddo

   end subroutine problem_initial_conditions
!-----------------------------------------------------------------------------
   subroutine sedov_tsl(user_vars, tsl_names)

      use constants,   only: pSUM
      use diagnostics, only: pop_vector
      use mpisetup,    only: proc, master, piernik_MPI_Allreduce

      implicit none

      real, dimension(:), intent(inout), allocatable                       :: user_vars
      character(len=*), dimension(:), intent(inout), allocatable, optional :: tsl_names
      real :: output

      if (present(tsl_names)) then
         call pop_vector(tsl_names, len(tsl_names(1)), ["foobar_sedov"])    !   add to header
      else
         ! do mpi stuff here...
         output = real(proc,8)
         call piernik_MPI_Allreduce(output, pSUM)
         if (master) call pop_vector(user_vars,[output])                 !   pop value
      endif

   end subroutine sedov_tsl

!> \brief Find hov close it the shockwave to the external edges and call expansion routine if necessary

   subroutine sedov_dist_to_edge

      use cg_leaves,     only: leaves
      use cg_level_base, only: base
      use cg_list,       only: cg_list_element
      use constants,     only: xdim, ydim, zdim, LO, HI
      use domain,        only: dom
      use fluidindex,    only: iarr_all_dn

      implicit none

      type(cg_list_element),  pointer :: cgl
      real, dimension(xdim:zdim, LO:HI) :: ddist
      integer :: i
      integer, parameter :: iprox = 2

      if (dtrig < 0.) return

      ddist = huge(1.)
      cgl => leaves%first
      do while (associated(cgl))
         if (any(cgl%cg%ext_bnd)) then
            !> \todo roll it to a nested loop
            if (dom%has_dir(xdim)) then
               if (cgl%cg%ext_bnd(xdim, LO)) then
                  do i = cgl%cg%is, cgl%cg%ie
                     if (any(cgl%cg%u(iarr_all_dn, i, :, :) > d0*dtrig)) then
                        ddist(xdim, LO) = min(ddist(xdim, LO), (cgl%cg%x(i) - cgl%cg%fbnd(xdim, LO))/cgl%cg%dx)
                        exit
                     endif
                  enddo
               endif
               if (cgl%cg%ext_bnd(xdim, HI)) then
                  do i = cgl%cg%ie, cgl%cg%is, -1
                     if (any(cgl%cg%u(iarr_all_dn, i, :, :) > d0*dtrig)) then
                        ddist(xdim, HI) = min(ddist(xdim, HI), (cgl%cg%fbnd(xdim, HI) - cgl%cg%x(i))/cgl%cg%dx)
                        exit
                     endif
                  enddo
               endif
            endif

            if (dom%has_dir(ydim)) then
               if (cgl%cg%ext_bnd(ydim, LO)) then
                  do i = cgl%cg%js, cgl%cg%je
                     if (any(cgl%cg%u(iarr_all_dn, :, i, :) > d0*dtrig)) then
                        ddist(ydim, LO) = min(ddist(ydim, LO), (cgl%cg%y(i) - cgl%cg%fbnd(ydim, LO))/cgl%cg%dy)
                        exit
                     endif
                  enddo
               endif
               if (cgl%cg%ext_bnd(ydim, HI)) then
                  do i = cgl%cg%je, cgl%cg%js, -1
                     if (any(cgl%cg%u(iarr_all_dn, :, i, :) > d0*dtrig)) then
                        ddist(ydim, HI) = min(ddist(ydim, HI), (cgl%cg%fbnd(ydim, HI) - cgl%cg%y(i))/cgl%cg%dy)
                        exit
                     endif
                  enddo
               endif
            endif

            if (dom%has_dir(zdim)) then
               if (cgl%cg%ext_bnd(zdim, LO)) then
                  do i = cgl%cg%ks, cgl%cg%ke
                     if (any(cgl%cg%u(iarr_all_dn, :, :, i) > d0*dtrig)) then
                        ddist(zdim, LO) = min(ddist(zdim, LO), (cgl%cg%z(i) - cgl%cg%fbnd(zdim, LO))/cgl%cg%dz)
                        exit
                     endif
                  enddo
               endif
               if (cgl%cg%ext_bnd(zdim, HI)) then
                  do i = cgl%cg%ke, cgl%cg%ks, -1
                     if (any(cgl%cg%u(iarr_all_dn, :, :, i) > d0*dtrig)) then
                        ddist(zdim, HI) = min(ddist(zdim, HI), (cgl%cg%fbnd(zdim, HI) - cgl%cg%z(i))/cgl%cg%dz)
                        exit
                     endif
                  enddo
               endif
            endif
         endif
         cgl => cgl%nxt
      enddo

      call base%expand(ddist(:,:) < iprox)

   end subroutine sedov_dist_to_edge

!> \brief Performa late initialization of the cg added after domain expansion

   subroutine sedov_late_init

      use cg_list,        only: cg_list_element
      use cg_list_dataop, only: expanded_domain
      use constants,      only: xdim, ydim, zdim, ION
      use dataio_pub,     only: die
      use fluidindex,     only: flind
      use func,           only: ekin, emag

      implicit none

      type(cg_list_element),  pointer :: cgl

      integer :: p

      cgl => expanded_domain%first
      do while (associated(cgl))
         if (cgl%cg%is_old) call die("[initproblem:sedov_late_init] Old piece on a new list")
         do p = 1, flind%energ
            associate (fl => flind%all_fluids(p)%fl)
            cgl%cg%u(fl%idn, :, :, :) = d0
            cgl%cg%u(fl%imx, :, :, :) = 0.
            cgl%cg%u(fl%imy, :, :, :) = 0.
            cgl%cg%u(fl%imz, :, :, :) = 0.
            cgl%cg%u(fl%ien, :, :, :) = p0/(fl%gam_1) + ekin(cgl%cg%u(fl%imx, :, :, :), cgl%cg%u(fl%imy, :, :, :), cgl%cg%u(fl%imz, :, :, :), cgl%cg%u(fl%idn, :, :, :))
            if (fl%tag == ION) then
               cgl%cg%b(xdim, :, :, :) = bx0
               cgl%cg%b(ydim, :, :, :) = by0
               cgl%cg%b(zdim, :, :, :) = bz0
               cgl%cg%u(fl%ien, :, :, :) = cgl%cg%u(fl%ien, :, :, :) + emag(cgl%cg%b(xdim, :, :, :), cgl%cg%b(ydim, :, :, :), cgl%cg%b(zdim, :, :, :)**2)
            endif
            end associate
         enddo
         cgl => cgl%nxt
      enddo

   end subroutine sedov_late_init

!> \brief mark the wave

   subroutine mark_overdensity

      use all_boundaries, only: all_bnd
      use cg_leaves,      only: leaves
      use cg_list,        only: cg_list_element
      use domain,         only: dom
      use fluidindex,     only: iarr_all_dn, flind

      implicit none

      type(cg_list_element), pointer :: cgl
      real :: dmax, diff
      integer :: id, i, j, k

      call all_bnd ! pretty likely overkill. \todo find a way to minimize calling this - perhaps manage a flag that says whether the boundaries are up to date or not

      cgl => leaves%first
      do while (associated(cgl))
         if (any(cgl%cg%leafmap)) then
            dmax = -huge(1.)
            do id = lbound(iarr_all_dn, dim=1), ubound(iarr_all_dn, dim=1)
               do i = cgl%cg%is, cgl%cg%ie
                  do j = cgl%cg%js, cgl%cg%je
                     do k = cgl%cg%ks, cgl%cg%ke
                        diff = maxval(abs(cgl%cg%u(id, i, j, k) - [ &
                             &            cgl%cg%u(id, i+dom%D_x, j, k), &
                             &            cgl%cg%u(id, i-dom%D_x, j, k), &
                             &            cgl%cg%u(id, i, j+dom%D_y, k), &
                             &            cgl%cg%u(id, i, j-dom%D_y, k), &
                             &            cgl%cg%u(id, i, j, k+dom%D_z), &
                             &            cgl%cg%u(id, i, j, k-dom%D_z) ] ) )
                        cgl%cg%refinemap(i, j, k) = (diff > ref_thr * d0)
                        dmax = max(dmax, diff)
                     enddo
                  enddo
               enddo
            enddo
            do id = 1, flind%energ
               do i = cgl%cg%is, cgl%cg%ie
                  do j = cgl%cg%js, cgl%cg%je
                     do k = cgl%cg%ks, cgl%cg%ke
                        diff = maxval(abs(cgl%cg%u(flind%all_fluids(id)%fl%ien, i, j, k) / [ &
                             &            cgl%cg%u(flind%all_fluids(id)%fl%ien, i+dom%D_x, j, k), &
                             &            cgl%cg%u(flind%all_fluids(id)%fl%ien, i-dom%D_x, j, k), &
                             &            cgl%cg%u(flind%all_fluids(id)%fl%ien, i, j+dom%D_y, k), &
                             &            cgl%cg%u(flind%all_fluids(id)%fl%ien, i, j-dom%D_y, k), &
                             &            cgl%cg%u(flind%all_fluids(id)%fl%ien, i, j, k+dom%D_z), &
                             &            cgl%cg%u(flind%all_fluids(id)%fl%ien, i, j, k-dom%D_z) ] ) )
                        cgl%cg%refinemap(i, j, k) = cgl%cg%refinemap(i, j, k) .or. (diff > ref_thr)
                        dmax = max(dmax, diff)
                     enddo
                  enddo
               enddo
            enddo
            cgl%cg%refine_flags%derefine = (dmax <  deref_thr*d0)
         endif
         cgl => cgl%nxt
      enddo

   end subroutine mark_overdensity

end module initproblem
