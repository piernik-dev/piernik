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

! Initial condition for Sedov-Taylor explosion
! Written by: M. Hanasz, March 2006

   implicit none

   private
   public  :: read_problem_par, problem_initial_conditions, problem_pointers

   integer(kind=4) :: n_sn
   real            :: d0, p0, bx0, by0, bz0, Eexpl, x0, y0, z0, r0, dt_sn, r, t_sn
   real :: ref_thr !< refinement threshold
   real :: deref_thr !< derefinement threshold

   namelist /PROBLEM_CONTROL/ d0, p0, bx0, by0, bz0, Eexpl, x0, y0, z0, r0, n_sn, dt_sn, ref_thr, deref_thr

contains

!-----------------------------------------------------------------------------

   subroutine problem_pointers

      use dataio_user, only: user_tsl
      use user_hooks,  only: problem_refine_derefine
#ifdef HDF5
      use dataio_user, only: user_vars_hdf5
#endif /* HDF5 */

      implicit none

#ifdef HDF5
      user_vars_hdf5 => sedov_vars_hdf5
#endif /* HDF5 */
      user_tsl       => sedov_tsl
      problem_refine_derefine => mark_overdensity

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

         diff_nml(PROBLEM_CONTROL)

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
   subroutine sedov_vars_hdf5(var, tab, ierrh, cg)

      use grid_cont, only: grid_container

      implicit none

      character(len=*),               intent(in)    :: var
      real(kind=4), dimension(:,:,:), intent(inout) :: tab
      integer,                        intent(inout) :: ierrh
      type(grid_container), pointer,  intent(in)    :: cg

      ierrh = 0
      select case (trim(var))
         case ("gid")  ! Grid_id
            tab(:,:,:) = real(cg%grid_id, kind=4)
         case default
            ierrh = -1
      end select

      if (.true. .or. cg%grid_id > 0) return ! suppress compiler warnings

   end subroutine sedov_vars_hdf5
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

!> \brief mark the wave

   subroutine mark_overdensity

      use cg_leaves,        only: leaves
      use cg_list,          only: cg_list_element
      use fluidindex,       only: iarr_all_dn, flind
!      use named_array_list, only: wna
      use refinement,       only: ref_flag

      implicit none

      type(cg_list_element), pointer :: cgl
      real :: dmax
      integer :: id

!      call leaves%internal_boundaries_4d(wna%fi) !< enable it as soon as c2f and f2c routines will work

      cgl => leaves%first
      do while (associated(cgl))
         if (any(cgl%cg%leafmap)) then
            dmax = -huge(1.)
            do id = lbound(iarr_all_dn, dim=1), ubound(iarr_all_dn, dim=1)
               dmax = max(dmax, maxval(cgl%cg%u(id, cgl%cg%is:cgl%cg%ie, cgl%cg%js:cgl%cg%je, cgl%cg%ks:cgl%cg%ke), mask=cgl%cg%leafmap))
            enddo
            do id = 1, flind%energ
               if ( maxval(cgl%cg%u(flind%all_fluids(id)%fl%ien, cgl%cg%is:cgl%cg%ie, cgl%cg%js:cgl%cg%je, cgl%cg%ks:cgl%cg%ke), mask=cgl%cg%leafmap) / &
                    minval(cgl%cg%u(flind%all_fluids(id)%fl%ien, cgl%cg%is:cgl%cg%ie, cgl%cg%js:cgl%cg%je, cgl%cg%ks:cgl%cg%ke), mask=cgl%cg%leafmap) > ref_thr) &
                    dmax = 2.*ref_thr*d0 !trick
            enddo
            cgl%cg%refine_flags = ref_flag( dmax >= ref_thr*d0, dmax < deref_thr*d0 )
         endif
         cgl => cgl%nxt
      enddo

   end subroutine mark_overdensity

end module initproblem
