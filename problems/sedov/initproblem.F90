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
#define RNG is:ie, js:je, ks:ke
module initproblem

! Initial condition for Sedov-Taylor explosion
! Written by: M. Hanasz, March 2006

   implicit none

   private
   public  :: read_problem_par, init_prob, problem_pointers

   integer :: n_sn
   real    :: d0, p0, bx0, by0, bz0, Eexpl, x0, y0, z0, r0, dt_sn, r, t_sn

   namelist /PROBLEM_CONTROL/ d0, p0, bx0, by0, bz0, Eexpl, x0, y0, z0, r0, n_sn, dt_sn

contains

!-----------------------------------------------------------------------------

   subroutine problem_pointers

      use dataio_user, only: user_plt_hdf5, user_vars_hdf5, user_tsl

      implicit none

      user_plt_hdf5  => sedov_plt_hdf5
      user_vars_hdf5 => sedov_vars_hdf5
      user_tsl       => sedov_tsl

   end subroutine problem_pointers

!-----------------------------------------------------------------------------

   subroutine read_problem_par

      use dataio_pub,  only: ierrh, par_file, namelist_errh, compare_namelist, cmdl_nml      ! QA_WARN required for diff_nml
      use domain,      only: dom, has_dir
      use mpi,         only: MPI_DOUBLE_PRECISION, MPI_INTEGER
      use mpisetup,    only: ibuff, rbuff, buffer_dim, master, slave, comm, ierr

      implicit none

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
      r0      = minval([dom%Lx, dom%Ly, dom%Lz]/dom%n_d(:), mask=has_dir(:))/2.
      n_sn    = 1
      dt_sn   = 0.0

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

         ibuff(1) = n_sn

      endif

      call MPI_Bcast(ibuff, buffer_dim, MPI_INTEGER,          0, comm, ierr)
      call MPI_Bcast(rbuff, buffer_dim, MPI_DOUBLE_PRECISION, 0, comm, ierr)

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

         n_sn         = ibuff(1)

      endif

   end subroutine read_problem_par
!-----------------------------------------------------------------------------
   subroutine init_prob

      use constants,  only: ION, DST
      use dataio_pub, only: msg, die, printinfo
      use fluidindex, only: flind, ibx, iby, ibz
      use fluidtypes, only: component_fluid
      use grid,       only: cga
      use grid_cont,  only: cg_list_element, grid_container
      use mpisetup,   only: master

      implicit none

      integer :: i, j, k, p
      type(component_fluid), pointer :: fl
      type(cg_list_element), pointer :: cgl
      type(grid_container), pointer :: cg

      if (flind%energ < flind%fluids) call die("[initproblem:init_prob] Not all fluids are adiabatic!")

      ! BEWARE:
      !  3 triple loop are completely unnecessary here, but this problem serves
      !  as an educational example

      do p = 1, flind%energ
         fl => flind%all_fluids(p)
         if (fl%tag == DST) call die("[initproblem:init_prob] This setup is not suitable for dust!")
         write(msg,*) "Working with ", fl%tag, " fluid."
         if (master) call printinfo(msg)

! Uniform equilibrium state

         cgl => cga%cg_leafs%cg_l(1)
         do while (associated(cgl))
            cg => cgl%cg

            do k = 1, cg%nz
               do j = 1, cg%ny
                  do i = 1, cg%nx
                     cg%u%arr(fl%idn,i,j,k) = d0
                     cg%u%arr(fl%imx,i,j,k) = 0.0
                     cg%u%arr(fl%imy,i,j,k) = 0.0
                     cg%u%arr(fl%imz,i,j,k) = 0.0
                     cg%u%arr(fl%ien,i,j,k) = p0/(fl%gam_1)
                     cg%u%arr(fl%ien,i,j,k) = cg%u%arr(fl%ien,i,j,k) + 0.5*(cg%u%arr(fl%imx,i,j,k)**2 +cg%u%arr(fl%imy,i,j,k)**2 + cg%u%arr(fl%imz,i,j,k)**2)/cg%u%arr(fl%idn,i,j,k)
                  enddo
               enddo
            enddo

! Explosion

            do k = 1, cg%nz
               do j = 1, cg%ny
                  do i = 1, cg%nx
                     r = sqrt( (cg%x(i)-x0)**2 + (cg%y(j)-y0)**2 + (cg%z(k)-z0)**2 )
                     if ( r**2 < r0**2) cg%u%arr(fl%ien,i,j,k)   = cg%u%arr(fl%ien,i,j,k) + Eexpl
                  enddo
               enddo
            enddo

            if (fl%tag == ION) then
               do k = 1, cg%nz
                  do j = 1, cg%ny
                     do i = 1, cg%nx
                        cg%b%arr(ibx,i,j,k) = bx0
                        cg%b%arr(iby,i,j,k) = by0
                        cg%b%arr(ibz,i,j,k) = bz0
                        cg%u%arr(fl%ien,i,j,k) = cg%u%arr(fl%ien,i,j,k) + 0.5*(cg%b%arr(ibx,i,j,k)**2 + cg%b%arr(iby,i,j,k)**2 + cg%b%arr(ibz,i,j,k)**2)
                     enddo
                  enddo
               enddo
            endif

            cgl => cgl%nxt
         enddo

      enddo

   end subroutine init_prob
!-----------------------------------------------------------------------------
   subroutine sedov_plt_hdf5(var, ij, xn, tab, ierrh, cg)

      use constants, only: xdim, ydim, zdim
      use grid_cont, only: grid_container

      implicit none

      character(len=*), intent(in)        :: var   !< quantity to be plotted
      integer, intent(in)                 :: ij    !< plane of plot
      integer(kind=8), intent(in)         :: xn    !< no. of cell at which we are slicing the local block
      integer, intent(inout)              :: ierrh !< error handling
      real, dimension(:,:), intent(inout) :: tab   !< array  containing given quantity
      type(grid_container), pointer, intent(in) :: cg

      ierrh = 0
      select case (var)
         case ("fooo")   ! Totally bogus quantity, just to check user_plt_hdf5 works
            if (ij==xdim) tab(:,:) = cg%u%arr(2, xn, cg%js:cg%je, cg%ks:cg%ke)*cg%u%arr(3, xn, cg%js:cg%je, cg%ks:cg%ke)* .123456789
            if (ij==ydim) tab(:,:) = cg%u%arr(2, cg%is:cg%ie, xn, cg%ks:cg%ke)*cg%u%arr(3, cg%is:cg%ie, xn, cg%ks:cg%ke)* .123456789
            if (ij==zdim) tab(:,:) = cg%u%arr(2, cg%is:cg%ie, cg%js:cg%je, xn)*cg%u%arr(3, cg%is:cg%ie, cg%js:cg%je, xn)* .123456789
         case default
            ierrh = -1
      end select

   end subroutine sedov_plt_hdf5
!-----------------------------------------------------------------------------
   subroutine sedov_vars_hdf5(var, tab, ierrh, cg)

      use grid_cont,  only: grid_container

      implicit none

      character(len=*), intent(in)                    :: var
      real(kind=4), dimension(:,:,:), intent(inout)   :: tab
      integer, intent(inout)                          :: ierrh
      type(grid_container), pointer, intent(in)       :: cg

      ierrh = 0
      select case (trim(var))
         case ("fooo")  ! Totally bogus quantity, just to check user_vars_hdf5 works
            tab(:,:,:) = .123456789
         case default
            ierrh = -1
      end select

      if (.true. .or. cg%empty) return ! suppress compiler warnings

   end subroutine sedov_vars_hdf5
!-----------------------------------------------------------------------------
   subroutine sedov_tsl(user_vars, tsl_names)

      use diagnostics, only: pop_vector
      use mpi,         only: MPI_DOUBLE_PRECISION, MPI_SUM
      use mpisetup,    only: proc, master, comm, ierr

      implicit none

      real, dimension(:), intent(inout), allocatable                       :: user_vars
      character(len=*), dimension(:), intent(inout), allocatable, optional :: tsl_names
      real :: output

      if (present(tsl_names)) then
         call pop_vector(tsl_names, len(tsl_names(1)), ["foobar_sedov"])    !   add to header
      else
         ! do mpi stuff here...
         call MPI_Allreduce(real(proc,8), output, 1, MPI_DOUBLE_PRECISION, MPI_SUM, comm, ierr)
         if (master) call pop_vector(user_vars,[output])                 !   pop value
      endif

   end subroutine sedov_tsl
!-----------------------------------------------------------------------------
end module initproblem
