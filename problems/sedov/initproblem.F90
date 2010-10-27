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
#include "piernik.def"
#include "macros.h"
#define RNG is:ie, js:je, ks:ke
module initproblem

! Initial condition for Sedov-Taylor explosion
! Written by: M. Hanasz, March 2006

   use problem_pub, only: problem_name, run_id

   integer :: n_sn
   real    :: d0, p0, bx0, by0, bz0, Eexpl, x0, y0, z0, r0, dt_sn, r, t_sn

   namelist /PROBLEM_CONTROL/  problem_name, run_id, &
                               d0,p0, bx0,by0,bz0, Eexpl,  x0,y0,z0, r0, &
                               n_sn, dt_sn
contains
!-----------------------------------------------------------------------------
   subroutine read_problem_par
      use dataio_public, only: ierrh, msg, par_file, namelist_errh, compare_namelist
      use grid,          only: dxmn
      use mpisetup,      only: cbuff_len, cbuff, ibuff, rbuff, buffer_dim, proc, comm, ierr, &
                               MPI_CHARACTER, MPI_DOUBLE_PRECISION, MPI_INTEGER
      use types,         only: idlen

      implicit none


      t_sn = 0.0

      problem_name = 'aaa'
      run_id  = 'aaa'
      d0      = 1.0
      p0      = 1.e-3
      bx0     =   0.
      by0     =   0.
      bz0     =   0.
      Eexpl   = 1.e5
      x0      = 0.0
      y0      = 0.0
      z0      = 0.0
      r0      = dxmn/2.
      n_sn    = 1
      dt_sn   = 0.0

      if (proc .eq. 0) then

         diff_nml(PROBLEM_CONTROL)

         cbuff(1) =  problem_name
         cbuff(2) =  run_id

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

      call MPI_Bcast(cbuff, cbuff_len*buffer_dim, MPI_CHARACTER,        0, comm, ierr)
      call MPI_Bcast(ibuff,    buffer_dim, MPI_INTEGER,          0, comm, ierr)
      call MPI_Bcast(rbuff,    buffer_dim, MPI_DOUBLE_PRECISION, 0, comm, ierr)

      if (proc /= 0) then

         problem_name = cbuff(1)
         run_id       = cbuff(2)(1:idlen)

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
      use arrays,         only: u, b
      use dataio_public,  only: msg, die, printinfo, user_plt_hdf5, user_vars_hdf5, user_tsl
      use fluidindex,     only: nvar, ibx, iby, ibz
      use grid,           only: x, y, z, nx, ny, nz
      use types,          only: component_fluid

      implicit none
      integer :: i, j, k, p

      type(component_fluid), pointer :: fl

      if (nvar%adiab < nvar%fluids) call die("[initproblem:init_prob] Not all fluids are adiabatic!")

      ! BEWARE:
      !  3 triple loop are completely unnecessary here, but this problem serves
      !  as an educational example

      do p = 1, nvar%adiab
         fl => nvar%all_fluids(p)
         if (fl%tag=="DST") call die("[initproblem:init_prob] This setup is not suitable for dust!")
         write(msg,*) "Working with ", fl%tag, " fluid."
         call printinfo(msg)


! Uniform equilibrium state
         do k = 1,nz
            do j = 1,ny
               do i = 1,nx
                  u(fl%idn,i,j,k) = d0
                  u(fl%imx,i,j,k) = 0.0
                  u(fl%imy,i,j,k) = 0.0
                  u(fl%imz,i,j,k) = 0.0
                  u(fl%ien,i,j,k) = p0/(fl%gam_1)
                  u(fl%ien,i,j,k) = u(fl%ien,i,j,k) + 0.5*(u(fl%imx,i,j,k)**2 +u(fl%imy,i,j,k)**2 + u(fl%imz,i,j,k)**2)/u(fl%idn,i,j,k)
               enddo
            enddo
         enddo

! Explosion

         do k = 1,nz
            do j = 1,ny
               do i = 1,nx
                  r = dsqrt( (x(i)-x0)**2 + (y(j)-y0)**2 + (z(k)-z0)**2 )
                  if ( r**2 < r0**2) u(fl%ien,i,j,k)   = u(fl%ien,i,j,k) + Eexpl
               enddo
            enddo
         enddo

         if (fl%tag=="ION") then
            do k = 1,nz
               do j = 1,ny
                  do i = 1,nx
                     b(ibx,i,j,k) = bx0
                     b(iby,i,j,k) = by0
                     b(ibz,i,j,k) = bz0
                     u(fl%ien,i,j,k) = u(fl%ien,i,j,k) + 0.5*(b(ibx,i,j,k)**2 + b(iby,i,j,k)**2 + b(ibz,i,j,k)**2)
                  enddo
               enddo
            enddo
         endif
      enddo
      user_plt_hdf5 => sedov_plt_hdf5
      user_vars_hdf5 => sedov_vars_hdf5
      user_tsl => sedov_tsl
      return
   end subroutine init_prob
!-----------------------------------------------------------------------------
   subroutine sedov_plt_hdf5(var,ij,xn,tab,ierrh)
      use arrays,        only: u
      use dataio_public, only: varlen
      use grid,          only: nb, nxb, nyb, nzb
      implicit none
      character(LEN=*), intent(in)        :: var   !< quantity to be plotted
      character(LEN=*), intent(in)        :: ij    !< plane of plot
      integer, intent(in)                 :: xn    !< no. of cell at which we are slicing the local block
      integer, intent(inout)              :: ierrh !< error handling
      real, dimension(:,:), intent(inout) :: tab   !< array  containing given quantity

      ierrh = 0
      select case (var)
         case ("fooo")   ! Totally bogus quantity, just to check user_plt_hdf5 works
            if (ij=="yz") tab(:,:) = u(2,xn,nb+1:nyb+nb,nb+1:nzb+nb)*u(3,xn,nb+1:nyb+nb,nb+1:nzb+nb)* .123456789
            if (ij=="xz") tab(:,:) = u(2,nb+1:nxb+nb,xn,nb+1:nzb+nb)*u(3,nb+1:nxb+nb,xn,nb+1:nzb+nb)* .123456789
            if (ij=="xy") tab(:,:) = u(2,nb+1:nxb+nb,nb+1:nyb+nb,xn)*u(3,nb+1:nxb+nb,nb+1:nyb+nb,xn)* .123456789
         case default
            ierrh = -1
      end select

   end subroutine sedov_plt_hdf5
!-----------------------------------------------------------------------------
   subroutine sedov_vars_hdf5(var,tab, ierrh)
      use arrays, only: u
      use grid,   only: is, ie, js, je, ks, ke
      implicit none
      character(len=*), intent(in)                    :: var
      real(kind=4), dimension(:,:,:), intent(inout)   :: tab
      integer, intent(inout)                          :: ierrh

      ierrh = 0
      select case (trim(var))
         case ("fooo")  ! Totally bogus quantity, just to check user_vars_hdf5 works
            tab(:,:,:) = .123456789
         case default
            ierrh = -1
      end select
      return
   end subroutine sedov_vars_hdf5
!-----------------------------------------------------------------------------
      subroutine sedov_tsl(user_vars, tsl_names)
         use diagnostics,     only: pop_vector
         implicit none
         real, dimension(:), intent(inout), allocatable                       :: user_vars
         character(len=*), dimension(:), intent(inout), allocatable, optional :: tsl_names

         if (present(tsl_names)) then
            call pop_vector(tsl_names, len(tsl_names(1)), ["foobar_sedov"])    !   add to header
         else
            call pop_vector(user_vars,[12345678.9])                            !   pop value
         endif
         return
      end subroutine sedov_tsl
!-----------------------------------------------------------------------------
end module initproblem
