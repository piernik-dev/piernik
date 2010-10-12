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
!    Initial implemetation of PIERNIK code was based on TVD split MHD code by
!    Ue-Li Pen
!        see: Pen, Arras & Wong (2003) for algorithm and
!             http://www.cita.utoronto.ca/~pen/MHD
!             for original source code "mhd.f90"
!
!    For full list of developers see $PIERNIK_HOME/license/pdt.txt
!

#include "piernik.def"
#include "macros.h"

module initproblem
! ----------------------------------------- !
! Initial condition for a 1D MHD shock tube !
! See: Ryu, Jones, ApJ 442:228-258, (1995)  !
!    and reference therein                  !
! ----------------------------------------- !

   use problem_pub, only: problem_name, run_id

   real             :: dl,vxl,vyl,vzl,bxl,byl,bzl,el
   real             :: dr,vxr,vyr,vzr,bxr,byr,bzr,er
   character(len=1) :: dir

   namelist /PROBLEM_CONTROL/  problem_name, &
      run_id,dl,vxl,vyl,vzl,bxl,byl,bzl,el,  &
      dr,vxr,vyr,vzr,bxr,byr,bzr,er

   contains

!-----------------------------------------------------------------------------

   subroutine read_problem_par
      use mpisetup, only: cbuff_len, cbuff, rbuff, buffer_dim, proc, comm, ierr, &
                           MPI_CHARACTER, MPI_DOUBLE_PRECISION
      use dataio_public, only: ierrh, msg, par_file, namelist_errh, compare_namelist

      implicit none

      problem_name = 'shock'
      run_id  = 'tst'

      if (proc .eq. 0) then

         diff_nml(PROBLEM_CONTROL)

         cbuff(1) =  problem_name
         cbuff(2) =  run_id

         rbuff(1) = dl
         rbuff(2) = vxl
         rbuff(3) = vyl
         rbuff(4) = vzl
         rbuff(5) = bxl
         rbuff(6) = byl
         rbuff(7) = bzl
         rbuff(8) = el
         rbuff(9) = dr
         rbuff(10) = vxr
         rbuff(11) = vyr
         rbuff(12) = vzr
         rbuff(13) = bxr
         rbuff(14) = byr
         rbuff(15) = bzr
         rbuff(16) = er

      endif

      call MPI_Bcast(cbuff, cbuff_len*buffer_dim, MPI_CHARACTER,        0, comm, ierr)
      call MPI_Bcast(rbuff,    buffer_dim, MPI_DOUBLE_PRECISION, 0, comm, ierr)

      if (proc /= 0) then

         problem_name = cbuff(1)
         run_id       = cbuff(2)(1:3)

         dl  = rbuff(1)
         vxl = rbuff(2)
         vyl = rbuff(3)
         vzl = rbuff(4)
         bxl = rbuff(5)
         byl = rbuff(6)
         bzl = rbuff(7)
         el  = rbuff(8)
         dr  = rbuff(9)
         vxr = rbuff(10)
         vyr = rbuff(11)
         vzr = rbuff(12)
         bxr = rbuff(13)
         byr = rbuff(14)
         bzr = rbuff(15)
         er  = rbuff(16)

      endif

   end subroutine read_problem_par

!-----------------------------------------------------------------------------

   subroutine init_prob

      use arrays,       only: u,b
      use grid,         only: x,y,z,nx,ny,nz
      use initionized,  only: idni,imxi,imyi,imzi
#ifndef ISO
      use initionized,  only: ieni, gamma_ion
      use mpisetup,     only: smallei
#endif /* ISO */
      implicit none

      integer  :: i,j,k
      real     :: xi,yj,zk
      real     :: vx,vy,vz,rho,pre,bx,by,bz

      call read_problem_par

!   Secondary parameters

      do j = 1,ny
         yj = y(j)
         do i = 1,nx
            xi = x(i)
            do k = 1,nz
               zk = z(k)

               if ((xi <= 0.5)) then
                  rho = dl
                  pre = el
                  vx  = vxl
                  vy  = vyl
                  vz  = vzl
                  bx  = bxl
                  by  = byl
                  bz  = bzl
               else
                  rho = dr
                  pre = er
                  vx  = vxr
                  vy  = vyr
                  vz  = vzr
                  bx  = bxr
                  by  = byr
                  bz  = bzr
               endif

               u(idni,i,j,k) = rho
               u(imxi,i,j,k) = vx*u(idni,i,j,k)
               u(imyi,i,j,k) = vy*u(idni,i,j,k)
               u(imzi,i,j,k) = vz*u(idni,i,j,k)
#ifndef ISO
               u(ieni,i,j,k) = pre ! pre here means eint
               u(ieni,i,j,k) = max(u(ieni,i,j,k), smallei)
               u(ieni,i,j,k) = u(ieni,i,j,k) +0.5*(vx**2+vy**2+vz**2)*u(idni,i,j,k)
#endif /* !ISO */
               b(1,i,j,k)   =  bx
               b(2,i,j,k)   =  by
               b(3,i,j,k)   =  bz

#ifndef ISO
               u(ieni,i,j,k)   = u(ieni,i,j,k) +0.5*sum(b(:,i,j,k)**2,1)
#endif /* !ISO */
            enddo
         enddo
      enddo
      write(*,*) maxval(b(3,:,:,:)), minval(b(3,:,:,:))
      return
   end subroutine init_prob

   subroutine user_plt(var,ij,xn,tab,ierrh)
      use arrays,         only: u,b
      use grid,           only: nb,nxb,nyb,nzb
      implicit none
      character(LEN=4)     :: var
      character(LEN=2)     :: ij
      integer              :: xn,ierrh
      real, dimension(:,:) :: tab

      ierrh = 0
      select case (var)
         case default
            ierrh = -1
      end select

   end subroutine user_plt

   subroutine user_hdf5(var,tab,ierrh)
!      use arrays,          only: u,b
!      use grid,            only: nb,nx,ny,nz
      implicit none
      character(LEN=4)     :: var
      real(kind=4), dimension(:,:,:) :: tab

      ierrh = 0
      select case (var)
         case default
            ierrh = -1
      end select

   end subroutine user_hdf5

end module initproblem
