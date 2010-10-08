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

! Initial condition for Keplerian disk
! Written by: M. Hanasz, March 2006

   use problem_pub, only: problem_name, run_id
   use mpisetup,    only: cbuff_len

   real    :: sigma0, Rin, R0, HtoR, eps, amp
   character(len=cbuff_len) :: sigma_model

   namelist /PROBLEM_CONTROL/  problem_name, run_id, sigma0, amp, &
                               Rin, R0, HtoR, sigma_model, eps

   contains

!-----------------------------------------------------------------------------

   subroutine read_problem_par

      use errh,     only: namelist_errh
      use mpisetup, only: cbuff_len, cbuff, rbuff, buffer_dim, proc, comm, ierr, &
                           MPI_CHARACTER, MPI_DOUBLE_PRECISION
      use dataio_public, only: ierrh, msg, par_file
      use func,          only: compare_namelist

      implicit none


      problem_name = 'stream_3D'
      run_id  = 'ts1'
      Sigma0  = 1.0
      Rin     = 1.0e-4
      R0      = 1.0
      HtoR    = 1.0
      sigma_model = 'hayashi'
      eps     = 1.0
      amp     = 1.e-5

      if(proc .eq. 0) then

         diff_nml(PROBLEM_CONTROL)

         cbuff(1) =  problem_name
         cbuff(2) =  run_id
         cbuff(3) =  sigma_model

         rbuff(1) = sigma0
         rbuff(2) = Rin
         rbuff(3) = R0
         rbuff(4) = HtoR
         rbuff(5) = eps
         rbuff(6) = amp

      end if

      call MPI_Bcast(cbuff, cbuff_len*buffer_dim, MPI_CHARACTER,        0, comm, ierr)
      call MPI_Bcast(rbuff,    buffer_dim, MPI_DOUBLE_PRECISION, 0, comm, ierr)

      if (proc / =0) then

         problem_name = cbuff(1)
         run_id       = cbuff(2)(1:3)
         sigma_model  = cbuff(3)

         sigma0       = rbuff(1)
         Rin          = rbuff(2)
         R0           = rbuff(3)
         HtoR         = rbuff(4)
         eps          = rbuff(5)
         amp          = rbuff(6)

      endif

   end subroutine read_problem_par

!-----------------------------------------------------------------------------

   real function dens_Rdistr(R,Rin,n)
      real, intent(in) :: R,Rin,n
      real :: ninv
      ninv = 1./n
      dens_Rdistr = (R - Rin)**ninv / R**(2.+ninv)
   end function dens_Rdistr

   subroutine init_prob
      use arrays,      only: u,gp
      use grid,        only: x,y,z,nx,ny,nz,nzd,dx,dy,dz,xmin,nb
      use constants,   only: newtong,pi,dpi
      use gravity,     only: ptmass, r_smooth, grav_pot2accel
      use initneutral, only: idnn, imxn, imyn, imzn
      use initdust,    only: idnd, imxd, imyd, imzd
      use initfluids,  only: cs_iso_neu, cs_iso_neu2
      use mpisetup,    only: smalld
#ifndef ISO
      use initneutral, only: ienn, gamma_neu
      use mpisetup,    only: smallei
#endif /* !ISO */
      implicit none

      integer :: i,j,k
      real :: xi,yj,zk, rc, H0, sqr_gm, rho0
      real :: n,norm, H, ninv
      real :: gradP, iOmega, ilook, gradgp
      real, dimension(3,nz) :: noise
      real, dimension(:), allocatable :: omega,omegad
#ifndef ISO
      real :: vx, vy, vz
#endif /* !ISO */

      allocate(omega(nx),omegad(nx))
!   Secondary parameters

      sqr_gm = sqrt(newtong*ptmass)

      write(*,*) 'dx = ',dx

      n = 0.5*Rin / (R0 - Rin)
      ninv = 1./n
      if(sigma_model == 'hayashi') then
         sigma0 = 0.2* R0**(-1.5)
      endif
      H0 = R0 * HtoR
      cs_iso_neu2 = H0 * (pi*newtong) * sigma0
      cs_iso_neu = sqrt(cs_iso_neu2)

      rho0 = sigma0 / (sqrt(dpi)*H0)

      norm = 1. / dens_Rdistr(R0,Rin,n)

      do k = 1,nz
        zk = z(k)
        do j = 1,ny
           yj = y(j)
           do i = 1,nx
              xi = x(i)
              rc = sqrt(xi**2+yj**2)
              H = HtoR * rc
              u(idnn,i,j,k) = max(rho0 * norm * dens_RdistR(rc,Rin,n) * exp(- 0.25 * zk**2 / H**2 ), smalld)
              u(idnd,i,j,k) = eps*u(idnn,i,j,k)
           enddo
        enddo
      enddo

      do i = 2,nx-1   ! 2d
         rc= x(i)*sqrt(2.0)
         gradgp=  0.5*(gp(i+1,i+1,max(nz/2,1))-gp(i-1,i-1,max(nz/2,1)))/dx/sqrt(2.)
         gradp = -0.5*(u(idnn,i+1,i+1,max(nz/2,1))-u(idnn,i-1,i-1,max(nz/2,1)))/dx &
                            /sqrt(2.)*cs_iso_neu2
         omega(i)  = sqrt( abs( (gradgp-gradp)/rc ) )
         omegad(i) = sqrt( abs(    gradgp/rc      ) )
      enddo
      omega(1)  = omega(2);  omega(nx)  = omega(nx-1)
      omegad(1) = omegad(2); omegad(nx) = omegad(nx-1)

      call random_seed()

      do j = 1,ny
         yj = y(j)
         do i = 1,nx
            xi = x(i)
            rc = sqrt(xi**2+yj**2)
            call random_number(noise)

            ilook = (rc-xmin)/dx/sqrt(2.) + 0.5 + nb
            iOmega = omega(int(ilook))+(rc-x(int(ilook))*sqrt(2.))*(omega(int(ilook)+1)-omega(int(ilook))) &
                         /(x(int(ilook)+1)-x(int(ilook)))/sqrt(2.)
!
!
            u(imxn,i,j,:) = -yj*iOmega*u(idnn,i,j,:)
            u(imyn,i,j,:) =  xi*iOmega*u(idnn,i,j,:)
            u(imzn,i,j,:) = 0.0
#ifndef ISO
            u(ienn,i,j,:) = cs_iso_neu2/(gamma_ion-1.0)*u(idnn,i,j,:)
            u(ienn,i,j,:) = max(u(ienn,i,j,:), smallei)
            u(ienn,i,j,:) = u(ienn,i,j,:) +0.5*(vx**2+vy**2+vz**2)*u(idnn,i,j,:)
#endif /* !ISO */

            iOmega = omegad(int(ilook))+(rc-x(int(ilook))*sqrt(2.))*(omegad(int(ilook)+1)-omegad(int(ilook))) &
                         /(x(int(ilook)+1)-x(int(ilook)))/sqrt(2.)
            u(imxd,i,j,:) = -yj*iOmega*u(idnd,i,j,:) + amp*(noise(1,:)-0.5)
            u(imyd,i,j,:) =  xi*iOmega*u(idnd,i,j,:) + amp*(noise(2,:)-0.5)
            u(imzd,i,j,:) = 0.0 + amp*(noise(3,:)-0.5)

         enddo
      enddo
      write(*,*) maxval(u(idnn,:,:,:)), minval(u(idnn,:,:,:))

      deallocate(omega,omegad)

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
      select case(var)
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
      select case(var)
         case default
            ierrh = -1
      end select

   end subroutine user_hdf5

end module initproblem

