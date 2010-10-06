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

   use problem_pub, only: problem_name, run_id

   character(len=32) :: fnoise
   real :: rhog, eps, amp, kx, kz
   real, dimension(8), save :: vec
   logical :: linear

   namelist /PROBLEM_CONTROL/  problem_name, run_id, &
                               rhog, eps, amp, fnoise, kx,kz, linear

   contains

!-----------------------------------------------------------------------------

   subroutine read_problem_par

      use errh,     only : namelist_errh
      use mpisetup, only : ierr, rbuff, cbuff, ibuff, lbuff, proc, buffer_dim, comm, &
           &               MPI_CHARACTER, MPI_DOUBLE_PRECISION, MPI_INTEGER, MPI_LOGICAL
      use dataio_public, only : cwd, msg, par_file
      use func,          only : compare_namelist

      implicit none

      integer :: ierrh

      linear = .false.
      problem_name = 'aaa'
      run_id  = 'aa'
      rhog    = 10.0
      eps     =  1.0
      amp    =  0.0
      kx     =  30.0
      kz     =  30.0

      if(proc .eq. 0) then

         diff_nml(PROBLEM_CONTROL)

         cbuff(1) =  problem_name
         cbuff(2) =  run_id
         cbuff(3) =  fnoise

         rbuff(1) = rhog
         rbuff(2) = eps
         rbuff(3) = amp
         rbuff(4) = kx
         rbuff(5) = kz

         lbuff(1) = linear

         call MPI_BCAST(cbuff, 32*buffer_dim, MPI_CHARACTER,        0, comm, ierr)
         call MPI_BCAST(ibuff,    buffer_dim, MPI_INTEGER,          0, comm, ierr)
         call MPI_BCAST(rbuff,    buffer_dim, MPI_DOUBLE_PRECISION, 0, comm, ierr)
         call MPI_BCAST(lbuff,    buffer_dim, MPI_LOGICAL,          0, comm, ierr)

      else

         call MPI_BCAST(cbuff, 32*buffer_dim, MPI_CHARACTER,        0, comm, ierr)
         call MPI_BCAST(ibuff,    buffer_dim, MPI_INTEGER,          0, comm, ierr)
         call MPI_BCAST(rbuff,    buffer_dim, MPI_DOUBLE_PRECISION, 0, comm, ierr)
         call MPI_BCAST(lbuff,    buffer_dim, MPI_LOGICAL,          0, comm, ierr)

         problem_name = cbuff(1)
         run_id       = cbuff(2)(1:3)
         fnoise       = cbuff(3)

         rhog         = rbuff(1)
         eps          = rbuff(2)
         amp          = rbuff(3)
         kx           = rbuff(4)
         kz           = rbuff(5)

         linear       = lbuff(1)

      endif

   end subroutine read_problem_par

!-----------------------------------------------------------------------------

   subroutine init_prob
      use mpisetup,    only : proc,pcoords
      use arrays,       only : u
      use grid,         only : x,y,z,nx,ny,nz,nzd,ymin,ymax,Lx,Lz
#ifdef NEUTRAL
      use initneutral,  only : gamma_neu,idnn,imxn,imyn,imzn, cs_iso_neu, eta_gas_neu, csvk
#endif /* NEUTRAL */
#ifdef DUST
      use initdust,     only : idnd,imxd,imyd,imzd, dragc_gas_dust
#endif /* DUST */
      use constants,    only : pi,dpi
      use shear,        only : omega, qshear
#ifndef ISO
      use initneutral,  only : ienn
#endif /* !ISO */
      implicit none

      real :: rcx, rcy
      real :: ux,uy,wx,wy,taus,eta,vk,beta !, inv
      integer :: i, j, k,n, clock
      real(kind=4), dimension(3,nx,ny,nz) :: noise
      integer, dimension(:), allocatable :: seed
      complex(kind=8), dimension(7) :: coeff
!      character(len=32) :: ala

      if(run_id == 'lnA') then
         coeff(4) = (-0.1691398, 0.0361553 ) ! u_x
         coeff(5) = ( 0.1336704, 0.0591695 ) ! u_y
         coeff(6) = ( 0.1691389,-0.0361555 ) ! u_z
         coeff(7) = ( 0.0000224,+0.0000212 ) ! gas dens
         coeff(1) = (-0.1398623, 0.0372951 ) ! w_x
         coeff(2) = ( 0.1305628, 0.0640574 ) ! w_y
         coeff(3) = ( 0.1639549,-0.0233277 ) ! w_z
         kx = dpi/Lx
         kz = dpi/Lz
      else if(run_id == 'lnB') then
         write(*,*) 'Lin B'
         coeff(4) = (-0.0174121,-0.2770347 ) ! u_x
         coeff(5) = ( 0.2767976,-0.0187568 ) ! u_y
         coeff(6) = ( 0.0174130, 0.2770423 ) ! u_z
         coeff(7) = (-0.0000067,-0.0000691 ) ! gas dens
         coeff(1) = ( 0.0462916,-0.2743072 ) ! w_x
         coeff(2) = ( 0.2739304, 0.0039293 ) ! w_y
         coeff(3) = ( 0.0083263, 0.2768866 ) ! w_z
         kx = dpi/Lx
         kz = dpi/Lz
      else
         kx = dpi/Lx
         kz = dpi/Lz
         coeff(:) = ( 0.0, 0.0 )
      endif

      call random_seed(size=n)
      allocate(seed(n))
      call system_clock(count=clock)
      seed = clock*proc + 37 * (/ (i-1, i = 1, n) /)
      call random_seed(put=seed)
      deallocate(seed)

!     write(ala,'(A8,I1,I1)') trim(fnoise),pcoords(1),pcoords(3)
!     write(*,*) trim(ala)
!     open(1,file=trim(fnoise),form='unformatted')
!     read(1) noise
!     close(1)


      taus = 1./dragc_gas_dust
      vk   = cs_iso_neu/csvk
      eta  = eta_gas_neu
!      beta = 2.0*omega*eta*vk
!      inv  = 1./(4.0*omega**2*taus**2 + (eps+1.0)**2)
!      ux =  beta*eps*taus * inv
!      wx = -beta*    taus * inv
!      uy = -0.5*beta/omega * (4.0*omega**2*taus**2 + eps + 1.0) * inv
!      wy = -0.5*beta/omega * (                       eps + 1.0) * inv

      beta = 2.0*taus*eta*vk
      wx  = -beta / ( (1.0+eps)**2 + taus**2 )
      ux  = -eps*wx
      wy  = (1+eps)/(2.0*taus) * wx
      uy  = (1.0 + eps + taus**2)/ (2.0*taus) * wx


      write(*,*) kx,kz
      write(*,*) ux,uy,wx,wy
      write(*,*) eta_gas_neu * cs_iso_neu / csvk / omega

      do i = 1,nx
         rcx = x(i)
         do j = 1,ny
            rcy = y(j)
            do k = 1,nz
#ifdef NEUTRAL
               u(idnn,i,j,k) = rhog
               u(imxn,i,j,k) = ux * rhog
               u(imyn,i,j,k) = uy * rhog
               u(imzn,i,j,k) = 0.0
#ifndef ISO
               u(ienn,i,j,:) = 1.0/(gamma_neu-1.0)
#endif /* !ISO */
#endif /* NEUTRAL */
#ifdef DUST
               u(idnd,i,j,k) = eps*rhog
               u(imxd,i,j,k) = wx * eps*rhog
               u(imyd,i,j,k) = wy * eps*rhog
               u(imzd,i,j,k) = 0.0

! Linear test
            if(linear) then
!               u(idnd,i,j,k) =  u(idnd,i,j,k) + amp*eps*dsin(kz*z(k))*dcos(kx*x(i))
               u(idnd,i,j,k) =  u(idnd,i,j,k) + amp*eps*dcos(kz*z(k))*dcos(kx*x(i))
! ...                u(idnd,i,j,k) =  u(idnd,i,j,k) + amp*eps*dcos(kx*x(i))*dcos(kz*z(k))
! B              u(idnd,i,j,k) =  u(idnd,i,j,k) + amp * eps *&
! B                 ( real(coeff(7))*dcos(kx*x(i)) - &
! B                  aimag(coeff(7))*dsin(kx*x(i))) * dcos(kz*z(k))
               u(imxd,i,j,k) =  u(imxd,i,j,k) + eta*vk*amp * &
                  ( real(coeff(1))*dcos(kx*x(i)) - &
                   aimag(coeff(1))*dsin(kx*x(i))) * dcos(kz*z(k))
               u(imyd,i,j,k) =  u(imyd,i,j,k) + eta*vk*amp * &
                  ( real(coeff(2))*dcos(kx*x(i)) - &
                   aimag(coeff(2))*dsin(kx*x(i))) * dcos(kz*z(k))
               u(imzd,i,j,k) =  u(imzd,i,j,k) + eta*vk*(-amp) * &
                  (aimag(coeff(3))*dcos(kx*x(i)) + &
                    real(coeff(3))*dsin(kx*x(i))) * dsin(kz*z(k))
               u(imxn,i,j,k) =  u(imxn,i,j,k) + eta*vk*amp * &
                  ( real(coeff(4))*dcos(kx*x(i)) - &
                   aimag(coeff(4))*dsin(kx*x(i))) * dcos(kz*z(k))
               u(imyn,i,j,k) =  u(imyn,i,j,k) + eta*vk*amp * &
                  ( real(coeff(5))*dcos(kx*x(i)) - &
                   aimag(coeff(5))*dsin(kx*x(i))) * dcos(kz*z(k))
               u(imzn,i,j,k) =  u(imzn,i,j,k) + eta*vk*(-amp) * &
                  (aimag(coeff(6))*dcos(kx*x(i)) + &
                    real(coeff(6))*dsin(kx*x(i))) * dsin(kz*z(k))
!               u(idnn,i,j,k) =  u(idnn,i,j,k) + amp * &
!                  ( real(coeff(7))*dcos(kx*x(i)) - &
!                   aimag(coeff(7))*dsin(kx*x(i))) * dcos(kz*z(k))
               u(idnn,i,j,k) =  u(idnn,i,j,k) + (eta*vk)**2 * amp * &
                  ( real(coeff(7))*dcos(kx*x(i)) - &
                   aimag(coeff(7))*dsin(kx*x(i))) * dcos(kz*z(k))
             endif
!-------

#endif /* DUST */
            enddo
         enddo
      enddo
      if(.not.linear) then
        call random_number(noise)
        u(imxd,:,:,:) = u(imxd,:,:,:) +amp -2.0*amp*noise(1,:,:,:) * u(idnd,:,:,:)
        u(imyd,:,:,:) = u(imyd,:,:,:) +amp -2.0*amp*noise(2,:,:,:) * u(idnd,:,:,:)
        u(imzd,:,:,:) = u(imzd,:,:,:) +amp -2.0*amp*noise(3,:,:,:) * u(idnd,:,:,:)
      endif

      write(*,*) 'linear = ',linear

      return
   end subroutine init_prob

!-----------------------------------------------------------------------------------------------------------------------------------

   subroutine compare(t)
      use constants,    only : dpi
      use arrays,       only : u
      use grid,         only : x,y,z,Lx,Lz
      implicit none
      real, intent(in) :: t
      complex(kind=8), dimension(8) :: coeff

      if(run_id == 'lnA') then
         coeff(4) = (-0.1691398, 0.0361553 ) ! u_x
         coeff(5) = ( 0.1336704, 0.0591695 ) ! u_y
         coeff(6) = ( 0.1691389,-0.0361555 ) ! u_z
         coeff(7) = ( 0.0000224,+0.0000212 ) ! gas dens
         coeff(1) = (-0.1398623, 0.0372951 ) ! w_x
         coeff(2) = ( 0.1305628, 0.0640574 ) ! w_y
         coeff(3) = ( 0.1639549,-0.0233277 ) ! w_z

         coeff(8) = (-0.3480127, 0.4190204 ) ! omega
      else if(run_id == 'lnB') then
         write(*,*) 'Lin B'
         coeff(4) = (-0.0174121,-0.2770347 ) ! u_x
         coeff(5) = ( 0.2767976,-0.0187568 ) ! u_y
         coeff(6) = ( 0.0174130, 0.2770423 ) ! u_z
         coeff(7) = (-0.0000067,-0.0000691 ) ! gas dens
         coeff(1) = ( 0.0462916,-0.2743072 ) ! w_x
         coeff(2) = ( 0.2739304, 0.0039293 ) ! w_y
         coeff(3) = ( 0.0083263, 0.2768866 ) ! w_z

         coeff(8) = ( 0.4998786, 0.0154764 ) ! omega
      endif
      kx = dpi/Lx
      kz = dpi/Lz

   end subroutine compare

end module initproblem
