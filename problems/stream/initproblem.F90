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

   use mpisetup,    only: cbuff_len
   implicit none

   private
   public :: read_problem_par, init_prob

   character(len=cbuff_len) :: fnoise
   real :: rhog, eps, amp, kx, kz
   real, dimension(8), save :: vec
   logical :: linear

   namelist /PROBLEM_CONTROL/  &
                               rhog, eps, amp, fnoise, kx,kz, linear

   contains

!-----------------------------------------------------------------------------

   subroutine read_problem_par

      use dataio_pub,    only: ierrh, par_file, namelist_errh, compare_namelist, cmdl_nml      ! QA_WARN required for diff_nml
      use mpisetup,      only: rbuff, cbuff, lbuff, buffer_dim, master, slave, comm, ierr
      use mpi,           only: MPI_CHARACTER, MPI_DOUBLE_PRECISION, MPI_LOGICAL

      implicit none

      linear = .false.
      rhog    = 10.0
      eps     =  1.0
      amp    =  0.0
      kx     =  30.0
      kz     =  30.0

      if (master) then

         diff_nml(PROBLEM_CONTROL)

         cbuff(1) =  fnoise

         rbuff(1) = rhog
         rbuff(2) = eps
         rbuff(3) = amp
         rbuff(4) = kx
         rbuff(5) = kz

         lbuff(1) = linear

      endif

      call MPI_Bcast(cbuff, cbuff_len*buffer_dim, MPI_CHARACTER,        0, comm, ierr)
      call MPI_Bcast(rbuff,           buffer_dim, MPI_DOUBLE_PRECISION, 0, comm, ierr)
      call MPI_Bcast(lbuff,           buffer_dim, MPI_LOGICAL,          0, comm, ierr)

      if (slave) then

         fnoise       = cbuff(1)

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
      use arrays,        only: u
      use constants,     only: pi, dpi
      use dataio_pub,    only: msg, printinfo, run_id
      use fluidindex,    only: flind
      use grid,          only: cg
      use mpisetup,      only: proc, dom
      use shear,         only: omega
      use types,         only: component_fluid
      use interactions,  only: dragc_gas_dust
#ifdef SHEAR
      use shear,         only: eta_gas, csvk
#endif /* SHEAR */
      implicit none

      real                                :: rcx, rcy, ux, uy, wx, wy, taus, eta, vk, beta !, inv
      integer                             :: i, j, k, n, clock
      real(kind=4), dimension(3, cg%nx, cg%ny, cg%nz) :: noise
      integer, dimension(:), allocatable  :: seed
      complex(kind=8), dimension(7)       :: coeff
      type(component_fluid), pointer      :: dst, neu

#ifdef DUST
      dst => flind%dst
#else /* !DUST */
      call warn("[initproblem]: Dust fluid not initialized. I hope you know what you are doing!"
#endif /* !DUST */
#ifdef NEUTRAL
      neu => flind%neu
#else /* !NEUTRAL */
      call warn("[initproblem]: Neutral fluid not initialized. I hope you know what you are doing!"
#endif /* !NEUTRAL */

      if (run_id == 'lnA') then
         call printinfo("Lin A")
         coeff(4) = (-0.1691398, 0.0361553 ) ! u_x
         coeff(5) = ( 0.1336704, 0.0591695 ) ! u_y
         coeff(6) = ( 0.1691389,-0.0361555 ) ! u_z
         coeff(7) = ( 0.0000224,+0.0000212 ) ! gas dens
         coeff(1) = (-0.1398623, 0.0372951 ) ! w_x
         coeff(2) = ( 0.1305628, 0.0640574 ) ! w_y
         coeff(3) = ( 0.1639549,-0.0233277 ) ! w_z
         kx = dpi/dom%Lx
         kz = dpi/dom%Lz
      else if (run_id == 'lnB') then
         call printinfo("Lin B")
         coeff(4) = (-0.0174121,-0.2770347 ) ! u_x
         coeff(5) = ( 0.2767976,-0.0187568 ) ! u_y
         coeff(6) = ( 0.0174130, 0.2770423 ) ! u_z
         coeff(7) = (-0.0000067,-0.0000691 ) ! gas dens
         coeff(1) = ( 0.0462916,-0.2743072 ) ! w_x
         coeff(2) = ( 0.2739304, 0.0039293 ) ! w_y
         coeff(3) = ( 0.0083263, 0.2768866 ) ! w_z
         kx = dpi/dom%Lx
         kz = dpi/dom%Lz
      else
         kx = dpi/dom%Lx
         kz = dpi/dom%Lz
         coeff(:) = ( 0.0, 0.0 )
      endif

      call random_seed(size=n)
      allocate(seed(n))
      call system_clock(count=clock)
      seed = clock*proc + 37 * (/ (i-1, i = 1, n) /)
      call random_seed(put=seed)
      deallocate(seed)

      taus = 1./dragc_gas_dust
      vk   = neu%cs/csvk
      eta  = eta_gas
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

      write(msg,*) 'kx = ',kx,' kz = ',kz
      call printinfo(msg)
      write(msg,*) 'ux = ',ux,' uy = ',uy
      call printinfo(msg)
      write(msg,*) 'wx = ',wx,' wy = ',wy
      call printinfo(msg)
      write(msg,*) '\eta vk / \Omega = ', eta_gas * neu%cs / csvk / omega
      call printinfo(msg)

      do i = 1, cg%nx
         rcx = cg%x(i)
         do j = 1, cg%ny
            rcy = cg%y(j)
            do k = 1, cg%nz
#ifdef NEUTRAL
               u(neu%idn,i,j,k) = rhog
               u(neu%imx,i,j,k) = ux * rhog
               u(neu%imy,i,j,k) = uy * rhog
               u(neu%imz,i,j,k) = 0.0
#ifndef ISO
               u(neu%ien,i,j,:) = 1.0/(gamma_neu-1.0)
#endif /* !ISO */
#endif /* NEUTRAL */
#ifdef DUST
               u(dst%idn,i,j,k) = eps*rhog
               u(dst%imx,i,j,k) = wx * eps*rhog
               u(dst%imy,i,j,k) = wy * eps*rhog
               u(dst%imz,i,j,k) = 0.0

! Linear test
            if (linear) then
!               u(dst%idn,i,j,k) =  u(dst%idn,i,j,k) + amp*eps*dsin(kz*cg%z(k))*dcos(kx*cg%x(i))
               u(dst%idn,i,j,k) =  u(dst%idn,i,j,k) + amp*eps*dcos(kz*cg%z(k))*dcos(kx*cg%x(i))
! ...                u(dst%idn,i,j,k) =  u(dst%idn,i,j,k) + amp*eps*dcos(kx*x(i))*dcos(kz*cg%z(k))
! B              u(dst%idn,i,j,k) =  u(dst%idn,i,j,k) + amp * eps *&
! B                 ( real(coeff(7))*dcos(kx*cg%x(i)) - &
! B                  aimag(coeff(7))*dsin(kx*cg%x(i))) * dcos(kz*z(k))
               u(dst%imx,i,j,k) =  u(dst%imx,i,j,k) + eta*vk*amp * &
                  ( real(coeff(1))*dcos(kx*cg%x(i)) - &
                   aimag(coeff(1))*dsin(kx*cg%x(i))) * dcos(kz*cg%z(k))
               u(dst%imy,i,j,k) =  u(dst%imy,i,j,k) + eta*vk*amp * &
                  ( real(coeff(2))*dcos(kx*cg%x(i)) - &
                   aimag(coeff(2))*dsin(kx*cg%x(i))) * dcos(kz*cg%z(k))
               u(dst%imz,i,j,k) =  u(dst%imz,i,j,k) + eta*vk*(-amp) * &
                  (aimag(coeff(3))*dcos(kx*cg%x(i)) + &
                    real(coeff(3))*dsin(kx*cg%x(i))) * dsin(kz*cg%z(k))
               u(neu%imx,i,j,k) =  u(neu%imx,i,j,k) + eta*vk*amp * &
                  ( real(coeff(4))*dcos(kx*cg%x(i)) - &
                   aimag(coeff(4))*dsin(kx*cg%x(i))) * dcos(kz*cg%z(k))
               u(neu%imy,i,j,k) =  u(neu%imy,i,j,k) + eta*vk*amp * &
                  ( real(coeff(5))*dcos(kx*cg%x(i)) - &
                   aimag(coeff(5))*dsin(kx*cg%x(i))) * dcos(kz*cg%z(k))
               u(neu%imz,i,j,k) =  u(neu%imz,i,j,k) + eta*vk*(-amp) * &
                  (aimag(coeff(6))*dcos(kx*cg%x(i)) + &
                    real(coeff(6))*dsin(kx*cg%x(i))) * dsin(kz*cg%z(k))
!               u(neu%idn,i,j,k) =  u(neu%idn,i,j,k) + amp * &
!                  ( real(coeff(7))*dcos(kx*cg%x(i)) - &
!                   aimag(coeff(7))*dsin(kx*cg%x(i))) * dcos(kz*cg%z(k))
               u(neu%idn,i,j,k) =  u(neu%idn,i,j,k) + (eta*vk)**2 * amp * &
                  ( real(coeff(7))*dcos(kx*cg%x(i)) - &
                   aimag(coeff(7))*dsin(kx*cg%x(i))) * dcos(kz*cg%z(k))
             endif
!-------

#endif /* DUST */
            enddo
         enddo
      enddo
      if (.not.linear) then
        call random_number(noise)
        u(dst%imx,:,:,:) = u(dst%imx,:,:,:) +amp -2.0*amp*noise(1,:,:,:) * u(dst%idn,:,:,:)
        u(dst%imy,:,:,:) = u(dst%imy,:,:,:) +amp -2.0*amp*noise(2,:,:,:) * u(dst%idn,:,:,:)
        u(dst%imz,:,:,:) = u(dst%imz,:,:,:) +amp -2.0*amp*noise(3,:,:,:) * u(dst%idn,:,:,:)
      endif

      write(msg,*) 'linear = ',linear
      call printinfo(msg)

      return
   end subroutine init_prob

!-----------------------------------------------------------------------------------------------------------------------------------

   subroutine compare(t)

      use constants,    only: dpi
      use arrays,       only: u
      use dataio_pub,   only: run_id
      use mpisetup,     only: dom

      implicit none

      real, intent(in) :: t
      complex(kind=8), dimension(8) :: coeff

      if (t==0) coeff(4) = (0,0) ! suppress compiler warnings

      if (run_id == 'lnA') then
         coeff(4) = (-0.1691398, 0.0361553 ) ! u_x
         coeff(5) = ( 0.1336704, 0.0591695 ) ! u_y
         coeff(6) = ( 0.1691389,-0.0361555 ) ! u_z
         coeff(7) = ( 0.0000224,+0.0000212 ) ! gas dens
         coeff(1) = (-0.1398623, 0.0372951 ) ! w_x
         coeff(2) = ( 0.1305628, 0.0640574 ) ! w_y
         coeff(3) = ( 0.1639549,-0.0233277 ) ! w_z

         coeff(8) = (-0.3480127, 0.4190204 ) ! omega
      else if (run_id == 'lnB') then
         coeff(4) = (-0.0174121,-0.2770347 ) ! u_x
         coeff(5) = ( 0.2767976,-0.0187568 ) ! u_y
         coeff(6) = ( 0.0174130, 0.2770423 ) ! u_z
         coeff(7) = (-0.0000067,-0.0000691 ) ! gas dens
         coeff(1) = ( 0.0462916,-0.2743072 ) ! w_x
         coeff(2) = ( 0.2739304, 0.0039293 ) ! w_y
         coeff(3) = ( 0.0083263, 0.2768866 ) ! w_z

         coeff(8) = ( 0.4998786, 0.0154764 ) ! omega
      endif

      kx = dpi/dom%Lx
      kz = dpi/dom%Lz

   end subroutine compare

end module initproblem
