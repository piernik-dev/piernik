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
!>
!! \brief ()
!!
!! In this module following namelist of parameters is specified:
!! \copydetails resistivity::init_resistivity
!<
module resistivity
! pulled by RESISTIVE
   implicit none
   private
   public  :: init_resistivity, timestep_resist, cleanup_resistivity, dt_resist, eta_max, &
      diffuseby_x, diffusebz_x, diffusebx_y, diffusebz_y, diffusebx_z, diffuseby_z
   real    :: cfl_resist, eta_0, eta_1, j_crit, jc2, deint_max
   integer :: eta_scale
   double precision :: d_eta_factor

   real :: eta_max, dt_resist, dt_eint
   integer, dimension(3) :: loc_eta_max
   real, dimension(:,:,:), allocatable, target :: wb, etahelp, eta
   logical, save :: inactive = .false.

   contains

   subroutine cleanup_resistivity
      implicit none
      if (allocated(wb) ) deallocate(wb)
      if (allocated(eta)) deallocate(eta)
      if (allocated(etahelp)) deallocate(etahelp)

   end subroutine cleanup_resistivity

!>
!! \brief Routine to set parameters values from namelist RESISTIVITY
!!
!! \n \n
!! @b RESISTIVITY
!! \n \n
!! <table border="+1">
!! <tr><td width="150pt"><b>parameter</b></td><td width="135pt"><b>default value</b></td><td width="200pt"><b>possible values</b></td><td width="315pt"> <b>description</b></td></tr>
!! <tr><td>cfl_resist</td><td>0.4  </td><td>real value   </td><td>\copydoc resistivity::cfl_resist</td></tr>
!! <tr><td>eta_0     </td><td>0.0  </td><td>real value   </td><td>\copydoc resistivity::eta_0    </td></tr>
!! <tr><td>eta_1     </td><td>0.0  </td><td>real value   </td><td>\copydoc resistivity::eta_1    </td></tr>
!! <tr><td>eta_scale </td><td>4    </td><td>integer value</td><td>\copydoc resistivity::eta_scale</td></tr>
!! <tr><td>j_crit    </td><td>1.0e6</td><td>real value   </td><td>\copydoc resistivity::j_crit   </td></tr>
!! <tr><td>deint_max </td><td>0.01 </td><td>real value   </td><td>\copydoc resistivity::deint_max</td></tr>
!! </table>
!! \n \n
!<
   subroutine init_resistivity

      use dataio_pub,    only: par_file, ierrh, namelist_errh, compare_namelist  ! QA_WARN required for diff_nml
      use dataio_pub,    only: warn, die
      use grid,          only: nx, ny, nz, has_dir, zdim, xdim, ydim
      use mpisetup,      only: rbuff, ibuff, ierr, buffer_dim, comm, master, slave
      use mpi,           only: MPI_INTEGER, MPI_DOUBLE_PRECISION

      implicit none

      namelist /RESISTIVITY/ cfl_resist, eta_0, eta_1, eta_scale, j_crit, deint_max

      cfl_resist  =  0.4
      eta_0       =  0.0
      eta_1       =  0.0
      eta_scale   =  4
      j_crit      =  1.0e6
      deint_max   =  0.01

      if (master) then

         diff_nml(RESISTIVITY)

         ibuff(1) = eta_scale

         rbuff(1) = cfl_resist
         rbuff(2) = eta_0
         rbuff(3) = eta_1
         rbuff(4) = j_crit
         rbuff(5) = deint_max

      endif

      call MPI_Bcast(ibuff,    buffer_dim, MPI_INTEGER,          0, comm, ierr)
      call MPI_Bcast(rbuff,    buffer_dim, MPI_DOUBLE_PRECISION, 0, comm, ierr)

      if (slave) then

         eta_scale          = ibuff(1)

         cfl_resist         = rbuff(1)
         eta_0              = rbuff(2)
         eta_1              = rbuff(3)
         j_crit             = rbuff(4)
         deint_max          = rbuff(5)

      endif

      if (eta_scale < 0) call die("eta_scale must be greater or equal 0")

      if (.not.allocated(wb) ) allocate(wb(nx,ny,nz) )
      if (.not.allocated(eta)) allocate(eta(nx,ny,nz))
      if (.not.allocated(etahelp)) allocate(etahelp(nx,ny,nz))

      jc2 = j_crit**2
      if (has_dir(zdim)) then
         d_eta_factor = 1./(6.+dble(eta_scale))
      else
         d_eta_factor = 1./(4.+dble(eta_scale))
      endif

      if (.not. all(has_dir(xdim:ydim))) then
         if (master) call warn("[resistivity:init_resistivity] Resistivity module needs both x and y dimension. Switching off.")
         inactive = .true.
      endif

   end subroutine init_resistivity

   subroutine compute_resist

      use arrays,       only: b, u
      use constants,    only: small
      use fluidindex,   only: ibx, iby, ibz
      use grid,         only: idl, xdim, ydim, zdim, nx, ny, nz, is, ie, js, je, ks, ke, has_dir
      use mpisetup,     only: comm, ierr
      use mpi,          only: MPI_DOUBLE_PRECISION, MPI_MAX, MPI_IN_PLACE
#ifndef ISO
      use fluidindex,   only: nvar
#endif /* !ISO */

      implicit none

!      real, dimension(nx,ny,nz), intent(inout) :: eta

      if (inactive) return
!BEWARE: uninitialized values are poisoning the wb(:,:,:) array
!BEWARE: significant differences between single-CPU run and multi-CPU run (due to uninits?)
!--- square current computing in cell corner step by step

!--- current_z
      wb(2:nx,2:ny,:) = (b(iby,2:nx,2:ny,:)-b(iby,1:nx-1,2:ny,:))*idl(xdim) - (b(ibx,2:nx,2:ny,:)-b(ibx,2:nx,1:ny-1,:))*idl(ydim)
      wb(1,:,:) = wb(2,:,:) ; wb(:,1,:) = wb(:,2,:)

      if (has_dir(zdim)) then
         eta(:,:,2:nz) = 0.25*( wb(:,:,2:nz) + wb(:,:,1:nz-1) )**2 ; eta(:,:,1) = eta(:,:,2)
      else
         eta = wb**2    ! BEWARE: is it correct?
      endif

      if (has_dir(zdim)) then
!--- current_x
         wb(:,2:ny,2:nz) = (b(ibz,:,2:ny,2:nz)-b(ibz,:,1:ny-1,2:nz))*idl(ydim) - (b(iby,:,2:ny,2:nz)-b(iby,:,2:ny,1:nz-1))*idl(zdim)

         eta(2:nx,:,:) = eta(2:nx,:,:) + 0.25*(wb(2:nx,:,:)+wb(1:nx-1,:,:))**2; eta(1,:,:) = eta(2,:,:)
!--- current_y
         wb(2:nx,:,2:nz) = (b(ibx,2:nx,:,2:nz)-b(ibx,2:nx,:,1:nz-1))*idl(zdim) - (b(ibz,2:nx,:,2:nz)-b(ibz,1:nx-1,:,2:nz))*idl(xdim)

         eta(:,2:ny,:) = eta(:,2:ny,:) + 0.25*(wb(:,2:ny,:)+wb(:,1:ny-1,:))**2; eta(:,1,:) = eta(:,2,:)
      endif

!--- wb = current**2
      wb(:,:,:) = eta(:,:,:)

      eta(:,:,:) = eta_0 + eta_1*sqrt(max(0.0,eta(:,:,:)- jc2 ))
!>
!! \todo Following lines are split into separate lines because of intel and gnu dbgs
!! shoud that be so? Is there any other solution instead splitting?
!<
      etahelp(2:nx-1,2:ny-1,:) = eta(1:nx-2,2:ny-1,:) + eta(3:nx,2:ny-1,:) + eta(2:nx-1,1:ny-2,:) + eta(2:nx-1,3:ny,:)
      etahelp(1,:,:) = etahelp(2,:,:) ; etahelp(nx,:,:) = etahelp(nx-1,:,:) ; etahelp(:,1,:) = etahelp(:,2,:) ; etahelp(:,ny,:) = etahelp(:,ny-1,:)
      if (has_dir(zdim)) then
         etahelp(:,:,2:nz-1) = etahelp(:,:,2:nz-1) + eta(:,:,1:nz-2) + eta(:,:,3:nz)
         etahelp(:,:,1) = etahelp(:,:,2) ; etahelp(:,:,nz) = etahelp(:,:,nz-1)
      endif
      etahelp = (etahelp + dble(eta_scale)*eta)*d_eta_factor

      where (eta > eta_0)
         eta = etahelp
      endwhere

      eta_max           = maxval(eta(is:ie,js:je,ks:ke))
      loc_eta_max       = maxloc(eta(is:ie,js:je,ks:ke))

      call MPI_Allreduce(MPI_IN_PLACE, eta_max, 1, MPI_DOUBLE_PRECISION, MPI_MAX, 0, comm, ierr)

#ifndef ISO
      dt_eint = deint_max * abs(minval(               &
                ( u(nvar%ion%ien,is:ie,js:je,ks:ke)           &
                - 0.5*( u(nvar%ion%imx,is:ie,js:je,ks:ke)**2  &
                      + u(nvar%ion%imy,is:ie,js:je,ks:ke)**2  &
                      + u(nvar%ion%imz,is:ie,js:je,ks:ke)**2 )&
                      /u(nvar%ion%idn,is:ie,js:je,ks:ke)      &
                - 0.5*( b(ibx,is:ie,js:je,ks:ke)**2   &
                      + b(iby,is:ie,js:je,ks:ke)**2   &
                      + b(ibz,is:ie,js:je,ks:ke)**2)) &
                      /( eta(is:ie,js:je,ks:ke)       &
                      *wb(is:ie,js:je,ks:ke)+small) ))
#endif /* !ISO */

      return

   end subroutine compute_resist

!-----------------------------------------------------------------------

   subroutine timestep_resist

      use constants, only: big
      use grid,      only: dxmn
      use mpisetup,  only: comm, ierr
      use mpi,       only: MPI_DOUBLE_PRECISION, MPI_MIN, MPI_IN_PLACE

      implicit none

      if (eta_max .ne. 0. .and. .not. inactive) then
         dt_resist = cfl_resist * dxmn**2 / (2. * eta_max)
#ifndef ISO
         dt_resist = min(dt_resist,dt_eint)
#endif /* !ISO */
      else
         dt_resist = big
      endif

      call MPI_Allreduce(MPI_IN_PLACE, dt_resist, 1, MPI_DOUBLE_PRECISION, MPI_MIN, 0, comm, ierr)

   end subroutine timestep_resist

!-----------------------------------------------------------------------------

   subroutine vanleer_limiter(f,a,b) ! ToDo: overload me or use class(*) if you dare

      implicit none

      real, dimension(:), intent(in)      :: a !< second order correction of left- or right- moving waves flux on the left cell boundary
      real, dimension(:), intent(in)      :: b !< second order correction of left- or right- moving waves flux on the right cell boundary
      real, dimension(:), intent(inout)   :: f !< second order flux correction for left- or right- moving waves
      ! locals
      real, dimension(size(a,1)) :: c !< a*b

      c = a*b                                                                    ! ToDO: OPTIMIZE ME
      where (c > 0.0)
         f = f+2.0*c/(a+b)
      endwhere

   end subroutine vanleer_limiter

   subroutine tvdd_1d(b1d,eta1d,idi,dt,wcu1d)
!      use fluxes, only: flimiter

      implicit none

      real, dimension(:), pointer, intent(in)    :: eta1d, b1d
      real, dimension(:), intent(out)            :: wcu1d
      real, intent(in)                           :: idi,dt

      real, dimension(size(b1d))                 :: w, wp, wm, b1
      integer                                    :: n

      if (inactive) return

      n = size(b1d)
      w(2:n)    = eta1d(2:n) * ( b1d(2:n) - b1d(1:n-1) )*idi ; w(1)  = w(2)
      b1(1:n-1) = b1d(1:n-1) + 0.5*(w(2:n) - w(1:n-1))*dt*idi; b1(n) = b1(n-1)

      w(2:n)    = eta1d(2:n) * ( b1(2:n) - b1(1:n-1) )*idi   ; w(1)  = w(2)
      wp(1:n-1) = 0.5*(w(2:n) - w(1:n-1))                    ; wp(n) = wp(n-1)
      wm(2:n)   = wp(1:n-1)                                  ; wm(1) = wm(2)

      call vanleer_limiter(w,wm,wp)
      wcu1d     = w*dt

   end subroutine tvdd_1d

!-------------------------------------------------------------------------------

   subroutine diffuseby_x

      use arrays,        only: wcu, b
      use mpisetup,      only: dt
      use fluidindex,    only: iby
      use grid,          only: has_dir, xdim, ydim, zdim, nx, ny, nz, idl
      use magboundaries, only: bnd_emf

      implicit none

      real, dimension(:), pointer :: b1d, eta1d
      real, dimension(nx)         :: wcu1d
      integer                     :: j, k

      if (inactive) return

      call compute_resist
      eta(:,:,1:nz-1) = 0.5*(eta(:,:,1:nz-1)+eta(:,:,2:nz))

      do j = 1, ny
         do k = 1, nz
            b1d    => b(iby,:,j,k)
            eta1d  => eta(:,j,k)
            call tvdd_1d(b1d, eta1d, idl(xdim), dt, wcu1d)
            wcu(:,j,k) = wcu1d
         enddo
      enddo

      if (has_dir(xdim)) call bnd_emf(wcu,'emfz','xdim')
      if (has_dir(ydim)) call bnd_emf(wcu,'emfz','ydim')
      if (has_dir(zdim)) call bnd_emf(wcu,'emfz','zdim')

   end subroutine diffuseby_x

   subroutine diffusebz_x

      use arrays,        only: wcu, b
      use mpisetup,      only: dt
      use fluidindex,    only: ibz
      use grid,          only: has_dir, xdim, ydim, zdim, nx, ny, nz, idl
      use magboundaries, only: bnd_emf

      implicit none

      real, dimension(:), pointer :: b1d, eta1d
      real, dimension(nx)         :: wcu1d
      integer                     :: j, k

      if (inactive) return

      call compute_resist
      eta(:,1:ny-1,:) = 0.5*(eta(:,1:ny-1,:)+eta(:,2:ny,:))

      do j = 1, ny
         do k = 1, nz
            b1d    => b(ibz,:,j,k)
            eta1d  => eta(:,j,k)
            call tvdd_1d(b1d, eta1d, idl(xdim), dt, wcu1d)
            wcu(:,j,k) = wcu1d
         enddo
      enddo
      if (has_dir(xdim)) call bnd_emf(wcu,'emfy','xdim')
      if (has_dir(ydim)) call bnd_emf(wcu,'emfy','ydim')
      if (has_dir(zdim)) call bnd_emf(wcu,'emfy','zdim')

   end subroutine diffusebz_x

   subroutine diffusebz_y

      use arrays,        only: wcu, b
      use mpisetup,      only: dt
      use fluidindex,    only: ibz
      use grid,          only: has_dir, xdim, ydim, zdim, nx, ny, nz, idl
      use magboundaries, only: bnd_emf

      implicit none

      real, dimension(:), pointer :: b1d, eta1d
      real, dimension(ny)         :: wcu1d
      integer                     :: i, k

      if (inactive) return

      call compute_resist
      eta(1:nx-1,:,:) = 0.5*(eta(1:nx-1,:,:)+eta(2:nx,:,:))

      do i = 1, nx
         do k = 1, nz
            b1d    => b(ibz,i,:,k)
            eta1d  => eta(i,:,k)
            call tvdd_1d(b1d, eta1d, idl(ydim), dt, wcu1d)
            wcu(i,:,k) = wcu1d
         enddo
      enddo
      if (has_dir(ydim)) call bnd_emf(wcu,'emfx','ydim')
      if (has_dir(zdim)) call bnd_emf(wcu,'emfx','zdim')
      if (has_dir(xdim)) call bnd_emf(wcu,'emfx','xdim')

   end subroutine diffusebz_y

   subroutine diffusebx_y

      use arrays,        only: wcu, b
      use mpisetup,      only: dt
      use fluidindex,    only: ibx
      use grid,          only: has_dir, xdim, ydim, zdim, nx, ny, nz, idl
      use magboundaries, only: bnd_emf

      implicit none

      real, dimension(:), pointer :: b1d, eta1d
      real, dimension(ny)         :: wcu1d
      integer                     :: i, k

      if (inactive) return

      call compute_resist
      eta(:,:,1:nz-1) = 0.5*(eta(:,:,1:nz-1)+eta(:,:,2:nz))

      do i = 1, nx
         do k = 1, nz
            b1d    => b(ibx,i,:,k)
            eta1d  => eta(i,:,k)
            call tvdd_1d(b1d, eta1d, idl(ydim), dt, wcu1d)
            wcu(i,:,k) = wcu1d
         enddo
      enddo
      if (has_dir(ydim)) call bnd_emf(wcu, 'emfz', 'ydim')
      if (has_dir(zdim)) call bnd_emf(wcu, 'emfz', 'zdim')
      if (has_dir(xdim)) call bnd_emf(wcu, 'emfz', 'xdim')

   end subroutine diffusebx_y

   subroutine diffusebx_z

      use arrays,        only: wcu, b
      use mpisetup,      only: dt
      use fluidindex,    only: ibx
      use grid,          only: has_dir, xdim, ydim, zdim, nx, ny, nz, idl
      use magboundaries, only: bnd_emf

      implicit none

      real, dimension(:), pointer :: b1d, eta1d
      real, dimension(nz)         :: wcu1d
      integer                     :: i, j

      if (inactive) return

      call compute_resist
      eta(:,1:ny-1,:) = 0.5*(eta(:,1:ny-1,:)+eta(:,2:ny,:))

      do i = 1, nx
         do j = 1, ny
            b1d    => b(ibx,i,j,:)
            eta1d  => eta(i,j,:)
            call tvdd_1d(b1d, eta1d, idl(zdim), dt, wcu1d)
            wcu(i,j,:) = wcu1d
         enddo
      enddo
      if (has_dir(zdim)) call bnd_emf(wcu, 'emfy', 'zdim')
      if (has_dir(xdim)) call bnd_emf(wcu, 'emfy', 'xdim')
      if (has_dir(ydim)) call bnd_emf(wcu, 'emfy', 'ydim')

   end subroutine diffusebx_z

   subroutine diffuseby_z

      use arrays,        only: wcu, b
      use mpisetup,      only: dt
      use fluidindex,    only: iby
      use grid,          only: has_dir, xdim, ydim, zdim, nx, ny, nz, idl
      use magboundaries, only: bnd_emf

      implicit none

      real, dimension(:), pointer :: b1d, eta1d
      real, dimension(nz)         :: wcu1d
      integer                     :: i, j

      if (inactive) return

      call compute_resist
      eta(1:nx-1,:,:) = 0.5*(eta(1:nx-1,:,:)+eta(2:nx,:,:))

      do i = 1, nx
         do j = 1, ny
            b1d    => b(iby,i,j,:)
            eta1d  => eta(i,j,:)
            call tvdd_1d(b1d, eta1d, idl(zdim), dt, wcu1d)
            wcu(i,j,:) = wcu1d
         enddo
      enddo
      if (has_dir(zdim)) call bnd_emf(wcu, 'emfx', 'zdim')
      if (has_dir(xdim)) call bnd_emf(wcu, 'emfx', 'xdim')
      if (has_dir(ydim)) call bnd_emf(wcu, 'emfx', 'ydim')

   end subroutine diffuseby_z

end module resistivity
