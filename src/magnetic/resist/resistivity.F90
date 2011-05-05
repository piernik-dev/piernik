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
!! \brief Module of routines that correspond to resistivity
!!
!! In this module following namelist of parameters is specified:
!! \copydetails resistivity::init_resistivity
!<
module resistivity
! pulled by RESISTIVE
   use types, only: value
   implicit none

   private
   public  :: init_resistivity, timestep_resist, cleanup_resistivity, dt_resist, etamax,   &
        &     diffuseby_x, diffusebz_x, diffusebx_y, diffusebz_y, diffusebx_z, diffuseby_z, &
        &     cu2max, deimin, eta1_active

   real    :: cfl_resist                     !< CFL factor for resistivity effect
   real    :: eta_0                          !< uniform resistivity
   real    :: eta_1                          !< anomalous resistivity
   real    :: j_crit                         !< critical value of current density
   real    :: jc2                            !< squared critical value of current density
   real    :: deint_max                      !< COMMENT ME
   integer :: eta_scale                      !< COMMENT ME
   real    :: dt_resist, dt_eint
   double precision :: d_eta_factor
   type(value) :: etamax, cu2max, deimin
   real, dimension(:,:,:), allocatable, target :: wb, eh, eta
   real, dimension(:,:,:), allocatable         :: dbx, dby, dbz
   logical, save :: eta1_active = .true.       !< resistivity off-switcher while eta_1 == 0.0

contains

   subroutine cleanup_resistivity

      implicit none

      if (allocated(eta)) deallocate(eta)
      if (allocated(wb) ) deallocate(wb)
      if (allocated(eh) ) deallocate(eh)
      if (allocated(dbx)) deallocate(dbx)
      if (allocated(dby)) deallocate(dby)
      if (allocated(dbz)) deallocate(dbz)

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

      use dataio_pub,    only: par_file, ierrh, namelist_errh, compare_namelist, cmdl_nml  ! QA_WARN required for diff_nml
      use dataio_pub,    only: die, code_progress
      use constants,     only: PIERNIK_INIT_BASE, zdim, xdim, ydim
      use grid,          only: cg
      use mpisetup,      only: rbuff, ibuff, ierr, comm, master, slave, has_dir, buffer_dim
      use mpi,           only: MPI_INTEGER, MPI_DOUBLE_PRECISION

      implicit none
      real :: dims_twice

      namelist /RESISTIVITY/ cfl_resist, eta_0, eta_1, eta_scale, j_crit, deint_max

      if (code_progress < PIERNIK_INIT_BASE) call die("[arrays:init_arrays] grid not initialized.")

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

      if (.not.allocated(eta)) allocate(eta(cg%nx, cg%ny, cg%nz))
#ifdef ISO
      if (eta_1 == 0.) then
         eta = eta_0
         eta1_active = .false.
      endif
#endif /* ISO */

      if (eta1_active) then
         if (.not.allocated(wb) ) allocate( wb(cg%nx, cg%ny, cg%nz))
         if (.not.allocated(eh) ) allocate( eh(cg%nx, cg%ny, cg%nz))
         if (.not.allocated(dbx)) allocate(dbx(cg%nx, cg%ny, cg%nz))
         if (.not.allocated(dby)) allocate(dby(cg%nx, cg%ny, cg%nz))
         if (.not.allocated(dbz)) allocate(dbz(cg%nx, cg%ny, cg%nz))

         if (.not.has_dir(xdim)) dbx = 0.0
         if (.not.has_dir(ydim)) dby = 0.0
         if (.not.has_dir(zdim)) dbz = 0.0

         jc2 = j_crit**2
         dims_twice = 2. * count(has_dir)
         d_eta_factor = 1./(dims_twice+dble(eta_scale))
      endif

   end subroutine init_resistivity

   subroutine compute_resist

      use arrays,       only: b, u
      use constants,    only: small, xdim, ydim, zdim, MINL, MAXL
      use fluidindex,   only: ibx, iby, ibz
      use grid,         only: cg
      use mpisetup,     only: comm, ierr, has_dir, mpifind
      use mpi,          only: MPI_DOUBLE_PRECISION
#ifndef ISO
      use fluidindex,   only: flind
#endif /* !ISO */

      implicit none

      if (.not.eta1_active) return
!> \deprecated BEWARE: uninitialized values are poisoning the wb(:,:,:) array - should change  with rev. 3893
!> \deprecated BEWARE: significant differences between single-CPU run and multi-CPU run (due to uninits?)
!--- square current computing in cell corner step by step

      if (has_dir(xdim)) then
         dbx(2:cg%nx,:,:) = (b(iby,2:cg%nx,:,:)-b(iby,1:cg%nx-1,:,:))*cg%idl(xdim) ;      dbx(1,:,:) = dbx(2,:,:)
      endif
      if (has_dir(ydim)) then
         dby(:,2:cg%ny,:) = (b(ibx,:,2:cg%ny,:)-b(ibx,:,1:cg%ny-1,:))*cg%idl(ydim) ;      dby(:,1,:) = dby(:,2,:)
      endif
      if (has_dir(zdim)) then
         dbz(:,:,2:cg%nz) = (b(iby,:,:,2:cg%nz)-b(iby,:,:,1:cg%nz-1))*cg%idl(zdim) ;      dbz(:,:,1) = dbz(:,:,2)
      endif

!--- current_z **2
      eh = dbx - dby
      if (has_dir(zdim)) then
         wb(:,:,2:cg%nz) =                   0.25*(eh(:,:,2:cg%nz) + eh(:,:,1:cg%nz-1))**2 ; wb(:,:,1) = wb(:,:,2)
      else
         wb = eh**2
      endif
!--- current_x **2
      eh = dby - dbz
      if (has_dir(xdim)) then
         wb(2:cg%nx,:,:) = wb(2:cg%nx,:,:) + 0.25*(eh(2:cg%nx,:,:) + eh(1:cg%nx-1,:,:))**2 ; wb(1,:,:) = wb(2,:,:)
      else
         wb = wb + eh**2
      endif
!--- current_y **2
      eh = dbz - dbx
      if (has_dir(ydim)) then
         wb(:,2:cg%ny,:) = wb(:,2:cg%ny,:) + 0.25*(eh(:,2:cg%ny,:) + eh(:,1:cg%ny-1,:))**2 ; wb(:,1,:) = wb(:,2,:)
      else
         wb = wb + eh**2
      endif

      eta(:,:,:) = eta_0 + eta_1 * sqrt( max(0.0,wb(:,:,:)- jc2 ))

      eh = 0.0
      if (has_dir(xdim)) then
         eh(2:cg%nx-1,:,:) = eh(2:cg%nx-1,:,:) + eta(1:cg%nx-2,:,:) + eta(3:cg%nx,:,:) ;  eh(1,:,:) = eh(2,:,:) ; eh(cg%nx,:,:) = eh(cg%nx-1,:,:)
      endif
      if (has_dir(ydim)) then
         eh(:,2:cg%ny-1,:) = eh(:,2:cg%ny-1,:) + eta(:,1:cg%ny-2,:) + eta(:,3:cg%ny,:) ;  eh(:,1,:) = eh(:,2,:) ; eh(:,cg%ny,:) = eh(:,cg%ny-1,:)
      endif
      if (has_dir(zdim)) then
         eh(:,:,2:cg%nz-1) = eh(:,:,2:cg%nz-1) + eta(:,:,1:cg%nz-2) + eta(:,:,3:cg%nz) ;  eh(:,:,1) = eh(:,:,2) ; eh(:,:,cg%nz) = eh(:,:,cg%nz-1)
      endif
      eh = real((eh + eta_scale*eta)*d_eta_factor)

      where (eta > eta_0)
         eta = eh
      endwhere

      etamax%val       = maxval(eta(cg%is:cg%ie, cg%js:cg%je, cg%ks:cg%ke))
      etamax%loc       = maxloc(eta(cg%is:cg%ie, cg%js:cg%je, cg%ks:cg%ke)) + [cg%nb, cg%nb, cg%nb]
      call mpifind(etamax, MAXL)
      call MPI_Bcast(etamax%val, 1, MPI_DOUBLE_PRECISION, 0, comm, ierr)

      cu2max%val       = maxval(wb(cg%is:cg%ie, cg%js:cg%je, cg%ks:cg%ke))
      cu2max%loc       = maxloc(wb(cg%is:cg%ie, cg%js:cg%je, cg%ks:cg%ke)) + [cg%nb, cg%nb, cg%nb]
      call mpifind(cu2max, MAXL)

#ifndef ISO
      wb = ( u(flind%ion%ien,:,:,:) - 0.5*( u(flind%ion%imx,:,:,:)**2  + u(flind%ion%imy,:,:,:)**2  + u(flind%ion%imz,:,:,:)**2 ) &
           / u(flind%ion%idn,:,:,:) - 0.5 * ( b(ibx,:,:,:)**2  +   b(iby,:,:,:)**2  +   b(ibz,:,:,:)**2))/ ( eta * wb+small)
      dt_eint = deint_max * abs(minval(wb(cg%is:cg%ie, cg%js:cg%je, cg%ks:cg%ke)))

      deimin%val       = minval(wb(cg%is:cg%ie, cg%js:cg%je, cg%ks:cg%ke))
      deimin%loc       = minloc(wb(cg%is:cg%ie, cg%js:cg%je, cg%ks:cg%ke)) + [cg%nb, cg%nb, cg%nb]
      call mpifind(deimin, MINL)
#endif /* !ISO */

   end subroutine compute_resist

!-----------------------------------------------------------------------

   subroutine timestep_resist

      use constants, only: big
      use grid,      only: cg
      use mpisetup,  only: comm, ierr
      use mpi,       only: MPI_DOUBLE_PRECISION, MPI_MIN, MPI_IN_PLACE

      implicit none

      if (etamax%val /= 0.) then
         dt_resist = cfl_resist * cg%dxmn**2 / (2. * etamax%val)
#ifndef ISO
         dt_resist = min(dt_resist,dt_eint)
#endif /* !ISO */
      else
         dt_resist = big
      endif

      call MPI_Allreduce(MPI_IN_PLACE, dt_resist, 1, MPI_DOUBLE_PRECISION, MPI_MIN, comm, ierr)

   end subroutine timestep_resist

!-----------------------------------------------------------------------------
!>
!! \brief
!! \todo overload me or use class(*) if you dare
!<
   subroutine vanleer_limiter(f,a,b)

      implicit none

      real, dimension(:), intent(in)      :: a !< second order correction of left- or right- moving waves flux on the left cell boundary
      real, dimension(:), intent(in)      :: b !< second order correction of left- or right- moving waves flux on the right cell boundary
      real, dimension(:), intent(inout)   :: f !< second order flux correction for left- or right- moving waves
      ! locals
      real, dimension(size(a,1)) :: c !< a*b

      c = a*b                                                                    !> \todo OPTIMIZE ME
      where (c > 0.0)
         f = f+2.0*c/(a+b)
      endwhere

   end subroutine vanleer_limiter

   subroutine tvdd_1d(b1d,eta1d,idi,dt,wcu1d)

      implicit none

      real, dimension(:), pointer, intent(in)    :: eta1d, b1d
      real, dimension(:), intent(out)            :: wcu1d
      real, intent(in)                           :: idi,dt

      real, dimension(size(b1d))                 :: w, wp, wm, b1
      integer                                    :: n

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
      use fluidindex,    only: iby
      use grid,          only: cg
      use magboundaries, only: bnd_emf
      use mpisetup,      only: dt, has_dir
      use constants,     only: xdim, zdim

      implicit none

      real, dimension(:), pointer :: b1d, eta1d
      real, dimension(cg%nx)      :: wcu1d
      integer                     :: j, k

      call compute_resist
      eta(:,:,1:cg%nz-1) = 0.5*(eta(:,:,1:cg%nz-1)+eta(:,:,2:cg%nz))

      do j = 1, cg%ny
         do k = 1, cg%nz
            b1d    => b(iby,:,j,k)
            eta1d  => eta(:,j,k)
            call tvdd_1d(b1d, eta1d, cg%idl(xdim), dt, wcu1d)
            wcu(:,j,k) = wcu1d
         enddo
      enddo

      do j = xdim, zdim
         if (has_dir(j)) call bnd_emf(wcu,'emfz',j)
      enddo

   end subroutine diffuseby_x

   subroutine diffusebz_x

      use arrays,        only: wcu, b
      use fluidindex,    only: ibz
      use grid,          only: cg
      use magboundaries, only: bnd_emf
      use mpisetup,      only: dt, has_dir
      use constants,     only: xdim, zdim

      implicit none

      real, dimension(:), pointer :: b1d, eta1d
      real, dimension(cg%nx)      :: wcu1d
      integer                     :: j, k

      call compute_resist
      eta(:,1:cg%ny-1,:) = 0.5*(eta(:,1:cg%ny-1,:)+eta(:,2:cg%ny,:))

      do j = 1, cg%ny
         do k = 1, cg%nz
            b1d    => b(ibz,:,j,k)
            eta1d  => eta(:,j,k)
            call tvdd_1d(b1d, eta1d, cg%idl(xdim), dt, wcu1d)
            wcu(:,j,k) = wcu1d
         enddo
      enddo

      do j = xdim, zdim
         if (has_dir(j)) call bnd_emf(wcu,'emfy',j)
      enddo

   end subroutine diffusebz_x

   subroutine diffusebz_y

      use arrays,        only: wcu, b
      use fluidindex,    only: ibz
      use grid,          only: cg
      use magboundaries, only: bnd_emf
      use mpisetup,      only: dt, has_dir
      use constants,     only: xdim, ydim, zdim

      implicit none

      real, dimension(:), pointer :: b1d, eta1d
      real, dimension(cg%ny)      :: wcu1d
      integer                     :: i, k

      call compute_resist
      eta(1:cg%nx-1,:,:) = 0.5*(eta(1:cg%nx-1,:,:)+eta(2:cg%nx,:,:))

      do i = 1, cg%nx
         do k = 1, cg%nz
            b1d    => b(ibz,i,:,k)
            eta1d  => eta(i,:,k)
            call tvdd_1d(b1d, eta1d, cg%idl(ydim), dt, wcu1d)
            wcu(i,:,k) = wcu1d
         enddo
      enddo

      do i = xdim, zdim
         if (has_dir(i)) call bnd_emf(wcu,'emfx',i)
      enddo

   end subroutine diffusebz_y

   subroutine diffusebx_y

      use arrays,        only: wcu, b
      use fluidindex,    only: ibx
      use grid,          only: cg
      use magboundaries, only: bnd_emf
      use mpisetup,      only: dt, has_dir
      use constants,     only: xdim, ydim, zdim

      implicit none

      real, dimension(:), pointer :: b1d, eta1d
      real, dimension(cg%ny)      :: wcu1d
      integer                     :: i, k

      call compute_resist
      eta(:,:,1:cg%nz-1) = 0.5*(eta(:,:,1:cg%nz-1)+eta(:,:,2:cg%nz))

      do i = 1, cg%nx
         do k = 1, cg%nz
            b1d    => b(ibx,i,:,k)
            eta1d  => eta(i,:,k)
            call tvdd_1d(b1d, eta1d, cg%idl(ydim), dt, wcu1d)
            wcu(i,:,k) = wcu1d
         enddo
      enddo

      do i = xdim, zdim
         if (has_dir(i)) call bnd_emf(wcu,'emfz',i)
      enddo

   end subroutine diffusebx_y

   subroutine diffusebx_z

      use arrays,        only: wcu, b
      use fluidindex,    only: ibx
      use grid,          only: cg
      use magboundaries, only: bnd_emf
      use mpisetup,      only: dt, has_dir
      use constants,     only: xdim, zdim

      implicit none

      real, dimension(:), pointer :: b1d, eta1d
      real, dimension(cg%nz)      :: wcu1d
      integer                     :: i, j

      call compute_resist
      eta(:,1:cg%ny-1,:) = 0.5*(eta(:,1:cg%ny-1,:)+eta(:,2:cg%ny,:))

      do i = 1, cg%nx
         do j = 1, cg%ny
            b1d    => b(ibx,i,j,:)
            eta1d  => eta(i,j,:)
            call tvdd_1d(b1d, eta1d, cg%idl(zdim), dt, wcu1d)
            wcu(i,j,:) = wcu1d
         enddo
      enddo

      do i = xdim, zdim
         if (has_dir(i)) call bnd_emf(wcu,'emfy',i)
      enddo

   end subroutine diffusebx_z

   subroutine diffuseby_z

      use arrays,        only: wcu, b
      use fluidindex,    only: iby
      use grid,          only: cg
      use magboundaries, only: bnd_emf
      use mpisetup,      only: dt, has_dir
      use constants,     only: xdim, zdim

      implicit none

      real, dimension(:), pointer :: b1d, eta1d
      real, dimension(cg%nz)      :: wcu1d
      integer                     :: i, j

      call compute_resist
      eta(1:cg%nx-1,:,:) = 0.5*(eta(1:cg%nx-1,:,:)+eta(2:cg%nx,:,:))

      do i = 1, cg%nx
         do j = 1, cg%ny
            b1d    => b(iby,i,j,:)
            eta1d  => eta(i,j,:)
            call tvdd_1d(b1d, eta1d, cg%idl(zdim), dt, wcu1d)
            wcu(i,j,:) = wcu1d
         enddo
      enddo

      do j = xdim, zdim
         if (has_dir(j)) call bnd_emf(wcu,'emfx',j)
      enddo

   end subroutine diffuseby_z

end module resistivity
