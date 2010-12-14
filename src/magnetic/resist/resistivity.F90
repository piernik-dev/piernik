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
   real :: eta_max_proc, eta_max_all
   integer, dimension(3) :: loc_eta_max
!!!! temporary solution, one should _NOT_ allocate 7 arrays of size nx*ny*nz !!!!!!
   real, dimension(:,:,:), allocatable :: w,wm,wp,dw,eta,b1
   real, dimension(:,:,:), allocatable :: wb,etahelp

   contains

   subroutine cleanup_resistivity
      implicit none
      if (allocated(w)  ) deallocate(w)
      if (allocated(wb) ) deallocate(wb)
      if (allocated(wm) ) deallocate(wm)
      if (allocated(wp) ) deallocate(wp)
      if (allocated(dw) ) deallocate(dw)
      if (allocated(eta)) deallocate(eta)
      if (allocated(b1))  deallocate(b1)
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
      use dataio_pub,    only: die
      use grid,          only: nx, ny, nz, has_dir, zdim
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

      if (.not.allocated(w)  ) allocate(w(nx,ny,nz)  )
      if (.not.allocated(wb) ) allocate(wb(nx,ny,nz) )
      if (.not.allocated(wm) ) allocate(wm(nx,ny,nz) )
      if (.not.allocated(wp) ) allocate(wp(nx,ny,nz) )
      if (.not.allocated(dw) ) allocate(dw(nx,ny,nz) )
      if (.not.allocated(eta)) allocate(eta(nx,ny,nz))
      if (.not.allocated(etahelp)) allocate(etahelp(nx,ny,nz))
      if (.not.allocated(b1) ) allocate(b1(nx,ny,nz) )

      jc2 = j_crit**2
      if (has_dir(zdim)) then
         d_eta_factor = 1./(6.+dble(eta_scale))
      else
         d_eta_factor = 1./(4.+dble(eta_scale))
      endif

   end subroutine init_resistivity

   subroutine compute_resist(eta,ici)
      use arrays,       only: b, u
      use constants,    only: small
      use fluidindex,   only: ibx, iby, ibz
      use func,         only: mshift, pshift
      use grid,         only: idl, xdim, ydim, zdim, nx, ny, nz, is, ie, js, je, ks, ke, has_dir
      use mpisetup,     only: comm, ierr
      use mpi,          only: MPI_DOUBLE_PRECISION, MPI_MAX
#ifndef ISO
      use fluidindex,   only: nvar
#endif /* !ISO */

      implicit none
      integer,intent(in)                       :: ici
      real, dimension(nx,ny,nz), intent(inout) :: eta

!--- square current computing in cell corner step by step

!--- current_z
      wb(:,:,:) = (b(iby,:,:,:)-mshift(b(iby,:,:,:),xdim))*idl(xdim)
      wb = wb -   (b(ibx,:,:,:)-mshift(b(ibx,:,:,:),ydim))*idl(ydim)

      eta(:,:,:) = 0.25*( wb(:,:,:) + mshift(wb(:,:,:),zdim) )**2

      if (has_dir(zdim)) then
!--- current_x
         wb(:,:,:) = (b(ibz,:,:,:)-mshift(b(ibz,:,:,:),ydim))*idl(ydim)
         wb = wb -   (b(iby,:,:,:)-mshift(b(iby,:,:,:),zdim))*idl(zdim)

         eta(:,:,:) = eta(:,:,:) + 0.25*( wb(:,:,:)+mshift(wb(:,:,:),xdim) )**2
!--- current_y
         wb(:,:,:) = (b(ibx,:,:,:)-mshift(b(ibx,:,:,:),zdim))*idl(zdim)
         wb = wb -   (b(ibz,:,:,:)-mshift(b(ibz,:,:,:),xdim))*idl(xdim)

         eta(:,:,:) = eta(:,:,:) +0.25*( wb(:,:,:) + mshift(wb(:,:,:),ydim))**2
      endif

!--- wb = current**2
      wb(:,:,:) = eta(:,:,:)

      eta(:,:,:) = eta_0 + eta_1*sqrt(max(0.0,eta(:,:,:)- jc2 ))
!>
!! \todo Following lines are splitted into separate lines because of intel and gnu dbgs
!! shoud that be so? Is there any other solution instead splitting?
!<
      if (has_dir(zdim)) then
         etahelp(:,:,:)  =   mshift(eta(:,:,:),xdim)
         etahelp = etahelp + pshift(eta(:,:,:),xdim)
         etahelp = etahelp + mshift(eta(:,:,:),ydim)
         etahelp = etahelp + pshift(eta(:,:,:),ydim)
         etahelp = etahelp + mshift(eta(:,:,:),zdim)
         etahelp = etahelp + pshift(eta(:,:,:),zdim)
         etahelp = (etahelp+dble(eta_scale)*eta(:,:,:))*d_eta_factor
         where (eta > eta_0)
            eta = etahelp
         endwhere
      else
         where (eta > eta_0)
            eta(:,:,:) = (mshift(eta(:,:,:),xdim) + pshift(eta(:,:,:),xdim) &
                         +mshift(eta(:,:,:),ydim) + pshift(eta(:,:,:),ydim) &
                         +dble(eta_scale)*eta(:,:,:))*d_eta_factor
         endwhere
      endif

      eta_max_proc      = maxval(eta(is:ie,js:je,ks:ke))
      loc_eta_max       = maxloc(eta(is:ie,js:je,ks:ke))

      eta_max = eta_max_proc

      call MPI_Reduce(eta_max_proc, eta_max_all, 1, MPI_DOUBLE_PRECISION, MPI_MAX, 0, comm, ierr)
      call MPI_Bcast (eta_max_all, 1, MPI_DOUBLE_PRECISION, 0, comm, ierr)

      eta_max = eta_max_all

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

! icx = xdim = 1, icy = ydim = 2, icz = zdim = 3
      eta(:,:,:)=0.5*(eta(:,:,:)+pshift(eta(:,:,:),ici))

      return

   end subroutine compute_resist

!-----------------------------------------------------------------------

   subroutine timestep_resist
      use constants, only: big
      use grid,      only: dxmn
      use mpisetup,  only: comm, ierr
      use mpi,       only: MPI_DOUBLE_PRECISION, MPI_MIN

      implicit none
      real :: dx2,dt_resist_min

      if (eta_max .ne. 0.) then
         dx2 = dxmn**2
         dt_resist = cfl_resist*dx2/(2.*eta_max)
#ifndef ISO
         dt_resist = min(dt_resist,dt_eint)
#endif /* !ISO */
      else
         dt_resist = big
      endif

      call MPI_Reduce(dt_resist, dt_resist_min, 1, MPI_DOUBLE_PRECISION, MPI_MIN, 0, comm, ierr)
      call MPI_Bcast (dt_resist_min, 1, MPI_DOUBLE_PRECISION, 0, comm, ierr)

      dt_resist = dt_resist_min

   end subroutine timestep_resist

!-----------------------------------------------------------------------------

   subroutine tvdd(ibi,ici,n)
      use arrays,    only: b, wcu
      use func,      only: mshift, pshift
      use grid,      only: idl
      use mpisetup,  only: dt

      implicit none
      real    :: idi
      integer :: ibi,ici,n


      idi = idl(n)
      eta = 0.0

      call compute_resist(eta,ici)

! HALF STEP
      w(:,:,:) = (b(ibi,:,:,:)-mshift(b(ibi,:,:,:),n))*idi
      w  = eta*w
      b1 = b(ibi,:,:,:)+0.5*(pshift(w,n)-w)*dt*idi

! FULL STEP
      w(:,:,:) = (b1(:,:,:)-mshift(b1(:,:,:),n))*idi
      w  = eta*w
      wp = 0.5*(pshift(w,n)-w)
      wm = 0.5*(w-mshift(w,n))
      dw = 0.
      where (wm*wp > 0.) dw=2.*wm*wp/(wm+wp)
      wcu = (w+dw)*dt
   end subroutine tvdd

!-------------------------------------------------------------------------------

   subroutine diffuseby_x
      use arrays,        only: wcu
      use fluidindex,    only: iby, icz
      use grid,          only: has_dir, xdim, ydim, zdim
      use magboundaries, only: bnd_emf
      implicit none

      call tvdd(iby,icz,xdim)
      if (has_dir(xdim)) call bnd_emf(wcu,'emfz','xdim')
      if (has_dir(ydim)) call bnd_emf(wcu,'emfz','ydim')
      if (has_dir(zdim)) call bnd_emf(wcu,'emfz','zdim')

   end subroutine diffuseby_x

   subroutine diffusebz_x
      use arrays,        only: wcu
      use fluidindex,    only: ibz, icy
      use grid,          only: has_dir, xdim, ydim, zdim
      use magboundaries, only: bnd_emf
      implicit none

      call tvdd(ibz,icy,xdim)
      if (has_dir(xdim)) call bnd_emf(wcu,'emfy','xdim')
      if (has_dir(ydim)) call bnd_emf(wcu,'emfy','ydim')
      if (has_dir(zdim)) call bnd_emf(wcu,'emfy','zdim')

   end subroutine diffusebz_x

   subroutine diffusebz_y
      use arrays,        only: wcu
      use fluidindex,    only: ibz, icx
      use grid,          only: has_dir, xdim, ydim, zdim
      use magboundaries, only: bnd_emf
      implicit none

      call tvdd(ibz,icx,ydim)
      if (has_dir(ydim)) call bnd_emf(wcu,'emfx','ydim')
      if (has_dir(zdim)) call bnd_emf(wcu,'emfx','zdim')
      if (has_dir(xdim)) call bnd_emf(wcu,'emfx','xdim')
   end subroutine diffusebz_y

   subroutine diffusebx_y
      use arrays,        only: wcu
      use fluidindex,    only: ibx, icz
      use grid,          only: has_dir, xdim, ydim, zdim
      use magboundaries, only: bnd_emf
      implicit none

      call tvdd(ibx,icz,ydim)
      if (has_dir(ydim)) call bnd_emf(wcu, 'emfz', 'ydim')
      if (has_dir(zdim)) call bnd_emf(wcu, 'emfz', 'zdim')
      if (has_dir(xdim)) call bnd_emf(wcu, 'emfz', 'xdim')
   end subroutine diffusebx_y

   subroutine diffusebx_z
      use arrays,        only: wcu
      use fluidindex,    only: ibx, icy
      use grid,          only: has_dir, xdim, ydim, zdim
      use magboundaries, only: bnd_emf
      implicit none

      call tvdd(ibx,icy,zdim)
      if (has_dir(zdim)) call bnd_emf(wcu, 'emfy', 'zdim')
      if (has_dir(xdim)) call bnd_emf(wcu, 'emfy', 'xdim')
      if (has_dir(ydim)) call bnd_emf(wcu, 'emfy', 'ydim')
   end subroutine diffusebx_z

   subroutine diffuseby_z
      use arrays,        only: wcu
      use fluidindex,    only: iby, icx
      use grid,          only: has_dir, xdim, ydim, zdim
      use magboundaries, only: bnd_emf
      implicit none

      call tvdd(iby,icx,zdim)
      if (has_dir(zdim)) call bnd_emf(wcu, 'emfx', 'zdim')
      if (has_dir(xdim)) call bnd_emf(wcu, 'emfx', 'xdim')
      if (has_dir(ydim)) call bnd_emf(wcu, 'emfx', 'ydim')
   end subroutine diffuseby_z

end module resistivity
