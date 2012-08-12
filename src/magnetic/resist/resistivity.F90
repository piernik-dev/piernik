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
   use constants, only: dsetnamelen
   use types,     only: value

   implicit none

   private
   public  :: init_resistivity, timestep_resist, cleanup_resistivity, dt_resist, etamax,   &
        &     diffuseb,cu2max, deimin, eta1_active, compute_resist

   real    :: cfl_resist                     !< CFL factor for resistivity effect
   real    :: eta_0                          !< uniform resistivity
   real    :: eta_1                          !< anomalous resistivity
   real    :: j_crit                         !< critical value of current density
   real    :: jc2                            !< squared critical value of current density
   real    :: deint_max                      !< COMMENT ME
   integer(kind=4) :: eta_scale              !< COMMENT ME
   real    :: dt_resist, dt_eint
   real(kind=8)                            :: d_eta_factor
   type(value)                             :: etamax, cu2max, deimin
   logical, save                           :: eta1_active = .true.       !< resistivity off-switcher while eta_1 == 0.0
   character(len=dsetnamelen), parameter   :: eta_n = "eta", wb_n = "wb", eh_n = "eh", dbx_n = "dbx", dby_n = "dby", dbz_n = "dbz"

contains

   subroutine cleanup_resistivity
      implicit none
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

      use cg_list_global, only: all_cg
      use constants,      only: PIERNIK_INIT_GRID, zdim, xdim, ydim, wcu_n
      use dataio_pub,     only: par_file, ierrh, namelist_errh, compare_namelist, cmdl_nml, lun  ! QA_WARN required for diff_nml
      use dataio_pub,     only: die, code_progress
      use domain,         only: dom
      use cg_list,        only: cg_list_element
      use cg_list_bnd,    only: leaves
      use mpisetup,       only: rbuff, ibuff, master, slave, piernik_MPI_Bcast
      use named_array,    only: qna

      implicit none

      real :: dims_twice
      type(cg_list_element),  pointer :: cgl

      namelist /RESISTIVITY/ cfl_resist, eta_0, eta_1, eta_scale, j_crit, deint_max

      if (code_progress < PIERNIK_INIT_GRID) call die("[resistivity:init_resistivity] grid not initialized.")

      cfl_resist   = 0.4
      eta_0        = 0.0
      eta_1        = 0.0
      eta_scale    = 4
      j_crit       = 1.0e6
      deint_max    = 0.01

      if (master) then

         diff_nml(RESISTIVITY)

         ibuff(1) = eta_scale

         rbuff(1) = cfl_resist
         rbuff(2) = eta_0
         rbuff(3) = eta_1
         rbuff(4) = j_crit
         rbuff(5) = deint_max

      endif

      call piernik_MPI_Bcast(ibuff)
      call piernik_MPI_Bcast(rbuff)

      if (slave) then

         eta_scale          = ibuff(1)

         cfl_resist         = rbuff(1)
         eta_0              = rbuff(2)
         eta_1              = rbuff(3)
         j_crit             = rbuff(4)
         deint_max          = rbuff(5)

      endif

      if (eta_scale < 0) call die("eta_scale must be greater or equal 0")

      call all_cg%reg_var(wcu_n)
      call all_cg%reg_var(eta_n)
      call all_cg%reg_var(wb_n)
      call all_cg%reg_var(eh_n)
      call all_cg%reg_var(dbx_n)
      call all_cg%reg_var(dby_n)
      call all_cg%reg_var(dbz_n)
#ifdef ISO
      if (eta_1 == 0.) then
         cgl => leaves%first
         do while (associated(cgl))
            cgl%cg%q(qna%ind(eta_n))%arr = eta_0
            cgl => cgl%nxt
         enddo
         etamax%val  = eta_0
         eta1_active = .false.
      endif
#endif /* ISO */

      if (eta1_active) then

         cgl => leaves%first
         do while (associated(cgl))
            if (.not. dom%has_dir(xdim)) cgl%cg%q(qna%ind(dbx_n))%arr = 0.0
            if (.not. dom%has_dir(ydim)) cgl%cg%q(qna%ind(dby_n))%arr = 0.0
            if (.not. dom%has_dir(zdim)) cgl%cg%q(qna%ind(dbz_n))%arr = 0.0
            cgl => cgl%nxt
         enddo

         jc2 = j_crit**2
         dims_twice = 2. * dom%eff_dim
         d_eta_factor = 1./(dims_twice+real(eta_scale, kind=8))
      endif

   end subroutine init_resistivity

   subroutine compute_resist

      use constants,   only: xdim, ydim, zdim, MAXL, oneq
      use dataio_pub,  only: die
      use domain,      only: dom, is_multicg
      use func,        only: ekin, emag
      use cg_list,     only: cg_list_element
      use cg_list_bnd, only: leaves
      use grid_cont,   only: grid_container
      use mpisetup,    only: piernik_MPI_Bcast
      use named_array, only: qna
#ifndef ISO
      use constants,   only: small, MINL
      use fluidindex,  only: flind
#endif /* !ISO */

      implicit none

      real, dimension(:,:,:), pointer :: p
      type(cg_list_element),  pointer :: cgl
      type(grid_container),   pointer :: cg
      real, dimension(:,:,:), pointer :: eta, dbx, dby, dbz, wb, eh

      if (.not.eta1_active) return
!> \deprecated BEWARE: uninitialized values are poisoning the wb(:,:,:) array - should change  with rev. 3893
!> \deprecated BEWARE: significant differences between single-CPU run and multi-CPU run (due to uninits?)
!--- square current computing in cell corner step by step

      cgl => leaves%first
      do while (associated(cgl))
         cg => cgl%cg
         eta => cg%q(qna%ind(eta_n))%arr
         dbx => cg%q(qna%ind(dbx_n))%arr
         dby => cg%q(qna%ind(dby_n))%arr
         dbz => cg%q(qna%ind(dbz_n))%arr
         wb => cg%q(qna%ind(wb_n))%arr
         eh => cg%q(qna%ind(eh_n))%arr

         if (dom%has_dir(xdim)) then
            dbx(2:cg%n_(xdim),:,:) = (cg%b(ydim,2:cg%n_(xdim),:,:)-cg%b(ydim,1:cg%n_(xdim)-1,:,:))*cg%idl(xdim) ; dbx(1,:,:) = dbx(2,:,:)
         endif
         if (dom%has_dir(ydim)) then
            dby(:,2:cg%n_(ydim),:) = (cg%b(xdim,:,2:cg%n_(ydim),:)-cg%b(xdim,:,1:cg%n_(ydim)-1,:))*cg%idl(ydim) ; dby(:,1,:) = dby(:,2,:)
         endif
         if (dom%has_dir(zdim)) then
            dbz(:,:,2:cg%n_(zdim)) = (cg%b(ydim,:,:,2:cg%n_(zdim))-cg%b(ydim,:,:,1:cg%n_(zdim)-1))*cg%idl(zdim) ; dbz(:,:,1) = dbz(:,:,2)
         endif

!--- current_z **2
         eh = dbx - dby
         if (dom%has_dir(zdim)) then
            wb(:,:,2:cg%n_(zdim)) =                         oneq*(eh(:,:,2:cg%n_(zdim)) + eh(:,:,1:cg%n_(zdim)-1))**2 ; wb(:,:,1) = wb(:,:,2)
         else
            wb = eh**2
         endif
!--- current_x **2
         eh = dby - dbz
         if (dom%has_dir(xdim)) then
            wb(2:cg%n_(xdim),:,:) = wb(2:cg%n_(xdim),:,:) + oneq*(eh(2:cg%n_(xdim),:,:) + eh(1:cg%n_(xdim)-1,:,:))**2 ; wb(1,:,:) = wb(2,:,:)
         else
            wb = wb + eh**2
         endif
!--- current_y **2
         eh = dbz - dbx
         if (dom%has_dir(ydim)) then
            wb(:,2:cg%n_(ydim),:) = wb(:,2:cg%n_(ydim),:) + oneq*(eh(:,2:cg%n_(ydim),:) + eh(:,1:cg%n_(ydim)-1,:))**2 ; wb(:,1,:) = wb(:,2,:)
         else
            wb = wb + eh**2
         endif

         eta(:,:,:) = eta_0 + eta_1 * sqrt( max(0.0,wb(:,:,:)- jc2 ))

         eh = 0.0
         if (dom%has_dir(xdim)) then
            eh(2:cg%n_(xdim)-1,:,:) = eh(2:cg%n_(xdim)-1,:,:) + eta(1:cg%n_(xdim)-2,:,:) + eta(3:cg%n_(xdim),:,:)
            eh(1,:,:) = eh(2,:,:) ; eh(cg%n_(xdim),:,:) = eh(cg%n_(xdim)-1,:,:)
         endif
         if (dom%has_dir(ydim)) then
            eh(:,2:cg%n_(ydim)-1,:) = eh(:,2:cg%n_(ydim)-1,:) + eta(:,1:cg%n_(ydim)-2,:) + eta(:,3:cg%n_(ydim),:)
            eh(:,1,:) = eh(:,2,:) ; eh(:,cg%n_(ydim),:) = eh(:,cg%n_(ydim)-1,:)
         endif
         if (dom%has_dir(zdim)) then
            eh(:,:,2:cg%n_(zdim)-1) = eh(:,:,2:cg%n_(zdim)-1) + eta(:,:,1:cg%n_(zdim)-2) + eta(:,:,3:cg%n_(zdim))
            eh(:,:,1) = eh(:,:,2) ; eh(:,:,cg%n_(zdim)) = eh(:,:,cg%n_(zdim)-1)
         endif
         eh = real((eh + eta_scale*eta)*d_eta_factor)

         where (eta > eta_0) eta = eh

         cgl => cgl%nxt
      enddo

      cg => leaves%first%cg
      if (is_multicg) call die("[resistivity:compute_resist] multiple grid pieces per procesor not implemented yet") !nontrivial get_extremum, wb, eta

      call leaves%get_extremum(qna%ind(eta_n), MAXL, etamax)
      call piernik_MPI_Bcast(etamax%val)
      call leaves%get_extremum(qna%ind(wb_n), MAXL, cu2max)
      etamax%assoc = dt_resist ; cu2max%assoc = dt_resist

#ifndef ISO
      cgl => leaves%first
      do while (associated(cgl))
         cg => cgl%cg
         eta => cg%q(qna%ind(eta_n))%arr
         wb => cg%q(qna%ind(wb_n))%arr
         wb = (cg%u(flind%ion%ien,:,:,:) - ekin(cg%u(flind%ion%imx,:,:,:), cg%u(flind%ion%imy,:,:,:), cg%u(flind%ion%imz,:,:,:), cg%u(flind%ion%idn,:,:,:)) - &
              emag(cg%b(xdim,:,:,:), cg%b(ydim,:,:,:), cg%b(zdim,:,:,:)))/ (eta(:,:,:) * wb+small)
         dt_eint = deint_max * abs(minval(cg%q(qna%ind(wb_n))%span(cg%ijkse)))
         cgl => cgl%nxt
      enddo

      call leaves%get_extremum(qna%ind(wb_n), MINL, deimin)
      deimin%assoc = dt_eint
#endif /* !ISO */
      NULLIFY(p)

   end subroutine compute_resist

!-----------------------------------------------------------------------

   subroutine timestep_resist(cg)

      use constants, only: big, I_ONE
      use grid_cont, only: grid_container
      use mpisetup,  only: comm, mpi_err
      use mpi,       only: MPI_DOUBLE_PRECISION, MPI_MIN, MPI_IN_PLACE

      implicit none

      type(grid_container), pointer, intent(in) :: cg

      if (etamax%val /= 0.) then
         dt_resist = cfl_resist * cg%dxmn**2 / (2. * etamax%val)
#ifndef ISO
         dt_resist = min(dt_resist,dt_eint)
#endif /* !ISO */
      else
         dt_resist = big
      endif

      call MPI_Allreduce(MPI_IN_PLACE, dt_resist, I_ONE, MPI_DOUBLE_PRECISION, MPI_MIN, comm, mpi_err)

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

      use constants,     only: half
      implicit none

      real, dimension(:), pointer, intent(in)    :: eta1d, b1d
      real, dimension(:), pointer, intent(out)   :: wcu1d
      real, intent(in)                           :: idi,dt

      real, dimension(size(b1d))                 :: w, wp, wm, b1
      integer                                    :: n

      n = size(b1d)
      w(2:n)    = eta1d(2:n) * ( b1d(2:n) - b1d(1:n-1) )*idi ;  w(1)  = w(2)
      b1(1:n-1) = b1d(1:n-1) + half*(w(2:n) - w(1:n-1))*dt*idi; b1(n) = b1(n-1)

      w(2:n)    = eta1d(2:n) * ( b1(2:n) - b1(1:n-1) )*idi   ; w(1)  = w(2)
      wp(1:n-1) = half*(w(2:n) - w(1:n-1))                   ; wp(n) = wp(n-1)
      wm(2:n)   = wp(1:n-1)                                  ; wm(1) = wm(2)

      call vanleer_limiter(w,wm,wp)
      wcu1d     = w*dt

   end subroutine tvdd_1d

!-------------------------------------------------------------------------------
!
! 6 routines have been substituted by one with parameters:
!   diffuseby_x  --> diffuseb(ibdir = ydim, sdir = xdim, etadir = zdim, emf = 'emfz', n1 = ydim, n2 = zdim)
!   diffusebz_x  --> diffuseb(ibdir = zdim, sdir = xdim, etadir = ydim, emf = 'emfy', n1 = ydim, n2 = zdim)
!   diffusebz_y  --> diffuseb(ibdir = zdim, sdir = ydim, etadir = xdim, emf = 'emfx', n1 = zdim, n2 = xdim)
!   diffusebx_y  --> diffuseb(ibdir = xdim, sdir = ydim, etadir = zdim, emf = 'emfz', n1 = zdim, n2 = xdim)
!   diffusebx_z  --> diffuseb(ibdir = xdim, sdir = zdim, etadir = ydim, emf = 'emfy', n1 = xdim, n2 = ydim)
!   diffuseby_z  --> diffuseb(ibdir = ydim, sdir = zdim, etadir = xdim, emf = 'emfx', n1 = xdim, n2 = ydim)

   subroutine diffuseb(ibdir, sdir)

      use constants,      only: xdim, ydim, zdim, ndims, half, I_ONE, wcu_n, idm, uv, INT4
      use domain,         only: dom
      use cg_list,        only: cg_list_element
      use global,         only: dt
      use cg_list_bnd,    only: leaves
      use grid_cont,      only: grid_container
      use magboundaries,  only: bnd_emf
      use named_array,    only: qna, wna

      implicit none

      integer(kind=4),  intent(in)   :: ibdir, sdir
      integer                        :: i1, i2, wcu_i, eta_i
      integer(kind=4)                :: n1, n2, etadir, dir, emf
      integer(kind=4), dimension(ndims) :: idml, idmh
      real, dimension(:),    pointer :: b1d, eta1d, wcu1d
      type(cg_list_element), pointer :: cgl
      type(grid_container),  pointer :: cg

      n1 = I_ONE + mod(sdir    ,   ndims)
      n2 = I_ONE + mod(sdir+I_ONE, ndims)
      etadir = sum([xdim,ydim,zdim]) - ibdir - sdir

      call compute_resist

      cgl => leaves%first
      do while (associated(cgl))
         cg => cgl%cg
         wcu_i = qna%ind(wcu_n)
         eta_i = qna%ind(eta_n)

!         select case (etadir)
!            case (xdim)
!               cg%q(eta_i)%arr(1:cg%n_(xdim)-1,:,:) = half*(cg%q(eta_i)%arr(1:cg%n_(xdim)-1,:,:)+cg%q(eta_i)%arr(2:cg%n_(xdim),:,:))
!            case (ydim)
!               cg%q(eta_i)%arr(:,1:cg%n_(ydim)-1,:) = half*(cg%q(eta_i)%arr(:,1:cg%n_(ydim)-1,:)+cg%q(eta_i)%arr(:,2:cg%n_(ydim),:))
!            case (zdim)
!               cg%q(eta_i)%arr(:,:,1:cg%n_(zdim)-1) = half*(cg%q(eta_i)%arr(:,:,1:cg%n_(zdim)-1)+cg%q(eta_i)%arr(:,:,2:cg%n_(zdim)))
!         end select

! following solution seems to be a bit faster than former select case
         idmh(:) = cg%n_(:) - idm(:,etadir)
         idml(:) = I_ONE + idm(:,etadir)
         cg%q(eta_i)%arr(:idmh(xdim),:idmh(ydim),:idmh(zdim)) = half*(cg%q(eta_i)%span(uv,idmh) + cg%q(eta_i)%span(idml,cg%n_))

         do i1 = lbound(cg%q(wcu_i)%arr,n1), ubound(cg%q(wcu_i)%arr,n1)
            do i2 = lbound(cg%q(wcu_i)%arr,n2), ubound(cg%q(wcu_i)%arr,n2)
               b1d   => cg%w(wna%bi)%get_sweep(sdir,ibdir,i1,i2)
               eta1d => cg%q(eta_i    )%get_sweep(sdir,      i1,i2)
               wcu1d => cg%q(wcu_i    )%get_sweep(sdir,      i1,i2)
               call tvdd_1d(b1d, eta1d, cg%idl(sdir), dt, wcu1d)
            enddo
         enddo

         cgl => cgl%nxt
      enddo

      cgl => leaves%first
      do while (associated(cgl))
         do dir = xdim, zdim
            emf = idm(etadir,dir) + 2_INT4
            if (dom%has_dir(dir)) call bnd_emf(cg%q(wcu_i)%arr, emf, dir, cgl%cg)
         enddo
         cgl => cgl%nxt
      enddo

   end subroutine diffuseb

end module resistivity
