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

!>
!! \brief Constrained Transport for face-centered magnetic field
!!
!! \details CT routines for RTVD and Riemann solvers.
!!
!! Beware: there is no support for AMR (and none planned).
!<

module ct
! pulled by MAGNETIC

   implicit none

   private
   public  :: magfield

contains

!>
!! \brief Subroutine computes magnetic field evolution
!!
!! The original RTVD MHD scheme incorporates magnetic field evolution via the Constrained transport (CT)
!! algorithm by Evans & Hawley (1988). The idea behind the CT scheme is to integrate numerically the induction equation
!! \f{equation}
!! \frac{\partial\vec{B}}{\partial t} = \nabla\times (\vec{v}  \times \vec{B}),
!! \f}
!! in a manner ensuring that the condition
!! \f{equation}
!! \nabla \cdot \vec{B} = 0,
!! \f}
!! is fulfilled to the machine accuracy.
!!
!! The divergence-free evolution of magnetic-field on a discrete computational
!! grid can be realized if the last equation holds for the discrete representation of the initial
!! condition, and that subsequent updates of \f$B\f$ do not change the total magnetic
!! flux threading cell faces.
!!
!! Time variations of magnetic flux threading surface \f$S\f$ bounded by contour \f$C\f$, due to Stokes theorem, can be written as
!! \f{equation}
!! \frac{\partial\Phi_S}{\partial t} = \frac{\partial}{\partial t} \int_S \vec{B} \cdot \vec{d\sigma} = \int_S \nabla \times ( \vec{v} \times \vec{B})\cdot \vec{d \sigma} =  \oint_{C} (\vec{v} \times \vec{B}) \cdot \vec{dl},
!! \f}
!! where \f$\vec{E}= \vec{v} \times \vec{B}\f$ is electric field, named also electromotive force (EMF).
!!
!! In a discrete representation variations of magnetic %fluxes, in %timestep \f$\Delta t\f$, threading faces of cell \f$(i,j,k)\f$ are
!! given by
!! \f{eqnarray}
!! \frac{\Phi^{x,n+1}_{i+1/2,j,k} - \Phi^{x,n}_{i+1/2,j,k}}{\Delta t} &=&
!! {E}^y_{i+1/2,j,k-1/2} \Delta y + {E}^z_{i+1/2,j+1/2,k} \Delta z     \nonumber\\
!! &-&{E}^y_{i+1/2,j,k+1/2} \Delta y - {E}^z_{i+1/2,j-1/2,k} \Delta z,  \nonumber
!! \f}
!! \f{eqnarray}
!! \frac{\Phi^{y,n+1}_{i,j+1/2,k} - \Phi^{y,n}_{i,j+1/2,k}}{\Delta t} &=&
!! {E}^x_{i,j+1/2,k+1/2} \Delta x + {E}^z_{i-1/2,j+1/2,k} \Delta z      \nonumber\\
!! &-&{E}^x_{i,j+1/2,k-1/2} \Delta x - {E}^z_{i+1/2,j+1/2,k} \Delta z,  \nonumber
!! \f}
!! \f{eqnarray}
!! \frac{\Phi^{z,n+1}_{i,j,k+1/2} - \Phi^{z,n}_{i,j,k+1/2}}{\Delta t} &=&
!! {E}^x_{i,j-1/2,k+1/2} \Delta x + {E}^y_{i+1/2,j,k+1/2} \Delta y      \nonumber\\
!! &-&{E}^x_{i,j+1/2,k+1/2} \Delta x - {E}^y_{i-1/2,j,k+1/2} \Delta y,  \nonumber
!! \f}
!! Variations of magnetic flux threading remaining three faces of cell \f$(i,j,k)\f$ can be written in a similar manner. We note that each EMF
!! contribution appears twice with opposite sign. Thus, the total change of magnetic flux, piercing all cell-faces, vanishes to machine accuracy.
!!
!! To present the idea in a slightly different way, following Pen et al. (2003),  we consider the three components of the induction
!! equation written in the explicit form:
!! \f{eqnarray}
!! \partial_t B_x &=& \partial_y (\mathrm{v_x B_y}) + \partial_z (v_x B_z) - \partial_y (v_y B_x) - \partial_z (v_z B_x),\\
!! \partial_t B_y &=& \partial_x (v_y B_x) + \partial_z (v_y B_z) - \partial_x (\mathrm{v_x B_y}) - \partial_z (v_z B_y),\\
!! \partial_t B_z &=& \partial_x (v_z B_x) + \partial_y (v_z B_y) - \partial_x (v_x B_z) - \partial_y (v_y B_z).
!! \f}
!!
!!We note that each combination of \f$v_a B_b\f$ appears twice in these equations. Let us consider \f$v_x B_y\f$.
!!Once  \f$v_x B_y\f$ is computed for numerical integration of the second equation, it should be also used for
!!integration of the first equation, to ensure cancellation of electromotive forces contributing to the total change of magnetic flux threading
!!cell faces.
!!
!!The scheme proposed by Pen et al. (2003), consists of the following steps:
!!\n (1) Computation of the edge-centered EMF component \f$ v_x B_y\f$.
!!\n (2) Update of \f$B_y\f$, according to the equation \f$\partial_t B_y = \partial_x (v_x B_y)\f$.
!!\n (3) Update of \f$B_x\f$, according to the equation \f$\partial_t B_x = \partial_y (v_x B_y)\f$.
!!
!!Analogous procedure applies to remaining EMF components.
!<

   subroutine tvdb(vibj, b, vg, n, dt, idi)

      use constants, only: big, half

      implicit none

      integer(kind=4),             intent(in)    :: n       !< array size
      real,                        intent(in)    :: dt      !< time step
      real,                        intent(in)    :: idi     !< cell length, depends on direction x, y or z
      real, dimension(:), pointer, intent(inout) :: vibj    !< face-centered electromotive force components (b*vg)
      real, dimension(:), pointer, intent(in)    :: b       !< magnetic field
      real, dimension(n),          intent(in)    :: vg      !< velocity in the center of cell boundary
! locals
      real, dimension(n)                         :: b1      !< magnetic field
      real, dimension(n)                         :: vibj1   !< face-centered electromotive force (EMF) components (b*vg)
      real, dimension(n)                         :: vh      !< velocity interpolated to the cell edges
      real                                       :: dti     !< dt/di
      real                                       :: v       !< auxiliary variable to compute EMF
      real                                       :: w       !< EMF component
      real                                       :: dw      !< The second-order correction to EMF component
      real                                       :: dwm     !< face centered EMF interpolated to left cell-edge
      real                                       :: dwp     !< face centered EMF interpolated to right cell-edge
      integer                                    :: i       !< auxiliary array indicator
      integer                                    :: ip      !< i+1
      integer                                    :: ipp     !< i+2
      integer                                    :: im      !< i-1

! unlike the B field, the vibj lives on the right cell boundary
      vh = 0.0

! velocity interpolation to the cell boundaries

      vh(1:n-1) =(vg(1:n-1)+ vg(2:n))*half;     vh(n) = vh(n-1)

      dti = dt*idi

! face-centered EMF components computation, depending on the sign of vh, the components are upwinded  to cell edges, leading to 1st order EMF

      where (vh > 0.)
         vibj1=b*vg
      elsewhere
         vibj1=eoshift(b*vg,1,boundary=big)
      endwhere

! values of magnetic field computation in Runge-Kutta half step

      ! TODO: GEOFACTOR missing
      b1(2:n) = b(2:n) - (vibj1(2:n) - vibj1(1:n-1)) * dti * half;    b1(1) = b(2)

      do i = 3, n-3
         ip  = i  + 1
         ipp = ip + 1
         im  = i  - 1
         v   = vh(i)

! recomputation of EMF components (w) with b1 and face centered EMF interpolation to cell-edges (dwp, dwm), depending on the sign of v.

         if (v > 0.0) then
            w   = vg(i) * b1(i)
            dwp = (vg(ip) * b1(ip) - w) * half
            dwm = (w - vg(im) * b1(im)) * half
         else
            w   = vg(ip) * b1(ip)
            dwp = (w - vg(ipp) * b1(ipp)) * half
            dwm = (vg(i) * b1(i) - w) * half
         endif

! the second-order corrections to the EMF components computation with the aid of the van Leer monotonic interpolation and 2nd order EMF computation

         dw=0.0
         if (dwm * dwp > 0.0) dw = 2.0 * dwm * dwp / (dwm + dwp)
         vibj(i) = (w + dw) * dt
      enddo

   end subroutine tvdb

!-------------------------------------------------------------------------------------------------------------------

   subroutine magfield(dir)

      use constants,   only: ndims, I_ONE
      use dataio_pub,  only: die
      use global,      only: cc_mag
#ifdef RESISTIVE
      use resistivity, only: diffuseb
#endif /* RESISTIVE */

      implicit none

      integer(kind=4), intent(in) :: dir

      integer(kind=4)             :: bdir, dstep

      if (cc_mag) call die("[ct:magfield] cell-centered magnetic field is not allowed for constrained transport")

      do dstep = 0, 1
         bdir  = I_ONE + mod(dir+dstep,ndims)
#ifdef IONIZED
         call advectb(bdir, dir)
#endif /* IONIZED */
#ifdef RESISTIVE
         call diffuseb(bdir, dir)
#endif /* RESISTIVE */
#if defined(IONIZED) || defined(RESISTIVE)
         call mag_add(dir, bdir)
#endif /* IONIZED || RESISTIVE */
      enddo

   end subroutine magfield

#ifdef IONIZED
!------------------------------------------------------------------------------------------

!>
!!   advectby_x --> advectb(bdir=ydim, vdir=xdim, emf='vxby')
!!   advectbz_x --> advectb(bdir=zdim, vdir=xdim, emf='vxbz')
!!   advectbx_y --> advectb(bdir=xdim, vdir=ydim, emf='vybx')
!!   advectbz_y --> advectb(bdir=zdim, vdir=ydim, emf='vybz')
!!   advectbx_z --> advectb(bdir=xdim, vdir=zdim, emf='vzbx')
!!   advectby_z --> advectb(bdir=ydim, vdir=zdim, emf='vzby')
!<

!>
!! \todo remove workaround for http://gcc.gnu.org/bugzilla/show_bug.cgi?id=48955
!<

   subroutine advectb(bdir, vdir)

      use cg_cost_data,     only: I_MHD
      use cg_leaves,        only: leaves
      use cg_list,          only: cg_list_element
      use constants,        only: xdim, ydim, zdim, LO, HI, ndims, INT4
      use dataio_pub,       only: die
      use domain,           only: dom
      use fluidindex,       only: flind
      use global,           only: dt
      use grid_cont,        only: grid_container
      use magboundaries,    only: bnd_emf
      use named_array_list, only: qna, wna

      implicit none

      integer(kind=4),       intent(in) :: bdir, vdir
      integer, pointer                  :: i1, i2, i1m, i2m
      integer, dimension(ndims), target :: ii, im
      integer                           :: rdir, i, j
      integer(kind=4)                   :: imom                   !< index of vdir momentum
      integer(kind=4)                   :: dir
      integer(kind=4), dimension(ndims) :: emf
      real, dimension(:), allocatable   :: vv, vv0 !< \todo workaround for bug in gcc-4.6, REMOVE ME
      real, dimension(:),    pointer    :: pm1, pm2, pd1, pd2
      type(cg_list_element), pointer    :: cgl
      type(grid_container),  pointer    :: cg
      real, dimension(:), pointer       :: vibj => null()

      imom = flind%ion%idn + int(vdir, kind=4)
      rdir = sum([xdim,ydim,zdim]) - bdir - vdir
      emf(vdir) = 1_INT4 ; emf(bdir) = 2_INT4 ; emf(rdir) = 3_INT4

      if (mod(3+vdir-bdir,3) == 1) then     ! even permutation
         i1 => ii(rdir) ; i1m => ii(rdir) ; i2 => ii(bdir) ; i2m => im(bdir)
      elseif (mod(3+vdir-bdir,3) == 2) then !  odd permutation
         i1 => ii(bdir) ; i1m => im(bdir) ; i2 => ii(rdir) ; i2m => ii(rdir)
      else
         call die('[ct:advectb] neither even nor odd permutation.')
         i1 => ii(rdir) ; i1m => ii(rdir) ; i2 => ii(rdir) ; i2m => ii(rdir) ! suppress compiler warnings
      endif

      cgl => leaves%first
      do while (associated(cgl))
         cg => cgl%cg
         call cg%costs%start

         if (any([allocated(vv), allocated(vv0)])) call die("[ct:advectb] vv or vv0 already allocated")
         allocate(vv(cg%n_(vdir)), vv0(cg%n_(vdir)))

         im(bdir) = cg%lhn(bdir, LO)
         do i = cg%lhn(bdir, LO) + dom%D_(bdir), cg%lhn(bdir, HI)
            ii(bdir) = i
            do j = cg%ijkse(rdir,LO), cg%ijkse(rdir,HI)
               ii(rdir) = j
               im(rdir) = ii(rdir)
               vv=0.0
               pm1 => cg%w(wna%fi)%get_sweep(vdir,imom,i1m,i2m)
               pm2 => cg%w(wna%fi)%get_sweep(vdir,imom,i1 ,i2 )
               pd1 => cg%w(wna%fi)%get_sweep(vdir,flind%ion%idn,i1m,i2m)
               pd2 => cg%w(wna%fi)%get_sweep(vdir,flind%ion%idn,i1 ,i2 )
               vv0(:) = (pm1+pm2)/(pd1+pd2) !< \todo workaround for bug in gcc-4.6, REMOVE ME
               !vv(2:cg%n_(vdir)-1)=(vv(1:cg%n_(vdir)-2) + vv(3:cg%n_(vdir)) + 2.0*vv(2:cg%n_(vdir)-1))*0.25
               vv(2:cg%n_(vdir)-1)=(vv0(1:cg%n_(vdir)-2) + vv0(3:cg%n_(vdir)) + 2.0*vv0(2:cg%n_(vdir)-1))*0.25 !< \todo workaround for bug in gcc-4.6, REMOVE ME
               vv(1)  = vv(2)
               vv(cg%n_(vdir)) = vv(cg%n_(vdir)-1)

               vibj => cg%q(qna%wai)%get_sweep(vdir, i1, i2)
               call tvdb(vibj, cg%w(wna%bi)%get_sweep(vdir,bdir, i1, i2), vv, cg%n_(vdir),dt, cg%idl(vdir))
               NULLIFY(pm1, pm2, pd1, pd2)

            enddo
            im(bdir) = ii(bdir)
         enddo

         do dir = xdim, zdim
            if (dom%has_dir(dir)) call bnd_emf(qna%wai,emf(dir),dir, cg)
         enddo

         deallocate(vv, vv0)

         call cg%costs%stop(I_MHD)
         cgl => cgl%nxt
      enddo
      NULLIFY(i1, i1m, i2, i2m)

   end subroutine advectb
#endif /* IONIZED */

!-------------------------------------------------------------------------------------------------------------------

#if defined(IONIZED) || defined(RESISTIVE)
   subroutine mag_add(dim1, dim2)

      use cg_cost_data,     only: I_MHD
      use cg_leaves,        only: leaves
      use cg_list,          only: cg_list_element
      use dataio_pub,       only: die
      use global,           only: cc_mag
      use grid_cont,        only: grid_container
      use all_boundaries,   only: all_mag_boundaries
      use user_hooks,       only: custom_emf_bnd
#ifdef RESISTIVE
      use constants,        only: wcu_n
      use dataio_pub,       only: die
      use domain,           only: is_multicg
      use named_array_list, only: qna
#endif /* RESISTIVE */

      implicit none

      integer(kind=4), intent(in)    :: dim1, dim2

      type(cg_list_element), pointer :: cgl
      type(grid_container),  pointer :: cg
#ifdef RESISTIVE
      real, dimension(:,:,:), pointer :: wcu
#endif /* RESISTIVE */

      if (cc_mag) call die("[ct:mag_add] cell-centered magnetic field is not allowed for constrained transport")

      cgl => leaves%first
      do while (associated(cgl))
         cg => cgl%cg
         call cg%costs%start

#ifdef RESISTIVE
! DIFFUSION FULL STEP
         wcu => cg%q(qna%ind(wcu_n))%arr
         if (is_multicg) call die("[ct:mag_add] multiple grid pieces per processor not implemented yet") ! not tested custom_emf_bnd
         if (associated(custom_emf_bnd)) call custom_emf_bnd(wcu)
         cg%b(dim2,:,:,:) = cg%b(dim2,:,:,:) -              wcu*cg%idl(dim1)
         cg%b(dim2,:,:,:) = cg%b(dim2,:,:,:) + pshift(wcu,dim1)*cg%idl(dim1)
         cg%b(dim1,:,:,:) = cg%b(dim1,:,:,:) +              wcu*cg%idl(dim2)
         cg%b(dim1,:,:,:) = cg%b(dim1,:,:,:) - pshift(wcu,dim2)*cg%idl(dim2)
#endif /* RESISTIVE */
#ifdef IONIZED
! ADVECTION FULL STEP
         if (associated(custom_emf_bnd)) call custom_emf_bnd(cg%wa)
         cg%b(dim2,:,:,:) = cg%b(dim2,:,:,:) - cg%wa*cg%idl(dim1)
         cg%wa = mshift(cg%wa,dim1)
         cg%b(dim2,:,:,:) = cg%b(dim2,:,:,:) + cg%wa*cg%idl(dim1)
         cg%b(dim1,:,:,:) = cg%b(dim1,:,:,:) - cg%wa*cg%idl(dim2)
         cg%wa = pshift(cg%wa,dim2)
         cg%b(dim1,:,:,:) = cg%b(dim1,:,:,:) + cg%wa*cg%idl(dim2)
#endif /* IONIZED */

         call cg%costs%stop(I_MHD)
         cgl => cgl%nxt
      enddo

      call all_mag_boundaries

   end subroutine mag_add
#endif /* IONIZED || RESISTIVE */

!-------------------------------------------------------------------------------------------------------------------

!>
!! \brief Function pshift makes one-cell, forward circular shift of 3D array in any direction
!! \param tab input array
!! \param d direction of the shift, where 1,2,3 corresponds to \a x,\a y,\a z respectively
!! \return real, dimension(size(tab,1),size(tab,2),size(tab,3))
!!
!! The function was written in order to significantly improve
!! the performance at the cost of the flexibility of original \p CSHIFT.
!<
   function pshift(tab, d)

      use dataio_pub,    only: warn

      implicit none

      real, dimension(:,:,:), intent(inout) :: tab
      integer(kind=4),        intent(in)    :: d

      integer :: ll
      real, dimension(size(tab,1),size(tab,2),size(tab,3)) :: pshift

      ll = size(tab,d)

      if (ll==1) then
         pshift = tab
         return
      endif

      if (d==1) then
         pshift(1:ll-1,:,:) = tab(2:ll,:,:); pshift(ll,:,:) = tab(1,:,:)
      else if (d==2) then
         pshift(:,1:ll-1,:) = tab(:,2:ll,:); pshift(:,ll,:) = tab(:,1,:)
      else if (d==3) then
         pshift(:,:,1:ll-1) = tab(:,:,2:ll); pshift(:,:,ll) = tab(:,:,1)
      else
         call warn('[ct:pshift]: Dim ill defined in pshift!')
      endif

   end function pshift

!>
!! \brief Function mshift makes one-cell, backward circular shift of 3D array in any direction
!! \param tab input array
!! \param d direction of the shift, where 1,2,3 corresponds to \a x,\a y,\a z respectively
!! \return real, dimension(size(tab,1),size(tab,2),size(tab,3))
!!
!! The function was written in order to significantly improve
!! the performance at the cost of the flexibility of original \p CSHIFT.
!<
   function mshift(tab,d)

      use dataio_pub,    only: warn

      implicit none

      real, dimension(:,:,:), intent(inout) :: tab
      integer(kind=4),        intent(in)    :: d

      integer :: ll
      real, dimension(size(tab,1) , size(tab,2) , size(tab,3)) :: mshift

      ll = size(tab,d)

      if (ll==1) then
         mshift = tab
         return
      endif

      if (d==1) then
         mshift(2:ll,:,:) = tab(1:ll-1,:,:); mshift(1,:,:) = tab(ll,:,:)
      else if (d==2) then
         mshift(:,2:ll,:) = tab(:,1:ll-1,:); mshift(:,1,:) = tab(:,ll,:)
      else if (d==3) then
         mshift(:,:,2:ll) = tab(:,:,1:ll-1); mshift(:,:,1) = tab(:,:,ll)
      else
         call warn('[ct:mshift]: Dim ill defined in mshift!')
      endif

   end function mshift

end module ct
