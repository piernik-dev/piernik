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
!    Hyperbolic divergence cleaning
!    A. Dedner et al., Journal of Computational Physics 175, 645-673 (2002)
!    Dr. Varadarajan Parthasarathy, Dr. Artur Gawryszczak, CAMK (Warszawa)
!    November 2017
!
#include "piernik.h"

module hdc

! pulled by ANY

   use constants, only: dsetnamelen, INVALID

   implicit none

   real, protected :: chspeed
   integer, protected :: igp = INVALID !< index of grad(psi) array
   character(len=dsetnamelen), parameter :: gradpsi_n = "grad_psi"

   private
   public :: chspeed, update_chspeed, map_chspeed, init_psi, glmdamping, eglm

contains

!>
!! \brief Allocate extra space for grad psi to speed up calculations
!!
!! To create the auxiliary array just use the call
!!    if (igp == INVALID) call aux_var
!! It should safe to call it multiple times from here - no new storage is
!! created as long as the index igp is consistent with wna entry.
!!
!! The divergence of the B field should be stored in cg%q(idivB)%arr(:,:,:) (use div_B, only: idivB)
!! The gradient of the psi field should be stored in cg%w(igp)%arr(xdim:zdim,:,:,:)
!!
!! These arrays will be automagically freed with destruction of grid containers.
!<

   subroutine aux_var

      use cg_list_global,   only: all_cg
      use constants,        only: INVALID, ndims
      use dataio_pub,       only: die
      use named_array_list, only: wna

      implicit none

      if (wna%exists(gradpsi_n)) then
         if (igp /= wna%ind(gradpsi_n)) call die ("[hdc:aux_var] wna%exists(gradpsi_n) .and. igp /= wna%ind(gradpsi_n)")
      else
         if (igp /= INVALID) call die ("[hdc:aux_var] .not. wna%exists(gradpsi_n) .and. igp /= INVALID")
         call all_cg%reg_var(gradpsi_n, dim4=ndims)
         igp = wna%ind(gradpsi_n)
      endif

   end subroutine aux_var

!>
!! \brief recalcualte the speed of propagation of psi waves
!!
!! This can be perhaps evaluated somewhere in timestep and reused here.
!<

   subroutine update_chspeed()

      use allreduce,    only: piernik_MPI_Allreduce
      use cg_cost_data, only: I_MHD
      use cg_leaves,    only: leaves
      use cg_list,      only: cg_list_element
      use constants,    only: GEO_XYZ, pMAX, pMIN, small, DIVB_HDC, RIEMANN_SPLIT, RIEMANN_UNSPLIT
      use dataio_pub,   only: die
      use domain,       only: dom
      use global,       only: use_fargo, cfl_glm, ch_grid, dt, divB_0_method, which_solver

      implicit none

      type(cg_list_element), pointer  :: cgl
      integer                         :: i, j, k

      if (divB_0_method /= DIVB_HDC) return
      if (all(which_solver /= [RIEMANN_SPLIT, RIEMANN_UNSPLIT])) call die("[hdc:update_chspeed] Only Riemann solvers have DIVB_HDC implemented")

      if (use_fargo) call die("[hdc:update_chspeed] FARGO is not implemented here yet.")
      if (dom%geometry_type /= GEO_XYZ) call die("[hdc:update_chspeed] non-cartesian geometry not implemented yet.")

      chspeed = merge(huge(1.), small, ch_grid)

      cgl => leaves%first
      do while (associated(cgl))
         call cgl%cg%costs%start

         if (ch_grid) then
            ! Rely only on grid properties. Psi is an artificial field and psi waves have to propagate as fast as stability permits.
            ! It leads to very bad values when time step drops suddenly (like on last timestep)
            if (dom%eff_dim > 0) &
                 chspeed = min(chspeed, cfl_glm * minval(cgl%cg%dl, mask=dom%has_dir) / dt)
         else
            ! Bind chspeed to fastest possible gas waves. Beware: check whether this works well with AMR.
            do k = cgl%cg%ks, cgl%cg%ke
               do j = cgl%cg%js, cgl%cg%je
                  do i = cgl%cg%is, cgl%cg%ie
                     chspeed = max(chspeed, point_chspeed(cgl%cg, i, j, k))
                  enddo
               enddo
            enddo
         endif

         call cgl%cg%costs%stop(I_MHD)
         cgl => cgl%nxt
      enddo

      call piernik_MPI_Allreduce(chspeed, merge(pMIN, pMAX, ch_grid))

   end subroutine update_chspeed

!> \brief put chspeed to cg%wa to allow for precise logging

   subroutine map_chspeed

      use cg_cost_data, only: I_MHD
      use cg_leaves,    only: leaves
      use cg_list,      only: cg_list_element
      use constants,    only: GEO_XYZ, small
      use dataio_pub,   only: die
      use domain,       only: dom
      use global,       only: cfl_glm, ch_grid, dt

      implicit none

      type(cg_list_element), pointer  :: cgl
      integer                         :: i, j, k

      ! no need to be as strict as in update_chspeed with dying

      if (dom%geometry_type /= GEO_XYZ) call die("[hdc:update_chspeed] non-cartesian geometry not implemented yet.")
      cgl => leaves%first
      do while (associated(cgl))
         call cgl%cg%costs%start

         if (ch_grid) then
            ! Rely only on grid properties. Psi is an artificial field and psi waves have to propagate as fast as stability permits.
            ! It leads to very bad values when time step drops suddenly (like on last timestep)
            if (dt > 0.) then
               cgl%cg%wa = merge(small, cfl_glm * minval(cgl%cg%dl, mask=dom%has_dir) / dt, dom%eff_dim == 0)
            else
               cgl%cg%wa = small
            endif
         else
            ! Bind chspeed to fastest possible gas waves. Beware: check whether this works well with AMR.
            do k = cgl%cg%ks, cgl%cg%ke
               do j = cgl%cg%js, cgl%cg%je
                  do i = cgl%cg%is, cgl%cg%ie
                     cgl%cg%wa(i, j, k) = point_chspeed(cgl%cg, i, j, k)
                  enddo
               enddo
            enddo
         endif

         call cgl%cg%costs%stop(I_MHD)
         cgl => cgl%nxt
      enddo

   end subroutine map_chspeed

!>
!! \brief chspeed at a point
!!
!! \beware The sound speed and |v|+c_s are independently coded in various places independently
!<

   real function point_chspeed(cg, i, j, k) result(chspeed)

      use constants,  only: small
      use dataio_pub, only: die
      use func,       only: emag, ekin
      use global,     only: cfl_glm
      use grid_cont,  only: grid_container
#ifndef ISO
      use constants,  only: xdim, ydim, zdim, half
      use domain,     only: dom
      use fluidindex, only: flind
      use fluids_pub, only: has_ion, has_neu, has_dst
      use fluidtypes, only: component_fluid
#endif /* !ISO */

      implicit none

      type(grid_container), pointer, intent(in) :: cg
      integer,                       intent(in) :: i
      integer,                       intent(in) :: j
      integer,                       intent(in) :: k

#ifndef ISO
      integer                         :: d
      real                            :: pmag, pgam
      class(component_fluid), pointer :: fl
#endif /* !ISO */

      chspeed = small  ! suppress -Wmaybe-uninitialized on chspeed
#ifdef ISO
      chspeed = cfl_glm * cg%cs_iso2(i, j, k)  ! BUG? should be rather sqrt(cs_iso2)
#else /* !ISO */
      if (has_ion) then
         fl => flind%ion
         pmag = emag(cg%b(xdim, i, j, k), cg%b(ydim, i, j, k), cg%b(zdim, i, j, k))  ! 1/2 |B|**2
         pgam = half * fl%gam * fl%gam_1 * (cg%u(fl%ien, i, j, k) - ekin(cg%u(fl%imx, i, j, k), cg%u(fl%imy, i, j, k), cg%u(fl%imz, i, j, k), cg%u(fl%idn, i, j, k)) - pmag)
         ! pgam = 1/2 * gamma * p = 1/2 * gamma * (gamma - 1) * (e - 1/2 * rho * |v|**2 - 1/2 * |B|**2)
         do d = xdim, zdim
            if (dom%has_dir(d)) then
               chspeed = max(chspeed, cfl_glm * ( &
                    abs(cg%u(fl%imx + d - xdim, i, j, k) / cg%u(fl%idn, i, j, k)) + &
                    sqrt( (pgam + pmag + sqrt( (pgam + pmag)**2 - 2 * pgam * cg%b(d, i, j, k)**2) ) / cg%u(fl%idn, i, j, k)  ) ) )
               ! Eqs. (14) and (15) JCoPh 229 (2010) 2117-2138
               ! c_h = cfl * \frac{\Delta l_min}{\Delta t}, where cfl = cfl_glm (or)
               ! c_h = cfl * max (fastest signal in the domain).
               ! fl%get_cs(i, j, k, cg%u, cg%b, cg%cs_iso2) returns upper estimate of fast magnetosonic wave
            endif
         enddo
      else if (has_neu) then
         fl => flind%neu
         pgam = half * fl%gam * fl%gam_1 * (cg%u(fl%ien, i, j, k) - ekin(cg%u(fl%imx, i, j, k), cg%u(fl%imy, i, j, k), cg%u(fl%imz, i, j, k), cg%u(fl%idn, i, j, k)))
         ! pgam = 1/2 * gamma * p = 1/2 * gamma * (gamma - 1) * (e - 1/2 * rho * |v|**2)
         do d = xdim, zdim
            if (dom%has_dir(d)) then
               chspeed = max(chspeed, cfl_glm * ( &
                    abs(cg%u(fl%imx + d - xdim, i, j, k) / cg%u(fl%idn, i, j, k)) + &
                    sqrt( 2*pgam  / cg%u(fl%idn, i, j, k)  ) ) )
            endif
         enddo
      else if (has_dst) then
         fl => flind%dst
         do d = xdim, zdim
            if (dom%has_dir(d)) then
               chspeed = max(chspeed, cfl_glm * abs(cg%u(fl%imx + d - xdim, i, j, k) / cg%u(fl%idn, i, j, k)))
            endif
         enddo
      else
         call die("[hdc:update_chspeed] Don't know what to do with chspeed without ION, NEU and DST")
      endif
#endif /* ISO */

   end function point_chspeed

!--------------------------------------------------------------------------------------------------------------
   subroutine init_psi

      use cg_leaves,        only: leaves
      use constants,        only: psi_n
      use named_array_list, only: qna

      implicit none

      if (qna%exists(psi_n)) call leaves%set_q_value(qna%ind(psi_n), 0.)

   end subroutine init_psi
!-----------------------------------------------------------------------------------------------------------

!>
!! \brief Parabolic damping to psi
!!
!! dedner A. Mignone et al. / Journal of Computational Physics 229 (2010) 5896–5920, eq. 9
!<
   subroutine glmdamping(half)

      use allreduce,        only: piernik_MPI_Allreduce
      use cg_cost_data,     only: I_MHD
      use cg_leaves,        only: leaves
      use cg_list,          only: cg_list_element
      use constants,        only: psi_n, DIVB_HDC, pMIN, RIEMANN_SPLIT, RIEMANN_UNSPLIT
      use dataio_pub,       only: die
      use domain,           only: dom
      use global,           only: glm_alpha, dt, divB_0_method, which_solver
      use named_array_list, only: qna

      implicit none

      logical, optional, intent(in) :: half

      type(cg_list_element), pointer :: cgl
      real :: fac, dt_eff

      if (divB_0_method /= DIVB_HDC) return ! I think it is equivalent to if (.not. qna%exists(psi_n))
      if (all(which_solver /= [RIEMANN_SPLIT, RIEMANN_UNSPLIT])) call die("[hdc:glmdamping] Only Riemann solvers have DIVB_HDC implemented")

      if (qna%exists(psi_n)) then

         fac = 0.

         dt_eff = dt
         if (present(half)) then
            if (half) dt_eff = dt/2.0
         endif

         if (dom%eff_dim > 0) then
            cgl => leaves%first
            do while (associated(cgl))
               call cgl%cg%costs%start

               fac = max(fac, glm_alpha*chspeed/(minval(cgl%cg%dl, mask=dom%has_dir)/dt_eff))

               call cgl%cg%costs%stop(I_MHD)
               cgl => cgl%nxt
            enddo
         endif
         fac = exp(-fac)
         call piernik_MPI_Allreduce(fac, pMIN)

         cgl => leaves%first
         do while (associated(cgl))
            call cgl%cg%costs%start

            cgl%cg%q(qna%ind(psi_n))%arr =  cgl%cg%q(qna%ind(psi_n))%arr * fac

            call cgl%cg%costs%stop(I_MHD)
            cgl => cgl%nxt
         enddo
      endif

! can be simplified to
!    if (qna%exists(psi_n)) call leaves%q_lin_comb( [qna%ind(psi_n), fac], qna%ind(psi_n))
! but for AMR we may decide to use different factors on different levels

   end subroutine glmdamping

!--------------------------------------------------------------------------------------------
   !>
   !! Eq.(38) Dedner et al. to be implemented
   !<
   subroutine eglm

      use all_boundaries, only: all_fluid_boundaries
#ifdef MAGNETIC
      use all_boundaries, only: all_mag_boundaries
#endif /* MAGNETIC */
      use cg_leaves,  only: leaves
      use cg_list,    only: cg_list_element
      use constants,  only: xdim, zdim, psi_n, GEO_XYZ, half, RIEMANN_SPLIT
      use dataio_pub, only: die
      use div_B,      only: divB, idivB
      use domain,     only: dom
      use fluidindex, only: flind
      use fluids_pub, only: has_ion
      use fluidtypes, only: component_fluid
      use global,     only: use_eglm, dt, which_solver
      use grid_cont,  only: grid_container
      use named_array_list, only: qna

      implicit none

      class(component_fluid), pointer  :: fl
      type(cg_list_element),  pointer  :: cgl
      type(grid_container),   pointer  :: cg
      integer                          :: i, j, k
      integer(kind=4)                  :: ipsi

      if (.not. use_eglm) return
      if (which_solver /= RIEMANN_SPLIT) call die("[hdc:eglm] Only Riemann solver has DIVB_HDC implemented")

      if (igp == INVALID) call aux_var
      ipsi = qna%ind(psi_n)
      call divB

      if (has_ion) then
         if (dom%geometry_type /= GEO_XYZ) call die("[hdc:update_chspeed] non-cartesian geometry not implemented yet.")
         fl => flind%ion
         cgl => leaves%first
         do while (associated(cgl))
            cg => cgl%cg
            do k = cgl%cg%ks, cgl%cg%ke
               do j = cgl%cg%js, cgl%cg%je
                  do i = cgl%cg%is, cgl%cg%ie
                     associate (im1 => i - Dom%D_x, ip1 => i + Dom%D_x, &
                          &     jm1 => j - Dom%D_y, jp1 => j + Dom%D_y, &
                          &     km1 => k - Dom%D_z, kp1 => k + Dom%D_z)

                        ! Gradient of psi
                        cg%w(igp)%arr(:,i,j,k) = half * [ &
                             (cg%q(ipsi)%arr(ip1,j,k) - cg%q(ipsi)%arr(im1,j,k)), &
                             (cg%q(ipsi)%arr(i,jp1,k) - cg%q(ipsi)%arr(i,jm1,k)), &
                             (cg%q(ipsi)%arr(i,j,kp1) - cg%q(ipsi)%arr(i,j,km1)) &
                             ] / cg%dl

                     end associate
                  enddo
               enddo
            enddo
            cgl=>cgl%nxt
         enddo

         ! don't fuse these loops - w(igp) depends on psi and then modifies psi

         cgl => leaves%first
         do while (associated(cgl))
            cg => cgl%cg
            do k = cgl%cg%ks, cgl%cg%ke
               do j = cgl%cg%js, cgl%cg%je
                  do i = cgl%cg%is, cgl%cg%ie
                     !Sources

                     ! momentum = momentum - dt*divB*B
                     cgl%cg%u(fl%imx:fl%imz,i,j,k) = cgl%cg%u(fl%imx:fl%imz,i,j,k) - dt * cg%q(idivB)%arr(i,j,k) * cgl%cg%b(xdim:zdim,i,j,k)

                     ! B = B - dt*divB*u
                     cgl%cg%b(xdim:zdim,i,j,k) = cgl%cg%b(xdim:zdim,i,j,k) - dt * cg%q(idivB)%arr(i,j,k) * (cgl%cg%u(fl%imx:fl%imz,i,j,k) / cgl%cg%u(fl%idn,i,j,k))

                     ! e = e - dt* (divB*u.B - B.grad(psi))
                     cgl%cg%u(fl%ien,i,j,k) = cgl%cg%u(fl%ien,i,j,k) - dt * ( &
                          cg%q(idivB)%arr(i,j,k) * dot_product(cgl%cg%u(fl%imx:fl%imz,i,j,k) / cgl%cg%u(fl%idn,i,j,k), cgl%cg%b(xdim:zdim,i,j,k)) - &
                          dot_product(cgl%cg%b(xdim:zdim,i,j,k), cg%w(igp)%arr(xdim:zdim,i,j,k)) )

                     ! psi = psi - dt*u.grad(psi), other term is calculated in damping
                     cgl%cg%q(ipsi)%arr(i,j,k) =  cgl%cg%q(ipsi)%arr(i,j,k) - &
                          dt * dot_product(cgl%cg%u(fl%imx:fl%imz,i,j,k) / cgl%cg%u(fl%idn,i,j,k), cg%w(igp)%arr(xdim:zdim,i,j,k))

                  enddo
               enddo
            enddo
            cgl=>cgl%nxt
         enddo
      endif

      ! OPT: to avoid these boundary exchanges one must provide div(B) and grad(psi) on bigger area and alter whole blocks.
      ! Thi may require extra guardcells
      ! Beware: higher orders of div(B) and grad(psi) may require boundary update at the beginning too
      call all_fluid_boundaries
      call leaves%leaf_arr3d_boundaries(ipsi)
#ifdef MAGNETIC
      call all_mag_boundaries
#endif /* MAGNETIC */

   end subroutine eglm

end module hdc
