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
! pulled by RIEMANN

   use constants, only: dsetnamelen, INVALID

  implicit none

  real, protected :: chspeed
  integer, protected :: igp = INVALID !< index of grad(psi) array
  character(len=dsetnamelen), parameter :: gradpsi_n = "grad_psi"

  private
  public :: chspeed, update_chspeed, init_psi, glmdamping, eglm

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

!>
  !! Acceptable values
  !!----------------------------------------------------------------
  !! Otvortex
  !! WENO3
  !! rk2: cfl = cfl_glm = 0.4 (check this: [divB:print_divB_norm])
  !! muscl_ : cfl = cfl_glm = 0.3
  !! LINEAR
  !! rk2: cfl = cfl_glm = 0.4 (check this: [divB:print_divB_norm])
  !! muscl_ : cfl = cfl_glm = 0.3
  !!---------------------------------------------------------------
  !! Field loop advection
  !! WENO3
  !! rk2: cfl = cfl_glm = 0.4
  !! muscl_ : cfl = cfl_glm = 0.3
  !! LINEAR
  !! rk2 : cfl = cfl_glm = 0.4
  !! muscl_ : cfl = cfl_glm = 0.3
  !!---------------------------------------------------------------
  !! c_h := cfl * \Delta l_min / \Delta t^n
  !! cfl = cfl_glm makes sense in lines 147 and 153
  !!---------------------------------------------------------------
  !! Values 0.4 and 0.3 are chosen from experience and sanity choice.
  !! Yet to test RJ tube and other tests.
  !!---------------------------------------------------------------
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

      use cg_leaves,  only: leaves
      use cg_list,    only: cg_list_element
      use constants,  only: GEO_XYZ, pMAX, small, xdim, ydim, zdim, half
      use dataio_pub, only: die
      use domain,     only: dom
      use fluidindex, only: flind
      use fluids_pub, only: has_ion
      use fluidtypes, only: component_fluid
      use func,       only: emag, ekin
      use global,     only: use_fargo, cfl_glm, ch_grid, dt
      use mpisetup,   only: piernik_MPI_Allreduce

      implicit none

      type(cg_list_element), pointer  :: cgl
      class(component_fluid), pointer :: fl
      integer                         :: i, j, k, d
      real                            :: pmag, pgam

      chspeed = huge(1.)
      if (has_ion) then
         if (use_fargo) call die("[hdc:update_chspeed] FARGO is not implemented here yet.")
         if (dom%geometry_type /= GEO_XYZ) call die("[hdc:update_chspeed] non-cartesian geometry not implemented yet.")
         chspeed = small
         fl => flind%ion
         cgl => leaves%first
         do while (associated(cgl))
            if (ch_grid) then
               ! Rely only on grid properties. Psi is an artificial field and psi waves have to propagate as fast as stability permits.
               ! It leads to very bad values when time step drops suddenly (like on last timestep)
               chspeed = max(chspeed, cfl_glm * minval(cgl%cg%dl, mask=dom%has_dir) / dt)
            else
               ! Bind chspeed to fastest possible gas waves. Beware: this may not always work well with AMR.
               do k = cgl%cg%ks, cgl%cg%ke
                  do j = cgl%cg%js, cgl%cg%je
                     do i = cgl%cg%is, cgl%cg%ie
                        pmag = emag(cgl%cg%b(xdim, i, j, k), cgl%cg%b(ydim, i, j, k), cgl%cg%b(zdim, i, j, k))  ! 1/2 |B|**2
                        pgam = half * fl%gam * fl%gam_1 * (cgl%cg%u(fl%ien, i, j, k) - ekin(cgl%cg%u(fl%imx, i, j, k), cgl%cg%u(fl%imy, i, j, k), cgl%cg%u(fl%imz, i, j, k), cgl%cg%u(fl%idn, i, j, k)) - pmag)
                        ! pgam = 1/2 * gamma * p = 1/2 * gamma * (gamma - 1) * (e - 1/2 * rho * |v|**2 - 1/2 * |B|**2)
                        do d = xdim, zdim
                           if (dom%has_dir(d)) then
                              chspeed = max(chspeed, cfl_glm * ( &
                                   abs(cgl%cg%u(fl%imx + d - xdim, i, j, k) / cgl%cg%u(fl%idn, i, j, k)) + &
                                   sqrt( (pgam + pmag + sqrt( (pgam + pmag)**2 - 2 * pgam * cgl%cg%b(d, i, j, k)**2) ) / cgl%cg%u(fl%idn, i, j, k)  ) &
                                   ) )
                              ! Eqs. (14) and (15) JCoPh 229 (2010) 2117-2138
                              ! c_h = cfl * \frac{\Delta l_min}{\Delta t}, where cfl = cfl_glm (or)
                              ! c_h = cfl * max (fastest signal in the domain).
                              ! fl%get_cs(i, j, k, cgl%cg%u, cgl%cg%b, cgl%cg%cs_iso2) returns upper estimate of fast magnetosonic wave
                           end if
                        enddo
                     enddo
                  enddo
               enddo
            endif
            cgl => cgl%nxt
         enddo
      endif

      call piernik_MPI_Allreduce(chspeed, pMAX)
  end subroutine update_chspeed
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
  !! Parabolic damping to psi
!<
  subroutine glmdamping

     use global,           only: glm_alpha, dt !, cfl
     use cg_list,          only: cg_list_element
     use grid_cont,        only: grid_container
     use cg_leaves,        only: leaves
     use constants,        only: psi_n
     use named_array_list, only: qna
     use domain,           only: dom

     implicit none

     type(cg_list_element), pointer :: cgl
     type(grid_container),  pointer :: cg

     if (qna%exists(psi_n)) then
        cgl => leaves%first
        do while (associated(cgl))
           cg => cgl%cg
           !cgl%cg%q(qna%ind(psi_n))%arr =  cgl%cg%q(qna%ind(psi_n))%arr * exp(-glm_alpha*cfl)
           ! minval(cg%dl,mask=dom%has_dir)
           cgl%cg%q(qna%ind(psi_n))%arr =  cgl%cg%q(qna%ind(psi_n))%arr * exp(-glm_alpha*chspeed/(minval(cg%dl,mask=dom%has_dir)/dt))
           cgl => cgl%nxt
        enddo
     endif

! can be simplified to
!    if (qna%exists(psi_n)) call leaves%q_lin_comb( [qna%ind(psi_n), exp(-glm_alpha*cfl)], qna%ind(psi_n))
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
     use constants,  only: xdim, zdim, psi_n, GEO_XYZ, half
     use dataio_pub, only: die
     use div_B,      only: divB, idivB
     use domain,     only: dom
     use fluidindex, only: flind
     use fluids_pub, only: has_ion
     use fluidtypes, only: component_fluid
     use global,     only: use_eglm, dt
     use grid_cont,  only: grid_container
     use named_array_list, only: qna

     implicit none

     class(component_fluid), pointer  :: fl
     type(cg_list_element),  pointer  :: cgl
     type(grid_container),   pointer  :: cg
     integer                          :: i, j, k, ipsi

     if (.not. use_eglm) return

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
     ! Beware: highre orders of div(B) and grad(psi) may require boundary update at the beginning too
     call all_fluid_boundaries
     call leaves%leaf_arr3d_boundaries(ipsi)
#ifdef MAGNETIC
     call all_mag_boundaries
#endif /* MAGNETIC */

   end subroutine eglm

end module hdc

!--------------------------------------------------------------------------------------------

