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
  public :: chspeed, update_chspeed, init_psi, glm_mhd, glmdamping, glm_3D, eglm

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

      use cg_leaves,  only: leaves
      use cg_list,    only: cg_list_element
      use constants,  only: GEO_XYZ, small
      use dataio_pub, only: die
      use domain,     only: dom
      use fluidindex, only: flind
      use fluids_pub, only: has_ion
      use fluidtypes, only: component_fluid
      use global,     only: use_fargo

      implicit none

      type(cg_list_element), pointer  :: cgl
      class(component_fluid), pointer :: fl
      integer                         :: i, j, k

      chspeed = huge(1.)
      if (has_ion) then
         if (use_fargo) call die("[hdc:update_chspeed] FARGO is not implemented here yet.")
         if (dom%geometry_type /= GEO_XYZ) call die("[hdc:update_chspeed] non-cartesian geometry not implemented yet.")
         chspeed = small
         fl => flind%ion
         cgl => leaves%first
         do while (associated(cgl))
            do k = cgl%cg%ks, cgl%cg%ke
               do j = cgl%cg%js, cgl%cg%je
                  do i = cgl%cg%is, cgl%cg%ie
                     chspeed = max(chspeed, maxval(abs(cgl%cg%u(fl%imx:fl%imz, i, j, k) / cgl%cg%u(fl%idn, i, j, k)) + fl%get_cs(i, j, k, cgl%cg%u, cgl%cg%b, cgl%cg%cs_iso2)))
                     ! for AMR or POLAR it may be better to have explicit dependence on cg%dl(:) here
                     ! OPT: it will also be cheaper
                  enddo
               enddo
            enddo
            cgl => cgl%nxt
         enddo
      endif

    return

  end subroutine update_chspeed
!--------------------------------------------------------------------------------------------------------------
  subroutine init_psi

     use cg_leaves,        only: leaves
     use constants,        only: psi_n
     use global,           only: psi_0
     use named_array_list, only: qna

     implicit none

     if (qna%exists(psi_n)) call leaves%set_q_value(qna%ind(psi_n), psi_0)

  end subroutine init_psi
  !---------------------------------------------------------------------------------------------------------------

!>
!! \brief calculate vector of fluxes for B_x and psi
!!
!! \details This follows Dedner et al. (2002), Eq. 42.
!! The result is already simplified: flux of B_x = psi_m and flux of psi = c_h^2 * B_x
!<

  subroutine glm_mhd(psi_l, psi_r, b_ccl, b_ccr, b_cc, psif)

    ! external procedures

    use constants,  only: half

    ! arguments

    implicit none

    real, dimension(:), intent(out) :: psif
    real, dimension(:), intent(in) :: psi_l
    real, dimension(:), intent(in) :: psi_r
    real, dimension(:), intent(out) :: b_cc
    real, dimension(:), intent(in) :: b_ccl
    real, dimension(:), intent(in) :: b_ccr

    b_cc = half * ((psi_r + psi_l) - chspeed * (b_ccr - b_ccl))
    psif = half * chspeed * (chspeed * (b_ccl + b_ccr) - (psi_r - psi_l))

  end subroutine glm_mhd

!-----------------------------------------------------------------------------------------------------------

!>
  !! Parabolic damping to psi
!<
  subroutine glmdamping

     use global,           only: cfl, glm_alpha
     use cg_list,          only: cg_list_element
     use cg_leaves,        only: leaves
     use constants,        only: psi_n
     use named_array_list, only: qna

     implicit none

     type(cg_list_element), pointer :: cgl

     if (qna%exists(psi_n)) then
        cgl => leaves%first
        do while (associated(cgl))
           cgl%cg%q(qna%ind(psi_n))%arr =  cgl%cg%q(qna%ind(psi_n))%arr * exp(-glm_alpha*cfl)
           ! for AMR or POLAR it may be better to have explicit dependence on cg%dl(:) here
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

     use cg_leaves,  only: leaves
     use cg_list,    only: cg_list_element
     use constants,  only: xdim, zdim, psi_n, GEO_XYZ, half
     use dataio_pub, only: die
     use div_B,      only: divB, idivB
     use domain,     only: dom
     use fluidindex, only: flind
     use fluids_pub, only: has_ion
     use fluidtypes, only: component_fluid
     use global,     only: use_eglm
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


                       !Sources

                       ! momentum = momentum -divB*B
                       cgl%cg%u(fl%imx:fl%imz,i,j,k) = cgl%cg%u(fl%imx:fl%imz,i,j,k) - cg%q(idivB)%arr(i,j,k)*(cgl%cg%b(xdim:zdim,i,j,k))

                       ! B = B -divB*u
                       cgl%cg%b(xdim:zdim,i,j,k) = cgl%cg%b(xdim:zdim,i,j,k) - cg%q(idivB)%arr(i,j,k)*(cgl%cg%u(fl%imx:fl%imz,i,j,k)/cgl%cg%u(fl%idn,i,j,k))

                       ! e = e - divB*u.B - B.grad(psi)

                       cgl%cg%u(fl%ien,i,j,k) = cgl%cg%u(fl%ien,i,j,k) - &
                            cg%q(idivB)%arr(i,j,k)*dot_product(cgl%cg%u(fl%imx:fl%imz,i,j,k)/cgl%cg%u(fl%idn,i,j,k),cgl%cg%b(xdim:zdim,i,j,k)) - &
                            dot_product(cgl%cg%b(xdim:zdim,i,j,k),cg%w(igp)%arr(xdim:zdim,i,j,k))

                       ! psi = psi - u.grad(psi), other term is calculated in damping
                       cgl%cg%q(qna%ind(psi_n))%arr(i,j,k) =  cgl%cg%q(qna%ind(psi_n))%arr(i,j,k) - &
                            dot_product(cgl%cg%u(fl%imx:fl%imz,i,j,k)/cgl%cg%u(fl%idn,i,j,k),cg%w(igp)%arr(xdim:zdim,i,j,k))

                    end associate
                 enddo
              enddo
           enddo
           cgl=>cgl%nxt
        enddo
     endif

     return
   end subroutine eglm

  subroutine glm_3D

     use all_boundaries,   only: all_mag_boundaries
     use cg_leaves,        only: leaves
     use cg_list,          only: cg_list_element
     use constants,        only: xdim, ydim, zdim, LO, HI, psi_n
     use domain,           only: dom
     use grid_cont,        only: grid_container
     use named_array_list, only: qna
     use global,           only: dt, glm_iter, use_hdc_3D

     implicit none

     integer :: ig, psii
     type(cg_list_element), pointer :: cgl
     type(grid_container),  pointer :: cg
     real, dimension(:,:,:), allocatable :: bbx, bby, bbz, pp

     if (.not. use_hdc_3D) then
        call glmdamping
        return
     endif

     psii = qna%ind(psi_n)

     do ig = 1, glm_iter
        call glmdamping
        call leaves%leaf_arr3d_boundaries(psii)
        call all_mag_boundaries
        cgl => leaves%first
        do while (associated(cgl))
           cg => cgl%cg
           allocate(bbx(cg%lhn(xdim, LO):cg%lhn(xdim, HI), cg%lhn(ydim, LO):cg%lhn(ydim, HI), cg%lhn(zdim, LO):cg%lhn(zdim, HI)), &
                &   bby(cg%lhn(xdim, LO):cg%lhn(xdim, HI), cg%lhn(ydim, LO):cg%lhn(ydim, HI), cg%lhn(zdim, LO):cg%lhn(zdim, HI)), &
                &   bbz(cg%lhn(xdim, LO):cg%lhn(xdim, HI), cg%lhn(ydim, LO):cg%lhn(ydim, HI), cg%lhn(zdim, LO):cg%lhn(zdim, HI)), &
                &    pp(cg%lhn(xdim, LO):cg%lhn(xdim, HI), cg%lhn(ydim, LO):cg%lhn(ydim, HI), cg%lhn(zdim, LO):cg%lhn(zdim, HI)))
           associate ( xs => cg%lhn(xdim, LO)+1, xe => cg%lhn(xdim, HI)-1, &
                &      ys => cg%lhn(ydim, LO)+1, ye => cg%lhn(ydim, HI)-1, &
                &      zs => cg%lhn(zdim, LO)+1, ze => cg%lhn(zdim, HI)-1 )

              pp(:,:,:) = cg%q(psii)%arr(:,:,:)
              if (dom%has_dir(xdim)) then
                 bbx(                  xs  :xe,   :, :) = &
                      &     cg%b(xdim, xs  :xe,   :, :) + dt/cg%dl(xdim) * 0.5 * ( &
                      ( cg%q(psii)%arr(xs-1:xe-1, :, :) - &
                      & cg%q(psii)%arr(xs+1:xe+1, :, :)) - &
                      ( 2*  cg%b(xdim, xs  :xe,   :, :) - &
                      &     cg%b(xdim, xs-1:xe-1, :, :) - &
                      &     cg%b(xdim, xs+1:xe+1, :, :)) * chspeed)

                 bbx(xs-1, :, :) = bbx(xs, :, :)
                 bbx(xe+1, :, :) = bbx(xe, :, :)
              endif
              if (dom%has_dir(xdim)) then
                 pp(                      xs  :xe,   :, :) = &
                      &      pp(          xs  :xe,   :, :) + dt/cg%dl(xdim) * 0.5 * chspeed * ( &
                      (        cg%b(xdim, xs-1:xe-1, :, :) - &
                      &        cg%b(xdim, xs+1:xe+1, :, :)) * chspeed - &
                      ( 2* cg%q(psii)%arr(xs  :xe,   :, :) - &
                      &    cg%q(psii)%arr(xs-1:xe-1, :, :) - &
                      &    cg%q(psii)%arr(xs+1:xe+1, :, :)))

                 pp(xs-1, :, :) = pp(xs, :, :)
                 pp(xe+1, :, :) = pp(xe, :, :)
              endif
              if (dom%has_dir(xdim)) cg%q(psii)%arr(:,:,:) = pp (:,:,:)
              if (dom%has_dir(xdim)) cg%b(xdim, :, :, :) = bbx(:,:,:)

              if (dom%has_dir(ydim)) then
                 bby(                  :, ys  :ye,   :) = &
                      &     cg%b(ydim, :, ys  :ye,   :) + dt/cg%dl(ydim) * 0.5 * ( &
                      ( cg%q(psii)%arr(:, ys-1:ye-1, :) - &
                      & cg%q(psii)%arr(:, ys+1:ye+1, :)) - &
                      ( 2*  cg%b(ydim, :, ys  :ye,   :) - &
                      &     cg%b(ydim, :, ys-1:ye-1, :) - &
                      &     cg%b(ydim, :, ys+1:ye+1, :)) * chspeed)

                 bby(:, ys-1, :) = bby(:, ys, :)
                 bby(:, ye+1, :) = bby(:, ye, :)
              endif
              if (dom%has_dir(ydim)) then
                 pp(                      :, ys  :ye,   :) = &
                      &      pp(          :, ys  :ye,   :) + dt/cg%dl(ydim) * 0.5 * chspeed * ( &
                      (        cg%b(ydim, :, ys-1:ye-1, :) - &
                      &        cg%b(ydim, :, ys+1:ye+1, :)) * chspeed - &
                      ( 2* cg%q(psii)%arr(:, ys  :ye,   :) - &
                      &    cg%q(psii)%arr(:, ys-1:ye-1, :) - &
                      &    cg%q(psii)%arr(:, ys+1:ye+1, :)))

                 pp(:, ys-1, :) = pp(:, ys, :)
                 pp(:, ye+1, :) = pp(:, ye, :)
              endif
              if (dom%has_dir(ydim)) cg%q(psii)%arr(:,:,:) = pp (:,:,:)
              if (dom%has_dir(ydim)) cg%b(ydim, :, :, :) = bby(:,:,:)

              if (dom%has_dir(zdim)) then
                 bbz(                  :, :, zs  :ze  ) = &
                      &     cg%b(zdim, :, :, zs  :ze  ) + dt/cg%dl(zdim) * 0.5 * ( &
                      ( cg%q(psii)%arr(:, :, zs-1:ze-1) - &
                      & cg%q(psii)%arr(:, :, zs+1:ze+1)) - &
                      ( 2*  cg%b(zdim, :, :, zs  :ze  ) - &
                      &     cg%b(zdim, :, :, zs-1:ze-1) - &
                      &     cg%b(zdim, :, :, zs+1:ze+1)) * chspeed)

                 bbz(:, :, zs-1) = bbz(:, :, zs)
                 bbz(:, :, ze+1) = bbz(:, :, ze)
              endif
              if (dom%has_dir(zdim)) then
                 pp(                      :, :, zs  :ze  ) = &
                      pp(                 :, :, zs  :ze  ) + dt/cg%dl(zdim) * 0.5 * chspeed * ( &
                      (        cg%b(zdim, :, :, zs-1:ze-1) - &
                      &        cg%b(zdim, :, :, zs+1:ze+1)) * chspeed - &
                      ( 2* cg%q(psii)%arr(:, :, zs  :ze  ) - &
                      &    cg%q(psii)%arr(:, :, zs-1:ze-1) - &
                      &    cg%q(psii)%arr(:, :, zs+1:ze+1)))

                 pp(:, :, zs-1) = pp(:, :, zs)
                 pp(:, :, ze+1) = pp(:, :, ze)
              endif
              if (dom%has_dir(zdim)) cg%q(psii)%arr(:,:,:) = pp (:,:,:)
              if (dom%has_dir(zdim)) cg%b(zdim, :, :, :) = bbz(:,:,:)

           end associate
           deallocate(bbx, bby, bbz, pp)

           cgl => cgl%nxt
        enddo
     enddo
     call leaves%leaf_arr3d_boundaries(psii)
     call all_mag_boundaries

  end subroutine glm_3D

end module hdc

!--------------------------------------------------------------------------------------------

