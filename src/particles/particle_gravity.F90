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

!>  \brief Auxiliary functions used by particle integrators

module particle_gravity
! pulled by NBODY

   implicit none

   private
   public :: update_particle_gravpot_and_acc, is_setacc_cic, is_setacc_int, is_setacc_tsc, mask_gpot1b, phi_pm_part, get_acc_model

   logical :: is_setacc_cic, is_setacc_int, is_setacc_tsc, mask_gpot1b


contains

   subroutine update_particle_gravpot_and_acc

      use constants,      only: ndims
      use cg_leaves,      only: leaves
      use cg_list,        only: cg_list_element
      use dataio_pub,     only: die
      use domain,         only: is_refined, is_multicg
      use gravity,        only: source_terms_grav
      use grid_cont,      only: grid_container
      use particle_types, only: pset
#ifdef NBODY_GRIDDIRECT
      use constants,      only: zero
#endif /* NBODY_GRIDDIRECT */
#ifdef VERBOSE
      use dataio_pub,     only: printinfo
#endif /* VERBOSE */

      implicit none

      type(grid_container),  pointer       :: cg
      type(cg_list_element), pointer       :: cgl

      integer                              :: n_part
      real,    dimension(:,:), allocatable :: dist
      integer, dimension(:,:), allocatable :: cells

#ifdef VERBOSE
      call printinfo('[particle_gravity:update_particle_gravpot_and_acc] Commencing update of particle gravpot & acceleration')
#endif /* VERBOSE */

      n_part = size(pset%p, dim=1)

      if (is_refined) call die("[particle_gravity:update_particle_gravpot_and_acc] AMR not implemented for particles yet")
      if (is_multicg) call die("[particle_gravity:update_particle_gravpot_and_acc] multi_cg not implemented for particles yet")
      cgl => leaves%first
      cg  => cgl%cg

#ifdef NBODY_GRIDDIRECT
      call update_gravpot_from_particles(n_part, cg, zero)
#endif /* NBODY_GRIDDIRECT */

      call source_terms_grav

      allocate(cells(n_part, ndims), dist(n_part, ndims))
      call locate_particles_in_cells(n_part, cg, cells, dist)

      call update_particle_density_array(n_part, cg, cells)

      call update_particle_potential_energy(n_part, cg, cells, dist)

      if (is_setacc_int) then
         call update_particle_acc_int(n_part, cg, cells, dist)
      elseif (is_setacc_cic) then
         call update_particle_acc_cic(n_part, cg, cells)
      elseif (is_setacc_tsc) then
         call update_particle_acc_tsc(cg)
      endif
      deallocate(cells, dist)

#ifdef VERBOSE
      call printinfo('[particle_gravity:update_particle_gravpot_and_acc] Finish update of particle gravpot & acceleration')
#endif /* VERBOSE */

   end subroutine update_particle_gravpot_and_acc

   subroutine update_particle_density_array(n_part, cg, cells)

      use constants,        only: nbdn_n, prth_n, ndims, one, zero, xdim, ydim, zdim
      use grid_cont,        only: grid_container
      use named_array_list, only: qna
      use particle_types,   only: pset

      implicit none

      integer,                          intent(in)    :: n_part
      type(grid_container), pointer,    intent(inout) :: cg
      integer, dimension(n_part,ndims), intent(in)    :: cells
      integer(kind=4)                                 :: ig
      integer                                         :: p
      real                                            :: factor

      factor = one
      ig = qna%ind(nbdn_n)
      cg%nbdn = zero
      call pset%map(ig,factor)

      ig = qna%ind(prth_n)
      cg%prth = zero
      do p = 1, n_part
         cg%q(ig)%arr(cells(p,xdim),cells(p,ydim),cells(p,zdim)) = cg%q(ig)%point(cells(p,:)) + one
      enddo

   end subroutine update_particle_density_array

   function phi_pm_part(pos, eps, mass)

      use constants, only: ndims, xdim, ydim, zdim
      use units,     only: newtong

      implicit none

      real, dimension(ndims), intent(in) :: pos
      real,                   intent(in) :: eps, mass
      real                               :: r, phi_pm_part

      r = sqrt(pos(xdim)**2 + pos(ydim)**2 + pos(zdim)**2 + eps**2)
      phi_pm_part = -newtong*mass / r

   end function phi_pm_part

   subroutine gravpot1b(p, cg, ig, eps2)

      use constants,      only: CENTER, LO, HI, xdim, ydim, zdim
      use grid_cont,      only: grid_container
      use particle_types, only: pset

      implicit none

      integer,                       intent(in)    :: p
      type(grid_container), pointer, intent(inout) :: cg
      integer(kind=4),               intent(in)    :: ig
      real,                          intent(in)    :: eps2
      integer                                      :: i, j, k

      do i = cg%lhn(xdim, LO), cg%lhn(xdim, HI)
         do j = cg%lhn(ydim, LO), cg%lhn(ydim, HI)
            do k = cg%lhn(zdim, LO), cg%lhn(zdim, HI)
               cg%q(ig)%arr(i,j,k) = cg%q(ig)%arr(i,j,k) + phi_pm_part([cg%coord(CENTER,xdim)%r(i) - pset%p(p)%pos(xdim), &
                                                                        cg%coord(CENTER,ydim)%r(j) - pset%p(p)%pos(ydim), &
                                                                        cg%coord(CENTER,zdim)%r(k) - pset%p(p)%pos(zdim)], eps2, pset%p(p)%mass)
            enddo
         enddo
      enddo

   end subroutine gravpot1b

#ifdef NBODY_GRIDDIRECT
   subroutine update_gravpot_from_particles(n_part, cg, eps2)

      use constants,        only: nbgp_n, zero
      use grid_cont,        only: grid_container
      use named_array_list, only: qna

      implicit none

      integer,                       intent(in)    :: n_part
      type(grid_container), pointer, intent(inout) :: cg
      real,                          intent(in)    :: eps2
      integer                                      :: p

      cg%nbgp = zero

      do p = 1, n_part
         call gravpot1b(p, cg, qna%ind(nbgp_n), eps2)
      enddo

   end subroutine update_gravpot_from_particles
#endif /* NBODY_GRIDDIRECT */

   subroutine locate_particles_in_cells(n_part, cg, cells, dist)

      use constants,      only: half, ndims, xdim, zdim, CENTER, LO, HI
      use domain,         only: dom
      use grid_cont,      only: grid_container
      use particle_types, only: pset

      implicit none

      integer,                          intent(in)  :: n_part
      type(grid_container),             intent(in)  :: cg
      integer, dimension(n_part,ndims), intent(out) :: cells
      real,    dimension(n_part,ndims), intent(out) :: dist
      integer                                       :: i, cdim

      do i = 1, n_part
         do cdim = xdim, zdim
            if ((pset%p(i)%pos(cdim) >= dom%edge(cdim, LO)) .and. (pset%p(i)%pos(cdim) <= dom%edge(cdim, HI))) then
               pset%p(i)%outside = .false.
            else
               pset%p(i)%outside = .true.
            endif

            cells(i, cdim) = int( half + (pset%p(i)%pos(cdim) - cg%coord(CENTER,cdim)%r(cg%ijkse(cdim, LO))) * cg%idl(cdim) )

            dist(i, cdim)  = pset%p(i)%pos(cdim) - ( cg%coord(CENTER, cdim)%r(cg%ijkse(cdim, LO)) + cells(i,cdim) * cg%dl(cdim) )
         enddo
      enddo

   end subroutine locate_particles_in_cells

!> \brief Determine potential energy in particle positions
    subroutine update_particle_potential_energy(n_part, cg, cells, dist)

      use constants,        only: gpot_n, ndims, half, two, xdim, ydim, zdim
      use grid_cont,        only: grid_container
      use named_array_list, only: qna
      use particle_func,    only: df_d_o2, d2f_d2_o2, d2f_dd_o2
      use particle_types,   only: pset

      implicit none

      integer,                          intent(in) :: n_part
      type(grid_container), pointer,    intent(in) :: cg
      integer, dimension(n_part,ndims), intent(in) :: cells
      real,    dimension(n_part,ndims), intent(in) :: dist
      integer                                      :: p
      integer(kind=4)                              :: ig
      real, dimension(n_part)                      :: dpot, d2pot

      ig = qna%ind(gpot_n)
      do p = 1, n_part
         dpot(p) = df_d_o2([cells(p, :)], cg, ig, xdim) * dist(p, xdim) + &
                   df_d_o2([cells(p, :)], cg, ig, ydim) * dist(p, ydim) + &
                   df_d_o2([cells(p, :)], cg, ig, zdim) * dist(p, zdim)

         d2pot(p) = d2f_d2_o2([cells(p, :)], cg, ig, xdim) * dist(p, xdim)**2 + &
                    d2f_d2_o2([cells(p, :)], cg, ig, ydim) * dist(p, ydim)**2 + &
                    d2f_d2_o2([cells(p, :)], cg, ig, zdim) * dist(p, zdim)**2 + &
                two*d2f_dd_o2([cells(p, :)], cg, ig, xdim, ydim) * dist(p, xdim)*dist(p, ydim) + &
                two*d2f_dd_o2([cells(p, :)], cg, ig, xdim, zdim) * dist(p, xdim)*dist(p, zdim)
         pset%p(p)%energy = cg%q(ig)%point(cells(p,:)) + dpot(p) + half * d2pot(p)
      enddo

   end subroutine update_particle_potential_energy

   subroutine update_particle_acc_int(n_part, cg, cells, dist)

      use constants,        only: gpot_n, ndims, xdim, ydim, zdim
      use grid_cont,        only: grid_container
      use named_array_list, only: qna
      use particle_func,    only: df_d_p, d2f_d2_p, d2f_dd_p
      use particle_types,   only: pset

      implicit none

      integer,                          intent(in) :: n_part
      type(grid_container), pointer,    intent(in) :: cg
      integer, dimension(n_part,ndims), intent(in) :: cells
      real, dimension(n_part, ndims),   intent(in) :: dist
      integer                                      :: i
      integer(kind=4)                              :: ig

      ig = qna%ind(gpot_n)
      do i = 1, n_part
         if ( (pset%p(i)%outside) .eqv. .false.) then
            pset%p(i)%acc(xdim) = - (df_d_p([cells(i, :)], cg, ig, xdim) + &
                                   d2f_d2_p([cells(i, :)], cg, ig, xdim)       * dist(i, xdim) + &
                                   d2f_dd_p([cells(i, :)], cg, ig, xdim, ydim) * dist(i, ydim) + &
                                   d2f_dd_p([cells(i, :)], cg, ig, xdim, zdim) * dist(i, zdim))

            pset%p(i)%acc(ydim) = -( df_d_p([cells(i, :)], cg, ig, ydim) + &
                                   d2f_d2_p([cells(i, :)], cg, ig, ydim)       * dist(i, ydim) + &
                                   d2f_dd_p([cells(i, :)], cg, ig, xdim, ydim) * dist(i, xdim) + &
                                   d2f_dd_p([cells(i, :)], cg, ig, ydim, zdim) * dist(i, zdim))

            pset%p(i)%acc(zdim) = -( df_d_p([cells(i, :)], cg, ig, zdim) + &
                                   d2f_d2_p([cells(i, :)], cg, ig, zdim)       * dist(i, zdim) + &
                                   d2f_dd_p([cells(i, :)], cg, ig, xdim, zdim) * dist(i, xdim) + &
                                   d2f_dd_p([cells(i, :)], cg, ig, ydim, zdim) * dist(i, ydim))
         !else
         !   call !funkcja liczaca pochodne z potencjalu policzonego z rozwiniecia multipolowego
         endif
      enddo

   end subroutine update_particle_acc_int

   subroutine update_particle_acc_cic(n_part, cg, cells)

      use constants,        only: idm, ndims, CENTER, xdim, ydim, zdim, half, zero, gp1b_n, gpot_n
      use grid_cont,        only: grid_container
      use named_array_list, only: qna
      use particle_types,   only: pset

      implicit none

      integer,                          intent(in)    :: n_part
      type(grid_container), pointer,    intent(inout) :: cg
      integer, dimension(n_part,ndims), intent(in)    :: cells

      integer                                         :: p
      integer(kind=4)                                 :: c, ig, cdim
      integer(kind=4), dimension(ndims,8), parameter  :: cijk = reshape([[0,0,0], [1,0,0], [0,1,0], [0,0,1], [1,1,0], [0,1,1], [1,0,1], [1,1,1]], [ndims,8])
      integer,         dimension(ndims)               :: cic_cells
      real,            dimension(ndims)               :: dxyz, axyz
      real(kind=8),    dimension(ndims,8)             :: fxyz
      real(kind=8),    dimension(8)                   :: wijk

      if (mask_gpot1b) then
         ig = qna%ind(gp1b_n)
      else
         ig = qna%ind(gpot_n)
      endif

      do p = 1, n_part
         if ((pset%p(p)%outside) .eqv. .false.) then
            do cdim = xdim, ndims
               if (pset%p(p)%pos(cdim) < cg%coord(CENTER, cdim)%r(cells(p,cdim))) then
                  cic_cells(cdim) = cells(p, cdim) - 1
               else
                  cic_cells(cdim) = cells(p, cdim)
               endif
               dxyz(cdim) = abs(pset%p(p)%pos(cdim) - cg%coord(CENTER, cdim)%r(cic_cells(cdim)))

            enddo

            wijk(1) = (cg%dx - dxyz(xdim))*(cg%dy - dxyz(ydim))*(cg%dz - dxyz(zdim))  !a(i  ,j  ,k  )
            wijk(2) =          dxyz(xdim) *(cg%dy - dxyz(ydim))*(cg%dz - dxyz(zdim))  !a(i+1,j  ,k  )
            wijk(3) = (cg%dx - dxyz(xdim))*         dxyz(ydim) *(cg%dz - dxyz(zdim))  !a(i  ,j+1,k  )
            wijk(4) = (cg%dx - dxyz(xdim))*(cg%dy - dxyz(ydim))*         dxyz(zdim)   !a(i  ,j  ,k+1)
            wijk(5) =          dxyz(xdim) *         dxyz(ydim) *(cg%dz - dxyz(zdim))  !a(i+1,j+1,k  )
            wijk(6) = (cg%dx - dxyz(xdim))*         dxyz(ydim) *         dxyz(zdim)   !a(i  ,j+1,k+1)
            wijk(7) =          dxyz(xdim) *(cg%dy - dxyz(ydim))*         dxyz(zdim)   !a(i+1,j  ,k+1)
            wijk(8) =          dxyz(xdim) *         dxyz(ydim) *         dxyz(zdim)   !a(i+1,j+1,k+1)
         !else multipole expansion for particles outside domain
         endif

         wijk = wijk/cg%dvol

         if (mask_gpot1b) then
            cg%gp1b = zero
            call gravpot1b(p, cg, ig, zero)
            cg%gp1b = -cg%gp1b + cg%gpot
         endif

         do cdim = xdim, zdim
            do c = 1, 8
               fxyz(cdim,c) = -(cg%q(ig)%point(cic_cells(:)+idm(cdim,:)+cijk(:,c)) - cg%q(ig)%point(cic_cells(:)-idm(cdim,:)+cijk(:,c)))
            enddo
            axyz(cdim) = sum(fxyz(cdim,:)*wijk(:))
         enddo
         pset%p(p)%acc(:) = half*axyz(:)*cg%idl(:)
      enddo

   end subroutine update_particle_acc_cic

   subroutine update_particle_acc_tsc(cg)

      use constants,        only: xdim, ydim, zdim, ndims, LO, HI, IM, I0, IP, CENTER, gp1b_n, gpot_n, idm, half, zero
      use domain,           only: dom
      use grid_cont,        only: grid_container
      use named_array_list, only: qna
      use particle_types,   only: pset

      implicit none


      type(grid_container), pointer, intent(inout) :: cg
      integer                                      :: p, cdim
      integer(kind=4)                              :: ig, i, j, k
      integer(kind=4), dimension(ndims, IM:IP)     :: ijkp
      integer(kind=4), dimension(ndims)            :: cur_ind, ind1, ind2
      real(kind=8),    dimension(ndims)            :: fxyz, axyz
      real                                         :: weight, delta_x, weight_tmp

      if (mask_gpot1b) then
         ig = qna%ind(gp1b_n)
      else
         ig = qna%ind(gpot_n)
      endif

      axyz(:) = 0.0

      do p = lbound(pset%p, dim=1), ubound(pset%p, dim=1)
         associate( part  => pset%p(p), &
                    idl   => cg%idl )

            if (any(part%pos < cg%fbnd(:,LO)) .or. any(part%pos > cg%fbnd(:,HI))) cycle

            if (mask_gpot1b) then
               cg%gp1b = zero
               call gravpot1b(p, cg, ig, zero)
               cg%gp1b = -cg%gp1b + cg%gpot
            endif

            do cdim = xdim, zdim
               if (dom%has_dir(cdim)) then
                  ijkp(cdim, I0) = nint((part%pos(cdim) - cg%coord(CENTER, cdim)%r(1))*cg%idl(cdim)) + 1   !!! BEWARE hardcoded magic
                  ijkp(cdim, IM) = max(ijkp(cdim, I0) - 1, cg%lhn(cdim, LO))
                  ijkp(cdim, IP) = min(ijkp(cdim, I0) + 1, cg%lhn(cdim, HI))
               else
                  ijkp(cdim, IM) = cg%ijkse(cdim, LO)
                  ijkp(cdim, I0) = cg%ijkse(cdim, LO)
                  ijkp(cdim, IP) = cg%ijkse(cdim, HI)
               endif
            enddo
            do i = ijkp(xdim, IM), ijkp(xdim, IP)
               do j = ijkp(ydim, IM), ijkp(ydim, IP)
                  do k = ijkp(zdim, IM), ijkp(zdim, IP)

                     cur_ind(:) = [i, j, k]
                     weight = 1.0

                     do cdim = xdim, zdim
                        if (.not.dom%has_dir(cdim)) cycle
                        delta_x = ( part%pos(cdim) - cg%coord(CENTER, cdim)%r(cur_ind(cdim)) ) * idl(cdim)

                        if (cur_ind(cdim) /= ijkp(cdim, 0)) then   !!! BEWARE hardcoded magic
                           weight_tmp = 1.125 - 1.5 * abs(delta_x) + half * delta_x**2
                        else
                           weight_tmp = 0.75 - delta_x**2
                        endif

                        weight = weight_tmp * weight
                        ind1 = cur_ind(:) + idm(cdim,:)
                        ind2 = cur_ind(:) - idm(cdim,:)
                        fxyz(cdim) = -half*(cg%q(ig)%point(ind1(:)) - cg%q(ig)%point(ind2(:)))*idl(cdim)
                     enddo
                     axyz(:) = axyz(:) + fxyz(:)*weight

                  enddo
               enddo
            enddo
         end associate

         pset%p(p)%acc(:) = axyz(:)

         axyz(:) = 0.0

      enddo

   end subroutine update_particle_acc_tsc

   subroutine get_acc_model(p, eps, acc2)

      use constants,      only: ndims, xdim, zdim
      use particle_types, only: pset

      implicit none

      integer,                intent(in)  :: p
      real,                   intent(in)  :: eps
      real, dimension(ndims), intent(out) :: acc2
      integer(kind=4)                     :: dir

      do dir = xdim, zdim
         acc2(dir) = -der_xyz(pset%p(p)%pos, 1.0e-8, eps, dir, pset%p(p)%mass)
      enddo

   end subroutine get_acc_model

   function der_xyz(pos, d, eps, dir, mass)

      use constants, only: idm, ndims, two

      implicit none

      real(kind=8),dimension(1,ndims), intent(in) :: pos
      real(kind=8),                    intent(in) :: d, eps
      integer(kind=4),                 intent(in) :: dir
      real,                            intent(in) :: mass
      real(kind=8)                                :: der_xyz

      der_xyz = ( phi_pm_part(pos(1,:)+real(idm(dir,:))*d, eps, mass) - phi_pm_part(pos(1,:)-real(idm(dir,:))*d, eps, mass) ) / (two*d)

   end function der_xyz

   subroutine direct_nbody_acc(mass, pos, acc, n, epot)

      use constants, only: ndims, zero
      use units,     only: newtong

      implicit none

      integer,                  intent(in)  :: n
      real, dimension(n),       intent(in)  :: mass
      real, dimension(n,ndims), intent(in)  :: pos
      real, dimension(n,ndims), intent(out) :: acc
      real,                     intent(out) :: epot

      integer                               :: i, j
      real, dimension(ndims)                :: rji, da
      real                                  :: r   ! | rji |
      real                                  :: r2  ! | rji |^2
      real                                  :: r3  ! | rji |^3

      acc(:,:) = zero
      epot     = zero

      do i = 1, n
         do j = i+1, n
            rji(:) = pos(j, :) - pos(i, :)

            r2 = sum(rji**2)
            r = sqrt(r2)
            r3 = r * r2

            ! add the {i,j} contribution to the total potential energy for the system
            epot = epot - newtong * mass(i) * mass(j) / r

            da(:) = newtong * rji(:) / r3

            acc(i,:) = acc(i,:) + mass(j) * da(:)
            acc(j,:) = acc(j,:) - mass(i) * da(:)
         enddo
      enddo

   end subroutine direct_nbody_acc

end module particle_gravity
