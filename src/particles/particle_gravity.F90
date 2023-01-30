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
   public :: update_particle_gravpot_and_acc, phi_pm_part, get_acc_model

contains

   subroutine update_particle_gravpot_and_acc

      use cg_leaves,        only: leaves
      use cg_list,          only: cg_list_element
      use cg_list_dataop,   only: ind_val
      use constants,        only: ndims, gp_n, gpot_n, gp1b_n, sgp_n, nbdn_n, prth_n, one, zero
      use dataio_pub,       only: die
      use domain,           only: is_refined
      use gravity,          only: source_terms_grav
      use grid_cont,        only: grid_container
      use named_array_list, only: qna
      use particle_maps,    only: map_particles
      use particle_pub,     only: mask_gpot1b
      use particle_types,   only: particle
      use particle_utils,   only: global_count_all_particles
#ifdef DROP_OUTSIDE_PART
      use particle_utils,   only: detach_particle
#endif /* DROP_OUTSIDE_PART */
#ifdef VERBOSE
      use dataio_pub,       only: printinfo
#endif /* VERBOSE */

      implicit none

      type(grid_container),  pointer       :: cg
      type(cg_list_element), pointer       :: cgl
      type(particle),        pointer       :: pset

      integer, dimension(ndims)            :: cell
      real,    dimension(ndims)            :: dist
      real                                 :: Mtot
      integer(kind=4)                      :: ig, ip, ib

#ifdef VERBOSE
      call printinfo('[particle_gravity:update_particle_gravpot_and_acc] Commencing update of particle gravpot & acceleration')
#endif /* VERBOSE */

#ifdef NBODY_GRIDDIRECT
      call update_gravpot_from_particles
#endif /* NBODY_GRIDDIRECT */
      call source_terms_grav
      call leaves%q_lin_comb([ ind_val(qna%ind(gp_n), 1.), ind_val(qna%ind(sgp_n), one)   ], qna%ind(gpot_n))

      if (global_count_all_particles() == 0) return

      Mtot = find_Mtot()

      cgl => leaves%first
      do while (associated(cgl))
         cg => cgl%cg
         cg%nbdn = zero
         cg%prth = zero
         cgl => cgl%nxt
      enddo

      if (is_refined) call die("[particle_gravity:update_particle_gravpot_and_acc] AMR not implemented yet")
      ! map_tsc contains a loop over cg and a call to update boundaries
      ! it gives O(#cg^2) cost and funny MPI errors when the number of cg differ from thread to thread
      ig = qna%ind(nbdn_n)
      call map_particles(ig, one)

      ip = qna%ind(prth_n)
      ig = qna%ind(gpot_n)
      ib = ig
      if (mask_gpot1b) ib = qna%ind(gp1b_n)

      cgl => leaves%first
      do while (associated(cgl))
         cg => cgl%cg

         pset => cg%pset%first
         do while (associated(pset))

            call locate_particle_in_cell(cg, pset%pdata%pos, cell, dist)

            call update_particle_density_array(ip, cg, cell, pset%pdata%outside)

            call update_particle_potential_energy(ig, cg, pset, cell, dist, Mtot)

#ifdef DROP_OUTSIDE_PART
            !Delete particles escaping the domain
            if (abs(pset%pdata%energy) < tiny(1.)) then
               call detach_particle(cg, pset)
               cycle
            endif
#endif /* DROP_OUTSIDE_PART */

            call update_particle_acc(ib, cg, pset, cell, dist)

            pset => pset%nxt
         enddo

         cgl => cgl%nxt
      enddo

#ifdef VERBOSE
      call printinfo('[particle_gravity:update_particle_gravpot_and_acc] Finish update of particle gravpot & acceleration')
#endif /* VERBOSE */

   end subroutine update_particle_gravpot_and_acc

   subroutine locate_particle_in_cell(cg, pos, cell, dist)

      use constants, only: ndims, xdim, zdim, LO
      use domain,    only: dom
      use grid_cont, only: grid_container

      implicit none

      type(grid_container),      intent(inout) :: cg
      real,    dimension(ndims), intent(in)    :: pos
      integer, dimension(ndims), intent(out)   :: cell
      real,    dimension(ndims), intent(out)   :: dist
      integer                                  :: cdim

      do cdim = xdim, zdim
         cell(cdim) = floor((pos(cdim) - dom%edge(cdim,LO)) * cg%idl(cdim), kind=4)
         dist(cdim) = pos(cdim) - ( dom%edge(cdim,LO) + cell(cdim) * cg%dl(cdim) )
      enddo

   end subroutine locate_particle_in_cell

   subroutine update_particle_density_array(ig, cg, cell, poutside)

      use constants, only: ndims, one, xdim, ydim, zdim
      use grid_cont, only: grid_container

      implicit none

      integer(kind=4),               intent(in)    :: ig
      type(grid_container), pointer, intent(inout) :: cg
      integer, dimension(ndims),     intent(in)    :: cell
      logical,                       intent(in)    :: poutside

      if (.not. poutside) cg%q(ig)%arr(cell(xdim),cell(ydim),cell(zdim)) = cg%q(ig)%point(cell) + one

   end subroutine update_particle_density_array

   function phi_pm_part(pos, mass)

      use constants,    only: ndims, xdim, ydim, zdim
      use particle_pub, only: r_soft
      use units,        only: newtong

      implicit none

      real, dimension(ndims), intent(in) :: pos
      real,                   intent(in) :: mass
      real                               :: r, phi_pm_part

      r = sqrt(pos(xdim)**2 + pos(ydim)**2 + pos(zdim)**2 + r_soft**2)
      phi_pm_part = -newtong*mass / r

   end function phi_pm_part

   subroutine gravpot1b(pset, cg, ig)

      use constants,      only: CENTER, LO, HI, xdim, ydim, zdim
      use grid_cont,      only: grid_container
      use particle_types, only: particle

      implicit none

      type(particle),       pointer, intent(in)    :: pset
      type(grid_container), pointer, intent(inout) :: cg
      integer(kind=4),               intent(in)    :: ig
      integer                                      :: i, j, k

      do i = cg%lhn(xdim, LO), cg%lhn(xdim, HI)
         do j = cg%lhn(ydim, LO), cg%lhn(ydim, HI)
            do k = cg%lhn(zdim, LO), cg%lhn(zdim, HI)
               cg%q(ig)%arr(i,j,k) = cg%q(ig)%arr(i,j,k) + phi_pm_part([cg%coord(CENTER,xdim)%r(i) - pset%pdata%pos(xdim), &
                                                                        cg%coord(CENTER,ydim)%r(j) - pset%pdata%pos(ydim), &
                                                                        cg%coord(CENTER,zdim)%r(k) - pset%pdata%pos(zdim)], pset%pdata%mass)
            enddo
         enddo
      enddo

   end subroutine gravpot1b

#ifdef NBODY_GRIDDIRECT
   subroutine update_gravpot_from_particles

      use cg_leaves,        only: leaves
      use cg_list,          only: cg_list_element
      use constants,        only: nbgp_n, zero
      use grid_cont,        only: grid_container
      use named_array_list, only: qna
      use particle_types,   only: particle

      implicit none

      type(grid_container),  pointer :: cg
      type(cg_list_element), pointer :: cgl
      type(particle),        pointer :: pset
      integer(kind=4)                :: ig

      ig = qna%ind(nbgp_n)
      cgl => leaves%first
      do while (associated(cgl))
         cg => cgl%cg

            cg%nbgp = zero

            pset => cg%pset%first
            do while (associated(pset))
               call gravpot1b(pset, cg, ig)
               pset => pset%nxt
            enddo

         cgl => cgl%nxt
      enddo

   end subroutine update_gravpot_from_particles
#endif /* NBODY_GRIDDIRECT */

!> \brief Find the total mass of particles

   real function find_Mtot() result(Mtot)

      use cg_leaves,      only: leaves
      use cg_list,        only: cg_list_element
      use constants,      only: pSUM
      use grid_cont,      only: grid_container
      use mpisetup,       only: piernik_MPI_Allreduce
      use particle_types, only: particle

      implicit none

      type(grid_container),  pointer :: cg
      type(cg_list_element), pointer :: cgl
      type(particle),        pointer :: pset

      Mtot = 0

      cgl => leaves%first
      do while (associated(cgl))
         cg => cgl%cg
         pset => cg%pset%first
         do while (associated(pset))
            if (.not. pset%pdata%outside .and. pset%pdata%phy) Mtot = Mtot + pset%pdata%mass
            !TO DO: include gas
            pset => pset%nxt
         enddo
         cgl => cgl%nxt
      enddo

      call piernik_MPI_Allreduce(Mtot, pSUM)

   end function find_Mtot

!> \brief Determine potential energy in particle positions
   subroutine update_particle_potential_energy(ig, cg, pset, cell, dist, Mtot)

      use constants,      only: ndims, half, two, xdim, ydim, zdim
      use grid_cont,      only: grid_container
      use particle_func,  only: df_d_o2, d2f_d2_o2, d2f_dd_o2
      use particle_types, only: particle
      use units,          only: newtong

      implicit none

      integer(kind=4),               intent(in)    :: ig
      type(grid_container), pointer, intent(in)    :: cg
      type(particle), pointer,       intent(inout) :: pset
      integer, dimension(ndims),     intent(in)    :: cell
      real,    dimension(ndims),     intent(in)    :: dist
      real,                          intent(in)    :: Mtot
      real                                         :: dpot, d2pot

      if (.not. pset%pdata%outside) then
         dpot = df_d_o2(cell, cg, ig, xdim) * dist(xdim) + &
                df_d_o2(cell, cg, ig, ydim) * dist(ydim) + &
                df_d_o2(cell, cg, ig, zdim) * dist(zdim)

         d2pot = d2f_d2_o2(cell, cg, ig, xdim) * dist(xdim)**2 + &
                 d2f_d2_o2(cell, cg, ig, ydim) * dist(ydim)**2 + &
                 d2f_d2_o2(cell, cg, ig, zdim) * dist(zdim)**2 + &
           two * d2f_dd_o2(cell, cg, ig, xdim, ydim) * dist(xdim)*dist(ydim) + &
           two * d2f_dd_o2(cell, cg, ig, xdim, zdim) * dist(xdim)*dist(zdim)
         pset%pdata%energy = pset%pdata%mass * (cg%q(ig)%point(cell) + dpot + half * d2pot)
      else
         pset%pdata%energy = -newtong * pset%pdata%mass * Mtot / norm2(pset%pdata%pos(:))
         if ((abs(pset%pdata%energy) < half * pset%pdata%mass * norm2(pset%pdata%vel(:)) **2)) pset%pdata%energy = 0.
      endif

   end subroutine update_particle_potential_energy

   subroutine update_particle_acc(ig, cg, pset, cell, dist)

      use constants,      only: ndims, zero
      use grid_cont,      only: grid_container
      use particle_pub,   only: is_setacc_cic, is_setacc_int, is_setacc_tsc, mask_gpot1b
      use particle_types, only: particle

      implicit none

      integer(kind=4),               intent(in)    :: ig
      type(grid_container), pointer, intent(inout) :: cg
      type(particle),       pointer, intent(in)    :: pset
      integer, dimension(ndims),     intent(in)    :: cell
      real,    dimension(ndims),     intent(in)    :: dist

      if (pset%pdata%outside) then
         pset%pdata%acc = update_particle_acc_exp(cg, pset%pdata%pos)
         return
      endif

      if (mask_gpot1b) then
         cg%gp1b = zero
         call gravpot1b(pset, cg, ig)
         cg%gp1b = -cg%gp1b + cg%gpot
      endif

      if (is_setacc_int) then
         pset%pdata%acc = update_particle_acc_int(ig, cg, cell, dist)
      elseif (is_setacc_cic) then
         pset%pdata%acc = update_particle_acc_cic(ig, cg, pset%pdata%pos, cell)
      elseif (is_setacc_tsc) then
         pset%pdata%acc = update_particle_acc_tsc(ig, cg, pset%pdata%pos)
      endif

   end subroutine update_particle_acc

   function update_particle_acc_int(ig, cg, cell, dist) result (acc)

      use constants,     only: ndims, xdim, ydim, zdim
      use grid_cont,     only: grid_container
      use particle_func, only: df_d_p, d2f_d2_p, d2f_dd_p

      implicit none

      integer(kind=4),               intent(in) :: ig
      type(grid_container), pointer, intent(in) :: cg
      integer, dimension(ndims),     intent(in) :: cell
      real,    dimension(ndims),     intent(in) :: dist
      real,    dimension(ndims)                 :: acc

      acc(xdim) = - (df_d_p(cell, cg, ig, xdim) + &
                    d2f_d2_p(cell, cg, ig, xdim)       * dist(xdim) + &
                    d2f_dd_p(cell, cg, ig, xdim, ydim) * dist(ydim) + &
                    d2f_dd_p(cell, cg, ig, xdim, zdim) * dist(zdim))


      acc(ydim) = -( df_d_p(cell, cg, ig, ydim) + &
                   d2f_d2_p(cell, cg, ig, ydim)       * dist(ydim) + &
                   d2f_dd_p(cell, cg, ig, xdim, ydim) * dist(xdim) + &
                   d2f_dd_p(cell, cg, ig, ydim, zdim) * dist(zdim))

      acc(zdim) = -( df_d_p(cell, cg, ig, zdim) + &
                   d2f_d2_p(cell, cg, ig, zdim)       * dist(zdim) + &
                   d2f_dd_p(cell, cg, ig, xdim, zdim) * dist(xdim) + &
                   d2f_dd_p(cell, cg, ig, ydim, zdim) * dist(ydim))

   end function update_particle_acc_int

   function update_particle_acc_cic(ig, cg, pos, cell) result (acc)

      use constants, only: idm, ndims, CENTER, xdim, ydim, zdim, half, I_ZERO, I_ONE
      use grid_cont, only: grid_container

      implicit none

      integer(kind=4),                    intent(in) :: ig
      type(grid_container), pointer,      intent(in) :: cg
      real,    dimension(ndims),          intent(in) :: pos
      integer, dimension(ndims),          intent(in) :: cell
      real,    dimension(ndims)                      :: acc

      integer(kind=4)                                :: c, cdim
      integer(kind=4), dimension(ndims,8), parameter :: cijk = reshape([[I_ZERO, I_ZERO, I_ZERO], &
           &                                                            [I_ONE,  I_ZERO, I_ZERO], [I_ZERO, I_ONE, I_ZERO], [I_ZERO, I_ZERO, I_ONE], &
           &                                                            [I_ONE,  I_ONE,  I_ZERO], [I_ZERO, I_ONE, I_ONE],  [I_ONE,  I_ZERO, I_ONE], &
           &                                                            [I_ONE,  I_ONE,  I_ONE]], [ndims,8_4])
      integer,         dimension(ndims)              :: cic_cells
      real,            dimension(ndims)              :: dxyz, axyz
      real(kind=8),    dimension(ndims,8)            :: fxyz
      real(kind=8),    dimension(8)                  :: wijk

      do cdim = xdim, zdim
         if (pos(cdim) < cg%coord(CENTER, cdim)%r(cell(cdim))) then
            cic_cells(cdim) = cell(cdim) - 1
         else
            cic_cells(cdim) = cell(cdim)
         endif
         dxyz(cdim) = abs(pos(cdim) - cg%coord(CENTER, cdim)%r(cic_cells(cdim)))
      enddo

      wijk(1) = (cg%dx - dxyz(xdim))*(cg%dy - dxyz(ydim))*(cg%dz - dxyz(zdim))  !a(i  ,j  ,k  )
      wijk(2) =          dxyz(xdim) *(cg%dy - dxyz(ydim))*(cg%dz - dxyz(zdim))  !a(i+1,j  ,k  )
      wijk(3) = (cg%dx - dxyz(xdim))*         dxyz(ydim) *(cg%dz - dxyz(zdim))  !a(i  ,j+1,k  )
      wijk(4) = (cg%dx - dxyz(xdim))*(cg%dy - dxyz(ydim))*         dxyz(zdim)   !a(i  ,j  ,k+1)
      wijk(5) =          dxyz(xdim) *         dxyz(ydim) *(cg%dz - dxyz(zdim))  !a(i+1,j+1,k  )
      wijk(6) = (cg%dx - dxyz(xdim))*         dxyz(ydim) *         dxyz(zdim)   !a(i  ,j+1,k+1)
      wijk(7) =          dxyz(xdim) *(cg%dy - dxyz(ydim))*         dxyz(zdim)   !a(i+1,j  ,k+1)
      wijk(8) =          dxyz(xdim) *         dxyz(ydim) *         dxyz(zdim)   !a(i+1,j+1,k+1)

      wijk = wijk / cg%dvol

      do cdim = xdim, zdim
         do c = 1, 8
            fxyz(cdim,c) = -(cg%q(ig)%point(cic_cells(:)+idm(cdim,:)+cijk(:,c)) - cg%q(ig)%point(cic_cells(:)-idm(cdim,:)+cijk(:,c)))
         enddo
         axyz(cdim) = sum(fxyz(cdim,:)*wijk(:))
      enddo

      acc(:) = half * axyz(:) * cg%idl(:)

   end function update_particle_acc_cic

   function update_particle_acc_tsc(ig, cg, pos) result (acc)

      use constants, only: xdim, ydim, zdim, ndims, LO, HI, IM, I0, IP, CENTER, idm, half, zero
      use domain,    only: dom
      use grid_cont, only: grid_container

      implicit none

      integer(kind=4),               intent(in) :: ig
      type(grid_container), pointer, intent(in) :: cg
      real, dimension(ndims),        intent(in) :: pos
      real, dimension(ndims)                    :: acc
      integer                                   :: i, j, k
      integer(kind=4)                           :: cdim
      integer, dimension(ndims, IM:IP)          :: ijkp
      integer, dimension(ndims)                 :: cur_ind
      real,    dimension(ndims)                 :: fxyz, axyz
      real                                      :: weight, delta_x, weight_tmp

      fxyz = zero ! suppress -Wmaybe-uninitialized
      axyz = zero

      do cdim = xdim, zdim
         if (dom%has_dir(cdim)) then
            ijkp(cdim, I0) = nint((pos(cdim) - cg%fbnd(cdim,LO)-cg%dl(cdim)/2.) * cg%idl(cdim) + int(cg%lhn(cdim, LO)) + dom%nb, kind=4)
            ijkp(cdim, IM) = max(ijkp(cdim, I0) - 1, int(cg%lhn(cdim, LO)))
            ijkp(cdim, IP) = min(ijkp(cdim, I0) + 1, int(cg%lhn(cdim, HI)))
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
                  delta_x = ( pos(cdim) - cg%coord(CENTER, cdim)%r(cur_ind(cdim)) ) * cg%idl(cdim)

                  if (cur_ind(cdim) /= ijkp(cdim, 0)) then   !!! BEWARE hardcoded magic
                     weight_tmp = 1.125 - 1.5 * abs(delta_x) + half * delta_x**2
                  else
                     weight_tmp = 0.75 - delta_x**2
                  endif

                  weight = weight_tmp * weight
                  fxyz(cdim) = -(cg%q(ig)%point(cur_ind(:) + idm(cdim,:)) - cg%q(ig)%point(cur_ind(:) - idm(cdim,:)))
               enddo
               axyz(:) = axyz(:) + fxyz(:) * weight
            enddo
         enddo
      enddo
      acc(:) = half * axyz(:) * cg%idl(:)

   end function update_particle_acc_tsc

!>
! \brief multipole expansion for particles outside domain
!<
   function update_particle_acc_exp(cg, pos) result (acc)

      use constants, only: ndims, xdim, ydim, zdim, half
      use grid_cont, only: grid_container
      use multipole, only: moments2pot

      implicit none

      type(grid_container), pointer, intent(in) :: cg
      real, dimension(ndims),        intent(in) :: pos
      real, dimension(ndims)                    :: acc
      integer(kind=4)                           :: cdim
      real                                      :: axyz
      real, dimension(5)                        :: tmp

      tmp = [0, 0, 1, 0, 0]

      do cdim = xdim, zdim
         axyz = - ( moments2pot( pos(xdim) + cg%dl(xdim) * tmp(cdim+2), pos(ydim) + cg%dl(ydim) * tmp(cdim+1), pos(zdim) + cg%dl(zdim) * tmp(cdim) ) - &
                    moments2pot( pos(xdim) - cg%dl(xdim) * tmp(cdim+2), pos(ydim) - cg%dl(ydim) * tmp(cdim+1), pos(zdim) - cg%dl(zdim) * tmp(cdim) ) )
         acc(cdim) = half * axyz * cg%idl(cdim)
      enddo

   end function update_particle_acc_exp

   subroutine get_acc_model(pos, mass, acc2)

      use constants, only: idm, ndims, two, xdim, zdim

      implicit none

      real, dimension(ndims), intent(in)  :: pos
      real,                   intent(in)  :: mass
      real, dimension(ndims), intent(out) :: acc2
      real                                :: d
      integer(kind=4)                     :: dir

      d = 1.0e-8
      do dir = xdim, zdim
         acc2(dir) = -( phi_pm_part(pos + real(idm(dir,:)) * d, mass) - phi_pm_part(pos - real(idm(dir,:)) * d, mass) ) / (two * d)
      enddo

   end subroutine get_acc_model

#if 0
   ! An unused function
   ! It may be useful for debuging

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

#endif /* 0 */

end module particle_gravity
