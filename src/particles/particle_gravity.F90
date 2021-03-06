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
   public :: update_particle_gravpot_and_acc, is_setacc_cic, is_setacc_int, is_setacc_tsc, mask_gpot1b, phi_pm_part, get_acc_model, eps

   logical :: is_setacc_cic, is_setacc_int, is_setacc_tsc, mask_gpot1b
   real    :: eps

contains

   subroutine update_particle_gravpot_and_acc

      use cg_leaves,        only: leaves
      use cg_list,          only: cg_list_element
      use cg_list_dataop,   only: ind_val
      use constants,        only: ndims, gp_n, gpot_n, sgp_n, one, pSUM
      use gravity,          only: source_terms_grav
      use grid_cont,        only: grid_container
      use named_array_list, only: qna
      use mpisetup,         only: piernik_MPI_Allreduce
      use particle_types,   only: particle

#ifdef VERBOSE
      use dataio_pub,       only: printinfo
#endif /* VERBOSE */

      implicit none

      type(grid_container),  pointer       :: cg
      type(cg_list_element), pointer       :: cgl
      type(particle), pointer              :: pset

      integer                              :: n_part, g_np !, k
      real,    dimension(:,:), allocatable :: dist
      integer, dimension(:,:), allocatable :: cells
      !integer, dimension(:),   allocatable :: pdel
      real :: Mtot

#ifdef VERBOSE
      call printinfo('[particle_gravity:update_particle_gravpot_and_acc] Commencing update of particle gravpot & acceleration')
#endif /* VERBOSE */

#ifdef NBODY_GRIDDIRECT
      call update_gravpot_from_particles
#endif /* NBODY_GRIDDIRECT */
      call source_terms_grav
      call leaves%q_lin_comb([ ind_val(qna%ind(gp_n), 1.), ind_val(qna%ind(sgp_n), one)   ], qna%ind(gpot_n))

      n_part = 0
      cgl => leaves%first
      do while (associated(cgl))
         pset => cgl%cg%pset%first
         do while (associated(pset))
            n_part = n_part + 1
            pset => pset%nxt
         enddo
         cgl => cgl%nxt
      enddo

      g_np = n_part
      call piernik_MPI_Allreduce(g_np, pSUM)

      if (g_np == 0) return

      Mtot = find_Mtot()

      cgl => leaves%first
      do while (associated(cgl))
         cg => cgl%cg

            allocate(cells(n_part, ndims), dist(n_part, ndims))
            call locate_particles_in_cells(n_part, cg, cells, dist)

            call update_particle_density_array(n_part, cg, cells)

            call update_particle_potential_energy(n_part, cg, cells, dist, Mtot)!, pdel)

            call update_particle_acc(n_part, cg, cells, dist)
            !Delete particles elsewhere > in update_particle_acc?
            !pset => cg%pset%first
            !do while (associated(pset))
            !   if (pdel(p)==1) then
            !      call pset%remove(k)
            !   endif
            !enddo

            deallocate(cells, dist)

         cgl => cgl%nxt
      enddo

#ifdef VERBOSE
      call printinfo('[particle_gravity:update_particle_gravpot_and_acc] Finish update of particle gravpot & acceleration')
#endif /* VERBOSE */

   end subroutine update_particle_gravpot_and_acc

   subroutine update_particle_density_array(n_part, cg, cells)

      use constants,        only: nbdn_n, prth_n, ndims, one, zero, xdim, ydim, zdim
      use dataio_pub,       only: die
      use domain,           only: is_multicg
      use grid_cont,        only: grid_container
      use named_array_list, only: qna
      use particle_maps,    only: map_particles
      use particle_types,   only: particle

      implicit none

      integer,                          intent(in)    :: n_part
      type(grid_container), pointer,    intent(inout) :: cg
      type(particle), pointer                         :: pset
      integer, dimension(n_part,ndims), intent(in)    :: cells
      integer(kind=4)                                 :: ig
      integer                                         :: p
      real                                            :: factor

      factor = one
      ig = qna%ind(nbdn_n)
      cg%nbdn = zero

      if (is_multicg) call die("[particle_gravity:update_particle_density_array] multicg not implemented yet")
      ! map_tsc contains a loop over cg and a call to update boundaries
      ! it gives O(#cg^2) cost and funny MPI errors when the number of cg differ from thread to thread
      call map_particles(ig,factor)

      ig = qna%ind(prth_n)
      cg%prth = zero
      p=1
      pset => cg%pset%first
      do while (associated(pset))
         if (.not. pset%pdata%outside) cg%q(ig)%arr(cells(p,xdim),cells(p,ydim),cells(p,zdim)) = cg%q(ig)%point(cells(p,:)) + one
         pset => pset%nxt
         p = p + 1
      enddo

   end subroutine update_particle_density_array

   function phi_pm_part(pos, mass)

      use constants, only: ndims, xdim, ydim, zdim
      use units,     only: newtong

      implicit none

      real, dimension(ndims), intent(in) :: pos
      real,                   intent(in) :: mass
      real                               :: r, phi_pm_part

      r = sqrt(pos(xdim)**2 + pos(ydim)**2 + pos(zdim)**2 + eps**2)
      phi_pm_part = -newtong*mass / r

   end function phi_pm_part

   subroutine gravpot1b(pset, cg, ig)

      use constants, only: CENTER, LO, HI, xdim, ydim, zdim
      use grid_cont, only: grid_container
      use particle_types, only: particle

      implicit none

      !integer,                       intent(in)    :: p
      type(particle), pointer, intent(in)          :: pset
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
      use particle_types, only: particle

      implicit none

      type(grid_container), pointer  :: cg
      type(cg_list_element), pointer :: cgl
      type(particle), pointer        :: pset
      integer                        :: p

      cgl => leaves%first
      do while (associated(cgl))
         cg => cgl%cg

            cg%nbgp = zero

            pset => cg%pset%first
            do while (associated(pset))
               call gravpot1b(pset, cg, qna%ind(nbgp_n))
               pset => pset%nxt
            enddo

         cgl => cgl%nxt
      enddo

   end subroutine update_gravpot_from_particles
#endif /* NBODY_GRIDDIRECT */

   subroutine locate_particles_in_cells(n_part, cg, cells, dist)

      use constants,      only: ndims, xdim, zdim, LO
      use domain,         only: dom
      use grid_cont,      only: grid_container
      use particle_types, only: particle

      implicit none

      integer,                          intent(in)    :: n_part
      type(grid_container),             intent(inout) :: cg
      integer, dimension(n_part,ndims), intent(out)   :: cells
      real,    dimension(n_part,ndims), intent(out)   :: dist

      type(particle), pointer                         :: pset
      integer                                         :: i, cdim

      i=1
      pset => cg%pset%first
      do while (associated(pset))
         call pset%pdata%is_outside()
         do cdim = xdim, zdim
            cells(i, cdim) = floor((pset%pdata%pos(cdim) - dom%edge(cdim,LO)) * cg%idl(cdim), kind=4)
            dist(i, cdim)  = pset%pdata%pos(cdim) - ( dom%edge(cdim,LO) + cells(i,cdim) * cg%dl(cdim) )
         enddo
         pset => pset%nxt
         i = i + 1
      enddo

   end subroutine locate_particles_in_cells

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
            if ((.not. pset%pdata%outside) .and. (pset%pdata%phy)) &
                 Mtot = Mtot + pset%pdata%mass
            !TO DO: include gas
            pset => pset%nxt
         enddo
         cgl => cgl%nxt
      enddo

      call piernik_MPI_Allreduce(Mtot, pSUM)

   end function find_Mtot

!> \brief Determine potential energy in particle positions
   subroutine update_particle_potential_energy(n_part, cg, cells, dist, Mtot)!, pdel)

      use constants,        only: gpot_n, ndims, half, two, xdim, ydim, zdim
      use grid_cont,        only: grid_container
      use named_array_list, only: qna
      use particle_func,    only: df_d_o2, d2f_d2_o2, d2f_dd_o2
      use particle_types,   only: particle
      use units,            only: newtong

      implicit none

      integer,                          intent(in) :: n_part
      type(grid_container), pointer,    intent(in) :: cg
      integer, dimension(n_part,ndims), intent(in) :: cells
      real,    dimension(n_part,ndims), intent(in) :: dist
      real,                             intent(in) :: Mtot

      type(particle), pointer                      :: pset
      integer                                      :: p
      !integer, dimension(n_part),       intent(out):: pdel
      integer(kind=4)                              :: ig
      real, dimension(n_part)                      :: dpot, d2pot

      ig = qna%ind(gpot_n)

      pset => cg%pset%first
      p=1
      do while (associated(pset))
         !pdel(p) = 0
         if (.not. pset%pdata%outside) then
            dpot(p) = df_d_o2([cells(p, :)], cg, ig, xdim) * dist(p, xdim) + &
                      df_d_o2([cells(p, :)], cg, ig, ydim) * dist(p, ydim) + &
                      df_d_o2([cells(p, :)], cg, ig, zdim) * dist(p, zdim)

            d2pot(p) = d2f_d2_o2([cells(p, :)], cg, ig, xdim) * dist(p, xdim)**2 + &
                       d2f_d2_o2([cells(p, :)], cg, ig, ydim) * dist(p, ydim)**2 + &
                       d2f_d2_o2([cells(p, :)], cg, ig, zdim) * dist(p, zdim)**2 + &
                 two * d2f_dd_o2([cells(p, :)], cg, ig, xdim, ydim) * dist(p, xdim)*dist(p, ydim) + &
                 two * d2f_dd_o2([cells(p, :)], cg, ig, xdim, zdim) * dist(p, xdim)*dist(p, zdim)
            pset%pdata%energy = pset%pdata%mass * (cg%q(ig)%point(cells(p,:)) + dpot(p) + half * d2pot(p))
         else
            pset%pdata%energy = -newtong * pset%pdata%mass * Mtot / norm2(pset%pdata%pos(:))
            if ((abs(pset%pdata%energy) < half * pset%pdata%mass * norm2(pset%pdata%vel(:)) **2)) then
               !pdel(p) = 1
               pset%pdata%energy = 0.
            endif

         endif
         pset => pset%nxt
         p = p + 1
      enddo

   end subroutine update_particle_potential_energy

   subroutine update_particle_acc(n_part, cg, cells, dist)

      use constants, only: ndims
      use grid_cont, only: grid_container
      implicit none

      integer,                          intent(in)    :: n_part
      type(grid_container), pointer,    intent(inout) :: cg
      integer, dimension(n_part,ndims), intent(in)    :: cells
      real, dimension(n_part, ndims),   intent(in)    :: dist

      if (is_setacc_int) then
         call update_particle_acc_int(n_part, cg, cells, dist)
      elseif (is_setacc_cic) then
         call update_particle_acc_cic(n_part, cg, cells)
      elseif (is_setacc_tsc) then
         call update_particle_acc_tsc(cg)
      endif

   end subroutine update_particle_acc

   subroutine update_particle_acc_int(n_part, cg, cells, dist)

      use constants,        only: gpot_n, ndims, xdim, ydim, zdim, zero
      use dataio_pub,       only: warn
      use grid_cont,        only: grid_container
      use named_array_list, only: qna
      use particle_func,    only: df_d_p, d2f_d2_p, d2f_dd_p
      use particle_types, only: particle

      implicit none

      integer,                          intent(in) :: n_part
      type(grid_container), pointer,    intent(in) :: cg
      type(particle), pointer                      :: pset
      integer, dimension(n_part,ndims), intent(in) :: cells
      real, dimension(n_part, ndims),   intent(in) :: dist
      integer                                      :: i
      integer(kind=4)                              :: ig

      ig = qna%ind(gpot_n)
      pset => cg%pset%first
      i=1
      do while (associated(pset))
         if (pset%pdata%outside) then
         ! multipole expansion for particles outside domain
            call warn('[particle_gravity:update_particle_acc_cic] particle is outside the domain - we need multipole expansion!!!')
            pset%pdata%acc = zero
            cycle
         endif

         pset%pdata%acc(xdim) = - (df_d_p([cells(i, :)], cg, ig, xdim) + &
                                   d2f_d2_p([cells(i, :)], cg, ig, xdim)       * dist(i, xdim) + &
                                   d2f_dd_p([cells(i, :)], cg, ig, xdim, ydim) * dist(i, ydim) + &
                                   d2f_dd_p([cells(i, :)], cg, ig, xdim, zdim) * dist(i, zdim))


         pset%pdata%acc(ydim) = -( df_d_p([cells(i, :)], cg, ig, ydim) + &
                                   d2f_d2_p([cells(i, :)], cg, ig, ydim)       * dist(i, ydim) + &
                                   d2f_dd_p([cells(i, :)], cg, ig, xdim, ydim) * dist(i, xdim) + &
                                   d2f_dd_p([cells(i, :)], cg, ig, ydim, zdim) * dist(i, zdim))

         pset%pdata%acc(zdim) = -( df_d_p([cells(i, :)], cg, ig, zdim) + &
                                   d2f_d2_p([cells(i, :)], cg, ig, zdim)       * dist(i, zdim) + &
                                   d2f_dd_p([cells(i, :)], cg, ig, xdim, zdim) * dist(i, xdim) + &
                                   d2f_dd_p([cells(i, :)], cg, ig, ydim, zdim) * dist(i, ydim))
         pset => pset%nxt
         i = i + 1
      enddo

   end subroutine update_particle_acc_int

   subroutine update_particle_acc_cic(n_part, cg, cells)

      use constants,        only: idm, ndims, CENTER, xdim, ydim, zdim, half, zero, gp1b_n, gpot_n, I_ZERO, I_ONE
      use dataio_pub,       only: warn
      use grid_cont,        only: grid_container
      use named_array_list, only: qna
      use particle_types,   only: particle

      implicit none

      integer,                          intent(in)    :: n_part
      type(grid_container), pointer,    intent(inout) :: cg
      type(particle), pointer                         :: pset
      integer, dimension(n_part,ndims), intent(in)    :: cells

      integer                                         :: p
      integer(kind=4)                                 :: c, ig, cdim
      integer(kind=4), dimension(ndims,8), parameter  :: cijk = reshape([[I_ZERO, I_ZERO, I_ZERO], &
           &                                                             [I_ONE,  I_ZERO, I_ZERO], [I_ZERO, I_ONE, I_ZERO], [I_ZERO, I_ZERO, I_ONE], &
           &                                                             [I_ONE,  I_ONE,  I_ZERO], [I_ZERO, I_ONE, I_ONE],  [I_ONE,  I_ZERO, I_ONE], &
           &                                                             [I_ONE,  I_ONE,  I_ONE]], [ndims,8_4])
      integer,         dimension(ndims)               :: cic_cells
      real,            dimension(ndims)               :: dxyz, axyz
      real(kind=8),    dimension(ndims,8)             :: fxyz
      real(kind=8),    dimension(8)                   :: wijk

      if (mask_gpot1b) then
         ig = qna%ind(gp1b_n)
      else
         ig = qna%ind(gpot_n)
      endif

      pset => cg%pset%first
      p=1
      do while (associated(pset))
         associate( part => pset%pdata )

            if (part%outside) then
            ! multipole expansion for particles outside domain
               call warn('[particle_gravity:update_particle_acc_cic] particle is outside the domain - we need multipole expansion!!!')
               part%acc = zero
               cycle
            endif

            if (mask_gpot1b) then
               cg%gp1b = zero
               call gravpot1b(pset, cg, ig)
               cg%gp1b = -cg%gp1b + cg%gpot
            endif

            do cdim = xdim, zdim
               if (part%pos(cdim) < cg%coord(CENTER, cdim)%r(cells(p,cdim))) then
                  cic_cells(cdim) = cells(p, cdim) - 1
               else
                  cic_cells(cdim) = cells(p, cdim)
               endif
               dxyz(cdim) = abs(part%pos(cdim) - cg%coord(CENTER, cdim)%r(cic_cells(cdim)))
            enddo

            wijk(1) = (cg%dx - dxyz(xdim))*(cg%dy - dxyz(ydim))*(cg%dz - dxyz(zdim))  !a(i  ,j  ,k  )
            wijk(2) =          dxyz(xdim) *(cg%dy - dxyz(ydim))*(cg%dz - dxyz(zdim))  !a(i+1,j  ,k  )
            wijk(3) = (cg%dx - dxyz(xdim))*         dxyz(ydim) *(cg%dz - dxyz(zdim))  !a(i  ,j+1,k  )
            wijk(4) = (cg%dx - dxyz(xdim))*(cg%dy - dxyz(ydim))*         dxyz(zdim)   !a(i  ,j  ,k+1)
            wijk(5) =          dxyz(xdim) *         dxyz(ydim) *(cg%dz - dxyz(zdim))  !a(i+1,j+1,k  )
            wijk(6) = (cg%dx - dxyz(xdim))*         dxyz(ydim) *         dxyz(zdim)   !a(i  ,j+1,k+1)
            wijk(7) =          dxyz(xdim) *(cg%dy - dxyz(ydim))*         dxyz(zdim)   !a(i+1,j  ,k+1)
            wijk(8) =          dxyz(xdim) *         dxyz(ydim) *         dxyz(zdim)   !a(i+1,j+1,k+1)

            wijk = wijk/cg%dvol

            do cdim = xdim, zdim
               do c = 1, 8
                  fxyz(cdim,c) = -(cg%q(ig)%point(cic_cells(:)+idm(cdim,:)+cijk(:,c)) - cg%q(ig)%point(cic_cells(:)-idm(cdim,:)+cijk(:,c)))
               enddo
               axyz(cdim) = sum(fxyz(cdim,:)*wijk(:))
            enddo

            part%acc(:) = half * axyz(:) * cg%idl(:)

         end associate
         pset => pset%nxt
         p = p + 1
      enddo

   end subroutine update_particle_acc_cic

   subroutine update_particle_acc_tsc(cg)

      use constants,        only: xdim, ydim, zdim, ndims, LO, HI, IM, I0, IP, CENTER, gp1b_n, gpot_n, idm, half, zero
      use domain,           only: dom
      use func,             only: operator(.equals.)
      use grid_cont,        only: grid_container
      use multipole,        only: moments2pot
      use named_array_list, only: qna
      use particle_types,   only: particle
      use domain,           only: dom

      implicit none


      type(grid_container), pointer, intent(inout) :: cg
      type(particle), pointer                      :: pset, pset2
      integer                                      :: i, j, k
      integer(kind=4)                              :: ig, cdim
      integer, dimension(ndims, IM:IP)             :: ijkp
      integer, dimension(ndims)                    :: cur_ind
      real,    dimension(ndims)                    :: fxyz, axyz
      real,    dimension(5)                        :: tmp
      real                                         :: weight, delta_x, weight_tmp

      if (mask_gpot1b) then
         ig = qna%ind(gp1b_n)
      else
         ig = qna%ind(gpot_n)
      endif

      fxyz = zero ! suppress -Wmaybe-uninitialized
      pset => cg%pset%first
      do while (associated(pset))

         !Delete particles escaping the domain
         if (pset%pdata%energy .equals. 0.0) then
            pset2 => pset%nxt
            call cg%pset%remove(pset)
            pset => pset2
            cycle
         endif

         associate( part => pset%pdata )

            axyz(:) = zero

            if (part%outside) then
               ! multipole expansion for particles outside domain
               tmp = [0, 0, 1, 0, 0]
               do cdim = xdim, zdim
                  fxyz(cdim) =   - ( moments2pot( part%pos(xdim) + cg%dl(xdim) * tmp(cdim+2), &
                                 part%pos(ydim) + cg%dl(ydim) * tmp(cdim+1), part%pos(zdim) + cg%dl(zdim) * tmp(cdim) ) - &
                                 moments2pot( part%pos(xdim) - cg%dl(xdim) * tmp(cdim+2), &
                                 part%pos(ydim) - cg%dl(ydim) * tmp(cdim+1), part%pos(zdim) - cg%dl(zdim) * tmp(cdim) ) )
               enddo
               axyz(:) = fxyz(:)
               part%acc(:) = half * axyz(:) * cg%idl(:)

               pset => pset%nxt
               cycle
            endif

            if (mask_gpot1b) then
               cg%gp1b = zero
               call gravpot1b(pset, cg, ig)
               cg%gp1b = -cg%gp1b + cg%gpot
            endif

            do cdim = xdim, zdim
               if (dom%has_dir(cdim)) then
                  ijkp(cdim, I0) = nint((part%pos(cdim) - cg%fbnd(cdim,LO)-cg%dl(cdim)/2.) * cg%idl(cdim) + int(cg%lhn(cdim, LO)) + dom%nb, kind=4)
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
                        delta_x = ( part%pos(cdim) - cg%coord(CENTER, cdim)%r(cur_ind(cdim)) ) * cg%idl(cdim)

                        if (cur_ind(cdim) /= ijkp(cdim, 0)) then   !!! BEWARE hardcoded magic
                           weight_tmp = 1.125 - 1.5 * abs(delta_x) + half * delta_x**2
                        else
                           weight_tmp = 0.75 - delta_x**2
                        endif

                        weight = weight_tmp * weight
                        fxyz(cdim) = -(cg%q(ig)%point(cur_ind(:) + idm(cdim,:)) - cg%q(ig)%point(cur_ind(:) - idm(cdim,:)))
                     enddo
                     axyz(:) = axyz(:) + fxyz(:)*weight
                  enddo
               enddo
            enddo
            part%acc(:) = half * axyz(:) * cg%idl(:)
         end associate
         pset => pset%nxt
      enddo

   end subroutine update_particle_acc_tsc

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
         acc2(dir) = -( phi_pm_part(pos+real(idm(dir,:))*d, mass) - phi_pm_part(pos-real(idm(dir,:))*d, mass) ) / (two*d)
      enddo

   end subroutine get_acc_model

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
