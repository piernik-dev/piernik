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
   public :: update_particle_gravpot_and_acc, is_setacc_cic, is_setacc_int, phi_pm_part, get_acc_model

   logical :: is_setacc_cic, is_setacc_int

contains

   subroutine update_particle_gravpot_and_acc(pset)

      use constants,      only: ndims, zero
      use cg_leaves,      only: leaves
      use cg_list,        only: cg_list_element
      use dataio_pub,     only: die
      use domain,         only: is_refined, is_multicg
      use grid_cont,      only: grid_container
      use particle_types, only: particle_set
#ifdef VERBOSE
      use dataio_pub,     only: printinfo
#endif /* VERBOSE */

      implicit none

      class(particle_set),   intent(inout) :: pset
      type(grid_container),  pointer       :: cg
      type(cg_list_element), pointer       :: cgl

      integer                              :: n_part
      real,    dimension(:,:), allocatable :: dist
      integer, dimension(:,:), allocatable :: cells

#ifdef VERBOSE
      call printinfo('[particle_gravity:update_particle_gravpot_and_acc] Commencing update of particle gravpot & acceleration')
#endif /* VERBOSE */

      n_part = size(pset%p, dim=1)
      allocate(cells(n_part, ndims), dist(n_part, ndims))

      if (is_refined) call die("[particle_gravity:update_particle_gravpot_and_acc] AMR not implemented for particles yet")
      if (is_multicg) call die("[particle_gravity:update_particle_gravpot_and_acc] multi_cg not implemented for particles yet")
      cgl => leaves%first
      cg  => cgl%cg

      call update_gravpot_from_particles(n_part, pset, cg, zero)

      call save_pot_pset(pset, cg)

      call locate_particles_in_cells(n_part, pset, cg, cells, dist)

      call update_particle_potential_energy(n_part, pset, cg, cells, dist)    !szukanie energii potencjalnej w punktach-polozeniach czastek

      if (is_setacc_int) then
         call update_particle_acc_int(n_part, pset, cg, cells, dist)
      elseif (is_setacc_cic) then
         call update_particle_acc_cic(n_part, pset, cg, cells)
      endif

#ifdef VERBOSE
      call printinfo('[particle_gravity:update_particle_gravpot_and_acc] Finish update of particle gravpot & acceleration')
#endif /* VERBOSE */

   end subroutine update_particle_gravpot_and_acc

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

   subroutine update_gravpot_from_particles(n_part, pset, cg, eps2)

      use constants,      only: CENTER, LO, HI, xdim, ydim, zdim, zero
      use grid_cont,      only: grid_container
      use particle_types, only: particle_set

      implicit none

      integer,                       intent(in)    :: n_part
      class(particle_set),           intent(in)    :: pset  !< particle list
      type(grid_container), pointer, intent(inout) :: cg
      real,                          intent(in)    :: eps2
      integer                                      :: i, j, k, p

      cg%nbgp = zero

      open(unit=999, file='pset.dat', action='write', position='append')
         do p = 1, n_part
            do i = cg%lhn(xdim, LO), cg%lhn(xdim, HI)
               do j = cg%lhn(ydim, LO), cg%lhn(ydim, HI)
                  do k = cg%lhn(zdim, LO), cg%lhn(zdim, HI)
                     cg%nbgp(i,j,k) = cg%nbgp(i,j,k) + phi_pm_part([cg%coord(CENTER,xdim)%r(i) - pset%p(p)%pos(xdim), &
                                                                    cg%coord(CENTER,ydim)%r(j) - pset%p(p)%pos(ydim), &
                                                                    cg%coord(CENTER,zdim)%r(k) - pset%p(p)%pos(zdim)], eps2, pset%p(p)%mass)
                  enddo
               enddo
            enddo
            write(999,*) p, pset%p(p)%pos
         enddo
      close(999)

   end subroutine update_gravpot_from_particles

   subroutine locate_particles_in_cells(n_part, pset, cg, cells, dist)

      use constants,      only: ndims, xdim, CENTER, LO, HI
      use grid_cont,      only: grid_container
      use particle_types, only: particle_set

      implicit none

      integer,                          intent(in)    :: n_part
      class(particle_set),              intent(inout) :: pset  !< particle list
      type(grid_container),             intent(in)    :: cg
      integer, dimension(n_part,ndims), intent(out)   :: cells
      real,    dimension(n_part,ndims), intent(out)   :: dist
      integer                                         :: i, cdim

      do i = 1, n_part
         do cdim = xdim, ndims
            if ((pset%p(i)%pos(cdim) >= cg%ijkse(cdim, LO)) .or. (pset%p(i)%pos(cdim) <= cg%ijkse(cdim, HI))) then
               pset%p(i)%outside = .false.
            else
               pset%p(i)%outside = .true.
            endif

            cells(i, cdim) = int( 0.5 + (pset%p(i)%pos(cdim) - cg%coord(CENTER,cdim)%r(0)) * cg%idl(cdim) )

            dist(i, cdim)  = pset%p(i)%pos(cdim) - ( cg%coord(CENTER, cdim)%r(0) + cells(i,cdim) * cg%dl(cdim) )
         enddo
      enddo

   end subroutine locate_particles_in_cells

   subroutine update_particle_potential_energy(n_part, pset, cg, cells, dist)

      use constants,        only: gpot_n, ndims, half, xdim, ydim, zdim
      use grid_cont,        only: grid_container
      use named_array_list, only: qna
      use particle_func,    only: df_d_o2, d2f_d2_o2, d2f_dd_o2
      use particle_types,   only: particle_set

      implicit none

      integer,                          intent(in)    :: n_part
      class(particle_set),              intent(inout) :: pset  !< particle list
      type(grid_container), pointer,    intent(in)    :: cg
      integer, dimension(n_part,ndims), intent(in)    :: cells
      real,    dimension(n_part,ndims), intent(in)    :: dist
      integer                                         :: p
      integer(kind=4)                                 :: ig
      real, dimension(n_part)                         :: dpot, d2pot

      ig = qna%ind(gpot_n)
      do p = 1, n_part
         dpot(p) = df_d_o2([cells(p, :)], cg, ig, xdim) * dist(p, xdim) + &
                   df_d_o2([cells(p, :)], cg, ig, ydim) * dist(p, ydim) + &
                   df_d_o2([cells(p, :)], cg, ig, zdim) * dist(p, zdim)

         d2pot(p) = d2f_d2_o2([cells(p, :)], cg, ig, xdim) * dist(p, xdim)**2 + &
                    d2f_d2_o2([cells(p, :)], cg, ig, ydim) * dist(p, ydim)**2 + &
                    d2f_d2_o2([cells(p, :)], cg, ig, zdim) * dist(p, zdim)**2 + &
                    2.0*d2f_dd_o2([cells(p, :)], cg, ig, xdim, ydim) * dist(p, xdim)*dist(p, ydim) + &
                    2.0*d2f_dd_o2([cells(p, :)], cg, ig, xdim, zdim) * dist(p, xdim)*dist(p, zdim)
         pset%p(p)%energy = cg%q(ig)%point(cells(p,:)) + dpot(p) + half * d2pot(p)
      enddo

   end subroutine update_particle_potential_energy

   subroutine update_particle_acc_int(n_part, pset, cg, cells, dist)

      use constants,        only: gpot_n, ndims, xdim, ydim, zdim
      use grid_cont,        only: grid_container
      use named_array_list, only: qna
      use particle_func,    only: df_d_p, d2f_d2_p, d2f_dd_p
      use particle_types,   only: particle_set

      implicit none

      integer,                          intent(in)    :: n_part
      class(particle_set),              intent(inout) :: pset
      type(grid_container), pointer,    intent(in)    :: cg
      integer, dimension(n_part,ndims), intent(in)    :: cells
      real, dimension(n_part, ndims),   intent(in)    :: dist
      integer                                         :: i
      integer(kind=4)                                 :: ig

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

   subroutine update_particle_acc_cic(n_part, pset, cg, cells)

      use constants,      only: ndims, CENTER, xdim, ydim, zdim, half, zero
      use grid_cont,      only: grid_container
      use particle_types, only: particle_set

      implicit none

      integer,                          intent(in)    :: n_part
      class(particle_set),              intent(inout) :: pset
      type(grid_container), pointer,    intent(in)    :: cg
      integer, dimension(n_part,ndims), intent(in)    :: cells

      integer                                         :: i, j, k, c, cdim
      integer                                         :: p
      integer(kind=8), dimension(n_part,ndims)        :: cic_cells
      real,            dimension(n_part,ndims)        :: dxyz
      real(kind=8),    dimension(n_part,8)            :: wijk, fx, fy, fz

      do i = 1, n_part
         pset%p(i)%acc = zero
         if ((pset%p(i)%outside) .eqv. .false.) then
            do cdim = xdim, ndims
               if (pset%p(i)%pos(cdim) < cg%coord(CENTER, cdim)%r(cells(i,cdim))) then
                  cic_cells(i, cdim) = cells(i, cdim) - 1
               else
                  cic_cells(i, cdim) = cells(i, cdim)
               endif
               dxyz(i, cdim) = abs(pset%p(i)%pos(cdim) - cg%coord(CENTER, cdim)%r(cic_cells(i,cdim)))

            enddo

            wijk(i, 1) = (cg%dx - dxyz(i, xdim))*(cg%dy - dxyz(i, ydim))*(cg%dz - dxyz(i, zdim)) !a(i  ,j  ,k  )
            wijk(i, 2) = (cg%dx - dxyz(i, xdim))*(cg%dy - dxyz(i, ydim))*         dxyz(i, zdim)  !a(i+1,j  ,k  )
            wijk(i, 3) = (cg%dx - dxyz(i, xdim))*         dxyz(i, ydim) *(cg%dz - dxyz(i, zdim)) !a(i  ,j+1,k  )
            wijk(i, 4) = (cg%dx - dxyz(i, xdim))*         dxyz(i, ydim) *         dxyz(i, zdim)  !a(i  ,j  ,k+1)
            wijk(i, 5) =          dxyz(i, xdim) *(cg%dy - dxyz(i, ydim))*(cg%dz - dxyz(i, zdim)) !a(i+1,j+1,k  )
            wijk(i, 6) =          dxyz(i, xdim) *(cg%dy - dxyz(i, ydim))*         dxyz(i, zdim)  !a(i  ,j+1,k+1)
            wijk(i, 7) =          dxyz(i, xdim) *         dxyz(i, ydim) *(cg%dz - dxyz(i, zdim)) !a(i+1,j  ,k+1)
            wijk(i, 8) =          dxyz(i, xdim) *         dxyz(i, ydim) *         dxyz(i, zdim)  !a(i+1,j+1,k+1)
         !else multipole expansion for particles outside domain
         endif

      enddo

      wijk = wijk/cg%dvol

      do p = 1, n_part
         c = 1
         do i = 0, 1
            do j = 0, 1
               do k = 0, 1
                  fx(p, c) = -(cg%gpot(cic_cells(p, xdim)+1+i, cic_cells(p, ydim)  +j, cic_cells(p, zdim)  +k) - cg%gpot(cic_cells(p, xdim)-1+i, cic_cells(p, ydim)  +j, cic_cells(p, zdim)  +k))
                  fy(p, c) = -(cg%gpot(cic_cells(p, xdim)  +i, cic_cells(p, ydim)+1+j, cic_cells(p, zdim)  +k) - cg%gpot(cic_cells(p, xdim)  +i, cic_cells(p, ydim)-1+j, cic_cells(p, zdim)  +k))
                  fz(p, c) = -(cg%gpot(cic_cells(p, xdim)  +i, cic_cells(p, ydim)  +j, cic_cells(p, zdim)+1+k) - cg%gpot(cic_cells(p, xdim)  +i, cic_cells(p, ydim)  +j, cic_cells(p, zdim)-1+k))
                  c = c + 1
               enddo
            enddo
         enddo
      enddo

      fx = half*fx*cg%idx
      fy = half*fy*cg%idy
      fz = half*fz*cg%idz

      do p = 1, n_part
         do c = 1, 8
            pset%p(p)%acc(xdim) = pset%p(p)%acc(xdim) + wijk(p, c) * fx(p, c)
            pset%p(p)%acc(ydim) = pset%p(p)%acc(ydim) + wijk(p, c) * fy(p, c)
            pset%p(p)%acc(zdim) = pset%p(p)%acc(zdim) + wijk(p, c) * fz(p, c)
         enddo
      enddo

   end subroutine update_particle_acc_cic

   subroutine get_acc_model(pset, p, eps, acc2)

      use constants,      only: ndims, xdim, zdim
      use grid_cont,      only: grid_container
      use particle_types, only: particle_set

      implicit none

      class(particle_set),    intent(in)  :: pset  !< particle list
      integer,                intent(in)  :: p
      real,                   intent(in)  :: eps
      real, dimension(ndims), intent(out) :: acc2
      integer(kind=4)                     :: dir

      do dir = xdim, zdim
         acc2(dir) = -der_xyz(pset%p(p)%pos, 1.0e-8, eps, dir)
      enddo

   end subroutine get_acc_model

   function der_xyz(pos, d, eps, dir)

      use constants, only: idm, ndims

      implicit none

      real(kind=8),dimension(1,ndims), intent(in) :: pos
      real(kind=8),                    intent(in) :: d, eps
      integer(kind=4),                 intent(in) :: dir
      real(kind=8)                                :: der_xyz

      der_xyz = ( phi_pm_part(pos(1,:)+real(idm(dir,:))*d, eps, 10.0) - phi_pm_part(pos(1,:)-real(idm(dir,:))*d, eps, 10.0) ) / (2.0*d)

   end function der_xyz

   subroutine save_pot_pset(pset, cg)

      use constants,      only: CENTER, LO, HI, xdim, ydim, zdim
      use dataio_pub,     only: printinfo, printio
      use grid_cont,      only: grid_container
      use particle_types, only: particle_set

      implicit none

      class(particle_set),           intent(in) :: pset  !< particle list
      type(grid_container), pointer, intent(in) :: cg

      integer                                   :: numer1, numer2, i, j, p
      logical                                   :: finish, save_potential

      save_potential = .true.
      finish         = .false.
      numer1 = 64
      numer2 = 20

      open(unit=90, file='pset.dat', action='write', position='append')
         do p=1, ubound(pset%p, dim=1)
            write(90,*) p, pset%p(p)%pos
         enddo
      close(90)
      if (save_potential) then
         call printio('[particle_integrators:save_pot_pset] Writing potential to files: potencjal1.dat, potencjal2.dat')
         open(unit=88, file='potencjal1.dat')
         open(unit=89, file='potencjal2.dat')
            do i = cg%lhn(xdim, LO), cg%lhn(xdim, HI)
               do j = cg%lhn(ydim, LO), cg%lhn(ydim, HI)
                  !do k = cg%lhn(zdim, LO), cg%lhn(zdim, HI)
                  write(88,*) i, j, numer1, cg%coord(CENTER, xdim)%r(i), cg%coord(CENTER, ydim)%r(j), cg%coord(CENTER, zdim)%r(numer1), cg%nbgp(i,j,numer1)
                  write(89,*) i, j, numer2, cg%coord(CENTER, xdim)%r(i), cg%coord(CENTER, ydim)%r(j), cg%coord(CENTER, zdim)%r(numer2), cg%nbgp(i,j,numer2)
                  !enddo
               enddo
               write(88,*)
               write(89,*)
            enddo
         close(88)
         close(89)

         if (finish) then
            call printinfo('[particle_integrators:save_pot_pset] Condition to end - finishing.')
            stop
         endif
      endif

   end subroutine save_pot_pset

end module particle_gravity
