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
!! \brief Star formation feedback + stellar particle formation
!!
!<
module star_formation
! pulled by NBODY

   use constants, only: dsetnamelen

   implicit none

   private
   public :: SF, initialize_id, attribute_id, pid_gen, dmass_stars

   integer(kind=4), parameter            :: giga = 1000000000
   integer(kind=4)                       :: pid_gen, maxpid, dpid
   real                                  :: dmass_stars
   character(len=dsetnamelen), parameter :: sfr_n   = "SFR_n"
contains

   subroutine SF(forward)

      use cg_leaves,        only: leaves
      use cg_list,          only: cg_list_element
      use constants,        only: ndims, xdim, ydim, zdim, LO, HI, CENTER, pi, nbdn_n
#ifdef COSM_RAYS
      use cr_data,          only: icr_H1, cr_table
      use initcosmicrays,   only: iarr_crn, cr_active
#endif /* COSM_RAYS */
      use fluidindex,       only: flind
      use fluidtypes,       only: component_fluid
      use global,           only: t, dt
      use grid_cont,        only: grid_container
      use named_array_list, only: qna
      use particle_types,   only: particle
      use particle_utils,   only: is_part_in_cg
      use units,            only: newtong, cm, sek, gram, erg
#ifdef THERM
      use thermal,          only: itemp
#endif /* THERM */

      implicit none

      logical, intent(in)             :: forward
      type(cg_list_element), pointer  :: cgl
      type(grid_container),  pointer  :: cg
      type(particle), pointer         :: pset
      class(component_fluid), pointer :: pfl
      integer                         :: ifl, i, j, k, i1, j1, k1, aijk1
      integer(kind=4)                 :: pid, ig, ir, n_SN
      real, dimension(ndims)          :: pos, vel, acc
      real                            :: dens_thr, sf_dens2dt, c_tau_ff, sfdf, eps_sf, frac, mass_SN, mass, ener, tdyn, tbirth, padd, t1, fact, stage, en_SN, en_SN01, en_SN09, mfdv, tinj, fpadd
      logical                         :: in, phy, out, fed, kick

      if (.not. forward) return

      tinj     = 6.5
      fpadd    = 1.8e40 * gram * cm /sek * 2.**0.38 * 2 * dt / tinj / 26  ! see Agertz+2013
      dens_thr = 0.035
      eps_sf   = 0.1
      n_SN     = 1000
      kick     = .false.

      mass_SN  = 100.0 * n_SN
      en_SN    = n_SN * 10.0**51 * erg
      en_SN01  = 0.1 * en_SN
      en_SN09  = (1 - 0.1 * cr_active) * en_SN
      c_tau_ff = sqrt(3.*pi/(32.*newtong))
      sfdf     = eps_sf / c_tau_ff * 2 * dt

      dmass_stars = 0.0
      ig = qna%ind(nbdn_n)
      ir = qna%ind(sfr_n)
      cgl => leaves%first
      do while (associated(cgl))
         cg => cgl%cg
         do ifl = 1, flind%fluids
            pfl => flind%all_fluids(ifl)%fl
            do i = cg%ijkse(xdim,LO), cg%ijkse(xdim,HI)
               do j = cg%ijkse(ydim,LO), cg%ijkse(ydim,HI)
                  do k = cg%ijkse(zdim,LO), cg%ijkse(zdim,HI)
                     if (cg%u(pfl%idn,i,j,k) > dens_thr) then
#ifdef THERM
                        if (cg%q(itemp)%arr(i,j,k) < 10**4) then
#endif /* THERM */
                           fed = .false.
                           sf_dens2dt = sfdf * cg%u(pfl%idn,i,j,k)**(3./2.)
                           mass    = sf_dens2dt * cg%dvol
                           pset => cg%pset%first
                           do while (associated(pset))
                              if (pos_in_1dim(pset%pdata%pos(xdim), cg%coord(LO,xdim)%r(i-1), cg%coord(HI,xdim)%r(i+1)) .and. &
                               &  pos_in_1dim(pset%pdata%pos(ydim), cg%coord(LO,ydim)%r(j-1), cg%coord(HI,ydim)%r(j+1)) .and. &
                               &  pos_in_1dim(pset%pdata%pos(zdim), cg%coord(LO,zdim)%r(k-1), cg%coord(HI,zdim)%r(k+1)) ) then

                                 if ((pset%pdata%tform >= -10.0) .and. (pset%pdata%mass < mass_SN)) then
                                    stage = aint(pset%pdata%mass/mass_SN)
                                    frac = sf_dens2dt / cg%u(pfl%idn,i,j,k)
                                    pset%pdata%vel      = (pset%pdata%mass * pset%pdata%vel + frac * cg%u(pfl%imx:pfl%imz,i,j,k) * cg%dvol) / (pset%pdata%mass + mass)
                                    pset%pdata%mass     =  pset%pdata%mass + mass
                                    dmass_stars         =  dmass_stars + mass
                                    cg%q(ir)%arr(i,j,k) = cg%q(ir)%arr(i,j,k) + mass
                                    cg%u(pfl%ien,i,j,k)         = (1 - frac) * cg%u(pfl%ien,i,j,k) !- frac * ekin(cg%u(pfl%imx,i,j,k), cg%u(pfl%imy,i,j,k), cg%u(pfl%imz,i,j,k), cg%u(pfl%idn,i,j,k))
                                    cg%u(pfl%idn,i,j,k)         = (1 - frac) * cg%u(pfl%idn,i,j,k)
                                    cg%u(pfl%imx:pfl%imz,i,j,k) = (1 - frac) * cg%u(pfl%imx:pfl%imz,i,j,k)
                                    if (aint(pset%pdata%mass/mass_SN) > stage) then
                                       if (.not. kick) then
                                          mfdv = (aint(pset%pdata%mass/mass_SN) - stage) / cg%dvol
#ifdef THERM
                                          cg%u(pfl%ien,i,j,k) = cg%u(pfl%ien,i,j,k)  + mfdv * en_SN09 ! adding SN energy
#endif /* THERM */
#ifdef COSM_RAYS
                                          if (cr_active > 0.0) cg%u(iarr_crn(cr_table(icr_H1)),i,j,k) = cg%u(iarr_crn(cr_table(icr_H1)),i,j,k) + mfdv * en_SN01  ! adding CR
#endif /* COSM_RAYS */
                                       endif
                                       pset%pdata%tform = t
                                    endif
                                    fed = .true.
                                    exit
                                 endif
                              endif
                              pset => pset%nxt
                           enddo
                           if (.not. fed) then
                              call attribute_id(pid)
                              pos = [cg%coord(CENTER, xdim)%r(i), cg%coord(CENTER, ydim)%r(j), cg%coord(CENTER, zdim)%r(k)]
                              vel = cg%u(pfl%imx:pfl%imz,i,j,k) / cg%u(pfl%idn,i,j,k)
                              frac = sf_dens2dt / cg%u(pfl%idn,i,j,k)
                              acc  = 0.0
                              ener = 0.0
                              tdyn = sqrt(3 * pi / (32 * newtong * cg%u(pfl%idn,i,j,k) + cg%q(ig)%arr(i,j,k)))
                              call is_part_in_cg(cg, pos, .true., in, phy, out)
                              dmass_stars = dmass_stars + mass
                              cg%q(ir)%arr(i,j,k) = cg%q(ir)%arr(i,j,k) + mass
                              cg%u(pfl%ien,i,j,k)          = (1 - frac) * cg%u(pfl%ien,i,j,k) !- frac * ekin(cg%u(pfl%imx,i,j,k), cg%u(pfl%imy,i,j,k), cg%u(pfl%imz,i,j,k), cg%u(pfl%idn,i,j,k))
                              cg%u(pfl%idn,i,j,k)          = (1 - frac) * cg%u(pfl%idn,i,j,k)
                              cg%u(pfl%imx:pfl%imz,i,j,k)  = (1 - frac) * cg%u(pfl%imx:pfl%imz,i,j,k)
                              tbirth = -10
                              if (mass > mass_SN) then
                                 if (.not. kick) then
                                    mfdv = aint(mass/mass_SN) / cg%dvol
#ifdef THERM
                                    cg%u(pfl%ien,i,j,k) = cg%u(pfl%ien,i,j,k)  + mfdv * en_SN09 ! adding SN energy
#endif /* THERM */
#ifdef COSM_RAYS
                                    if (cr_active > 0.0) cg%u(iarr_crn(cr_table(icr_H1)),i,j,k) = cg%u(iarr_crn(cr_table(icr_H1)),i,j,k) + mfdv * en_SN01  ! adding CR
#endif /* COSM_RAYS */
                                 endif
                                 tbirth = t
                              endif
                              call cg%pset%add(pid, mass, pos, vel, acc, ener, in, phy, out, tbirth, tdyn)
                           endif
#ifdef THERM
                        endif
#endif /* THERM */
                     endif
                  enddo
               enddo
            enddo
         enddo
! KICK
         if (kick) then
            pset => cg%pset%first
            do while (associated(pset))
               if (t < pset%pdata%tform + 10) then
                  t1 = t - pset%pdata%tform
                  do ifl = 1, flind%fluids
                     pfl => flind%all_fluids(ifl)%fl
                     do i = cg%ijkse(xdim,LO), cg%ijkse(xdim,HI)
                        if ( pos_in_1dim(pset%pdata%pos(xdim), cg%coord(LO,xdim)%r(i-1), cg%coord(HI,xdim)%r(i+1)) ) then
                           do j = cg%ijkse(ydim,LO), cg%ijkse(ydim,HI)
                              if ( pos_in_1dim(pset%pdata%pos(ydim), cg%coord(LO,ydim)%r(j-1), cg%coord(HI,ydim)%r(j+1)) ) then
                                 do k = cg%ijkse(zdim,LO), cg%ijkse(zdim,HI)
                                    if ( pos_in_1dim(pset%pdata%pos(zdim), cg%coord(LO,zdim)%r(k-1), cg%coord(HI,zdim)%r(k+1)) ) then
                                       i1 = nint((pset%pdata%pos(xdim) - cg%coord(CENTER,xdim)%r(i)) / cg%dx)
                                       j1 = nint((pset%pdata%pos(ydim) - cg%coord(CENTER,ydim)%r(j)) / cg%dy)
                                       k1 = nint((pset%pdata%pos(zdim) - cg%coord(CENTER,zdim)%r(k)) / cg%dz)
                                       aijk1 = abs(i1) + abs(j1) + abs(k1)
                                       if (t1 < tinj) then
                                          fact = 0.0
                                          if (aijk1 > 0.0) fact = 1.0 / sqrt(real(aijk1))
                                          padd = pset%pdata%mass * fpadd / cg%dvol
                                          !else
                                          !   padd = 3.6 * 10**4 * pset%pdata%mass/200 * 2*dt/40.0 / cg%dvol / 26    ! should use initial mass, not current mass
                                          !endif

                                          ! Momentum kick
                                          cg%u(pfl%ien,i,j,k) = cg%u(pfl%ien,i,j,k) - 0.5*(cg%u(pfl%imx,i,j,k)**2 +cg%u(pfl%imy,i,j,k)**2 + cg%u(pfl%imz,i,j,k)**2)/cg%u(pfl%idn,i,j,k)  ! remove ekin
                                          cg%u(pfl%imx,i,j,k) = cg%u(pfl%imx,i,j,k) + fact * i1 * padd
                                          cg%u(pfl%imy,i,j,k) = cg%u(pfl%imy,i,j,k) + fact * j1 * padd
                                          cg%u(pfl%imz,i,j,k) = cg%u(pfl%imz,i,j,k) + fact * k1 * padd
                                          cg%u(pfl%ien,i,j,k) = cg%u(pfl%ien,i,j,k) + 0.5*(cg%u(pfl%imx,i,j,k)**2 +cg%u(pfl%imy,i,j,k)**2 + cg%u(pfl%imz,i,j,k)**2)/cg%u(pfl%idn,i,j,k)  ! add new ekin
                                       endif
                                       if (aijk1 == 0) then
                                          if ((t1 - dt < tinj) .and. (t1 + dt > tinj)) then    ! Instantaneous injection Agertz
                                             mfdv = aint(pset%pdata%mass/mass_SN) / cg%dvol
#ifdef THERM
                                             cg%u(pfl%ien,i,j,k) = cg%u(pfl%ien,i,j,k)  + mfdv * en_SN09  ! adding SN energy
#endif /* THERM */
#ifdef COSM_RAYS
                                             if (cr_active > 0.0) cg%u(iarr_crn(cr_table(icr_H1)),i,j,k) = cg%u(iarr_crn(cr_table(icr_H1)),i,j,k) + mfdv * en_SN01  ! adding CR
#endif /* COSM_RAYS */
                                          endif
                                       endif
                                    endif
                                 enddo
                              endif
                           enddo
                        endif
                     enddo
                  enddo
               endif
               pset => pset%nxt
            enddo
         endif

         cgl => cgl%nxt
      enddo

   end subroutine SF

   logical function pos_in_1dim(pos, pl, pr) result(isin)

      implicit none

      real, intent(in) :: pos, pl, pr

      isin = (pl < pos .and. pr > pos)

   end function pos_in_1dim

   subroutine initialize_id()

      use mpisetup, only: proc, nproc

      implicit none

      dpid = int(giga/nproc, kind=4)
      pid_gen = proc * dpid
      !maxpid = (proc+1) * dpid
      !print *, proc, pid_gen, maxpid, dpid

   end subroutine initialize_id

   subroutine attribute_id(pid)

      use constants, only: I_ONE
      use mpisetup,  only: proc, nproc

      implicit none

      integer(kind=4), intent(out)    :: pid

      dpid = int(giga/nproc, kind=4)
      maxpid = (proc + I_ONE) * dpid
      if (pid_gen >= maxpid) maxpid = (proc+I_ONE) * dpid + giga

      pid_gen = pid_gen + I_ONE
      if (pid_gen >= maxpid) then
         print *, 'pool of pid full for proc ', proc, pid_gen, maxpid
         pid_gen = proc * dpid + giga
         maxpid = (proc + I_ONE) * dpid + giga
      endif
      pid = pid_gen

   end subroutine attribute_id

end module star_formation
