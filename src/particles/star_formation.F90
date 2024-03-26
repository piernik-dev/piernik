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
   real                                  :: dmass_stars, dens_thr
   character(len=dsetnamelen), parameter :: sfr_n   = "SFR_n"

contains

   subroutine SF(forward)

      use cg_leaves,        only: leaves
      use cg_list,          only: cg_list_element
      use constants,        only: ndims, xdim, ydim, zdim, LO, HI, CENTER, pi, nbdn_n
      use domain,           only: dom
      use fluidindex,       only: flind
      use fluidtypes,       only: component_fluid
      use func,             only: ekin
      use global,           only: t, dt
      use grid_cont,        only: grid_container
      use named_array_list, only: qna
      use particle_func,    only: particle_in_area, ijk_of_particle, l_neighb_part, r_neighb_part
      use particle_types,   only: particle
      use particle_utils,   only: is_part_in_cg
      use units,            only: newtong, cm, sek, gram, erg
#ifdef COSM_RAYS
      use initcosmicrays,   only: cr_active
#endif /* COSM_RAYS */

      implicit none

      logical, intent(in)               :: forward
      type(cg_list_element), pointer    :: cgl
      type(grid_container),  pointer    :: cg
      type(particle), pointer           :: pset
      class(component_fluid), pointer   :: pfl
      integer(kind=4)                   :: pid, ig, ir, n_SN, ifl, i, j, k, aijk1
      integer(kind=4), dimension(ndims) :: ijk1, ijkp, ijkl, ijkr
      real, dimension(ndims)            :: pos, vel, acc
      real, dimension(ndims,LO:HI)      :: sector
      real                              :: sf_dens2dt, c_tau_ff, sfdf, eps_sf, frac, mass_SN, mass, ener, tdyn, tbirth, padd, t1, tj, stage, en_SN, en_SN01, en_SN09, mfdv, tini, tinj, fpadd
      logical                           :: in, phy, out, fin, fed, kick, tcond1, tcond2

      if (.not. forward) return

      tini     = 10.0
      tinj     = 6.5
      fpadd    = 1.8e40 * gram * cm /sek * 2.**0.38 * 2 * dt / tinj / 26  ! see Agertz+2013
      dens_thr = 0.035
      eps_sf   = 0.1
      n_SN     = 1000
      kick     = .false.

      mass_SN  = 100.0 * n_SN
      en_SN    = n_SN * 10.0**51 * erg
      en_SN01  = 0.1 * en_SN
#ifdef COSM_RAYS
      en_SN09  = (1 - 0.1 * cr_active) * en_SN
#else /* !COSM_RAYS */
      en_SN09  = 0.0
#endif /* !COSM_RAYS */
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
               sector(xdim,:) = [cg%coord(LO,xdim)%r(i-1), cg%coord(HI,xdim)%r(i+1)]
               do j = cg%ijkse(ydim,LO), cg%ijkse(ydim,HI)
                  sector(ydim,:) = [cg%coord(LO,ydim)%r(j-1), cg%coord(HI,ydim)%r(j+1)]
                  do k = cg%ijkse(zdim,LO), cg%ijkse(zdim,HI)
                     sector(zdim,:) = [cg%coord(LO,zdim)%r(k-1), cg%coord(HI,zdim)%r(k+1)]
                     if (.not.check_threshold(cg, pfl%idn, i, j, k)) cycle
                     fed = .false.
                     sf_dens2dt = sfdf * cg%u(pfl%idn,i,j,k)**(3./2.)
                     mass       = sf_dens2dt * cg%dvol
                     pset => cg%pset%first
                     do while (associated(pset))
                        if ((pset%pdata%tform + tini >= 0.0) .and. (pset%pdata%mass < mass_SN)) then
                           if (particle_in_area(pset%pdata%pos, sector)) then
                              stage = aint(pset%pdata%mass / mass_SN)
                              frac = sf_dens2dt / cg%u(pfl%idn,i,j,k)
                              pset%pdata%vel      = (pset%pdata%mass * pset%pdata%vel + frac * cg%u(pfl%imx:pfl%imz,i,j,k) * cg%dvol) / (pset%pdata%mass + mass)
                              pset%pdata%mass     =  pset%pdata%mass + mass
                              call sf_fed(cg, pfl, i, j, k, ir, mass, 1 - frac)
                              if (aint(pset%pdata%mass / mass_SN) > stage) then
                                 if (.not. kick) then
                                    mfdv = (aint(pset%pdata%mass / mass_SN) - stage) / cg%dvol
                                    call sf_inject(cg, pfl%ien, i, j, k, mfdv * en_SN09, mfdv * en_SN01)
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
                        call is_part_in_cg(cg, pos, .true., in, phy, out, fin)
                        call sf_fed(cg, pfl, i, j, k, ir, mass, 1 - frac)
                        tbirth = -tini
                        if (mass > mass_SN) then
                           if (.not. kick) then
                              mfdv = aint(mass/mass_SN) / cg%dvol
                              call sf_inject(cg, pfl%ien, i, j, k, mfdv * en_SN09, mfdv * en_SN01)
                           endif
                           tbirth = t
                        endif
                        call cg%pset%add(pid, mass, pos, vel, acc, ener, in, phy, out, fin, tbirth, tdyn)
                     endif
                  enddo
               enddo
            enddo
         enddo
! KICK
         if (kick) then
            pset => cg%pset%first
            do while (associated(pset))
               t1 = t - pset%pdata%tform
               tj = t1 - tinj
               tcond1 = (tj < 0.0)
               tcond2 = (abs(tj) < dt)
               if (t1 < tini .and. (tcond1 .or. tcond2)) then
                  ijkp = ijk_of_particle(pset%pdata%pos, dom%edge(:,LO), cg%idl)
                  ijkl = l_neighb_part(ijkp, cg%ijkse(:,LO))
                  ijkr = r_neighb_part(ijkp, cg%ijkse(:,HI))
                  do ifl = 1, flind%fluids
                     pfl => flind%all_fluids(ifl)%fl
                     do i = ijkl(xdim), ijkr(xdim)
                        do j = ijkl(ydim), ijkr(ydim)
                           do k = ijkl(zdim), ijkr(zdim)
                              ijk1 = nint((pset%pdata%pos - [cg%coord(CENTER,xdim)%r(i), cg%coord(CENTER,ydim)%r(j), cg%coord(CENTER,zdim)%r(k)]) * cg%idl, kind=4)
                              aijk1 = sum(abs(ijk1))
                              if (aijk1 > 0.0 .and. tcond1) then
                                 padd = pset%pdata%mass * fpadd / cg%dvol / sqrt(real(aijk1))
                              !else
                              !   padd = 3.6 * 10**4 * pset%pdata%mass/200 * 2*dt/40.0 / cg%dvol / 26    ! should use initial mass, not current mass
                              !endif

                                 ! Momentum kick
                                 cg%u(pfl%ien,i,j,k) = cg%u(pfl%ien,i,j,k) - ekin(cg%u(pfl%imx,i,j,k), cg%u(pfl%imy,i,j,k), cg%u(pfl%imz,i,j,k), cg%u(pfl%idn,i,j,k))
                                 cg%u(pfl%imx:pfl%imz,i,j,k) = cg%u(pfl%imx:pfl%imz,i,j,k) + ijk1 * padd
                                 cg%u(pfl%ien,i,j,k) = cg%u(pfl%ien,i,j,k) + ekin(cg%u(pfl%imx,i,j,k), cg%u(pfl%imy,i,j,k), cg%u(pfl%imz,i,j,k), cg%u(pfl%idn,i,j,k))
                              else if (aijk1 == 0 .and. tcond2) then    ! Instantaneous injection Agertz
                                 mfdv = aint(pset%pdata%mass / mass_SN) / cg%dvol
                                 call sf_inject(cg, pfl%ien, i, j, k, mfdv * en_SN09, mfdv * en_SN01)
                              endif
                           enddo
                        enddo
                     enddo
                  enddo
               endif
               pset => pset%nxt
            enddo
         endif

         cgl => cgl%nxt
      enddo

   end subroutine SF

   subroutine sf_fed(cg, pfl, i, j, k, ir, mass, frac1)

      use fluidtypes, only: component_fluid
      use grid_cont,  only: grid_container

      implicit none

      type(grid_container),   pointer :: cg
      class(component_fluid), pointer :: pfl
      integer(kind=4),     intent(in) :: i, j, k, ir
      real,                intent(in) :: mass, frac1

      dmass_stars                 = dmass_stars         + mass
      cg%q(ir)%arr(i,j,k)         = cg%q(ir)%arr(i,j,k) + mass
      cg%u(pfl%ien,i,j,k)         = frac1 * cg%u(pfl%ien,i,j,k) !- frac * ekin(cg%u(pfl%imx,i,j,k), cg%u(pfl%imy,i,j,k), cg%u(pfl%imz,i,j,k), cg%u(pfl%idn,i,j,k))
      cg%u(pfl%idn,i,j,k)         = frac1 * cg%u(pfl%idn,i,j,k)
      cg%u(pfl%imx:pfl%imz,i,j,k) = frac1 * cg%u(pfl%imx:pfl%imz,i,j,k)

   end subroutine sf_fed

   subroutine sf_inject(cg, ien, i, j, k, mft, mfcr)

      use grid_cont,      only: grid_container
#ifdef COSM_RAYS
      use cr_data,        only: icr_H1, cr_table
      use initcosmicrays, only: iarr_crn, cr_active
#endif /* COSM_RAYS */

      implicit none

      type(grid_container), pointer :: cg
      integer(kind=4),   intent(in) :: ien, i, j, k
      real,              intent(in) :: mft, mfcr

#ifdef THERM
      cg%u(ien,i,j,k) = cg%u(ien,i,j,k)  + mft  ! adding SN energy
#endif /* THERM */
#ifdef COSM_RAYS
      if (cr_active > 0.0) cg%u(iarr_crn(cr_table(icr_H1)),i,j,k) = cg%u(iarr_crn(cr_table(icr_H1)),i,j,k) + mfcr  ! adding CR
#endif /* COSM_RAYS */

      return
      if (cg%u(ien,i,j,k) > mft * mfcr) return ! suppress compiler warnings on unused arguments

   end subroutine sf_inject

   logical function check_threshold(cg, idn, i, j, k) result(thres)

      use grid_cont, only: grid_container
#ifdef THERM
      use thermal,   only: itemp
#endif /* THERM */

      implicit none

      type(grid_container), pointer :: cg
      integer(kind=4),   intent(in) :: idn, i, j, k

      thres = (cg%u(idn,i,j,k) > dens_thr)
#ifdef THERM
      thres = thres .and. (cg%q(itemp)%arr(i,j,k) < 10**4)
#endif /* THERM */

end function check_threshold

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
