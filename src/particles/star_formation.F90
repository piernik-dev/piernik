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
  implicit none

  private
  public :: SF, initialize_id, attribute_id, pid_gen, dmass_stars

  integer(kind=4)       :: pid_gen, maxpid, dpid

contains

  subroutine SF(forward)

    use cg_leaves,             only: leaves
    use cg_list,               only: cg_list_element
    use constants,             only: ndims, xdim, ydim, zdim, LO, HI, CENTER, pi, nbdn_n
#ifdef COSM_RAYS
    use cr_data,               only: icr_H1, cr_table
    use initcosmicrays,        only: iarr_crn, cr_active
#endif /* COSM_RAYS */
    use fluidindex,            only: flind
    use fluidtypes,            only: component_fluid
    use func,                  only: operator(.equals.), ekin
    use global,                only: t, dt
    use grid_cont,             only: grid_container
    use named_array_list,      only: wna, qna
    use particle_types,        only: particle
    use particle_utils,        only: is_part_in_cg
#ifdef THERM
    use thermal,               only: itemp
#endif /* THERM */
    use units,                 only: newtong, erg, cm, sek, gram

    logical, intent(in)                                :: forward
    type(cg_list_element), pointer                     :: cgl
    type(grid_container),  pointer                     :: cg
    type(particle), pointer                            :: pset
    class(component_fluid), pointer                    :: pfl
    integer                                            :: ifl, i, j, k, i1, j1, k1
    real                                               :: dens_thr, sf_dens, c_tau_ff, eps_sf, frac, mass_SN
    logical                                            :: fed, kick
    integer(kind=4)                                    :: pid, ig, stage, n_SN
    real, dimension(ndims)                             :: pos, vel, acc
    real                                               :: dmass_stars, mass, ener, tdyn, tbirth, padd, t1, fact
    logical                                            :: in, phy, out

    if (.not. forward) return
    dens_thr = 0.035
    kick = .false.
    eps_sf = 0.1
    n_SN = 1000
    mass_SN = 100.0 * n_SN
    dmass_stars = 0.0
    ig = qna%ind(nbdn_n)
    cgl => leaves%first
    do while (associated(cgl))
       cg => cgl%cg
       do ifl = 1, flind%fluids
          pfl => flind%all_fluids(ifl)%fl
          do i = cg%ijkse(xdim,LO), cg%ijkse(xdim,HI)
             do j = cg%ijkse(ydim,LO), cg%ijkse(ydim,HI)
                do k = cg%ijkse(zdim,LO), cg%ijkse(zdim,HI)
                   if (cg%u(pfl%idn,i,j,k) .gt. dens_thr) then
#ifdef THERM
                      if (cg%q(itemp)%arr(i,j,k) .lt. 10**4)) then
#endif /* THERM */
                      fed = .false.
                      c_tau_ff = sqrt(3.*pi/(32.*newtong))
                      sf_dens = eps_sf / c_tau_ff * cg%u(pfl%idn,i,j,k)**(3./2.)
                      pset => cg%pset%first
                      do while (associated(pset))
                         if (cg%coord(LO,xdim)%r(i) .lt. pset%pdata%pos(1)+cg%dx .and. cg%coord(HI,xdim)%r(i) .gt. pset%pdata%pos(1)-cg%dx) then
                            if (cg%coord(LO,ydim)%r(j) .lt. pset%pdata%pos(2)+cg%dy .and. cg%coord(HI,ydim)%r(j) .gt. pset%pdata%pos(2)-cg%dy) then
                               if (cg%coord(LO,zdim)%r(k) .lt. pset%pdata%pos(3)+cg%dz .and. cg%coord(HI,zdim)%r(k) .gt. pset%pdata%pos(3)-cg%dz) then

                                  if ((pset%pdata%tform .ge. -10.0) .and. (pset%pdata%mass .lt. mass_SN)) then
                                     stage = aint(pset%pdata%mass/mass_SN)
                                     frac = sf_dens * 2*dt / cg%u(pfl%idn, i, j, k)
                                     pset%pdata%vel(1:3) = (pset%pdata%mass *pset%pdata%vel(1:3) + frac * cg%u(pfl%imx:pfl%imz,i,j,k) * cg%dvol) / (pset%pdata%mass + sf_dens * cg%dvol * 2*dt)
                                     pset%pdata%mass     =  pset%pdata%mass + sf_dens * cg%dvol * 2*dt
                                     dmass_stars = dmass_stars + sf_dens * cg%dvol * 2*dt
                                     cg%q(qna%ind("SFR_n"))%arr(i,j,k)  = cg%q(qna%ind("SFR_n"))%arr(i,j,k) + sf_dens * cg%dvol * 2*dt
                                     cg%u(pfl%ien, i, j, k)          = (1-frac) * cg%u(pfl%ien, i, j, k) !- frac * ekin(cg%u(pfl%imx,i,j,k), cg%u(pfl%imy,i,j,k), cg%u(pfl%imz,i,j,k), cg%u(pfl%idn, i, j, k))
                                     cg%w(wna%fi)%arr(pfl%idn,i,j,k) = (1 - frac) * cg%w(wna%fi)%arr(pfl%idn,i,j,k)
                                     cg%u(pfl%imx:pfl%imz, i, j, k)  = (1 - frac) * cg%u(pfl%imx:pfl%imz, i, j, k)
                                     if (aint(pset%pdata%mass/mass_SN) .gt. stage) then
                                        if (.not. kick) then
#ifdef THERM
                                           cg%u(pfl%ien,i,j,k) = cg%u(pfl%ien,i,j,k)  + (aint(pset%pdata%mass/mass_SN) - stage) * (1-0.1*cr_active) * n_SN * 10.0**51 * erg / cg%dvol ! adding SN energy
#endif /* THERM */
#ifdef COSM_RAYS
                                           if (cr_active > 0.0) cg%w(wna%fi)%arr(iarr_crn(cr_table(icr_H1)),i,j,k) = cg%w(wna%fi)%arr(iarr_crn(cr_table(icr_H1)),i,j,k) + (aint(pset%pdata%mass/mass_SN) - stage) * 0.1 * n_SN * 10.0**51 * erg / cg%dvol  ! adding CR
#endif /* COSM_RAYS */
                                        endif
                                        pset%pdata%tform = t
                                        !print *, 'particle filled', pset%pdata%pid, aint(pset%pdata%mass/mass_SN) - stage
                                     endif
                                     fed = .true.
                                     exit
                                  endif
                               endif
                            endif
                         endif
                         pset => pset%nxt
                      enddo
                      if (.not. fed) then
                         call attribute_id(pid)
                         pos(1) = cg%coord(CENTER, xdim)%r(i)
                         pos(2) = cg%coord(CENTER, ydim)%r(j)
                         pos(3) = cg%coord(CENTER, zdim)%r(k)
                         mass   = sf_dens * cg%dvol * 2*dt
                         frac   = mass / cg%u(pfl%idn, i, j, k) / cg%dvol
                         vel(1) = cg%u(pfl%imx, i, j, k) / cg%w(wna%fi)%arr(pfl%idn,i,j,k)
                         vel(2) = cg%u(pfl%imy, i, j, k) / cg%w(wna%fi)%arr(pfl%idn,i,j,k)
                         vel(3) = cg%u(pfl%imz, i, j, k) / cg%w(wna%fi)%arr(pfl%idn,i,j,k)
                         acc(:)=0.0
                         ener = 0.0
                         tdyn = sqrt(3*pi/(32*newtong*(cg%w(wna%fi)%arr(pfl%idn,i,j,k))+cgl%cg%q(ig)%arr(i,j,k)))
                         call is_part_in_cg(cg, pos, in, phy, out)
                         dmass_stars = dmass_stars + mass
                         cg%q(qna%ind("SFR_n"))%arr(i,j,k)  = cg%q(qna%ind("SFR_n"))%arr(i,j,k) + sf_dens* cg%dvol * 2*dt
                         cg%u(pfl%ien, i, j, k)          = (1-frac) * cg%u(pfl%ien, i, j, k) !- frac * ekin(cg%u(pfl%imx,i,j,k), cg%u(pfl%imy,i,j,k), cg%u(pfl%imz,i,j,k), cg%u(pfl%idn, i, j, k))
                         cg%w(wna%fi)%arr(pfl%idn,i,j,k) = (1 - frac) * cg%w(wna%fi)%arr(pfl%idn,i,j,k)
                         cg%u(pfl%imx:pfl%imz, i, j, k)  = (1 - frac) * cg%u(pfl%imx:pfl%imz, i, j, k)
                         tbirth = -10
                         if (mass .gt. mass_SN) then
                            if (.not. kick) then
#ifdef THERM
                               cg%u(pfl%ien,i,j,k) = cg%u(pfl%ien,i,j,k)  + aint(mass/mass_SN) * (1-0.1*cr_active) * n_SN * 10.0**51 * erg / cg%dvol ! adding SN energy
#endif /* THERM */
#ifdef COSM_RAYS
                               if (cr_active > 0.0) cg%w(wna%fi)%arr(iarr_crn(cr_table(icr_H1)),i,j,k) = cg%w(wna%fi)%arr(iarr_crn(cr_table(icr_H1)),i,j,k) + aint(mass/mass_SN) * 0.1 * n_SN * 10.0**51 * erg / cg%dvol  ! adding CR
#endif /* COSM_RAYS */
                            endif
                            tbirth = t
                         endif
                         call cg%pset%add(pid, mass, pos, vel, acc, ener, in, phy, out, tbirth, tdyn)
                      endif
                   endif
                enddo
             enddo
          enddo
       enddo
! KICK
       if (kick) then
          pset => cg%pset%first
          do while (associated(pset))
             if (t .lt. pset%pdata%tform + 10) then  
                t1 = t - pset%pdata%tform
                do ifl = 1, flind%fluids
                   pfl => flind%all_fluids(ifl)%fl
                   do i = cg%ijkse(xdim,LO), cg%ijkse(xdim,HI)
                      if (cg%coord(LO,xdim)%r(i) .lt. pset%pdata%pos(1)+cg%dx .and. cg%coord(HI,xdim)%r(i) .gt. pset%pdata%pos(1)-cg%dx) then
                         do j = cg%ijkse(ydim,LO), cg%ijkse(ydim,HI)
                            if (cg%coord(LO,ydim)%r(j) .lt. pset%pdata%pos(2)+cg%dy .and. cg%coord(HI,ydim)%r(j) .gt. pset%pdata%pos(2)-cg%dy) then
                               do k = cg%ijkse(zdim,LO), cg%ijkse(zdim,HI)
                                  if (cg%coord(LO,zdim)%r(k) .lt. pset%pdata%pos(3)+cg%dz .and. cg%coord(HI,zdim)%r(k) .gt. pset%pdata%pos(3)-cg%dz) then
                                     i1 = nint((pset%pdata%pos(1)-cg%coord(CENTER,xdim)%r(i)) / cg%dx)
                                     j1 = nint((pset%pdata%pos(2)-cg%coord(CENTER,ydim)%r(j)) / cg%dy)
                                     k1 = nint((pset%pdata%pos(3)-cg%coord(CENTER,zdim)%r(k)) / cg%dz)
                                     fact = 0.0
                                     if (abs(i1)+abs(j1)+abs(k1) .gt. 0.0) fact = 1.0 / sqrt(real(abs(i1) + abs(j1) + abs(k1)))
                                     if (t1 .lt. 6.5) then
                                        padd = pset%pdata%mass * 1.8 * 10.0**40 *gram * cm /sek * (1/0.5)**0.38 *2*dt/6.5 / cg%dvol / 26  ! see Agertz+2013
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
                                     if (abs(i1) + abs(j1) + abs(k1) == 0) then
                                        if ((t1-dt < 6.5) .and. ((t1+dt) .gt. 6.5)) then    ! Instantaneous injection Agertz
#ifdef THERM
                                           cg%u(pfl%ien,i,j,k) = cg%u(pfl%ien,i,j,k)  + aint(pset%pdata%mass/mass_SN) * (1-0.1*cr_active) * n_SN * 10.0**51 * erg / cg%dvol  ! adding SN energy
#endif /* THERM */
#ifdef COSM_RAYS
                                           if (cr_active > 0.0) cg%w(wna%fi)%arr(iarr_crn(cr_table(icr_H1)),i,j,k) = cg%w(wna%fi)%arr(iarr_crn(cr_table(icr_H1)),i,j,k) + aint(pset%pdata%mass/mass_SN) * 0.1 * n_SN * 10.0**51 * erg / cg%dvol  ! adding CR
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


  subroutine initialize_id()

    use mpisetup,             only: proc, nproc

    dpid = int(1000000000/nproc)
    pid_gen =  proc * dpid
    !maxpid = (proc+1) * dpid
    !print *, proc, pid_gen, maxpid, dpid
    
  end subroutine initialize_id

  subroutine attribute_id(pid)

    use mpisetup,             only: proc, nproc

    integer, intent(out)    :: pid

    dpid = int(1000000000/nproc)
    maxpid = (proc+1) * dpid
    if (pid_gen .gt. maxpid) then
       maxpid = (proc+1) * dpid + 1000000000
    endif

    pid_gen = pid_gen+1
    if (pid_gen .ge. maxpid) then
       print *, 'pool of pid full for proc ', proc, pid_gen, maxpid
       pid_gen = proc * dpid + 1000000000
       maxpid = (proc+1) * dpid + 1000000000
    endif
    pid     = pid_gen
    
    
  end subroutine attribute_id

end module star_formation
