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

  implicit none

  private
  public :: SF, initialize_id, attribute_id, pid_gen, dmass_stars

  integer(kind=4)       :: pid_gen, maxpid, dpid, dmass_stars

contains

  subroutine SF(forward)

    use cg_leaves,             only: leaves
    use cg_list,               only: cg_list_element
    use constants,             only: ndims, xdim, ydim, zdim, LO, HI, CENTER, pi, nbdn_n
    use cr_data,               only: icr_H1, cr_table
    use fluidindex,            only: flind
    use fluidtypes,            only: component_fluid
    use func,                  only: operator(.equals.), ekin
    use global,                only: nstep, t, dt
    use grid_cont,             only: grid_container
    use initcosmicrays,        only: iarr_crn, cr_active
    use mpisetup,              only: proc
    use named_array_list,      only: wna, qna
    use particle_types,        only: particle
    use particle_utils,        only: is_part_in_cg
    use thermal,               only: itemp
    use units,                 only: newtong, erg

    logical, intent(in)                                :: forward
    type(cg_list_element), pointer                     :: cgl
    type(grid_container),  pointer                     :: cg
    type(particle), pointer                            :: pset
    class(component_fluid), pointer                    :: pfl
    integer                                            :: ifl, i, j, k
    real                                               :: dens_thr, sf_dens, c_tau_ff, eps_sf, frac, mass_SN
    logical                                            :: fed
    integer(kind=4)                                    :: pid, ig, stage
    real, dimension(ndims)                             :: pos, vel, acc
    real                                               :: mass, ener, tdyn
    logical                                            :: in, phy, out, cond

    if (.not. forward) return
    dens_thr = 0.035
    eps_sf = 0.1
    mass_SN = 100.0
    ig = qna%ind(nbdn_n)
    cgl => leaves%first
    do while (associated(cgl))
       cg => cgl%cg
       do ifl = 1, flind%fluids
          pfl => flind%all_fluids(ifl)%fl
          do i = cg%ijkse(xdim,LO), cg%ijkse(xdim,HI)
             do j = cg%ijkse(ydim,LO), cg%ijkse(ydim,HI)
                do k = cg%ijkse(zdim,LO), cg%ijkse(zdim,HI)
                   if ((cg%u(pfl%idn,i,j,k) .gt. dens_thr) .and. (cg%q(itemp)%arr(i,j,k) .lt. 10**4)) then
                      fed = .false.
                      c_tau_ff = sqrt(3.*pi/(32.*newtong))
                      sf_dens = eps_sf / c_tau_ff * cg%u(pfl%idn,i,j,k)**(3./2.)
                      pset => cg%pset%first
                      do while (associated(pset))
                         if (cg%coord(LO,xdim)%r(i) .lt. pset%pdata%pos(1)+cg%dx .and. cg%coord(HI,xdim)%r(i) .gt. pset%pdata%pos(1)-cg%dx) then
                            if (cg%coord(LO,ydim)%r(j) .lt. pset%pdata%pos(2)+cg%dy .and. cg%coord(HI,ydim)%r(j) .gt. pset%pdata%pos(2)-cg%dy) then
                               if (cg%coord(LO,zdim)%r(k) .lt. pset%pdata%pos(3)+cg%dz .and. cg%coord(HI,zdim)%r(k) .gt. pset%pdata%pos(3)-cg%dz) then

                                  if ((pset%pdata%tform .ge. 0.0) .and. (pset%pdata%mass .lt. 10**5)) then
                                     stage = aint(pset%pdata%mass/mass_SN)
                                     pset%pdata%mass = pset%pdata%mass + sf_dens * cg%dvol * 2*dt
                                     frac = sf_dens * 2*dt / cg%u(pfl%idn, i, j, k)
                                     cg%u(pfl%ien, i, j, k)          = cg%u(pfl%ien, i, j, k) - frac * ekin(cg%u(pfl%imx,i,j,k), cg%u(pfl%imy,i,j,k), cg%u(pfl%imz,i,j,k), cg%u(pfl%idn, i, j, k))
                                     cg%w(wna%fi)%arr(pfl%idn,i,j,k) = (1 - frac) * cg%w(wna%fi)%arr(pfl%idn,i,j,k)
                                     cg%u(pfl%imx:pfl%imz, i, j, k)  = (1 - frac) * cg%u(pfl%imx:pfl%imz, i, j, k)
                                     if (aint(pset%pdata%mass/mass_SN) .gt. stage) then
                                        cg%u(pfl%ien,i,j,k) = cg%u(pfl%ien,i,j,k)  + (aint(pset%pdata%mass/mass_SN) - stage) * (1-0.1*cr_active) * 10.0**51 * erg / cg%dvol ! adding SN energy
#ifdef COSM_RAYS
                                        if (cr_active > 0.0) cg%w(wna%fi)%arr(iarr_crn(cr_table(icr_H1)),i,j,k) = cg%w(wna%fi)%arr(iarr_crn(cr_table(icr_H1)),i,j,k) + (aint(pset%pdata%mass/mass_SN) - stage) * 0.1 * 10.0**51 * erg / cg%dvol  ! adding CR
#endif /* COSM_RAYS */
                                     endif
                                     fed = .true.
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
                         call cg%pset%add(pid, mass, pos, vel, acc, ener, in, phy, out, t, tdyn)
                         cg%u(pfl%ien, i, j, k)          = cg%u(pfl%ien, i, j, k) - frac * ekin(cg%u(pfl%imx,i,j,k), cg%u(pfl%imy,i,j,k), cg%u(pfl%imz,i,j,k), cg%u(pfl%idn, i, j, k))
                         cg%w(wna%fi)%arr(pfl%idn,i,j,k) = (1 - frac) * cg%w(wna%fi)%arr(pfl%idn,i,j,k)
                         cg%u(pfl%imx:pfl%imz, i, j, k)  = (1 - frac) * cg%u(pfl%imx:pfl%imz, i, j, k)
                         if (mass .gt. mass_SN) then
                            cg%u(pfl%ien,i,j,k) = cg%u(pfl%ien,i,j,k)  + aint(mass/mass_SN) * (1-0.1*cr_active) * 10.0**51 * erg / cg%dvol ! adding SN energy
#ifdef COSM_RAYS
                            if (cr_active > 0.0) cg%w(wna%fi)%arr(iarr_crn(cr_table(icr_H1)),i,j,k) = cg%w(wna%fi)%arr(iarr_crn(cr_table(icr_H1)),i,j,k) + aint(mass/mass_SN) * 0.1 * 10.0**51 * erg / cg%dvol  ! adding CR
#endif /* COSM_RAYS */
                         endif
                      endif
                   endif
                enddo
             enddo
          enddo
       enddo
       cgl => cgl%nxt
    enddo
    
  end subroutine SF


  subroutine SF_crit(pfl, cg, i, j, k, tdyn, cond)

    use constants,             only: pi
    use crhelpers,             only: divv_i
    use fluidtypes,            only: component_fluid
    use grid_cont,             only: grid_container
    use named_array_list,      only: wna
    use units,                 only: fpiG, kboltz, mH
    use thermal,               only: calc_tcool, itemp

    logical, intent(out)                      :: cond
    real,    intent(in)                       :: tdyn
    real                                      :: density_thr, G, RJ, tcool, kbgmh, temp
    integer, intent(in)                       :: i, j, k
    type(grid_container), pointer, intent(in) :: cg
    class(component_fluid), pointer           :: pfl
    
    density_thr=0.035!004

    G = fpiG/(4*pi)
    cond = .false.
    if ((abs(cg%z(k)) > 7000) .or. ((cg%x(i)**2+cg%y(j)**2) > 20000**2)) return ! no SF in the stream
    if (cg%w(wna%fi)%arr(pfl%idn,i,j,k) .lt. density_thr) return   ! threshold density
    if (cg%q(divv_i)%arr(i,j,k) .ge. 0) return                     ! convergent flow
    if (cg%w(wna%fi)%arr(pfl%idn,i,j,k) * cg%dvol .lt. 1.2 * 10**6) return   ! part mass > 3 10^5
    temp = cg%q(itemp)%arr(i,j,k)
    RJ = 2.8 * sqrt(temp/1000) * sqrt(3*pi/(32*G*cg%w(wna%fi)%arr(pfl%idn,i,j,k)))
    !print *, 'Jeans mass', RJ!4*pi/3 * RJ**3 * cg%w(wna%fi)%arr(pfl%idn,i,j,k), pfl%cs, 2.8 * sqrt(temp/1000)!, pi/6.0 * pfl%cs**3 / G**(3.0/2) / cg%w(wna%fi)%arr(pfl%idn,i,j,k)**0.5
    if (cg%w(wna%fi)%arr(pfl%idn,i,j,k) * cg%dvol .lt. 4*pi/3 * RJ**3 * cg%w(wna%fi)%arr(pfl%idn,i,j,k)) return !, pi/6.0 * pfl%cs**3 / G**(3.0/2) / cg%w(wna%fi)%arr(pfl%idn,i,j,k)**0.5)) return  !pi/6.0 * pfl%cs**3 / G**(3.0/2) / cg%w(wna%fi)%arr(pfl%idn,i,j,k)**0.5 ) return    ! Jeans mass
    kbgmh  = kboltz / (pfl%gam_1 * mH)
    call calc_tcool(temp, cg%w(wna%fi)%arr(pfl%idn,i,j,k), kbgmh, tcool)
    !print *, tcool, tdyn
    if (tcool .gt. tdyn) return
    cond = .true.
    
  end subroutine SF_crit

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
