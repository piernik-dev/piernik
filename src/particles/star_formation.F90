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
    use fluidindex,            only: flind
    use fluidtypes,            only: component_fluid
    use func,                  only: ekin
    use global,                only: nstep, t
    use grid_cont,             only: grid_container
    use mpisetup,              only: proc
    use named_array_list,      only: wna, qna
    !use particle_types,        only: particle
    use particle_utils,        only: is_part_in_cg
    use units,                 only: fpiG

    logical, intent(in)                                :: forward
    type(cg_list_element), pointer                     :: cgl
    type(grid_container),  pointer                     :: cg
    !type(particle), pointer                            :: pset
    class(component_fluid), pointer                    :: pfl
    !real, dimension(:,:,:), pointer                    :: dens
    integer                                            :: ifl, i, j, k
    real                                               :: thresh, G
    integer(kind=4)                                    :: pid, ig
    real, dimension(ndims)                             :: pos, vel, acc
    real                                               :: mass, ener, tdyn, frac
    logical                                            :: in, phy, out, cond

    dmass_stars = 0.0
    if (.not. forward) return
    G = fpiG/(4*pi)
    ig = qna%ind(nbdn_n)
    !thresh = 0.005
    frac = 0.1
    cgl => leaves%first
    do while (associated(cgl))
       cg => cgl%cg

       do ifl = 1, flind%fluids
          pfl => flind%all_fluids(ifl)%fl
          !dens =>  cg%w(wna%fi)%span(pfl%idn,cg%ijkse)
          do i = cg%ijkse(xdim,LO), cg%ijkse(xdim,HI)
             do j = cg%ijkse(ydim,LO), cg%ijkse(ydim,HI)
                do k = cg%ijkse(zdim,LO), cg%ijkse(zdim,HI)
                   tdyn = sqrt(3*pi/(32*G*(cg%w(wna%fi)%arr(pfl%idn,i,j,k))+cgl%cg%q(ig)%arr(i,j,k)))
                   call SF_crit(pfl, cg, i, j, k, tdyn, cond)
                   if (cond) then
                      call attribute_id(pid)
                      pos(1) = cg%coord(CENTER, xdim)%r(i)
                      pos(2) = cg%coord(CENTER, ydim)%r(j)
                      pos(3) = cg%coord(CENTER, zdim)%r(k)
                      mass = frac * cg%w(wna%fi)%arr(pfl%idn,i,j,k) * cg%dvol !* 0.75
                      vel(1) = cg%u(pfl%imx, i, j, k) / cg%w(wna%fi)%arr(pfl%idn,i,j,k)
                      vel(2) = cg%u(pfl%imy, i, j, k) / cg%w(wna%fi)%arr(pfl%idn,i,j,k)
                      vel(3) = cg%u(pfl%imz, i, j, k) / cg%w(wna%fi)%arr(pfl%idn,i,j,k)
                      acc(:)=0.0
                      ener = 0.0
                      call is_part_in_cg(cg, pos, in, phy, out)
                      !print *, 'SF!', cg%u(pfl%ien, i, j, k), 0.1 * ekin(cg%u(pfl%imx,i,j,k), cg%u(pfl%imy,i,j,k), cg%u(pfl%imz,i,j,k), cg%u(pfl%idn,i,j,k))
                      !print *, proc, pid, tdyn
                      call cg%pset%add(pid, mass, pos, vel, acc, ener, in, phy, out, t, tdyn)
                      cg%u(pfl%ien, i, j, k)          = cg%u(pfl%ien, i, j, k) - frac * ekin(cg%u(pfl%imx,i,j,k), cg%u(pfl%imy,i,j,k), cg%u(pfl%imz,i,j,k), cg%u(pfl%idn,i,j,k))! * 0.75
                      cg%w(wna%fi)%arr(pfl%idn,i,j,k) = (1-frac) * cg%w(wna%fi)%arr(pfl%idn,i,j,k)
                      cg%u(pfl%imx:pfl%imz, i, j, k)  = (1-frac) * cg%u(pfl%imx:pfl%imz, i, j, k)
                      dmass_stars = dmass_stars + mass
                   endif
                enddo
             enddo
          enddo
       enddo
       cgl => cgl%nxt
    enddo
    call feedback()

  end subroutine SF

  subroutine feedback

    use cg_leaves,             only: leaves
    use cg_list,               only: cg_list_element
    use constants,             only: ndims, xdim, ydim, zdim, LO, HI, CENTER
    use cr_data,               only: icr_H1, cr_table
    use fluidindex,            only: flind
    use fluidtypes,            only: component_fluid
    use func,                  only: operator(.equals.)
    use global,                only: nstep, t, dt
    use grid_cont,             only: grid_container
    use initcosmicrays,        only: iarr_crn, cr_active
    use mpisetup,              only: proc
    use named_array_list,      only: wna
    use particle_types,        only: particle
    use units,                 only: clight, gram, cm, sek

    type(cg_list_element), pointer                     :: cgl
    type(grid_container),  pointer                     :: cg
    type(particle), pointer                            :: pset
    class(component_fluid), pointer                    :: pfl
    integer                                            :: ifl, i, j, k, i1, j1, k1
    real                                               :: thresh
    integer(kind=4)                                    :: pid
    real, dimension(ndims)                             :: pos, vel, acc
    real                                               :: mass, ener, t1, tdyn, msf, fact, padd


    cgl => leaves%first
    do while (associated(cgl))
       cg => cgl%cg

       pset => cg%pset%first
       do while (associated(pset))
          !print*, pset%pdata%pid,  pset%pdata%tform
          !cycle
          if (t .lt. pset%pdata%tform + 120.0) then ! + 12*pset%pdata%tdyn?
             t1 = t - pset%pdata%tform
             tdyn = pset%pdata%tdyn
             msf = 0.0
             if (t1 .gt. -1*6.5) then
                msf = pset%pdata%mass * ( (1+t1/tdyn) * exp(-t1/tdyn) - (1+(t1+2*dt)/tdyn) * exp(-(t1+2*dt)/tdyn))
                !msf = pset%pdata%mass
                pset%pdata%mass = pset%pdata%mass - 0.25*msf
             endif
             do ifl = 1, flind%fluids
                pfl => flind%all_fluids(ifl)%fl
                do i = cg%ijkse(xdim,LO), cg%ijkse(xdim,HI)
                   if (cg%coord(LO,xdim)%r(i) .lt. pset%pdata%pos(1)+cg%dx .and. cg%coord(HI,xdim)%r(i) .gt. pset%pdata%pos(1)-cg%dx) then
                      do j = cg%ijkse(ydim,LO), cg%ijkse(ydim,HI)
                         if (cg%coord(LO,ydim)%r(j) .lt. pset%pdata%pos(2)+cg%dy .and. cg%coord(HI,ydim)%r(j) .gt. pset%pdata%pos(2)-cg%dy) then
                            do k = cg%ijkse(zdim,LO), cg%ijkse(zdim,HI)
                               if (cg%coord(LO,zdim)%r(k) .lt. pset%pdata%pos(3)+cg%dz .and. cg%coord(HI,zdim)%r(k) .gt. pset%pdata%pos(3)-cg%dz) then
                                  ! Agertz
                                  i1 = (pset%pdata%pos(1)-cg%coord(CENTER,xdim)%r(i)) / cg%dx
                                  j1 = (pset%pdata%pos(2)-cg%coord(CENTER,ydim)%r(j)) / cg%dy
                                  k1 = (pset%pdata%pos(3)-cg%coord(CENTER,zdim)%r(k)) / cg%dz
                                  !print *, pset%pdata%pos(1), cg%coord(CENTER,xdim)%r(i), i, i1
                                  if (abs(i1)+abs(j1)+abs(k1) .gt. 0.0) fact = 1.0 / sqrt(real(abs(i1) + abs(j1) + abs(k1)))
                                  if (t1 .lt. 6.5) then
                                     padd = pset%pdata%mass * 2.0 * 10.0**40 *gram * cm /sek / cg%dvol / 26.0 * 2*dt/6.5
                                  else
                                     padd =  36000 * pset%pdata%mass/42 / cg%dvol / 26.0 * 2*dt/113.5     ! should use initial mass, not current mass
                                  endif
                                  
                                  ! Momentum kick
                                  cg%u(pfl%ien,i,j,k) = cg%u(pfl%ien,i,j,k) - 0.5*(cg%u(pfl%imx,i,j,k)**2 +cg%u(pfl%imy,i,j,k)**2 + cg%u(pfl%imz,i,j,k)**2)/cg%u(pfl%idn,i,j,k)  ! remove ekin
                                  cg%u(pfl%imx,i,j,k) = cg%u(pfl%imx,i,j,k) + fact * i1 * padd
                                  cg%u(pfl%imy,i,j,k) = cg%u(pfl%imy,i,j,k) + fact * j1 * padd
                                  cg%u(pfl%imz,i,j,k) = cg%u(pfl%imz,i,j,k) + fact * k1 * padd
                                  cg%u(pfl%ien,i,j,k) = cg%u(pfl%ien,i,j,k) + 0.5*(cg%u(pfl%imx,i,j,k)**2 +cg%u(pfl%imy,i,j,k)**2 + cg%u(pfl%imz,i,j,k)**2)/cg%u(pfl%idn,i,j,k)  ! add new ekin

                                  !Butsky
                                  if (abs(i1) + abs(j1) + abs(k1) == 0) then
                                     if (t1 .gt. 6.5) then
                                        cg%w(wna%fi)%arr(pfl%idn,i,j,k) = cg%w(wna%fi)%arr(pfl%idn,i,j,k) + 0.25 * msf / cg%dvol                               ! adding mass
                                        cg%u(pfl%imx:pfl%imz, i, j, k)  = cg%u(pfl%imx:pfl%imz, i, j, k)  + 0.25 * msf / cg%dvol * pset%pdata%vel(:)           ! adding momentum
                                        cg%u(pfl%ien,i,j,k)             = cg%u(pfl%ien,i,j,k)             + 0.25 * msf / cg%dvol * sum(pset%pdata%vel(:)**2)   ! adding kinetic energy
                                        cg%u(pfl%ien,i,j,k)             = cg%u(pfl%ien,i,j,k)  + 0.00001 * 0.25 * (1 - 0.9*cr_active) * msf / cg%dvol * clight**2 
                                        if (cr_active > 0.0) cg%w(wna%fi)%arr(iarr_crn(cr_table(icr_H1)),i,j,k) = cg%w(wna%fi)%arr(iarr_crn(cr_table(icr_H1)),i,j,k) + 0.1 * 0.00001 * 0.25 *  msf / cg%dvol * clight**2
                                     endif
                                     if ((t1-dt < 6.5) .and. ((t1+dt) .gt. 6.5)) then
                                        cg%u(pfl%ien,i,j,k) = cg%u(pfl%ien,i,j,k)  + 0.75 * 0.00001 * 0.25 * (pset%pdata%mass + 0.25*msf) / cg%dvol * clight**2   ! adding SN energy
                                        !print *, "SN ENERGY INJECTED!", 0.00001 * (pset%pdata%mass + 0.25*msf) / cg%dvol * clight**2, pset%pdata%mass + 0.25*msf
#ifdef COSM_RAYS
                                        if (cr_active > 0.0) cg%w(wna%fi)%arr(iarr_crn(cr_table(icr_H1)),i,j,k) = cg%w(wna%fi)%arr(iarr_crn(cr_table(icr_H1)),i,j,k) + 0.1 * 0.00001 * 0.25 * (pset%pdata%mass + 0.25*msf) / cg%dvol * clight**2   ! adding CR
#endif /* COSM_RAYS */
                                     endif
                                  endif
                                  !print *, 'Feedbacking!', pset%pdata%tform
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

       cgl => cgl%nxt
    enddo
    
  end subroutine feedback


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
    RJ = 2.8 * sqrt(temp/1000) * sqrt(3*pi/(32*G*cg%w(wna%fi)%arr(pfl%idn,i,j,k))) 
    !print *, 'Jeans mass', 4*pi/3 * RJ**3 * cg%w(wna%fi)%arr(pfl%idn,i,j,k), pfl%cs, 2.8 * sqrt(temp/1000)!, pi/6.0 * pfl%cs**3 / G**(3.0/2) / cg%w(wna%fi)%arr(pfl%idn,i,j,k)**0.5
    if (cg%w(wna%fi)%arr(pfl%idn,i,j,k) * cg%dvol .lt. 4*pi/3 * RJ**3 * cg%w(wna%fi)%arr(pfl%idn,i,j,k)) return !, pi/6.0 * pfl%cs**3 / G**(3.0/2) / cg%w(wna%fi)%arr(pfl%idn,i,j,k)**0.5)) return  !pi/6.0 * pfl%cs**3 / G**(3.0/2) / cg%w(wna%fi)%arr(pfl%idn,i,j,k)**0.5 ) return    ! Jeans mass
    kbgmh  = kboltz / (pfl%gam_1 * mH)
    temp = cg%q(itemp)%arr(i,j,k)
    call calc_tcool(temp, cg%w(wna%fi)%arr(pfl%idn,i,j,k), kbgmh, tcool)
    !print *, tcool, tdyn
    if (tcool .gt. tdyn) return
    !print *, 'yay', cg%w(wna%fi)%arr(pfl%idn,i,j,k) * cg%dvol, pi/6.0 * pfl%cs**3 / fpiG**(3.0/2) / cg%w(wna%fi)%arr(pfl%idn,i,j,k)**0.5 * (4*pi)**(3.0/2)
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
