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
  public :: SF, initialize_id, attribute_id, pid_gen

  integer(kind=4)       :: pid_gen, maxpid, dpid

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
    real                                               :: mass, ener, tdyn
    logical                                            :: in, phy, out, cond

    
    if (.not. forward) return
    G = fpiG/(4*pi)
    ig = qna%ind(nbdn_n)
    thresh = 0.005
    cgl => leaves%first
    do while (associated(cgl))
       cg => cgl%cg

       do ifl = 1, flind%fluids
          pfl => flind%all_fluids(ifl)%fl
          !dens =>  cg%w(wna%fi)%span(pfl%idn,cg%ijkse)
          do i = cg%ijkse(xdim,LO), cg%ijkse(xdim,HI)
             do j = cg%ijkse(ydim,LO), cg%ijkse(ydim,HI)
                do k = cg%ijkse(zdim,LO), cg%ijkse(zdim,HI)
                   call SF_crit(pfl, cg, i, j, k, cond)
                   if (cond) then
                      !if (cg%w(wna%fi)%arr(pfl%idn,i,j,k) .gt. thresh) then
                      call attribute_id(pid)
                      !pid = (nstep+1)  + i*10 + j*2 +k + proc    !  Do better
                      pos(1) = cg%coord(CENTER, xdim)%r(i)
                      pos(2) = cg%coord(CENTER, xdim)%r(j)
                      pos(3) = cg%coord(CENTER, xdim)%r(k)
                      mass = 0.1 * cg%w(wna%fi)%arr(pfl%idn,i,j,k) * cg%dvol
                      vel(1) = cg%u(pfl%imx, i, j, k) / cg%w(wna%fi)%arr(pfl%idn,i,j,k)
                      vel(2) = cg%u(pfl%imy, i, j, k) / cg%w(wna%fi)%arr(pfl%idn,i,j,k)
                      vel(3) = cg%u(pfl%imz, i, j, k) / cg%w(wna%fi)%arr(pfl%idn,i,j,k)
                      acc(:)=0.0
                      ener = 0.0
                      call is_part_in_cg(cg, pos, in, phy, out)
                      !print *, 'SF!', cg%u(pfl%ien, i, j, k), 0.1 * ekin(cg%u(pfl%imx,i,j,k), cg%u(pfl%imy,i,j,k), cg%u(pfl%imz,i,j,k), cg%u(pfl%idn,i,j,k))
                      tdyn = sqrt(3*pi/(32*G*(cg%w(wna%fi)%arr(pfl%idn,i,j,k))+cgl%cg%q(ig)%arr(i,j,k)))
                      print *, proc, pid, tdyn
                      call cg%pset%add(pid, mass, pos, vel, acc, ener, in, phy, out, t, tdyn)
                      cg%u(pfl%ien, i, j, k)          = cg%u(pfl%ien, i, j, k) - 0.1 * ekin(cg%u(pfl%imx,i,j,k), cg%u(pfl%imy,i,j,k), cg%u(pfl%imz,i,j,k), cg%u(pfl%idn,i,j,k))
                      cg%w(wna%fi)%arr(pfl%idn,i,j,k) = 0.9 * cg%w(wna%fi)%arr(pfl%idn,i,j,k)
                      cg%u(pfl%imx:pfl%imz, i, j, k)  = 0.9 * cg%u(pfl%imx:pfl%imz, i, j, k)
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
    use constants,             only: ndims, xdim, ydim, zdim, LO, HI
    use cr_data,               only: icr_H1, cr_table
    use fluidindex,            only: flind
    use fluidtypes,            only: component_fluid
    use global,                only: nstep, t, dt
    use grid_cont,             only: grid_container
    use initcosmicrays,        only: iarr_crn
    use mpisetup,              only: proc
    use named_array_list,      only: wna
    use particle_types,        only: particle
    use units,                 only: clight

    type(cg_list_element), pointer                     :: cgl
    type(grid_container),  pointer                     :: cg
    type(particle), pointer                            :: pset
    class(component_fluid), pointer                    :: pfl
    integer                                            :: ifl, i, j, k
    real                                               :: thresh
    integer(kind=4)                                    :: pid
    real, dimension(ndims)                             :: pos, vel, acc
    real                                               :: mass, ener, t1, tdyn, msf


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
             msf = pset%pdata%mass * ( (1+t1/tdyn) * exp(-t1/tdyn) - (1+(t1+2*dt)/tdyn) * exp(-(t1+2*dt)/tdyn))
             !pset%pdata%mass = pset%pdata%mass - 0.25*msf
             do ifl = 1, flind%fluids
                pfl => flind%all_fluids(ifl)%fl
                do i = cg%ijkse(xdim,LO), cg%ijkse(xdim,HI)
                   if (cg%coord(LO,xdim)%r(i) .lt. pset%pdata%pos(1) .and. cg%coord(HI,xdim)%r(i) .gt. pset%pdata%pos(1)) then
                      do j = cg%ijkse(ydim,LO), cg%ijkse(ydim,HI)
                         if (cg%coord(LO,ydim)%r(j) .lt. pset%pdata%pos(2) .and. cg%coord(HI,ydim)%r(j) .gt. pset%pdata%pos(2)) then
                            do k = cg%ijkse(zdim,LO), cg%ijkse(zdim,HI)
                               if (cg%coord(LO,zdim)%r(k) .lt. pset%pdata%pos(3) .and. cg%coord(HI,zdim)%r(k) .gt. pset%pdata%pos(3)) then
                                  cg%w(wna%fi)%arr(pfl%idn,i,j,k) = cg%w(wna%fi)%arr(pfl%idn,i,j,k) + 0.25 * msf / cg%dvol                               ! adding mass
                                  cg%u(pfl%imx:pfl%imz, i, j, k)  = cg%u(pfl%imx:pfl%imz, i, j, k)  + 0.25 * msf / cg%dvol * pset%pdata%vel(:)           ! adding momentum
                                  !print *, 'cr', cg%w(wna%fi)%arr(iarr_crn(cr_table(icr_H1)),i,j,k), 0.1 * 0.00001 * msf / cg%dvol * clight**2 
                                  cg%w(wna%fi)%arr(iarr_crn(cr_table(icr_H1)),i,j,k) = cg%w(wna%fi)%arr(iarr_crn(cr_table(icr_H1)),i,j,k) + 0.1 * 0.00001 * msf / cg%dvol * clight**2 
                                  cg%u(pfl%ien,i,j,k)             = cg%u(pfl%ien,i,j,k)             + 0.25 * msf / cg%dvol * sum(pset%pdata%vel(:)**2)   ! adding kinetic energy
                                  print *, pset%pdata%mass, pset%pdata%tform, msf, cg%u(pfl%ien,i,j,k), 0.00001 * msf / cg%dvol * clight**2
                                  cg%u(pfl%ien,i,j,k)              = cg%u(pfl%ien,i,j,k)              + 0.00001 * msf / cg%dvol * clight**2              ! adding SN energy
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


  subroutine SF_crit(pfl, cg, i, j, k, cond)

    use constants,             only: pi
    use crhelpers,             only: divv_i
    use fluidtypes,            only: component_fluid
    use grid_cont,             only: grid_container
    use named_array_list,      only: wna
    use units,                 only: fpiG

    logical, intent(out)                      :: cond
    real                                      :: density_thr, G
    integer, intent(in)                       :: i, j, k
    type(grid_container), pointer, intent(in) :: cg
    class(component_fluid), pointer           :: pfl
    
    density_thr=0.005

    G = fpiG/(4*pi)
    cond = .false.
    if (cg%w(wna%fi)%arr(pfl%idn,i,j,k) .lt. density_thr) return   ! threshold density
    if (cg%q(divv_i)%arr(i,j,k) .ge. 0) return                     ! convergent flow
    if (cg%w(wna%fi)%arr(pfl%idn,i,j,k) * cg%dvol .lt. 3.0 * 10**6) return   ! part mass > 3 10^5
    if (cg%w(wna%fi)%arr(pfl%idn,i,j,k) * cg%dvol .lt. pi/6.0 * pfl%cs**3 / G**(3.0/2) / cg%w(wna%fi)%arr(pfl%idn,i,j,k)**0.5 ) return    ! Jeans mass

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
