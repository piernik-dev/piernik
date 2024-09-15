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

module initproblem

! wind accretion problem

   implicit none

   private
   public :: read_problem_par, wind_profile, problem_initial_conditions, problem_pointers

   real :: mstar, mdot, rin, rdust, vel_scale, tburst, dburst, vburst, mburst

   namelist /PROBLEM_CONTROL/  mstar, mdot, rin, rdust, vel_scale, tburst, dburst, vburst, mburst

contains

!-----------------------------------------------------------------------------

   subroutine problem_pointers

      use user_hooks, only: problem_customize_solution

      implicit none

      problem_customize_solution => problem_customize_solution_wind

   end subroutine problem_pointers

!-----------------------------------------------------------------------------

   subroutine read_problem_par

      use bcast,      only: piernik_MPI_Bcast
      use dataio_pub, only: nh
      use mpisetup,   only: rbuff, master, slave

      implicit none

      mstar     = 0. !< stellar mass
      mdot      = 0. !< mass loss
      rin       = 0. !< inner boundary for wind solution
      rdust     = 0. !< outer boundary for customized wind solution
      vel_scale = 1. !< velocity factor
      tburst    = 0. !< burst initial time
      dburst    = 0. !< burst duration since tburst
      vburst    = 0. !< velocity factor during burst
      mburst    = 0. !< density factor during burst

      if (master) then

         if (.not.nh%initialized) call nh%init()
         open(newunit=nh%lun, file=nh%tmp1, status="unknown")
         write(nh%lun,nml=PROBLEM_CONTROL)
         close(nh%lun)
         open(newunit=nh%lun, file=nh%par_file)
         nh%errstr=""
         read(unit=nh%lun, nml=PROBLEM_CONTROL, iostat=nh%ierrh, iomsg=nh%errstr)
         close(nh%lun)
         call nh%namelist_errh(nh%ierrh, "PROBLEM_CONTROL")
         read(nh%cmdl_nml,nml=PROBLEM_CONTROL, iostat=nh%ierrh)
         call nh%namelist_errh(nh%ierrh, "PROBLEM_CONTROL", .true.)
         open(newunit=nh%lun, file=nh%tmp2, status="unknown")
         write(nh%lun,nml=PROBLEM_CONTROL)
         close(nh%lun)
         call nh%compare_namelist()

         rbuff(1)  = mstar
         rbuff(2)  = mdot
         rbuff(3)  = rin
         rbuff(4)  = rdust
         rbuff(5)  = vel_scale
         rbuff(6)  = tburst
         rbuff(7)  = dburst
         rbuff(8)  = vburst
         rbuff(9)  = mburst

      endif

      call piernik_MPI_Bcast(rbuff)

      if (slave) then

         mstar     = rbuff(1)
         mdot      = rbuff(2)
         rin       = rbuff(3)
         rdust     = rbuff(4)
         vel_scale = rbuff(5)
         tburst    = rbuff(6)
         dburst    = rbuff(7)
         vburst    = rbuff(8)
         mburst    = rbuff(9)

      endif

   end subroutine read_problem_par

!-----------------------------------------------------------------------------

   subroutine wind_profile(r, vel, dens)
   ! velocity and density profiles for isothermal wind

      use fluidindex, only: flind
      use constants,  only: fpi ! four*pi
      use units,      only: newtong

      implicit none
      real, intent(in)  :: r
      real, intent(out) :: vel, dens
      real              :: r_c, vel0

      ! r_c = G*M/2/c_s**2
      r_c = newtong*mstar/2./flind%ion%cs2
      vel = 2. * vel_scale/(1 + exp(2*(1 - r/r_c))) * flind%ion%cs
      if (r > rin) then
         dens = mdot/fpi/r**2/vel
      else
         vel0 = 2. * vel_scale/(1 + exp(2*(1 - rin/r_c))) * flind%ion%cs
         dens = mdot/fpi/rin**2/vel0
      endif

   end subroutine wind_profile

!-----------------------------------------------------------------------------

   subroutine problem_initial_conditions

      use cg_leaves,   only: leaves
      use cg_list,     only: cg_list_element
      use constants,   only: xdim, ydim, zdim, LO, HI
      use fluidindex,  only: flind
      use fluidtypes,  only: component_fluid
      use global,      only: smalld
      use grid_cont,   only: grid_container

      implicit none

      class(component_fluid), pointer :: fl
      integer                         :: i, j, k
      real                            :: xi, yj, vx, vy, vz, zk, rc, vel,&
                                       dens, phi, theta
      type(cg_list_element),  pointer :: cgl
      type(grid_container),   pointer :: cg

      fl => flind%ion

      cgl => leaves%first
      do while (associated(cgl))
         cg => cgl%cg

         do i = cg%lhn(xdim,LO), cg%lhn(xdim,HI)
            xi = cg%x(i)
            do j = cg%lhn(ydim,LO), cg%lhn(ydim,HI)
               yj = cg%y(j)
               do k = cg%lhn(zdim,LO), cg%lhn(zdim,HI)
                  zk = cg%z(k)
                  rc = sqrt(xi**2 + yj**2 + zk**2)

                  call wind_profile(rc, vel, dens)

                  cg%u(fl%idn,i,j,k) = max(dens, smalld)
                  if (fl%ien > 1) then
                     cg%u(fl%ien,i,j,k) = fl%cs2/(fl%gam_1)*cg%u(fl%idn,i,j,k)
                     cg%u(fl%ien,i,j,k) = cg%u(fl%ien,i,j,k) + 0.5*vel**2*cg%u(fl%idn,i,j,k)
                  endif


                  phi = atan2(yj, xi)
                  theta = acos(zk/rc)
                  vx = vel*sin(theta)*cos(phi)
                  vy = vel*sin(theta)*sin(phi)
                  vz = vel*cos(theta)

                  cg%u(fl%imx,i,j,k) = vx*cg%u(fl%idn,i,j,k)
                  cg%u(fl%imy,i,j,k) = vy*cg%u(fl%idn,i,j,k)
                  cg%u(fl%imz,i,j,k) = vz*cg%u(fl%idn,i,j,k)

               enddo
            enddo
         enddo

         cgl => cgl%nxt
      enddo

   end subroutine problem_initial_conditions

   subroutine problem_customize_solution_wind(forward)
      ! impose wind inflow after sweeps and modify the fluid data

      use cg_leaves,        only: leaves
      use cg_list,          only: cg_list_element
      use constants,        only: xdim, ydim, zdim, LO, HI
      use global,           only: smalld, t
      use grid_cont,        only: grid_container
      use all_boundaries,   only: all_fluid_boundaries
      use fluidindex,       only: flind
      use fluidtypes,       only: component_fluid

      implicit none

      class(component_fluid), pointer :: fl
      logical, intent(in)             :: forward
      integer                         :: i, j, k
      type(cg_list_element), pointer  :: cgl
      type(grid_container),  pointer  :: cg
      real                            :: xi, yj, vx, vy, vz, zk, rc, vel,&
                                       dens, phi, theta

      fl => flind%ion

      cgl => leaves%first
      do while (associated(cgl))
         cg => cgl%cg

         do i = cg%lhn(xdim,LO), cg%lhn(xdim,HI)
            xi = cg%x(i)
            do j = cg%lhn(ydim,LO), cg%lhn(ydim,HI)
               yj = cg%y(j)
               do k = cg%lhn(zdim,LO), cg%lhn(zdim,HI)
                  zk = cg%z(k)
                  rc = sqrt(xi**2 + yj**2 + zk**2)

                  if (rc < rdust) then
                     call wind_profile(rc, vel, dens)
                     if ((tburst > 0.) .and. (t > tburst) .and. (t < tburst + dburst)) then
                        vel = vel*vburst
                        dens = dens*mburst
                     endif

                     cg%u(fl%idn,i,j,k) = max(dens, smalld)

                     phi = atan2(yj, xi)
                     theta = acos(zk/rc)
                     vx = vel*sin(theta)*cos(phi)
                     vy = vel*sin(theta)*sin(phi)
                     vz = vel*cos(theta)

                     cg%u(fl%imx,i,j,k) = vx*cg%u(fl%idn,i,j,k)
                     cg%u(fl%imy,i,j,k) = vy*cg%u(fl%idn,i,j,k)
                     cg%u(fl%imz,i,j,k) = vz*cg%u(fl%idn,i,j,k)
                  endif

               enddo
            enddo
         enddo

         cgl => cgl%nxt
      enddo

      call all_fluid_boundaries

      return
      if (forward) i = j ! suppress compiler warnings on unused arguments

   end subroutine problem_customize_solution_wind

end module initproblem
