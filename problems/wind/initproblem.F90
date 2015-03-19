! wind problem
!
#include "piernik.h"

module initproblem

   implicit none

   private
   public :: read_problem_par, vel_profile, problem_initial_conditions, problem_pointers

   real :: mdot

   namelist /PROBLEM_CONTROL/  mdot

contains

!-----------------------------------------------------------------------------

   subroutine problem_pointers

      use user_hooks, only: problem_customize_solution

      implicit none

      problem_customize_solution => problem_customize_solution_wind

   end subroutine problem_pointers

!-----------------------------------------------------------------------------

   subroutine read_problem_par

      use dataio_pub, only: nh      ! QA_WARN required for diff_nml
      use mpisetup,   only: rbuff, master, slave, piernik_MPI_Bcast

      implicit none

      mdot = 1.

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

         rbuff(1)  = mdot

      endif

      call piernik_MPI_Bcast(rbuff)

      if (slave) then

         mdot = rbuff(1)

      endif

   end subroutine read_problem_par

!-----------------------------------------------------------------------------

   subroutine vel_profile(r, vel, dens)
   ! velocity profile for isothermal wind

      use fluidindex, only: flind
      use constants,  only: fpi

      implicit none
      real, intent(in)  :: r
      real, intent(out) :: vel, dens
      real, parameter   :: r_smooth = 1.

      vel = 2./(1 + exp(2*(1 - r))) * flind%ion%cs
      dens = mdot/fpi/(r**2 + r_smooth**2)/vel

   end subroutine vel_profile

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

                  call vel_profile(rc, vel, dens)
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
      use global,           only: smalld
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

                  if (rc < 10) then
                     call vel_profile(rc, vel, dens)
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
