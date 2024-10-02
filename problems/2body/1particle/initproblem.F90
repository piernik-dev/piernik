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

! Initial condition for Sedov-Taylor explosion
! Written by: M. Hanasz, March 2006

   use constants, only: dsetnamelen, ndims

   implicit none

   private
   public  :: read_problem_par, problem_initial_conditions, problem_initial_nbody, problem_pointers

   real, dimension(ndims)   :: vel                 !< initial velocity for single-particle test

   namelist /PROBLEM_CONTROL/ vel

   character(len=dsetnamelen), parameter :: apot_n = "apot" !< name of the analytical potential field

contains
!-----------------------------------------------------------------------------
   subroutine problem_pointers
      use user_hooks, only: problem_customize_solution
#ifdef HDF5
      use dataio_user, only: user_vars_hdf5
#endif /* HDF5 */

      implicit none

      problem_customize_solution => list_particles
#ifdef HDF5
      user_vars_hdf5   => particle_error_vars
#endif /* HDF5 */

   end subroutine problem_pointers
!-----------------------------------------------------------------------------
   subroutine read_problem_par

      use bcast,          only: piernik_MPI_Bcast
      use cg_list_global, only: all_cg
      use dataio_pub,     only: nh
      use mpisetup,       only: rbuff, master, slave

      implicit none

      ! namelist default parameter values
      vel         = [0., -1., 0.]

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

         rbuff(1:1+size(vel)-1) = vel

      endif

      call piernik_MPI_Bcast(rbuff)

      if (slave) then

         vel = rbuff(1:1+size(vel)-1)

      endif

      call all_cg%reg_var(apot_n)

   end subroutine read_problem_par
!-----------------------------------------------------------------------------
   subroutine problem_initial_conditions

      use cg_leaves,        only: leaves
      use cg_list,          only: cg_list_element
      use fluidindex,       only: flind

      implicit none

      integer                        :: p
      type(cg_list_element), pointer :: cgl

      cgl => leaves%first
      do while (associated(cgl))
         associate(cg => cgl%cg)
            do p = 1, flind%fluids
               associate(fl => flind%all_fluids(p)%fl)
                  cg%u(fl%idn,RNG) = 0.0
                  cg%u(fl%imx,RNG) = 0.0
                  cg%u(fl%imy,RNG) = 0.0
                  cg%u(fl%imz,RNG) = 0.0
               end associate
            enddo
         end associate
         cgl => cgl%nxt
      enddo

   end subroutine problem_initial_conditions

   subroutine problem_initial_nbody

      use constants,      only: I_ONE
      use domain,         only: dom
      use particle_utils, only: add_part_in_proper_cg

      implicit none

      call add_part_in_proper_cg(I_ONE, 1., dom%C_, vel, [0.0, 0.0, 0.0], 0.0)

   end subroutine problem_initial_nbody

   subroutine particle_error_vars(var, tab, ierrh, cg)

      use dataio_pub,       only: die
      use func,             only: operator(.notequals.)
      use grid_cont,        only: grid_container
      use named_array_list, only: qna

      implicit none

      character(len=*),               intent(in)    :: var
      real, dimension(:,:,:),         intent(inout) :: tab
      integer,                        intent(inout) :: ierrh
      type(grid_container), pointer,  intent(in)    :: cg

      if (.not. qna%exists(apot_n)) call die("[initproblem:maclaurin_error_vars] Cannot find apot_n")

      call compute_particle_potential

      ierrh = 0
      select case (trim(var))
         case ("errp")
            tab(:,:,:) = cg%q(qna%ind(apot_n))%span(cg%ijkse) - cg%sgp(RNG)
         case ("relerr")
            where (cg%q(qna%ind(apot_n))%span(cg%ijkse) .notequals. 0.)
               tab(:,:,:) = cg%sgp(RNG)/cg%q(qna%ind(apot_n))%span(cg%ijkse) -1.
            elsewhere
               tab(:,:,:) = 0.
            endwhere
         case default
            ierrh = -1
      end select

   end subroutine particle_error_vars

   subroutine compute_particle_potential

      use cg_leaves,        only: leaves
      use cg_list,          only: cg_list_element
      use grid_cont,        only: grid_container
      use named_array_list, only: qna
      use particle_types,   only: particle
      use units,            only: newtong

      implicit none

      integer                        :: i, j, k
      integer(kind=4)                :: apot_i
      type(cg_list_element), pointer :: cgl
      type(grid_container),  pointer :: cg
      type(particle), pointer    :: pset

      apot_i = qna%ind(apot_n)
      call leaves%set_q_value(apot_i, 0.)

      cgl => leaves%first
      do while (associated(cgl))
         cg => cgl%cg
         do k = cg%ks, cg%ke
            do j = cg%js, cg%je
               do i = cg%is, cg%ie
#ifdef NBODY
                  pset => cg%pset%first
                  do while (associated(pset))
                     cg%q(apot_i)%arr(i, j, k) = cg%q(apot_i)%arr(i, j, k) - newtong * pset%pdata%mass / &
                          sqrt(sum((pset%pdata%pos - [cg%x(i), cg%y(j), cg%z(k)])**2))
                     pset => pset%nxt
                  enddo
#endif /* NBODY */
               enddo
            enddo
         enddo

         cgl => cgl%nxt
      enddo

   end subroutine compute_particle_potential

   subroutine list_particles(fwd)

      use cg_leaves,  only: leaves
      use cg_list,    only: cg_list_element
      use dataio_pub, only: msg, printinfo
      use global,     only: t
      use particle_types, only: particle

      implicit none

      logical, intent(in) :: fwd

      type(cg_list_element), pointer :: cgl
      type(particle), pointer    :: pset

      cgl => leaves%first
      do while (associated(cgl))

#ifdef NBODY
         pset => cgl%cg%pset%first
         do while (associated(pset))
            write(msg, '(a,i7,a,g12.4,a,g12.4,2(a,3g12.4),a)')"particle(",  pset%pdata%pid, ") t=", t, " m=", pset%pdata%mass, " @[", pset%pdata%pos, "] ->[", pset%pdata%vel, " ]"
            if (pset%pdata%outside) msg = trim(msg) // " (is outside)"
            call printinfo(msg)
            pset => pset%nxt
         enddo
#endif /* NBODY */

         cgl => cgl%nxt
      enddo

      if (.false.) write(msg, '(l1)') fwd

   end subroutine list_particles

end module initproblem
