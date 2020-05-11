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
!! \brief Diagnostics of particle module
!!
!! This module contains particle diagnostics related things that are needed by I/O and debug routines.
!<

module particle_utils
! pulled by GRAV

   implicit none

   private
   public :: max_pvel_1d, add_part_in_proper_cg, print_all_particles, is_part_in_cg
   public :: max_pacc_3d, particle_diagnostics, twodtscheme, dump_diagnose, tot_energy, d_energy, tot_angmom, d_angmom, count_all_particles, part_leave_cg

   real    :: tot_angmom           !< angular momentum of set of the particles
   real    :: tot_energy           !< total energy of set of the particles
   real    :: d_energy             !< error of energy of set of the particles in succeeding timesteps
   real    :: d_angmom             !< error of angular momentum in succeeding timensteps
   logical :: twodtscheme
   logical :: dump_diagnose        !< dump diagnose for each particle to a seperate log file
   integer, parameter :: npf = 12  !< number of single particle fields

contains

   subroutine print_all_particles

      !use cg_leaves, only: leaves
      !use cg_list,   only: cg_list_element

      implicit none

      !type(cg_list_element), pointer :: cgl

      !!!!!!!!!Does nothing for now

      !cgl => leaves%first
      !do while (associated(cgl))
      !   call cgl%cg%pset%print
      !   cgl => cgl%nxt
      !enddo

   end subroutine print_all_particles

   subroutine max_pvel_1d(cg, max_v)

      use constants, only: ndims, xdim, zdim, zero
      use grid_cont, only: grid_container
      use particle_types, only: particle

      implicit none

      type(grid_container), pointer, intent(in)  :: cg
      type(particle), pointer                    :: pset
      real, dimension(ndims),        intent(out) :: max_v
      integer(kind=4)                            :: cdim
      real                                       :: v_tmp

      !Better way to do this? Was easier with arrays
      max_v = zero
      do cdim = xdim, zdim
         pset => cg%pset%first
         do while (associated(pset))
            v_tmp = abs(pset%pdata%vel(cdim))
            if (v_tmp > max_v(cdim)) max_v(cdim) = v_tmp
            pset => pset%nxt
         enddo
      enddo

   end subroutine max_pvel_1d

   subroutine max_pacc_3d(cg, max_pacc)

      use constants, only: big, CENTER, half, LO, xdim, zdim, zero
      use grid_cont, only: grid_container
      use mpisetup,  only: proc
      use types,     only: value
      use particle_types, only: particle

      implicit none

      type(grid_container), pointer, intent(in)  :: cg
      type(particle), pointer                    :: pset
      type(value),                   intent(out) :: max_pacc
      real                                       :: acc2, max_acc
      integer(kind=4) :: cdim

      max_pacc%assoc = big

      max_acc  = zero
      pset => cg%pset%first
      do while (associated(pset))
         acc2 = sum(pset%pdata%acc(:)**2)
         if (acc2 > max_acc) then
            max_acc = acc2
            max_pacc%coords(:) = pset%pdata%pos(:)
            !max_pacc%proc = i !> \todo it might be an information about extremum particle, but the scheme of log file is to print the process number
         endif
         pset => pset%nxt
      enddo
      max_pacc%val = sqrt(max_acc)
      do cdim = xdim, zdim
         max_pacc%loc(cdim) = int( half + (max_pacc%coords(cdim) - cg%coord(CENTER,cdim)%r(cg%ijkse(cdim, LO))) * cg%idl(cdim) )
      enddo
      max_pacc%proc = proc

   end subroutine max_pacc_3d

   subroutine particle_diagnostics

#ifdef VERBOSE
      use dataio_pub, only: msg, printinfo
#endif /* VERBOSE */

      implicit none

      real, save    :: init_energy          !< total initial energy of set of the particles
      real, save    :: init_angmom          !< initial angular momentum of set of the particles
      logical, save :: first_run_lf = .true.

      call get_angmom_totener(tot_angmom, tot_energy)

      if (first_run_lf) then
         init_energy = tot_energy
         init_angmom = tot_angmom
         first_run_lf = .false.
      endif
      d_energy = log_error(tot_energy, init_energy)
      d_angmom = log_error(tot_angmom, init_angmom)

#ifdef VERBOSE
      write(msg,'(a,3(1x,e12.5))') '[particle_utils:particle_diagnostics] Total energy: initial, current, error ', init_energy, tot_energy, d_energy
      call printinfo(msg)
      write(msg,'(a,3(1x,e12.5))') '[particle_utils:particle_diagnostics] ang_momentum: initial, current, error ', init_angmom, tot_angmom, d_angmom
      call printinfo(msg)
#endif /* VERBOSE */

      if (dump_diagnose) call dump_particles_to_textfile

   end subroutine particle_diagnostics

   real function log_error(vcurr,vinit)

      use constants, only: zero
      use func,      only: operator(.equals.)

      implicit none

      real, intent(in) :: vcurr, vinit

      if (vinit .equals. zero) then
         log_error = vcurr
         return
      endif
      log_error = (vcurr - vinit) / vinit
      if (log_error .equals. zero) then
         log_error = zero
      else
         log_error = log(abs(log_error))
      endif

   end function log_error

   subroutine get_angmom_totener(ang_momentum, total_energy)

      use cg_leaves, only: leaves
      use cg_list,   only: cg_list_element
      use constants, only: xdim, ydim, zdim, zero
      use particle_types, only: particle

      implicit none

      real, intent(out)              :: ang_momentum
      real, intent(out)              :: total_energy !< total energy of set of particles
      real                           :: L1, L2, L3
      type(cg_list_element), pointer :: cgl
      type(particle), pointer        :: pset

      ang_momentum = zero
      total_energy = zero

      cgl => leaves%first
      do while (associated(cgl))

         pset => cgl%cg%pset%first
         do while (associated(pset))
            associate( part => pset%pdata )
               L1 = part%pos(ydim) * part%vel(zdim) - part%pos(zdim) * part%vel(ydim)
               L2 = part%pos(zdim) * part%vel(xdim) - part%pos(xdim) * part%vel(zdim)
               L3 = part%pos(xdim) * part%vel(ydim) - part%pos(ydim) * part%vel(xdim)
               ang_momentum = ang_momentum + part%mass * sqrt(L1**2 + L2**2 + L3**2)
               total_energy = total_energy + part%energy
             end associate
             pset => pset%nxt
         enddo

         cgl => cgl%nxt
      enddo

   end subroutine get_angmom_totener

   subroutine is_part_in_cg(cg, pos, in, phy, out)

      use constants,     only: LO, HI, ndims, xdim, ydim, zdim, I_ONE, LEFT, RIGHT
      use domain,        only: dom
      use grid_cont,     only: grid_container
      use particle_func, only: particle_in_area

      implicit none

      type(grid_container), pointer, intent(in)  :: cg
      real, dimension(ndims),        intent(in)  :: pos
      logical,                       intent(out) :: in, phy, out
      real, dimension(ndims,2)                   :: bnd1, bnd2

      in  = .false.
      phy = .false.
      out = .false.

      !There is probably a better way to write this
      bnd1(:,1) = [cg%coord(LEFT, xdim)%r(cg%lh1(xdim,LO)), cg%coord(LEFT, ydim)%r(cg%lh1(ydim,LO)), cg%coord(LEFT, zdim)%r(cg%lh1(zdim,LO))]
      bnd1(:,2) = [cg%coord(RIGHT,xdim)%r(cg%lh1(xdim,HI)), cg%coord(RIGHT,ydim)%r(cg%lh1(ydim,HI)), cg%coord(RIGHT,zdim)%r(cg%lh1(zdim,HI))]

      bnd2(:,1) = [cg%coord(LEFT, xdim)%r(cg%ijkse(xdim,LO)+I_ONE), cg%coord(LEFT, ydim)%r(cg%ijkse(ydim,LO)+I_ONE), cg%coord(LEFT, zdim)%r(cg%ijkse(zdim,LO)+I_ONE)]
      bnd2(:,2) = [cg%coord(RIGHT,xdim)%r(cg%ijkse(xdim,HI)-I_ONE), cg%coord(RIGHT,ydim)%r(cg%ijkse(ydim,HI)-I_ONE), cg%coord(RIGHT,zdim)%r(cg%ijkse(zdim,HI)-I_ONE)]

      if (particle_in_area(pos, bnd2))    in  = .true.
      if (particle_in_area(pos, cg%fbnd)) phy = .true.
      !Ghost particle
      if (particle_in_area(pos,bnd1))     out = .true.
      if (.not. particle_in_area(pos, dom%edge)) then
         call cg_outside_dom(pos, cg%fbnd, phy)
         if (phy) then
            in  = .true.
            out = .true.
         endif
      endif

    end subroutine is_part_in_cg

   subroutine cg_outside_dom(pos, fbnd, phy)

      use constants, only: LO, HI, ndims, xdim, zdim
      use domain,    only: dom
      use func,      only: operator(.equals.), operator(.notequals.)

      implicit none

      real, dimension(ndims),       intent(in)  :: pos
      real, dimension(ndims,LO:HI), intent(in)  :: fbnd
      logical,                      intent(out) :: phy
      integer                                   :: cdim, k, count, count2, count3
      integer, dimension(ndims)                 :: tmp

      phy = .false.
      count2 = 0
      do cdim = xdim, zdim
         if (pos(cdim) < dom%edge(cdim,LO)) then
            count2 = count2 + 1
            tmp(cdim) = LO
         else if (pos(cdim) > dom%edge(cdim,HI)) then
            count2 = count2 + 1
            tmp(cdim) = HI
         else
            tmp(cdim) = 0
         endif
      enddo

      !CORNER
      if (count2 == 3) then
         count3 = 0
         do cdim = xdim, zdim
            if (fbnd(cdim,tmp(cdim)) .equals. dom%edge(cdim,tmp(cdim))) count3 = count3 + 1
         enddo
         if (count3 == 3) phy = .true.
         return

      !EDGE
      else if (count2 == 2) then
         do cdim = xdim, zdim
            if (tmp(cdim) /= 0) then
               if (fbnd(cdim,tmp(cdim)) .notequals. dom%edge(cdim,tmp(cdim))) return
            else
               if ( (pos(cdim) < fbnd(cdim,LO)) .or. (pos(cdim) > fbnd(cdim,HI)) ) return
            endif
         enddo
         phy = .true.
         return

      !FACE
      else
         do cdim = xdim, zdim
            if ( ((pos(cdim) < dom%edge(cdim,LO)) .and. (fbnd(cdim,LO) .equals. dom%edge(cdim,LO))) .or. (((pos(cdim) > dom%edge(cdim,HI)) .and. (fbnd(cdim,HI) .equals. dom%edge(cdim,HI)))) ) then
               count = 0
               do k = xdim, zdim
                  if (k /= cdim) then
                     if ( (pos(k) > fbnd(k,LO)) .and. (pos(k) < fbnd(k,HI)) ) then
                        count = count + 1
                     endif
                  endif
               enddo
               if (count == 2) then
                  phy = .true.
                  return
               endif
            endif

         enddo
      endif

    end subroutine cg_outside_dom

   subroutine add_part_in_proper_cg(pid, mass, pos, vel, acc, ener)

      use cg_leaves,  only: leaves
      use cg_list,    only: cg_list_element
      use constants,  only: ndims

      implicit none

      integer(kind=4),        intent(in) :: pid
      real, dimension(ndims), intent(in) :: pos, vel
      real, dimension(ndims), intent(in) :: acc
      real,                   intent(in) :: ener
      real,                   intent(in) :: mass
      type(cg_list_element), pointer     :: cgl
      logical                            :: in,phy,out

      cgl => leaves%first
      do while (associated(cgl))
         call is_part_in_cg(cgl%cg, pos, in, phy, out)
         if (phy .or. out) then
            call cgl%cg%pset%add(pid, mass, pos, vel, acc, ener, in, phy, out)
            return
         endif

         cgl => cgl%nxt
      enddo

   end subroutine add_part_in_proper_cg

   ! Sends leaving particles between processors, and creates ghosts
   subroutine part_leave_cg()

      use cg_leaves,     only: leaves
      use cg_level_base, only: base
      use cg_list,       only: cg_list_element
      use constants,     only: ndims, I_ONE, I_TWO, LO, HI
      use domain,        only: dom
      use grid_cont,     only: grid_container
      use mpi,           only: MPI_DOUBLE_PRECISION, MPI_INTEGER
      use mpisetup,      only: proc, comm, mpi_err, FIRST, LAST
      use ppp,           only: ppp_main
      use particle_func, only: particle_in_area
      use particle_types, only: particle

      implicit none

      integer, dimension(FIRST:LAST)     :: nsend, nrecv, counts, countr, disps, dispr
      integer                            :: i, j, ind, b
      integer(kind=4)                    :: pid
      real, dimension(ndims)             :: pos, vel, acc
      real, dimension(:), allocatable    :: part_info, part_info2
      real                               :: mass, ener
      type(cg_list_element), pointer     :: cgl
      type(grid_container),  pointer     :: cg
      type(particle), pointer            :: pset, pset2
      logical                            :: in, phy, out, phy_out
      character(len=*), parameter        :: ts_label = "leave_cg"

      call ppp_main%start(ts_label)

      nsend = 0
      nrecv = 0

      !Count number of particles to be sent
      cgl => leaves%first
      do while (associated(cgl))
         cg => cgl%cg
         do j = FIRST, LAST
            pset => cg%pset%first
            do while (associated(pset))
               if (j == proc) then
                  pset => pset%nxt
                  cycle
               endif
               if (.not. pset%pdata%in) then
                  ! TO CHECK: PARTICLES CHANGING CG OUTSIDE DOMAIN?
                  associate ( gsej => base%level%dot%gse(j) )
                    do b = lbound(gsej%c(:), dim=1), ubound(gsej%c(:), dim=1)
                       if (particle_in_area(pset%pdata%pos, [(gsej%c(b)%se(:,LO) - dom%n_d(:)/2. - I_ONE) * cg%dl(:), (gsej%c(b)%se(:,HI) - dom%n_d(:)/2. + I_TWO)*cg%dl(:)])) then
                          nsend(j) = nsend(j) + 1 ! WON'T WORK in AMR!!!
                       else if (pset%pdata%outside) then
                           call cg_outside_dom(pset%pdata%pos, [(gsej%c(b)%se(:,LO) - dom%n_d(:)/2.) * cg%dl(:), (gsej%c(b)%se(:,HI) - dom%n_d(:)/2. + I_ONE) * cg%dl(:)], phy_out)
                           if (phy_out) nsend(j) = nsend(j) + 1
                        endif
                     enddo
                  end associate
               endif
               pset => pset%nxt
            enddo
         enddo
         cgl => cgl%nxt
      enddo

      !Exchange information about particles numbers to be sent / received
      call MPI_Alltoall(nsend, I_ONE, MPI_INTEGER, nrecv, I_ONE, MPI_INTEGER, comm, mpi_err)

      !Store data of particles to be sent
      allocate(part_info(sum(nsend(:))*npf))
      ind = 1
      cgl => leaves%first
      do while (associated(cgl))
         associate( cg => cgl%cg )
            do j = FIRST, LAST
               pset => cg%pset%first
               do while (associated(pset))
                  if (j == proc) then
                     pset => pset%nxt
                     cycle ! TO DO IN AMR: ADD PARTICLES CHANGING CG INSIDE PROCESSOR
                  endif
                  if (.not. pset%pdata%in) then
                     associate ( gsej => base%level%dot%gse(j) )
                       do b = lbound(gsej%c(:), dim=1), ubound(gsej%c(:), dim=1)
                           if (particle_in_area(pset%pdata%pos, [(gsej%c(b)%se(:,LO) - dom%n_d(:)/2. - I_ONE) * cg%dl(:), (gsej%c(b)%se(:,HI) - dom%n_d(:)/2. + I_TWO)*cg%dl(:)])) then
                              part_info(ind:ind+npf-1) = collect_single_part_fields(ind, pset%pdata)
                           else if (pset%pdata%outside) then
                              call cg_outside_dom(pset%pdata%pos, [(gsej%c(b)%se(:,LO) - dom%n_d(:)/2.) * cg%dl(:), (gsej%c(b)%se(:,HI) - dom%n_d(:)/2. + I_ONE) * cg%dl(:)], phy_out)
                              if (phy_out) part_info(ind:ind+npf-1) = collect_single_part_fields(ind, pset%pdata)
                           endif
                        enddo
                     end associate
                  endif
                  pset => pset%nxt
               enddo
            enddo

            !Remove particles out of cg
            pset => cg%pset%first
            do while (associated(pset))
               if (.not. pset%pdata%out) then
                  pset2 => pset%nxt
                  call cg%pset%remove(pset)
                  pset => pset2
                  cycle
               endif
               pset => pset%nxt
            enddo

         end associate
         cgl => cgl%nxt
      enddo

      !Send / receive particle data
      counts = npf*nsend
      allocate(part_info2(sum(nrecv)*npf))
      countr = npf*nrecv
      disps(FIRST) = 0
      dispr(FIRST) = 0
      do j = FIRST+I_ONE,LAST
         disps(j) = disps(j-1) + counts(j-1)
         dispr(j) = dispr(j-1) + countr(j-1)
      enddo

      call MPI_Alltoallv(part_info, counts, disps, MPI_DOUBLE_PRECISION, part_info2, countr, dispr, MPI_DOUBLE_PRECISION, comm, mpi_err)

      !Add particles in cgs
      cgl => leaves%first
      do while (associated(cgl))
         associate(cg => cgl%cg)
            ind = 1
            do j = FIRST, LAST
               if (nrecv(j) /= 0) then
                  do i = 1, nrecv(j)
                     pos = part_info2(ind+2:ind+4)
                     call is_part_in_cg(cg, pos, in, phy, out) ! TO DO IN AMR USE GRID_ID TO CUT THE SEARCH SHORT
                     if (.not. out) then
                        print *, 'error, particle', part_info2(ind), 'cannot be attributed!' ! NON-AMR ONLY
                     endif
                     if (out) then
                        pid = nint(part_info2(ind), kind=4)
                        mass = part_info2(ind+1)
                        vel  = part_info2(ind+5:ind+7)
                        acc  = part_info2(ind+8:ind+10)
                        ener = part_info2(ind+11)
                        call cg%pset%add(pid, mass, pos, vel, acc, ener, in, phy, out)
                     endif
                     ind = ind + npf
                  enddo
               endif
            enddo
         end associate
         cgl => cgl%nxt
      enddo

      deallocate(part_info2)
      deallocate(part_info)

      call ppp_main%stop(ts_label)

    end subroutine part_leave_cg

   function collect_single_part_fields(ind, p) result(pinfo)

      use particle_types, only: particle_data

      implicit none

      real, dimension(npf)          :: pinfo
      integer,        intent(inout) :: ind
      type(particle_data), intent(in)    :: p

      pinfo(1)    = p%pid
      if (.not. (pinfo(1) >=1)) print *, 'error, id cannot be zero', ind, pinfo(1)
      pinfo(2)    = p%mass
      pinfo(3:5)  = p%pos
      pinfo(6:8)  = p%vel
      pinfo(9:11) = p%acc
      pinfo(12)   = p%energy

      ind = ind + npf

   end function collect_single_part_fields

   integer(kind=4) function count_all_particles() result(pcount)

      use cg_leaves,      only: leaves
      use cg_list,        only: cg_list_element
      use constants,      only: I_ONE
      use particle_types, only: particle

      implicit none

      type(cg_list_element), pointer :: cgl
      type(particle), pointer    :: pset

      pcount = 0
      cgl => leaves%first
      do while (associated(cgl))
         pset => cgl%cg%pset%first
         do while (associated(pset))
            if (pset%pdata%phy) pcount = pcount + I_ONE
            pset => pset%nxt
         enddo
         cgl => cgl%nxt
      enddo

   end function count_all_particles

   subroutine dump_particles_to_textfile

      use cg_leaves,        only: leaves
      use cg_list,          only: cg_list_element
      use constants,        only: half, ndims
      use global,           only: t, dt
      use particle_gravity, only: get_acc_model
      use particle_types,   only: particle

      implicit none

      real                           :: kdt
      real, dimension(ndims)         :: acc2
      integer                        :: i, lun_out
      type(cg_list_element), pointer :: cgl
      type(particle), pointer        :: pset

      kdt = dt ; if (.not.twodtscheme) kdt = half * dt

      open(newunit=lun_out, file='leapfrog_out.log', status='unknown',  position='append')
      cgl => leaves%first
      do while (associated(cgl))

         pset => cgl%cg%pset%first
         do while (associated(pset))
            associate( part => pset%pdata )
               call get_acc_model(part%pos, part%mass, acc2)
               write(lun_out, '(a,I3.3,1X,19(E13.6,1X))') 'particle', i, t, kdt, part%mass, part%pos, part%vel, part%acc, acc2(:), part%energy
             end associate
             pset => pset%nxt
         enddo

         cgl => cgl%nxt
      enddo
      close(lun_out)

   end subroutine dump_particles_to_textfile

end module particle_utils
