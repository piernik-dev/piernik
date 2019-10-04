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
! pulled by ANY

   implicit none

   private
   public :: max_pvel_1d, add_part_in_proper_cg, count_all_particles, print_all_particles, is_part_in_cg, part_leave_cg
#ifdef NBODY
   public :: max_pacc_3d, particle_diagnostics, twodtscheme, dump_diagnose, tot_energy, d_energy, tot_angmom, d_angmom

   real    :: tot_angmom           !< angular momentum of set of the particles
   real    :: tot_energy           !< total energy of set of the particles
   real    :: d_energy             !< error of energy of set of the particles in succeeding timesteps
   real    :: d_angmom             !< error of angular momentum in succeeding timensteps
   logical :: twodtscheme
   logical :: dump_diagnose        !< dump diagnose for each particle to a seperate log file
#endif /* NBODY */

contains

   subroutine max_pvel_1d(cg, max_v)

      use constants, only: ndims, xdim, zdim
      use grid_cont, only: grid_container

      implicit none

      type(grid_container), pointer, intent(in)  :: cg
      real, dimension(ndims),        intent(out) :: max_v
      integer(kind=4)                            :: cdim

      do cdim = xdim, zdim
         max_v(cdim)  = maxval(abs(cg%pset%p(:)%vel(cdim)))
      enddo

   end subroutine max_pvel_1d

#ifdef NBODY
   subroutine max_pacc_3d(cg, max_pacc)

      use constants, only: big, CENTER, half, LO, xdim, zdim, zero
      use grid_cont, only: grid_container
      use mpisetup,  only: proc
      use types,     only: value

      implicit none

      type(grid_container), pointer, intent(in)  :: cg
      type(value),                   intent(out) :: max_pacc
      real                                       :: acc2, max_acc
      integer                                    :: i, n_part
      integer(kind=4) :: cdim

      max_pacc%assoc = big

      max_acc  = zero
      n_part = size(cg%pset%p, dim=1)

      do i = 1, n_part
         acc2 = sum(cg%pset%p(i)%acc(:)**2)
         if (acc2 > max_acc) then
            max_acc = acc2
            max_pacc%coords(:) = cg%pset%p(i)%pos(:)
            !max_pacc%proc = i !> \todo it might be an information about extremum particle, but the scheme of log file is to print the process number
         endif
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

      implicit none

      real, intent(out)              :: ang_momentum
      real, intent(out)              :: total_energy !< total energy of set of particles
      integer                        :: i
      real                           :: L1, L2, L3
      type(cg_list_element), pointer :: cgl

      ang_momentum = zero
      total_energy = zero

      cgl => leaves%first
      do while (associated(cgl))

         do i = 1, size(cgl%cg%pset%p, dim=1)
            associate( part => cgl%cg%pset%p(i) )
               L1 = part%pos(ydim) * part%vel(zdim) - part%pos(zdim) * part%vel(ydim)
               L2 = part%pos(zdim) * part%vel(xdim) - part%pos(xdim) * part%vel(zdim)
               L3 = part%pos(xdim) * part%vel(ydim) - part%pos(ydim) * part%vel(xdim)
               ang_momentum = ang_momentum + part%mass * sqrt(L1**2 + L2**2 + L3**2)
               total_energy = total_energy + part%energy
            end associate
         enddo

         cgl => cgl%nxt
      enddo

   end subroutine get_angmom_totener
#endif /* NBODY */

   subroutine is_part_in_cg(cg, pos, in, phy, out)

      use grid_cont,      only: grid_container
      use constants,     only: LO, HI, ndims, xdim, ydim, zdim, I_ONE, LEFT, RIGHT
      use domain,        only: dom
      use particle_func, only: particle_in_area

      implicit none

      real, dimension(ndims), intent(in)  :: pos
      logical,                intent(out) :: in, phy, out
      integer                             :: cdim, k, count, count2, count3
      integer, dimension(ndims)           :: tmp
      real, dimension(ndims,2)            :: bnd1, bnd2
      type(grid_container), pointer, intent(in)      :: cg

      in  = .false.
      phy = .false.
      out = .false.

      !There is probably a better way to write this
      bnd1(:,1) = [cg%coord(LEFT, xdim)%r(cg%lh1(xdim,LO)), cg%coord(LEFT, ydim)%r(cg%lh1(ydim,LO)), cg%coord(LEFT, zdim)%r(cg%lh1(zdim,LO))]
      bnd1(:,2) = [cg%coord(RIGHT,xdim)%r(cg%lh1(xdim,HI)), cg%coord(RIGHT,ydim)%r(cg%lh1(ydim,HI)), cg%coord(RIGHT,zdim)%r(cg%lh1(zdim,HI))]

      bnd2(:,1) = [cg%coord(LEFT, xdim)%r(cg%ijkse(xdim,LO)+I_ONE), cg%coord(LEFT, ydim)%r(cg%ijkse(ydim,LO)+I_ONE), cg%coord(LEFT, zdim)%r(cg%ijkse(zdim,LO)+I_ONE)]
      bnd2(:,2) = [cg%coord(RIGHT,xdim)%r(cg%ijkse(xdim,HI)-I_ONE), cg%coord(RIGHT,ydim)%r(cg%ijkse(ydim,HI)-I_ONE), cg%coord(RIGHT,zdim)%r(cg%ijkse(zdim,HI)-I_ONE)]


      if (particle_in_area(pos, bnd2))        in  = .true.
      if (particle_in_area(pos, cg%fbnd)) phy = .true.
      !Ghost particle
      if (particle_in_area(pos,bnd1))         out=.true.

      if (.not. particle_in_area(pos, dom%edge)) then
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
               if (cg%fbnd(cdim,tmp(cdim)) == dom%edge(cdim,tmp(cdim))) count3 = count3 + 1
            enddo
            if (count3 == 3) phy = .true.
            return

            !EDGE
         else if (count2 == 2) then
            do cdim = xdim, zdim
               if (tmp(cdim) /= 0) then
                  if (cg%fbnd(cdim,tmp(cdim)) /= dom%edge(cdim,tmp(cdim))) return
               else
                  if ( (pos(cdim) < cg%fbnd(cdim,LO)) .or. (pos(cdim) > cg%fbnd(cdim,HI)) ) return
               endif
            enddo
            phy = .true.
            return

            !FACE
         else
            do cdim = xdim, zdim
               if ( ((pos(cdim) < dom%edge(cdim,LO)) .and. (cg%fbnd(cdim,LO) == dom%edge(cdim,LO))) .or. (((pos(cdim) > dom%edge(cdim,HI)) .and. (cg%fbnd(cdim,HI) == dom%edge(cdim,HI)))) ) then
                  count = 0
                  do k = xdim, zdim
                     if (k /= cdim) then
                        if ( (pos(k) > cg%fbnd(k,LO)) .and. (pos(k) < cg%fbnd(k,HI)) ) then
                           count = count + 1
                        endif
                     endif
                  enddo
                  if (count == 2) then
                     phy =.true.
                     return
                  endif
               endif

            enddo
         endif
      endif

   end subroutine is_part_in_cg

#ifdef NBODY
   subroutine add_part_in_proper_cg(pid, mass, pos, vel, acc, ener)
#else /* !NBODY */
   subroutine add_part_in_proper_cg(pid, mass, pos, vel)
#endif /* !NBODY */


      use cg_leaves,     only: leaves
      use cg_list,       only: cg_list_element
      use constants,     only: ndims
      use dataio_pub,    only: msg, printinfo
      implicit none

      integer,                intent(in) :: pid
      real, dimension(ndims), intent(in) :: pos, vel
#ifdef NBODY
      real, dimension(ndims), intent(in) :: acc
      real,                   intent(in) :: ener
#endif /* NBODY */
      real,                   intent(in) :: mass
      type(cg_list_element), pointer     :: cgl
      logical                            :: in,phy,out


      cgl => leaves%first
      do while (associated(cgl))
         call is_part_in_cg(cgl%cg, pos, in, phy, out)
         if (phy .or. out) then
            call printinfo(msg)
#ifdef NBODY
            call cgl%cg%pset%add(pid, mass, pos, vel, acc, ener, in, phy, out)
#else /* !NBODY */
            call cgl%cg%pset%add(pid, mass, pos, vel)
#endif /* !NBODY */
            return
         endif

         cgl => cgl%nxt
      enddo


   end subroutine add_part_in_proper_cg


   ! Sends leaving particles between processors, and creates ghosts
   subroutine part_leave_cg()

     use constants,      only: ndims, I_ONE, LO, HI
     use grid_cont,      only: grid_container
     use mpisetup,       only: proc, comm, mpi_err, FIRST, LAST
     use mpi,            only: MPI_REAL, MPI_INTEGER
     use cg_level_base,  only: base
     use domain,         only: dom
     use particle_func,  only: particle_in_area
     use cg_leaves,      only: leaves
     use cg_list,        only: cg_list_element

     implicit none

     integer                                  :: i,j, pid, npart,k, b
     real,                 dimension(ndims)   :: pos, vel, acc
     real                                     :: mass,ener
     integer, dimension(:), allocatable       :: pdel
     real(kind=4), dimension(:), allocatable  :: part_info, part_info2
     integer,   dimension(FIRST:LAST)         :: nsend, nrecv, counts, countr, disps, dispr
     integer                                  :: ind
     type(cg_list_element), pointer           :: cgl
     type(grid_container),  pointer           :: cg
     logical                                  :: already, in, phy, out

     nsend=0
     nrecv=0

     !Count number of particles to be sent
     cgl => leaves%first
     do while (associated(cgl))
        cg => cgl%cg
          npart=size(cg%pset%p,dim=1)
          do i=1,npart
             if (.not. cg%pset%p(i)%in) then
                do j=FIRST, LAST
                   if (j==proc) cycle
                   ! TO DO: ADD CONDITION FOR PARTICLES CHANGING CG OUTSIDE DOMAIN
                   do b = lbound(base%level%dot%gse(j)%c(:), dim=1), ubound(base%level%dot%gse(j)%c(:), dim=1)
                      if ( particle_in_area( cg%pset%p(i)%pos, [ (base%level%dot%gse(j)%c(b)%se(:,LO) - dom%n_d(:)/2. - I_ONE) * cg%dl(:), (base%level%dot%gse(j)%c(b)%se(:,HI) - dom%n_d(:)/2.+I_ONE)*cg%dl(:)])) nsend(j) = nsend(j) + 1 ! WON'T WORK in AMR!!!
                   enddo
                enddo
             endif
          enddo
          cgl => cgl%nxt
       enddo

       !Exchange information about particles numbers to be sent / received
       call MPI_Alltoall(nsend, I_ONE, MPI_INTEGER, nrecv, I_ONE, MPI_INTEGER, comm, mpi_err)

       !Store data of particles to be sent
       allocate(part_info(sum(nsend(:))*12))
       ind=1
       cgl => leaves%first
       do while (associated(cgl))
          associate(cg => cgl%cg)
            npart=size(cg%pset%p,dim=1)
            allocate(pdel(npart))
            pdel=0
            do j=FIRST, LAST
               do i=1,npart
                  if (.not. cg%pset%p(i)%out) pdel(i)=1
                  if (j==proc) cycle ! TO DO IN AMR: ADD PARTICLES CHANGING CG INSIDE PROCESSOR
                  if (.not. cg%pset%p(i)%in) then
                     do b = lbound(base%level%dot%gse(j)%c(:), dim=1), ubound(base%level%dot%gse(j)%c(:), dim=1)
                        if ( particle_in_area( cg%pset%p(i)%pos, [(base%level%dot%gse(j)%c(b)%se(:,LO) - dom%n_d(:)/2. - I_ONE) * cg%dl(:), (base%level%dot%gse(j)%c(b)%se(:,HI) - dom%n_d(:)/2. + I_ONE)*cg%dl(:)])) then
                           part_info(ind)=cg%pset%p(i)%pid
                           part_info(ind+1)=cg%pset%p(i)%mass
                           part_info(ind+2:ind+4)=cg%pset%p(i)%pos
                           part_info(ind+5:ind+7)=cg%pset%p(i)%vel
                           part_info(ind+8:ind+10)=cg%pset%p(i)%acc
                           part_info(ind+11)=cg%pset%p(i)%energy
                           ind=ind+12
                        endif
                     enddo
                  endif
               enddo
            enddo

            !Remove particles out of cg
            k = 1
            do i = 1, npart
               if (pdel(i)==1) then
                  call cg%pset%remove(k)
               else
                  k = k + 1
               endif
            enddo
            deallocate(pdel)

          end associate
          cgl => cgl%nxt
       enddo

       !Send / receive particle data
       counts=12*nsend
       allocate(part_info2(sum(nrecv)*12))
       countr=12*nrecv
       disps(FIRST) = 0
       dispr(FIRST) = 0
       do j=FIRST+I_ONE,LAST
          disps(j) = disps(j-1) + counts(j-1)
          dispr(j) = dispr(j-1) + countr(j-1)
       enddo

       call MPI_Alltoallv(part_info, counts, disps, MPI_REAL, part_info2, countr, dispr, MPI_REAL, comm, mpi_err)

       !Add particles in cgs
       cgl => leaves%first
       do while (associated(cgl))
          associate(cg => cgl%cg)
            ind=1
            do j=FIRST, LAST
               if (nrecv(j) /= 0) then
                  do i=1, nrecv(j)
                     pos=part_info2(ind+2:ind+4)
                     call is_part_in_cg(cg, pos, in, phy, out) ! TO DO IN AMR USE GRID_ID TO CUT THE SEARCH SHORT
                     if (.not. out) print *, 'error, particle', part_info2(ind), 'cannot be attributed!' ! NON-AMR ONLY
                     if (out) then
                        pid=nint(part_info2(ind))
                        already = .false.
                        do k = 1, size(cg%pset%p, dim=1)
                           if (pid == cg%pset%p(k)%pid) then
                              already = .true.
                              exit
                           endif
                        enddo
                        if (.not. already) then
                           mass=part_info2(ind+1)
                           vel=part_info2(ind+5:ind+7)
                           acc=part_info2(ind+8:ind+10)
                           ener=part_info2(ind+11)
                           call cg%pset%add(pid, mass, pos, vel, acc, ener, in, phy, out)
                        endif
                     endif
                     ind=ind+12
                  enddo
               endif
            enddo
          end associate
          cgl => cgl%nxt
       enddo

       deallocate(part_info2)
       deallocate(part_info)


   end subroutine part_leave_cg


   integer(kind=4) function count_all_particles() result(pcount)

      use cg_leaves, only: leaves
      use cg_list,   only: cg_list_element

      implicit none

      integer(kind=4)                :: i
      type(cg_list_element), pointer :: cgl

      pcount = 0
      cgl => leaves%first
      do while (associated(cgl))
         do i = 1, size(cgl%cg%pset%p, dim=1, kind=4)
            if (cgl%cg%pset%p(i)%phy) pcount = pcount + 1
         enddo
         cgl => cgl%nxt
      enddo

   end function count_all_particles

   subroutine print_all_particles

      use cg_leaves, only: leaves
      use cg_list,   only: cg_list_element

      implicit none

      type(cg_list_element), pointer :: cgl

      cgl => leaves%first
      do while (associated(cgl))
         call cgl%cg%pset%print
         cgl => cgl%nxt
      enddo

   end subroutine print_all_particles

#ifdef NBODY
   subroutine dump_particles_to_textfile

      use cg_leaves,        only: leaves
      use cg_list,          only: cg_list_element
      use constants,        only: half, ndims
      use global,           only: t, dt
      use particle_gravity, only: get_acc_model

      implicit none

      real                           :: kdt
      real, dimension(ndims)         :: acc2
      integer                        :: i, lun_out
      type(cg_list_element), pointer :: cgl

      kdt = dt ; if (.not.twodtscheme) kdt = half * dt

      open(newunit=lun_out, file='leapfrog_out.log', status='unknown',  position='append')
      cgl => leaves%first
      do while (associated(cgl))

         do i = 1, size(cgl%cg%pset%p, dim=1)
            associate( part => cgl%cg%pset%p(i) )
               call get_acc_model(part%pos, part%mass, acc2)
               write(lun_out, '(a,I3.3,1X,19(E13.6,1X))') 'particle', i, t, kdt, part%mass, part%pos, part%vel, part%acc, acc2(:), part%energy
            end associate
         enddo

         cgl => cgl%nxt
      enddo
      close(lun_out)

   end subroutine dump_particles_to_textfile
#endif /* NBODY */

end module particle_utils
