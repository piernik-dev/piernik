! $Id$
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
!! \brief Module of boundary conditions for magnetic fields, magnetic vector potential and electromotive forces
!<
module magboundaries
! pulled by MAGNETIC
   implicit none

   private
   public :: bnd_a, bnd_b, bnd_emf, all_mag_boundaries

contains

   subroutine bnd_a(A)

      use cg_leaves,  only: leaves
      use dataio_pub, only: die
      use domain,     only: is_mpi_noncart, is_multicg
      use grid_cont,  only: grid_container
      use mpisetup,   only: have_mpi

      implicit none

      real, dimension(:,:,:,:)      :: A  !< vector potential of magnetic field
      integer(kind=4)               :: i
      type(grid_container), pointer :: cg

      if (is_multicg) call die("[magboundaries:bnd_a] multiple grid pieces per procesor not implemented yet") !nontrivial MPI_Waitall

      if (have_mpi .and. is_mpi_noncart) call die("[magboundaries:bnd_a] is_mpi_noncart is not implemented") !procn, psize
      call die("[magboundaries:bnd_a] Unimplemented")

      cg => leaves%first%cg

      if (.false.) i=int(A(1,1,1,1), kind=4)

   end subroutine bnd_a

   subroutine bnd_b(dir, cg)

      use constants,     only: ndims, xdim, ydim, zdim, LO, HI, I_TWO, I_THREE, &
           &                   BND_MPI, BND_FC, BND_MPI_FC, BND_PER, BND_REF, BND_OUT, BND_OUTD, BND_OUTH, BND_OUTHD, BND_COR, BND_SHE
      use dataio_pub,    only: msg, warn, die
      use grid_cont,     only: grid_container
      use mpisetup,      only: master

      implicit none

      integer(kind=4),               intent(in)    :: dir
      type(grid_container), pointer, intent(inout) :: cg

      integer(kind=4)                              :: side
      logical, save                                :: frun = .true.
      logical, save,   dimension(ndims,LO:HI)      :: bnd_not_provided = .false.
      integer(kind=4), dimension(ndims,LO:HI)      :: l, r

! Non-MPI boundary conditions
      if (frun) then
         bnd_not_provided(:,         :) = (cg%bnd(:,:) == BND_REF)       .or. (cg%bnd(:,         :) == BND_MPI)
         bnd_not_provided(xdim:ydim, :) = bnd_not_provided(xdim:ydim, :) .or. (cg%bnd(xdim:ydim, :) == BND_COR)
         bnd_not_provided(xdim,      :) = bnd_not_provided(xdim,      :) .or. (cg%bnd(xdim,      :) == BND_SHE)
      endif

      if (bnd_not_provided(dir,LO) .and. bnd_not_provided(dir,HI)) return  ! avoid triple case

      l = cg%lhn ; r = l

      do side = LO, HI
         select case (cg%bnd(dir, side))
            case (BND_MPI, BND_REF)
               ! Do nothing
            case (BND_FC, BND_MPI_FC)
               call die("[magboundaries:bnd_b] fine-coarse interfaces not implemented yet")
            case (BND_COR)
               if (dir == zdim) then
                  write(msg,'(2(a,i3))') "[magboundaries:bnd_b]: Boundary condition ",cg%bnd(dir, side)," not implemented in ",dir
                  if (master) call warn(msg)
               endif
            case (BND_SHE)
               if (dir /= xdim) then
                  write(msg,'(2(a,i3))') "[magboundaries:bnd_b]: Boundary condition ",cg%bnd(dir, side)," not implemented in ",dir
                  if (master) call warn(msg)
               endif
            case (BND_PER)
            case (BND_OUT, BND_OUTD, BND_OUTH, BND_OUTHD)
                  l(dir,:) = cg%lhn(dir,side)
                  r(dir,:) = cg%lhn(dir,side) + I_THREE - I_TWO*side
                  cg%b(:,l(xdim,LO):l(xdim,HI),l(ydim,LO):l(ydim,HI),l(zdim,LO):l(zdim,HI)) = cg%b(:,r(xdim,LO):r(xdim,HI),r(ydim,LO):r(ydim,HI),r(zdim,LO):r(zdim,HI))
            case default
               write(msg,'(2(a,i3))') "[magboundaries:bnd_b]: Boundary condition ",cg%bnd(dir, side)," not implemented in ",dir
               if (master) call warn(msg)
         end select
      enddo

   end subroutine bnd_b

!=====================================================================================================

   subroutine bnd_emf(ivar, emfdir, dir, cg)

      use constants,  only: ndims, xdim, ydim, zdim, LO, HI, I_ZERO, I_ONE, &
                            BND_MPI, BND_FC, BND_MPI_FC, BND_PER, BND_REF, BND_OUT, BND_OUTD, BND_OUTH, BND_OUTHD, BND_COR, BND_SHE
      use dataio_pub, only: msg, warn, die
      use grid_cont,  only: grid_container
      use mpisetup,   only: master

      implicit none

      integer,                       intent(in)    :: ivar      !< index in cg%q array
      integer(kind=4),               intent(in)    :: emfdir
      integer(kind=4),               intent(in)    :: dir
      type(grid_container), pointer, intent(inout) :: cg

#ifndef ZERO_BND_EMF
      real, dimension(:,:,:), allocatable          :: dvar
#endif /* !ZERO_BND_EMF */
      real, dimension(:,:,:), pointer              :: p3, p3a
      real                                         :: bndsign
      logical,                         save        :: frun = .true.
      logical, dimension(ndims,LO:HI), save        :: bnd_not_provided = .false.
      logical                                      :: zndiff
      integer(kind=4), dimension(ndims,LO:HI)      :: l, r
      integer(kind=4), dimension(LO:HI)            :: sbase, edge, nbcells, sidebase
      integer(kind=4)                              :: ssign, side, ib

      if (frun) then
         bnd_not_provided(:, :)         = (cg%bnd(:,:) == BND_PER)       .or. (cg%bnd(:,         :) == BND_MPI)
         bnd_not_provided(xdim:ydim, :) = bnd_not_provided(xdim:ydim, :) .or. (cg%bnd(xdim:ydim, :) == BND_COR)
         bnd_not_provided(xdim, :)      = bnd_not_provided(xdim,      :) .or. (cg%bnd(xdim,      :) == BND_SHE)
         frun = .false.
      endif

      if (bnd_not_provided(dir,LO) .and. bnd_not_provided(dir,HI)) return  ! avoid triple case

      bndsign = huge(1.0); edge=huge(I_ONE); nbcells=huge(I_ONE); zndiff=.false.; sidebase=huge(I_ONE)
      ! the code below should not use these values and the compiler should not complain on possible use of uninitialized variables.

      if (any(cg%bnd(dir,LO) == [BND_REF, BND_OUT, BND_OUTD, BND_OUTH, BND_OUTHD]) .or. &
       &  any(cg%bnd(dir,HI) == [BND_REF, BND_OUT, BND_OUTD, BND_OUTH, BND_OUTHD])) then
         call compute_bnd_indxs(emfdir, cg%n_b(dir),edge,nbcells,sidebase,bndsign,zndiff)
      endif

      l = cg%lhn ; r = l

      do side = LO, HI
         select case (cg%bnd(dir, side))
            case (BND_MPI, BND_PER)
               ! Do nothing
            case (BND_FC, BND_MPI_FC)
               call die("[magboundaries:bnd_emf] fine-coarse interfaces not implemented yet")
            case (BND_COR)
               if (dir == zdim) then
                  write(msg,'(a,i3,a,i1,a,i3)') "[magboundaries:bnd_emf]: Boundary condition ",cg%bnd(dir, side)," not implemented for ",emfdir, " in ", dir
                  if (master) call warn(msg)
               endif
            case (BND_SHE)
               if (dir /= xdim) then
                  write(msg,'(a,i3,a,i1,a,i3)') "[magboundaries:bnd_emf]: Boundary condition ",cg%bnd(dir, side)," not implemented for ",emfdir, " in ", dir
                  if (master) call warn(msg)
               endif
            case (BND_REF)
               sbase(:)  = [nbcells(LO)+I_ONE, edge(HI)] ; ssign = int(2*side-3, kind=4)
               if (zndiff) then
                  l(dir,:) = edge(side) ; p3 => cg%q(ivar)%span(l) ; p3 = 0.0
               endif
               do ib=1,nbcells(side)
                  l(dir,:) = sbase(side)+ssign*ib    ; p3  => cg%q(ivar)%span(l)
                  r(dir,:) = sidebase(side)-ssign*ib ; p3a => cg%q(ivar)%span(r)
                  p3 = bndsign * p3a
               enddo
            case (BND_OUT, BND_OUTD, BND_OUTH, BND_OUTHD)
               sbase(:)  = [I_ZERO, edge(HI)]
#ifdef ZERO_BND_EMF
               l(dir,LO) = sbase(side)+I_ONE ; l(dir,HI) = sbase(side)+nbcells(side)
               p3 => cg%q(ivar)%span(l) ; p3 = 0.0
#else /* !ZERO_BND_EMF */
               l(dir,:) = 1 ; allocate(dvar(l(xdim,HI) ,l(ydim,HI), l(zdim,HI)))
               edge(side) = edge(side) + HI - side ; nbcells(side) = nbcells(side) + HI - side
!               l(dir,:) = sidebase(side)+HI-side ; r(dir,:) = l(dir,:)-1 original
               l(dir,:) = edge(side)+HI-side ; r(dir,:) = l(dir,:) - I_ONE
               dvar(:,:,:) = cg%q(ivar)%span(l) - cg%q(ivar)%span(r)
               r(dir,:) = edge(side) ; p3a => cg%q(ivar)%span(r)
               do ib=1,nbcells(side)
                  l(dir,:) = sbase(side) + ib ; p3 => cg%q(ivar)%span(l)
                  p3 = p3a + real(ib+sbase(side)-edge(side))*dvar
               enddo
               deallocate(dvar)
#endif /* ZERO_BND_EMF */
            case default
               write(msg,'(a,i3,a,i1,a,i3)') "[magboundaries:bnd_emf]: Boundary condition ",cg%bnd(dir, side)," not implemented for ",emfdir, " in ", dir
               if (master) call warn(msg)
         end select
      enddo

   end subroutine bnd_emf
!>
!! \brief Routine delivers common boundary cells indexes in cases of reflection or outflow boundary types
!!
!! \param bndcase
!!    1 - v component compatible with direction;
!!    2 - b component compatible with direction or emf component incompatible with direction;
!!    3 - other cases; BEWARE: magic integers
!<
   subroutine compute_bnd_indxs(bndcase, ndirb, edge, nbcells, rrbase, bndsign, zndiff)

      use constants, only: LO, HI, I_ONE, I_TWO
      use domain,    only: dom

      implicit none

      integer(kind=4),                   intent(in)  :: bndcase
      integer(kind=4),                   intent(in)  :: ndirb   !< cg%{nxb,nyb,nzb} depanding on the current direction
      integer(kind=4), dimension(LO:HI), intent(out) :: edge    !< index of the left and right edge of physical domain for emf
      integer(kind=4), dimension(LO:HI), intent(out) :: nbcells !< number of cells in a loop at left and right boundaries
      integer(kind=4), dimension(LO:HI), intent(out) :: rrbase  !< COMMENT ME
      real,                              intent(out) :: bndsign !< 1. or -1. to change the sign or not
      logical,                           intent(out) :: zndiff  !< COMMENT ME

      select case (bndcase)
         case (1)
            edge(LO)    = dom%nb
            nbcells(LO) = dom%nb - I_ONE
            bndsign     = -1.
         case (2)
            edge(LO)    = dom%nb + I_ONE
            nbcells(LO) = dom%nb
            bndsign     = -1.
         case (3)
            edge(LO)    = dom%nb
            nbcells(LO) = dom%nb
            bndsign     = 1.
      end select  ! (bndcase)

      zndiff      = (edge(LO) - nbcells(LO) == 1)
!     nbcells(HI) = dom%nb - edge(LO) + nbcells(LO)
      nbcells(HI) = I_TWO * dom%nb - edge(LO)
      edge(HI)    = ndirb + edge(LO)
      rrbase(LO)  = edge(LO)
      rrbase(HI)  = ndirb + nbcells(LO) + I_ONE  ! = edge(HI) + 1 - edge(LO) + nbcells(LO)

   end subroutine compute_bnd_indxs

   subroutine all_mag_boundaries

      use cg_leaves,        only: leaves
      use cg_list,          only: cg_list_element
      use cg_list_global,   only: all_cg
      use constants,        only: xdim, zdim
      use domain,           only: dom
      use named_array_list, only: wna

      implicit none

      type(cg_list_element), pointer :: cgl
      integer(kind=4) :: dir

      do dir = xdim, zdim
         if (dom%has_dir(dir)) call all_cg%internal_boundaries_4d(wna%bi, dim=dir) ! should be more selective (modified leaves?)
      enddo

      cgl => leaves%first
      do while (associated(cgl))
         do dir = xdim, zdim
            if (dom%has_dir(dir)) call bnd_b(dir, cgl%cg)
         enddo
         cgl => cgl%nxt
      enddo

   end subroutine all_mag_boundaries

end module magboundaries
