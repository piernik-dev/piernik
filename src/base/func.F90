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
!! \brief (KK) Module that contains unclassified functions
!!
!! This module should be empty. Every function or subroutine placed here belong
!! elsewhere. We are yet unsure where to put them.
!! \todo Move all structures elsewhere
!! \warning Procedures \a dipole and \a rn_angles were moved to sn_sources.F90
!<
module func
   implicit none
   private
   public :: pshift, mshift, fix_string, ekin, emag, L2norm, get_extremum, sanitize_smallx_checks
   integer, parameter :: one = 1
contains

!>
!! \brief Function pshift makes one-cell, forward circular shift of 3D array in any direction
!! \param tab input array
!! \param d direction of the shift, where 1,2,3 corresponds to \a x,\a y,\a z respectively
!! \return real, dimension(size(tab,1),size(tab,2),size(tab,3))
!!
!! The function was written in order to significantly improve
!! the performance at the cost of the flexibility of original \p CSHIFT.
!<
   function pshift(tab, d)

      use dataio_pub,    only: warn

      implicit none

      real, dimension(:,:,:), intent(inout) :: tab
      integer(kind=4), intent(in) :: d

      integer :: ll
      real, dimension(size(tab,1),size(tab,2),size(tab,3)) :: pshift

      ll = size(tab,d)

      if (ll==1) then
         pshift = tab
         return
      endif

      if (d==1) then
         pshift(1:ll-1,:,:) = tab(2:ll,:,:); pshift(ll,:,:) = tab(1,:,:)
      else if (d==2) then
         pshift(:,1:ll-1,:) = tab(:,2:ll,:); pshift(:,ll,:) = tab(:,1,:)
      else if (d==3) then
         pshift(:,:,1:ll-1) = tab(:,:,2:ll); pshift(:,:,ll) = tab(:,:,1)
      else
         call warn('[func:pshift]: Dim ill defined in pshift!')
      endif

      return
   end function pshift

!>
!! \brief Function mshift makes one-cell, backward circular shift of 3D array in any direction
!! \param tab input array
!! \param d direction of the shift, where 1,2,3 corresponds to \a x,\a y,\a z respectively
!! \return real, dimension(size(tab,1),size(tab,2),size(tab,3))
!!
!! The function was written in order to significantly improve
!! the performance at the cost of the flexibility of original \p CSHIFT.
!<
   function mshift(tab,d)

      use dataio_pub,    only: warn

      implicit none

      real, dimension(:,:,:), intent(inout) :: tab
      integer(kind=4) :: d

      integer :: ll
      real, dimension(size(tab,1) , size(tab,2) , size(tab,3)) :: mshift

      ll = size(tab,d)

      if (ll==1) then
         mshift = tab
         return
      endif

      if (d==1) then
         mshift(2:ll,:,:) = tab(1:ll-1,:,:); mshift(1,:,:) = tab(ll,:,:)
      else if (d==2) then
         mshift(:,2:ll,:) = tab(:,1:ll-1,:); mshift(:,1,:) = tab(:,ll,:)
      else if (d==3) then
         mshift(:,:,2:ll) = tab(:,:,1:ll-1); mshift(:,:,1) = tab(:,:,ll)
      else
         call warn('[func:mshift]: Dim ill defined in mshift!')
      endif

      return
   end function mshift

!> \todo move to a better place (possibly dataio_pub since stop types module stops using die routine)
   subroutine get_extremum(tab, minmax, prop, cg)

      use constants,  only: MINL, MAXL, ndims, xdim, ydim, zdim
      use dataio_pub, only: msg, warn, die
      use domain,     only: has_dir
      use grid,       only: all_cg
      use grid_cont,  only: grid_container
      use mpi,        only: MPI_DOUBLE_PRECISION, MPI_INTEGER, MPI_STATUS_IGNORE
      use mpisetup,   only: mpifind, comm, ierr, master, proc, FIRST
      use types,      only: value

      implicit none

      real, dimension(:,:,:), intent(in), pointer  :: tab
      integer(kind=4),        intent(in)  :: minmax
      type(value),            intent(out) :: prop
      type(grid_container), pointer, intent(in) :: cg

      integer, parameter :: tag1 = 11
      integer, parameter :: tag2 = 12

      if (all_cg%cnt > 1) call die("[func:get_extremum] multiple grid pieces per procesor not implemented yet") !nontrivial

      select case (minmax)
         case (MINL)
            prop%val = minval(tab)
            prop%loc = minloc(tab) + [cg%nb, cg%nb, cg%nb]
         case (MAXL)
            prop%val = maxval(tab)
            prop%loc = maxloc(tab) + [cg%nb, cg%nb, cg%nb]
         case default
            write(msg,*) "[func:get_extremum]: I don't know what to do with minmax = ", minmax
            call warn(msg)
      end select

      call mpifind(prop, minmax)

      if (proc == prop%proc) then
         where (.not. has_dir(:)) prop%coords(:) = 0.
         if (has_dir(xdim)) prop%coords(xdim) = cg%x(prop%loc(xdim))
         if (has_dir(ydim)) prop%coords(ydim) = cg%y(prop%loc(ydim))
         if (has_dir(zdim)) prop%coords(zdim) = cg%z(prop%loc(zdim))
      endif

      if (prop%proc /= 0) then
         if (proc == prop%proc) then ! slave
            call MPI_Send (prop%loc,    ndims, MPI_INTEGER,          FIRST, tag1, comm, ierr)
            call MPI_Send (prop%coords, ndims, MPI_DOUBLE_PRECISION, FIRST, tag2, comm, ierr)
         endif
         if (master) then
            call MPI_Recv (prop%loc,    ndims, MPI_INTEGER,          prop%proc, tag1, comm, MPI_STATUS_IGNORE, ierr)
            call MPI_Recv (prop%coords, ndims, MPI_DOUBLE_PRECISION, prop%proc, tag2, comm, MPI_STATUS_IGNORE, ierr)
         endif
      endif

   end subroutine get_extremum

!>
!! \brief Find sane values for smalld and smallp
!!
!! \deprecated This routine was moved here from initfluids just to resolve circular dependencies in some setups.
!! \todo Find a better file for it
!<

   subroutine sanitize_smallx_checks

      use constants,  only: big_float, DST, I_ONE
      use dataio_pub, only: warn, msg, die
      use fluidindex, only: flind, ibx, iby, ibz
      use fluidtypes, only: component_fluid
      use global,     only: smalld, smallp
      use gc_list,    only: cg_list_element
      use grid,       only: all_cg
      use grid_cont,  only: grid_container
      use mpi,        only: MPI_IN_PLACE, MPI_DOUBLE_PRECISION, MPI_MIN
      use mpisetup,   only: master, comm, ierr

      implicit none

      type(cg_list_element), pointer :: cgl
      type(grid_container), pointer :: cg
      type(component_fluid), pointer   :: fl
      integer                          :: i
      real, pointer, dimension(:,:,:)  :: dn, mx, my, mz, en, bx, by, bz
      real, parameter                  :: safety_factor = 1.e-4
      real, parameter                  :: max_dens_span = 5.0
      real                             :: maxdens, span

      maxdens = 0.0

      if (all_cg%cnt > 1) call die("[func:sanitize_smallx_checks] multiple grid pieces per procesor not fully implemented yet") !nontrivial maxval, minval

      cgl => all_cg%first
      do while (associated(cgl))
         cg => cgl%cg

         bx => cg%b%arr(ibx,:,:,:)
         by => cg%b%arr(iby,:,:,:)
         bz => cg%b%arr(ibz,:,:,:)

         if (smalld >= big_float) then
            do i = lbound(flind%all_fluids,1), ubound(flind%all_fluids,1)
               fl => flind%all_fluids(i)
               dn => cg%u%arr(fl%idn,:,:,:)

               maxdens = max(maxval(dn), maxdens)
               smalld  = min(minval(dn), smalld)
            enddo
            span   = log10(maxdens) - log10(smalld)
            smalld = smalld * safety_factor
            call MPI_Allreduce(MPI_IN_PLACE, smalld, I_ONE, MPI_DOUBLE_PRECISION, MPI_MIN, comm, ierr)
            if (master) then
               write(msg,'(A,ES11.4)') "[func:sanitize_smallx_checks] adjusted smalld to ", smalld
               call warn(msg)
               if (span > max_dens_span) then
                  write(msg,'(A,I3,A)') "[func:sanitize_smallx_checks] density spans over ", int(span), " orders of magnitude!"
                  call warn(msg)
               endif
            endif
         endif

         if (smallp >= big_float) then
            do i = lbound(flind%all_fluids,1), ubound(flind%all_fluids,1)
               fl => flind%all_fluids(i)
               if (fl%tag == DST) cycle
               dn => cg%u%arr(fl%idn,:,:,:)
               mx => cg%u%arr(fl%imx,:,:,:)
               my => cg%u%arr(fl%imy,:,:,:)
               mz => cg%u%arr(fl%imz,:,:,:)
               if (fl%has_energy) then
                  en => cg%u%arr(fl%ien,:,:,:)
                  if (fl%is_magnetized) then
                     smallp = min( minval( en - ekin(mx,my,mz,dn) - emag(bx,by,bz))/fl%gam_1, smallp)
                  else
                     smallp = min( minval( en - ekin(mx,my,mz,dn))/fl%gam_1, smallp )
                  endif
               else
                  smallp = min( minval( cg%cs_iso2*dn ), smallp )
               endif
            enddo
            smallp = smallp * safety_factor
            call MPI_Allreduce(MPI_IN_PLACE, smallp, I_ONE, MPI_DOUBLE_PRECISION, MPI_MIN, comm, ierr)
            if (smallp < 0.) then
               write(msg,'(A,ES11.4,A)') "[func:sanitize_smallx_checks] Negative smallp detected! smallp=",smallp," may indicate nonphysical initial conditions."
               if (master) call warn(msg)
               smallp = tiny(1.)
            endif
            if (master) then
               write(msg,'(A,ES11.4)') "[func:sanitize_smallx_checks] adjusted smallp to ", smallp
               call warn(msg)
            endif
         endif

         cgl => cgl%nxt
      enddo

      if (associated(dn)) nullify(dn)
      if (associated(mx)) nullify(mx)
      if (associated(my)) nullify(my)
      if (associated(mz)) nullify(mz)
      if (associated(en)) nullify(en)
      if (associated(bx)) nullify(bx)
      if (associated(by)) nullify(by)
      if (associated(bz)) nullify(bz)
      if (associated(fl)) nullify(fl)

   end subroutine sanitize_smallx_checks

   function fix_string(str) result (outstr)
      implicit none
      character(len=*), intent(in)  :: str
      character(len=len(str)) :: outstr

      integer            :: i
      character(len=one) :: c

      do i=1, len(str)
         outstr(i:i) = " "
      enddo

      do i=1, len(str)
         c = str(i:i)
         outstr(i:i) = ''
         if ( is_lowercase(c) .or. is_uppercase(c) .or. is_digit(c) .or. c=='_' .or. c=='-' ) then
            outstr(i:i) = c
         else
            return
         endif
      enddo
      return
   end function fix_string

   logical function is_lowercase(c)
      implicit none
      character(len=one),  intent(in) :: c

      is_lowercase = .false.

      if (ichar(c) >= 97 .and. ichar(c) <= 122) is_lowercase=.true.
   end function is_lowercase

   logical function is_uppercase(c)
      implicit none
      character(len=one), intent(in) :: c

      is_uppercase = .false.

      if (ichar(c) >= 65 .and. ichar(c) <= 90) is_uppercase=.true.
   end function is_uppercase

   logical function is_digit(c)
      implicit none
      character(len=one), intent(in) :: c

      is_digit = .false.

      if (ichar(c) >= 48 .and. ichar(c) <= 57) is_digit=.true.
   end function is_digit

   elemental real function L2norm(x1,x2,x3,y1,y2,y3)
      implicit none
      real, intent(in) :: x1, x2, x3
      real, intent(in) :: y1, y2, y3

      L2norm = sqrt( (x1 - y1)**2 + (x2 - y2)**2 + (x3 - y3)**2 )
   end function L2norm

   elemental real function emag(bx,by,bz)
      use constants,  only: half
      implicit none
      real, intent(in) :: bx, by, bz

      emag = half*(bx**2 + by**2 + bz**2)

   end function emag

   elemental real function ekin(mx,my,mz,dn)
      implicit none
      real, intent(in) :: mx, my, mz, dn

      ekin = emag(mx,my,mz)/dn
   end function ekin

end module func
