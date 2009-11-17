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
!    Initial implemetation of PIERNIK code was based on TVD split MHD code by
!    Ue-Li Pen
!        see: Pen, Arras & Wong (2003) for algorithm and
!             http://www.cita.utoronto.ca/~pen/MHD
!             for original source code "mhd.f90"
!
!    For full list of developers see $PIERNIK_HOME/license/pdt.txt
!
#include "piernik.def"
!>
!! \brief [DW] Module that contains all routines related to %interactions between fluids
!!
!!
!!
!! In this module a namelist of parameters is specified:
!!
!! \copydetails interactions::init_interactions
!<
module interactions

   use constants

   real, allocatable, dimension(:,:,:,:)  :: omx0        !< array to store x component of initial velocity
   real, allocatable, dimension(:,:,:,:)  :: omy0        !< array to store y component of initial velocity
   real, allocatable, dimension(:,:)      :: alfsup      !< xy-array of values between 0 and 1, 1 where boundary support is used
   real, allocatable, dimension(:,:)      :: collfaq     !< nfluid x nfluid array of collision factors
   real :: collision_factor                              !< collision factor
   real :: cfl_interact                                  !< Courant factor for %interactions


   contains
!>
!! \brief Routine that sets the initial values of %interactions parameters from namelist @c INTERACTIONS
!>
!!
!! \n \n
!! @b INTERACTIONS
!! \n \n
!! <table border="+1">
!! <tr><td width="150pt"><b>parameter</b></td><td width="135pt"><b>default value</b></td><td width="200pt"><b>possible values</b></td><td width="315pt"> <b>description</b></td></tr>
!! <tr><td>collision_factor</td><td>0.0</td><td>real value, between 0 and 1</td><td>\copydoc interactions::collision_factor</td></tr>
!! <tr><td>cfl_interact    </td><td>0.8</td><td>real value, between 0 and 1</td><td>\copydoc interactions::cfl_interact    </td></tr>
!! </table>
!! \n \n
!<
   subroutine init_interactions
      use mpisetup
      use fluidindex,   only : nfluid
#ifdef DUST
      use fluidindex,   only : i_dst
      use initdust,     only : dragc_gas_dust
#endif /* DUST */
      implicit none
      character(LEN=100) :: par_file, tmp_log_file

      namelist /INTERACTIONS/ collision_factor, cfl_interact

      par_file = trim(cwd)//'/problem.par'
      tmp_log_file = trim(cwd)//'/tmp.log'

      collision_factor  = 0.0
      cfl_interact      = 0.8

      if(proc .eq. 0) then
         open(1,file=par_file)
            read(unit=1,nml=INTERACTIONS)
         close(1)
         open(3, file=tmp_log_file, position='append')
            write(unit=3,nml=INTERACTIONS)
         close(3)

         rbuff(1)  = collision_factor
         rbuff(2)  = cfl_interact

         call MPI_BCAST(rbuff,    buffer_dim, MPI_DOUBLE_PRECISION, 0, comm, ierr)

     else

         call MPI_BCAST(rbuff,    buffer_dim, MPI_DOUBLE_PRECISION, 0, comm, ierr)

         collision_factor     = rbuff(1)
         cfl_interact         = rbuff(2)

      endif
#ifdef COLLISIONS
      allocate(collfaq(nfluid,nfluid))
      collfaq = collision_factor
      collfaq(i_dst,:) = dragc_gas_dust
      collfaq(:,i_dst) = dragc_gas_dust
#endif /* COLLISIONS */

   end subroutine init_interactions
!>
!! \brief Routine that governs the type of interaction
!! \param sweep string of characters that points out the current sweep direction
!! \param i1 integer, number of column in the first direction after one pointed out by sweep
!! \param i2 integer, number of column in the second direction after one pointed out by sweep
!! \param n number of elements in the spatial dimension of fluid sweep array
!! \param du sweep array of fluid corrections
!! \param uu sweep fluid array
!<
   subroutine fluid_interactions(sweep, i1, i2, n, du, uu)
      use fluidindex,   only : nvar
      implicit none
      integer               :: i1,i2,n
      real                  :: dt
      real, dimension(nvar,n) :: du,ddu,uu
      character sweep*6

      du=0.0
#ifdef COLLISIONS
      call dragforce(sweep,i1,i2, n, ddu, uu)
      du = du + ddu
#endif /* COLLISIONS */

#ifdef BND_MOTION_SUPPORT
      call support_bnd_rotation(sweep,i1,i2, n, ddu, uu)
      du = du + ddu
#endif /* BND_MOTION_SUPPORT */

   end subroutine fluid_interactions

#ifdef COLLISIONS
!>
!! \brief Routine that computes contribution from a drag force
!! \param sweep string of characters that points out the current sweep direction
!! \param i1 integer, number of column in the first direction after one pointed out by sweep
!! \param i2 integer, number of column in the second direction after one pointed out by sweep
!! \param n number of elements in the spatial dimension of fluid sweep array
!! \param du sweep array of fluid corrections
!! \param uu sweep fluid array
!<
   subroutine dragforce(sweep, i1, i2, n, du, uu)
      use fluidindex,   only : nvar,nfluid,iarr_all_dn,iarr_all_mx,iarr_all_my,nadiab
#ifndef ISO
      use fluidindex,   only : iarr_all_en
#endif /* !ISO */
      use grid,         only : maxxyz,x,y,z
      implicit none
      integer               :: i1,i2,n
      real                  :: dt
      real, dimension(nvar,n) :: du,ddu,uu
      character sweep*6
      integer ifl,jfl
      real, dimension(nfluid,nfluid,n) :: flch
      real, dimension(nfluid,n)        :: colls,velcoor
      real, dimension(maxxyz)          :: r1,r2
      real    :: a1
      integer :: rend

      du=0.0
!-----collisions between fluids------
      select case(sweep)
         case('xsweep')
            a1    = y(i1)
            rend  = size(x)
            r1(1:rend) = x(:)                            ! r1   max(size(x),size(y),size(z))
            r2(1:rend) = a1 / (r1(1:rend)*r1(1:rend) + a1 * a1)
         case('ysweep')
            a1    = x(i2)
            rend  = size(y)
            r1(1:rend) = y(1:rend)
            r2(1:rend) = a1 / (r1(1:rend)*r1(1:rend) + a1 * a1)
         case('zsweep')
            a1    = 1.0
            rend  = size(z)
            r1(1:rend) = 0.0
            r2(1:rend) = 1.0
      end select

      do ifl=1,nfluid
         do jfl=1,nfluid
            if(ifl .ne. jfl) then
               flch(ifl,jfl,:)= collfaq(ifl,jfl)*uu(iarr_all_dn(ifl),:) * uu(iarr_all_dn(jfl),:) &
                        *( uu(iarr_all_mx(jfl),:)/uu(iarr_all_dn(jfl),:) &
                         - uu(iarr_all_mx(ifl),:)/uu(iarr_all_dn(ifl),:))
            else
               flch(ifl,jfl,:)=0.0
            endif
         enddo
      enddo
      do ifl=1,nfluid
         colls(ifl,:)=sum(flch(ifl,:,:),1)
      enddo
#ifndef ISO
      du(iarr_all_en(1:nadiab),:)=uu(iarr_all_mx(1:nadiab),:)*colls(1:nadiab,:)
#endif /* ISO */
      du(iarr_all_mx,:)=colls

   end subroutine dragforce
#endif /* COLLISIONS */

#ifdef BND_MOTION_SUPPORT
!>
!! \brief Routine that computes a support for fluid rotation in area close to domain boundaries
!! \param sweep string of characters that points out the current sweep direction
!! \param i1 integer, number of column in the first direction after one pointed out by sweep
!! \param i2 integer, number of column in the second direction after one pointed out by sweep
!! \param n number of elements in the spatial dimension of fluid sweep array
!! \param du sweep array of fluid corrections
!! \param uu sweep fluid array
!<
   subroutine support_bnd_rotation(sweep, i1, i2, n, du, uu)
      use fluidindex,   only : nvar,nfluid,iarr_all_dn,iarr_all_mx,iarr_all_my,nadiab
#ifndef ISO
      use fluidindex,   only : iarr_all_en
#endif /* !ISO */
      use grid,         only : maxxyz,x,y,z
      implicit none
      integer               :: i1,i2,n
      real                  :: dt
      real, dimension(nvar,n) :: du,ddu,uu
      character sweep*6
      integer i,j,k,ifl,ni1,ni2
      real, dimension(nfluid,n) :: kplrsup,velcoor,vel0
      real, dimension(maxxyz)   :: r1,r2
      real :: a1
      integer :: rend,ii

      du=0.0
      select case(sweep)
         case('xsweep')
            a1    = y(i1)
            rend  = size(x)
            r1(1:rend) = x(:)                            ! r1   max(size(x),size(y),size(z))
            r2(1:rend) = a1 / (r1(1:rend)*r1(1:rend) + a1 * a1)
            vel0(:,:)  = omx0(:,:,i1,i2)
         case('ysweep')
            a1    = x(i2)
            rend  = size(y)
            r1(1:rend) = y(1:rend)
            r2(1:rend) = a1 / (r1(1:rend)*r1(1:rend) + a1 * a1)
            vel0(:,:)  = omy0(:,i2,:,i1)
         case('zsweep')
            a1    = 1.0
            rend  = size(z)
            r1(1:rend) = 0.0
            r2(1:rend) = 1.0
      end select

      if(sweep .ne. 'zsweep') then
         do ifl=1,nfluid
#ifdef KEPL_SUPP_SIMX
            velcoor(ifl,:)=uu(iarr_all_mx(ifl),:)/uu(iarr_all_dn(ifl),:)
            kplrsup(ifl,:)=-alfsup(:,i1)*(velcoor(ifl,:)-vel0(ifl,:)) !*uu(iarr_all_dn(ifl),:)
#else /* KEPL_SUPP_SIMX */
            velcoor(ifl,:)=(uu(iarr_all_mx(ifl),:)*a1-uu(iarr_all_my(ifl),:)*r1)/uu(iarr_all_dn(ifl),:)*r2
            kplrsup(ifl,:)=-alfsup(:,i1)*(velcoor(ifl,:)-vel0(ifl,:))!*uu(iarr_all_dn(ifl),:)
#endif /* KEPL_SUPP_SIMX */
         enddo
      else
         kplrsup(:,:)=0.0
      endif
#ifndef ISO
!      Duus(iena(fadiab),:)=(uu(imxa(fadiab),:)*kplrsup*uu(idna(fadiab),:)*dt &
!              +0.5*(kplrsup*uu(idna(fadiab),:)*dt)**2)/uu(idna(fadiab),:)
      du(iarr_all_en(1:nadiab),:)=kplrsup(1:nadiab,:)*uu(iarr_all_mx(1:nadiab),:)
#endif /* ISO */
      du(iarr_all_mx,:)=kplrsup(:,:)*uu(iarr_all_dn,:)

   end subroutine support_bnd_rotation
#endif /* BND_MOTION_SUPPORT */

end module interactions
