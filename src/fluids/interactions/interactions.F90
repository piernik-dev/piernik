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
#include "macros.h"
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

   real, allocatable, dimension(:,:,:,:)  :: omx0        !< array to store x component of initial velocity
   real, allocatable, dimension(:,:,:,:)  :: omy0        !< array to store y component of initial velocity
   real, allocatable, dimension(:,:)      :: alfsup      !< xy-array of values between 0 and 1, 1 where boundary support is used
   real, allocatable, dimension(:,:)      :: collfaq     !< nvar%fluids x nvar%fluids array of collision factors
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
      use mpisetup,      only: proc, rbuff, MPI_DOUBLE_PRECISION, buffer_dim, ierr, comm
      use fluidindex,    only: nvar
      use dataio_public, only: par_file, cwd, ierrh, namelist_errh
      use func,          only: compare_namelist

#ifdef DUST
      use initdust,      only: dragc_gas_dust
#endif /* DUST */

      implicit none

      namelist /INTERACTIONS/ collision_factor, cfl_interact

      collision_factor  = 0.0
      cfl_interact      = 0.8

      if (proc .eq. 0) then

         diff_nml(INTERACTIONS)

         rbuff(1)  = collision_factor
         rbuff(2)  = cfl_interact

      endif

      call MPI_Bcast(rbuff,    buffer_dim, MPI_DOUBLE_PRECISION, 0, comm, ierr)

      if (proc /= 0) then

         collision_factor     = rbuff(1)
         cfl_interact         = rbuff(2)

      endif

#ifdef COLLISIONS
      allocate(collfaq(nvar%fluids,nvar%fluids))
      collfaq = collision_factor
      collfaq(nvar%dst%pos,:) = dragc_gas_dust
      collfaq(:,nvar%dst%pos) = dragc_gas_dust
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
!!
!! This is the place for contributions from any fluid %interactions that can be defined in this module. \n
!! The manner to add another contribution is following:
!! \code
!!     call defined_interaction(sweep,i1,i2,n,ddu,uu)
!!     du = du + ddu
!! \endcode
!! where \c defined_interaction has to be specified as a subroutine in this module.
!<
   subroutine fluid_interactions(sweep, i1, i2, n, du, uu)
      use fluidindex,   only: nvar
      implicit none
      integer, intent(in)   :: i1,i2,n
      real, dimension(nvar%all,n)  :: du,uu
      character(len=6), intent(in) :: sweep
#ifdef COLLISIONS
      real, dimension(nvar%all,n)  :: ddu
#endif /* COLLISIONS */

      du=0.0
#ifdef COLLISIONS
      call dragforce(sweep,i1,i2, n, ddu, uu)
      du = du + ddu
#endif /* COLLISIONS */

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
      use fluidindex,   only: nvar,iarr_all_dn,iarr_all_mx,iarr_all_my
#ifndef ISO
      use fluidindex,   only: iarr_all_en
#endif /* !ISO */
      use grid,         only: maxxyz,x,y,z
      implicit none
      integer, intent(in)   :: i1,i2,n
      real, dimension(nvar%all,n)  :: du,uu
      character(len=6), intent(in) :: sweep
      integer  :: ifl,jfl
      real, dimension(nvar%fluids,nvar%fluids,n) :: flch
      real, dimension(nvar%fluids,n)             :: colls
      real, dimension(maxxyz)                    :: r1,r2
      real    :: a1
      integer :: rend

      du=0.0
!-----collisions between fluids------
      select case (sweep)
         case ('xsweep')
            a1    = y(i1)
            rend  = size(x)
            r1(1:rend) = x(:)                            ! r1   max(size(x),size(y),size(z))
            r2(1:rend) = a1 / (r1(1:rend)*r1(1:rend) + a1 * a1)
         case ('ysweep')
            a1    = x(i2)
            rend  = size(y)
            r1(1:rend) = y(1:rend)
            r2(1:rend) = a1 / (r1(1:rend)*r1(1:rend) + a1 * a1)
         case ('zsweep')
            a1    = 1.0
            rend  = size(z)
            r1(1:rend) = 0.0
            r2(1:rend) = 1.0
      end select

      do ifl=1,nvar%fluids
         do jfl=1,nvar%fluids
            if (ifl .ne. jfl) then
               flch(ifl,jfl,:)= collfaq(ifl,jfl)*uu(iarr_all_dn(ifl),:) * uu(iarr_all_dn(jfl),:) &
                        *( uu(iarr_all_mx(jfl),:)/uu(iarr_all_dn(jfl),:) &
                         - uu(iarr_all_mx(ifl),:)/uu(iarr_all_dn(ifl),:))
            else
               flch(ifl,jfl,:)=0.0
            endif
         enddo
      enddo
      do ifl=1,nvar%fluids
         colls(ifl,:)=sum(flch(ifl,:,:),1)
      enddo
#ifndef ISO
      du(iarr_all_en(1:nvar%adiab),:)=uu(iarr_all_mx(1:nvar%adiab),:)*colls(1:nvar%adiab,:)
#endif /* !ISO */
      du(iarr_all_mx,:)=colls

   end subroutine dragforce
#endif /* COLLISIONS */

end module interactions
