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
#include "macros.h"
!>
!! \brief (DW) Module that contains all routines related to %interactions between fluids
!!
!! In this module a namelist of parameters is specified:
!!
!! \copydetails interactions::init_interactions
!<
module interactions
! pulled by ANY

   implicit none

   private
   public :: init_interactions, fluid_interactions, collfaq, cfl_interact, dragc_gas_dust, has_interactions, &
      & interactions_grace_passed

   real, allocatable, dimension(:,:)      :: collfaq     !< nvar%fluids x nvar%fluids array of collision factors
   real :: collision_factor                              !< collision factor
   real :: cfl_interact                                  !< Courant factor for %interactions
   real :: dragc_gas_dust                                !< \deprecated remove me
   real :: taus                                          !< stopping time
   logical :: has_interactions

   interface
      subroutine fluid_interactions_iface(dens,vel,acc)
         implicit none
         real, dimension(:,:), pointer, intent(in)  :: dens
         real, dimension(:,:), pointer, intent(in)  :: vel
         real, dimension(size(dens,1),size(dens,2)), intent(out) :: acc
      end subroutine fluid_interactions_iface
   end interface

   procedure(fluid_interactions_iface), pointer :: fluid_interactions

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
!! <tr><td>dragc_gas_dust  </td><td>0.0</td><td>real value</td>                 <td>\copydoc interactions::dragc_gas_dust</td></tr>
!! <tr><td>cfl_interact    </td><td>0.8</td><td>real value, between 0 and 1</td><td>\copydoc interactions::cfl_interact</td></tr>
!! </table>
!! \n \n
!<
   subroutine init_interactions

      use dataio_pub,    only: par_file, ierrh, namelist_errh, compare_namelist, cmdl_nml      ! QA_WARN required for diff_nml
      use mpisetup,      only: master, slave, lbuff, rbuff, buffer_dim, ierr, comm!, grace_period_passed
      use mpi,           only: MPI_DOUBLE_PRECISION, MPI_LOGICAL

      implicit none

      namelist /INTERACTIONS/ collision_factor, cfl_interact, dragc_gas_dust, has_interactions

      collision_factor  = 0.0
      cfl_interact      = 0.8
      dragc_gas_dust    = 0.0

      has_interactions  = .false.

      if (master) then

         diff_nml(INTERACTIONS)

         rbuff(1)  = collision_factor
         rbuff(2)  = cfl_interact
         rbuff(3)  = dragc_gas_dust

         lbuff(1)  = has_interactions

      endif

      call MPI_Bcast(rbuff,    buffer_dim, MPI_DOUBLE_PRECISION, 0, comm, ierr)
      call MPI_Bcast(lbuff,    buffer_dim, MPI_LOGICAL,          0, comm, ierr)

      if (slave) then

         collision_factor     = rbuff(1)
         cfl_interact         = rbuff(2)
         dragc_gas_dust       = rbuff(3)

         has_interactions     = lbuff(1)

      endif

      fluid_interactions => fluid_interactions_dummy
!      if (grace_period_passed()) call interactions_grace_passed

   end subroutine init_interactions

   subroutine interactions_grace_passed

      use fluidindex,    only: nvar
#ifdef VERBOSE
      use dataio_pub,    only: printinfo
#endif /* VERBOSE */

      implicit none

      if (dragc_gas_dust > 0.0 .or. collision_factor > 0.0) then
#ifdef VERBOSE
         call printinfo("[interactions:interactions_grace_passed] Initializing aerodynamical drag")
#endif /* VERBOSE */
         allocate(collfaq(nvar%fluids,nvar%fluids))
         collfaq = collision_factor
         collfaq(nvar%dst%pos,:) = dragc_gas_dust
         collfaq(:,nvar%dst%pos) = dragc_gas_dust

         taus = 1. / dragc_gas_dust
         fluid_interactions => fluid_interactions_aero_drag
         has_interactions = .true.    !> \depracted BEWARE: temporary hack
      endif

   end subroutine interactions_grace_passed

   subroutine fluid_interactions_dummy(dens, velx, acc)

      implicit none

      real, dimension(:,:), intent(in)  :: dens
      real, dimension(:,:), intent(in)  :: velx
      real, dimension(size(dens,1),size(dens,2)), intent(out) :: acc

      acc = 0.0

      if (.false.) acc(1,1) = acc(1,1) + 0.*velx(1,1) ! suppress compiler warning on unused argument

   end subroutine fluid_interactions_dummy

   subroutine fluid_interactions_aero_drag(dens, velx, acc)

      use fluidindex,       only: nvar

      implicit none

      real, dimension(:,:), intent(in), pointer :: dens
      real, dimension(:,:), intent(in), pointer :: velx
      real, dimension(size(dens,1),size(dens,2)), intent(out) :: acc

      acc(nvar%dst%pos,:) = -dragc_gas_dust * (velx(nvar%dst%pos,:) - velx(nvar%neu%pos,:))
      acc(nvar%neu%pos,:) = -acc(nvar%dst%pos,:) * dens(nvar%dst%pos,:) / dens(nvar%neu%pos,:)
!      acc(nvar%neu%pos,:) = -dragc_gas_dust * dens(nvar%dst%pos,:) / dens(nvar%neu%pos,:) * ( velx(nvar%neu%pos,:) - velx(nvar%dst%pos,:) )

   end subroutine fluid_interactions_aero_drag

#ifdef FLUID_INTERACTIONS_DW
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
   subroutine fluid_interactions_dw(sweep, i1, i2, n, du, uu)
      use fluidindex,   only: nvar
      implicit none
      integer, intent(in)   :: i1,i2,n
      real, dimension(nvar%all,n)  :: du,uu
      character(len=*), intent(in) :: sweep
      real, dimension(nvar%all,n)  :: ddu

      du=0.0
      call dragforce(sweep,i1,i2, n, ddu, uu)
      du = du + ddu

   end subroutine fluid_interactions_dw

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
      use fluidindex,   only: nvar, iarr_all_dn, iarr_all_mx
      use grid,         only: cg
#ifndef ISO
      use fluidindex,   only: iarr_all_en
#endif /* !ISO */
      implicit none
      real                  :: a1
      integer, intent(in)   :: i1, i2, n
      real, dimension(nvar%all,n)  :: du,uu
      character(len=*), intent(in) :: sweep
      integer  :: ifl,jfl,rend
      real, dimension(nvar%fluids,nvar%fluids,n) :: flch
      real, dimension(nvar%fluids,n)             :: colls
      real, dimension(cg%maxxyz)                    :: r1,r2

      du=0.0
!-----collisions between fluids------
      select case (sweep)
         case ('xsweep')
            a1    = cg%y(i1)
            rend  = size(cg%x)
            r1(1:rend) = cg%x(:)                            ! r1   max(size(x),size(y),size(z))
            r2(1:rend) = a1 / (r1(1:rend)*r1(1:rend) + a1 * a1)
         case ('ysweep')
            a1    = cg%x(i2)
            rend  = size(cg%y)
            r1(1:rend) = cg%y(1:rend)
            r2(1:rend) = a1 / (r1(1:rend)*r1(1:rend) + a1 * a1)
         case ('zsweep')
            a1    = 1.0
            rend  = size(cg%z)
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
      du(iarr_all_en(1:nvar%energ),:)=uu(iarr_all_mx(1:nvar%energ),:)*colls(1:nvar%energ,:)
#endif /* !ISO */
      du(iarr_all_mx,:)=colls

   end subroutine dragforce
#endif /* FLUID_INTERACTIONS_DW */
end module interactions
