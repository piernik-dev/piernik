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
!! \brief Module that contains all routines related to %interactions between fluids
!!
!! In this module a namelist of parameters is specified:
!!
!! \copydetails interactions::init_interactions
!<
module interactions
! pulled by ANY

   use constants, only: cbuff_len

   implicit none

   private
   public :: init_interactions, cleanup_interactions, fluid_interactions, collfaq, cfl_interact, dragc_gas_dust, has_interactions, &
      & interactions_grace_passed, epstein_factor, balsara_implicit_interactions, update_grain_size

   real, allocatable, dimension(:,:)          :: collfaq            !< flind%fluids x flind%fluids array of collision factors
   real                                       :: collision_factor   !< collision factor
   real                                       :: cfl_interact       !< Courant factor for %interactions
   real                                       :: dragc_gas_dust     !< \deprecated remove me
!   real                                       :: taus               !< stopping time
   real                                       :: grain_size         !< size of dust grains in cm
   real                                       :: grain_dens         !< density of dust grains in g/cm^3
   real, dimension(:), allocatable, protected :: epstein_factor     !< grain_size * grain_dens / c_s for iso case
   character(len=cbuff_len)                   :: interactions_type  !< type of interactions between fluids
   logical                                    :: has_interactions

   interface
      function fluid_interactions_iface(dens,vel) result(acc)
         implicit none
         real, dimension(:,:), pointer, intent(in)  :: dens
         real, dimension(:,:), pointer, intent(in)  :: vel
         real, dimension(size(dens,1),size(dens,2)) :: acc
      end function fluid_interactions_iface
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
!! <tr><td>collision_factor</td><td>0.0    </td><td>real value, between 0 and 1</td><td>\copydoc interactions::collision_factor</td></tr>
!! <tr><td>cfl_interact    </td><td>0.8    </td><td>real value, between 0 and 1</td><td>\copydoc interactions::cfl_interact    </td></tr>
!! <tr><td>dragc_gas_dust  </td><td>0.0    </td><td>real value                 </td><td>\copydoc interactions::dragc_gas_dust  </td></tr>
!! <tr><td>has_interactions</td><td>.false.</td><td>logical value              </td><td>\copydoc interactions::has_interactions</td></tr>
!! <tr><td>grain_size      </td><td>10.0   </td><td>real value                 </td><td>\copydoc interactions::dragc_gas_dust  </td></tr>
!! <tr><td>grain_dens      </td><td>1.6    </td><td>real value                 </td><td>\copydoc interactions::dragc_gas_dust  </td></tr>
!! <tr><td>interactions_type</td><td>'none'</td><td>'none'/'aerodrag'/'aerodrag_ep'/'dragforce_dw'/'balsara'</td><td>\copydoc interactions::interactions_type</td></tr>
!! </table>
!! \n \n
!<
   subroutine init_interactions

      use constants,     only: PIERNIK_INIT_FLUIDS, cbuff_len
      use dataio_pub,    only: die, code_progress, par_file, ierrh, namelist_errh, compare_namelist, cmdl_nml, lun      ! QA_WARN required for diff_nml
      use fluidindex,    only: flind
      use mpisetup,      only: master, slave, cbuff, lbuff, rbuff, buffer_dim, mpi_err, comm, FIRST!, grace_period_passed
      use mpi,           only: MPI_DOUBLE_PRECISION, MPI_LOGICAL, MPI_CHARACTER
      use units,         only: cm, gram

      implicit none

      namelist /INTERACTIONS/ collision_factor, cfl_interact, dragc_gas_dust, has_interactions, grain_size, grain_dens

      if (code_progress < PIERNIK_INIT_FLUIDS) call die("[interactions:init_interactions] fluids not initialized.")

      collision_factor  = 0.0
      cfl_interact      = 0.8
      dragc_gas_dust    = 0.0
      grain_size        = 10.0
      grain_dens        = 1.6

      has_interactions  = .false.

      interactions_type = 'none'

      if (master) then

         diff_nml(INTERACTIONS)

         rbuff(1)  = collision_factor
         rbuff(2)  = cfl_interact
         rbuff(3)  = dragc_gas_dust
         rbuff(4)  = grain_size
         rbuff(5)  = grain_dens

         lbuff(1)  = has_interactions

         cbuff(1)  = interactions_type

      endif

      call MPI_Bcast(rbuff,           buffer_dim, MPI_DOUBLE_PRECISION, FIRST, comm, mpi_err)
      call MPI_Bcast(lbuff,           buffer_dim, MPI_LOGICAL,          FIRST, comm, mpi_err)
      call MPI_Bcast(cbuff, cbuff_len*buffer_dim, MPI_CHARACTER,        FIRST, comm, mpi_err)

      if (slave) then

         collision_factor     = rbuff(1)
         cfl_interact         = rbuff(2)
         dragc_gas_dust       = rbuff(3)
         grain_size           = rbuff(4)
         grain_dens           = rbuff(5)

         has_interactions     = lbuff(1)

         interactions_type    = cbuff(1)

      endif

      fluid_interactions => fluid_interactions_dummy
      grain_dens = grain_dens * gram * cm**(-2)  ! It is always grain_size * grain_dens, hence unit g/cm2 in grain_dens
      if (.not.allocated(epstein_factor)) allocate(epstein_factor(flind%fluids))

      epstein_factor = 0.0

!      if (grace_period_passed()) call interactions_grace_passed

   end subroutine init_interactions

   subroutine cleanup_interactions

      implicit none

      if (allocated(epstein_factor)) deallocate(epstein_factor)

   end subroutine cleanup_interactions

   subroutine interactions_grace_passed

      use constants,     only: DST
      use dataio_pub,    only: warn, printinfo, die, msg
      use fluidindex,    only: flind
      use mpisetup,      only: master

      implicit none

      logical, save :: warned = .false.
      integer       :: i

      if (associated(flind%dst)) then
         do i = 1, flind%fluids
            if (flind%all_fluids(i)%fl%tag /= DST) then
               epstein_factor(flind%all_fluids(i)%fl%pos) = grain_size * grain_dens / flind%all_fluids(i)%fl%cs !BEWARE iso assumed
               if (epstein_factor(flind%all_fluids(i)%fl%pos) <= 0.0) &
                  call warn("[interactions:interactions_grace_passed] epstein_factor <= 0.0, that's not good :/")
            else
               epstein_factor(flind%all_fluids(i)%fl%pos) = 0.0
            endif
         enddo
         ! has_interactions = .true.    !> \deprecated BEWARE: temporary hack,  switches on timestep_interactions, don't needed in implicit solver??
      else
         if (.not. warned .and. master) call warn("[interactions:interactions_grace_passed] Cannot initialize aerodynamical drag because dust does not exist.")
         warned = .true.
      endif

#ifndef BALSARA
      if (dragc_gas_dust > 0.0 .or. collision_factor > 0.0) then
         if (associated(flind%dst)) then
            if (master) call printinfo("[interactions:interactions_grace_passed] Initializing aerodynamical drag")
            allocate(collfaq(flind%fluids,flind%fluids))
            collfaq = collision_factor
            collfaq(flind%dst%pos,:) = dragc_gas_dust
            collfaq(:,flind%dst%pos) = dragc_gas_dust

            select case (interactions_type)
               case ('none')
                  ! already pointed to fluid_interactions_dummy
               case ('aerodrag')
                  fluid_interactions => fluid_interactions_aero_drag
               case ('aerodrag_ep')
                  fluid_interactions => fluid_interactions_aero_drag_ep
               case ('dragforce_dw')
                  fluid_interactions => dragforce
               case ('balsara')
                  call die("[interactions:interactions_grace_passed] precompiler flag BALSARA has to be defined")
               case default
                  write(msg,'(2a)') "[interactions:interactions_grace_passed] unknown value of interactions_type ", interactions_type
                  call die(msg)
            end select
            has_interactions = .true.    !> \deprecated BEWARE: temporary hack,  switches on timestep_interactions
         else
            if (.not. warned .and. master) call warn("[interactions:interactions_grace_passed] Cannot initialize aerodynamical drag because dust does not exist.")
            warned = .true.
         endif
      endif
#endif /* !BALSARA */

   end subroutine interactions_grace_passed

   function fluid_interactions_dummy(dens, velx) result(acc)

      implicit none

      real, dimension(:,:), pointer, intent(in)  :: dens
      real, dimension(:,:), pointer, intent(in)  :: velx
      real, dimension(size(dens,1),size(dens,2)) :: acc

      acc = 0.0

      if (.false.) acc(1,1) = acc(1,1) + 0.*velx(1,1) ! suppress compiler warning on unused argument

   end function fluid_interactions_dummy

   function fluid_interactions_aero_drag(dens, velx) result(acc)

      use fluidindex,       only: flind

      implicit none

      real, dimension(:,:), intent(in), pointer  :: dens
      real, dimension(:,:), intent(in), pointer  :: velx
      real, dimension(size(dens,1),size(dens,2)) :: acc

      acc(flind%dst%pos,:) = -dragc_gas_dust * (velx(flind%dst%pos,:) - velx(flind%neu%pos,:))
      acc(flind%neu%pos,:) = -acc(flind%dst%pos,:) * dens(flind%dst%pos,:) / dens(flind%neu%pos,:)
!      acc(flind%neu%pos,:) = -dragc_gas_dust * dens(flind%dst%pos,:) / dens(flind%neu%pos,:) * ( velx(flind%neu%pos,:) - velx(flind%dst%pos,:) )

   end function fluid_interactions_aero_drag

   function fluid_interactions_aero_drag_ep(dens, velx) result(acc)

      use fluidindex,       only: flind

      implicit none

      real, dimension(:,:), intent(in), pointer  :: dens
      real, dimension(:,:), intent(in), pointer  :: velx
      real, dimension(size(dens,1),size(dens,2)) :: acc

      acc(flind%dst%pos,:) = -dens(flind%neu%pos,:) / epstein_factor(flind%neu%pos) * (velx(flind%dst%pos,:) - velx(flind%neu%pos,:))
      acc(flind%neu%pos,:) = -acc(flind%dst%pos,:) * dens(flind%dst%pos,:) / dens(flind%neu%pos,:)

   end function fluid_interactions_aero_drag_ep

!>
!! \details
!! Balsara Dinshaw S., Tilley David A., Rettig Terrence, Brittain Sean D. MNRAS (2009) 397: 24.
!! Tilley, David A., Balsara, Dinshaw S. MNRAS (2008) 389: 1058.
!<
   subroutine balsara_implicit_interactions(u1,u0,vx,cs_iso2,istep)

      use constants,  only: half, one, zero
      use fluidindex, only: iarr_all_dn, iarr_all_mx, flind
      use global,     only: dt

      implicit none

      real, dimension(:,:), intent(inout)      :: u1
      real, dimension(:,:), intent(in)         :: u0
      real, dimension(:,:), intent(in)         :: vx
      real, dimension(:), pointer,  intent(in) :: cs_iso2
      integer, intent(in)                      :: istep

      real, dimension(flind%fluids,size(u1,2)) :: vprim
      real, dimension(size(u1,2))              :: delta, drag
      !>
      !! \deprecated BEWARE: this bit assumes that we have 2 fluids && u1 == u0 - \grad F
      !! \todo 2) half-time step should be \le \frac{1}{2}\frac{c_s}{drag * \rho\prim_? |v'_d - v'_g|}
      !! \todo 3) what if not isothermal?
      !! \todo 4) remove hardcoded integers
      !<
      if (epstein_factor(flind%neu%pos) <= zero) return
      drag(:) = dt*half / (grain_size * grain_dens) * sqrt( cs_iso2(:) + abs( vx(1,:) - vx(2,:) )**2)

      delta(:) = one + drag(:) * (u1(iarr_all_dn(1),:) + u1(iarr_all_dn(2),:))
      delta(:) = one/delta(:)

      if (istep == 1) then
         vprim(1,:) =  delta(:)*( (1./u1(iarr_all_dn(1),:) + drag(:))*u1(iarr_all_mx(1),:) + drag(:)*u1(iarr_all_mx(2),:) )
         vprim(2,:) =  delta(:)*( (1./u1(iarr_all_dn(2),:) + drag(:))*u1(iarr_all_mx(2),:) + drag(:)*u1(iarr_all_mx(1),:) )
      else
         vprim(2,:) =  delta(:)*(  &
            drag(:)*( u0(iarr_all_dn(2),:)/u1(iarr_all_dn(2),:)*u0(iarr_all_mx(1),:) - u0(iarr_all_dn(1),:)/u1(iarr_all_dn(2),:)*u0(iarr_all_mx(2),:) &
                    + u1(iarr_all_mx(1),:) ) + u1(iarr_all_mx(2),:) * ( 1./u1(iarr_all_dn(2),:) + drag(:) )  )
         vprim(1,:) =  delta(:)*(  &
            drag(:)*( u0(iarr_all_dn(1),:)/u1(iarr_all_dn(1),:)*u0(iarr_all_mx(2),:) - u0(iarr_all_dn(2),:)/u1(iarr_all_dn(1),:)*u0(iarr_all_mx(1),:) &
                    + u1(iarr_all_mx(2),:) ) + u1(iarr_all_mx(1),:) * ( 1./u1(iarr_all_dn(1),:) + drag(:) )  )
      endif
      u1(iarr_all_mx,:) = u1(iarr_all_dn,:) * vprim(:,:)
      return
   end subroutine balsara_implicit_interactions

   subroutine update_grain_size(new_size)
      use fluidindex, only: flind
      use constants,  only: DST
      implicit none
      real, intent(in) :: new_size
      integer          :: i

      grain_size = new_size

      do i = 1, flind%fluids
         if (flind%all_fluids(i)%fl%tag /= DST) then
            epstein_factor(flind%all_fluids(i)%fl%pos) = grain_size * grain_dens / flind%all_fluids(i)%fl%cs !BEWARE iso assumed
         else
            epstein_factor(flind%all_fluids(i)%fl%pos) = 0.0
         endif
      enddo

   end subroutine update_grain_size

!>
!! \brief Function that computes contribution from a drag force for all fluids
!! \param dens density flind%fluids x ncells array
!! \param velx velocity component of the current sweep direction
!! \result acc acceleration due to drag force for all fluids
!<
   function dragforce(dens,velx) result(acc)

      use fluidindex, only: flind

      implicit none

      real, dimension(:,:), intent(in), pointer  :: dens
      real, dimension(:,:), intent(in), pointer  :: velx
      real, dimension(size(dens,1),size(dens,2)) :: acc
      integer                                    :: ifl, jfl

      acc(:,:) = 0.0
!-----collisions between fluids------
      do ifl=1,flind%fluids
         do jfl=1,flind%fluids
            if (ifl /= jfl) acc(ifl,:) = acc(ifl,:) + collfaq(ifl,jfl) * dens(ifl,:) * dens(jfl,:) * (velx(jfl,:) - velx(ifl,:)) ! / dens(ifl,:) ?
         enddo
      enddo

   end function dragforce

end module interactions
