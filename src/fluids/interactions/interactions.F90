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
   public :: init_interactions, cleanup_interactions, fluid_interactions_exec, collfaq, cfl_interact, dragc_gas_dust, has_interactions, &
      & interactions_grace_passed, epstein_factor, balsara_implicit_interactions, update_grain_size

   real, allocatable, dimension(:,:)          :: collfaq            !< flind%fluids x flind%fluids array of collision factors
   real                                       :: collision_factor   !< collision factor
   real                                       :: cfl_interact       !< Courant factor for %interactions
   real                                       :: dragc_gas_dust     !< Drag coefficient \deprecated remove me
   real                                       :: grain_size         !< size of dust grains in cm
   real                                       :: grain_dens         !< density of dust grains in g/cm^3
   real                                       :: grain_dens_x_size  !< size times density of dust grains in g/cm^2
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
!! <tr><td>grain_size      </td><td>10.0   </td><td>real value                 </td><td>\copydoc interactions::grain_size      </td></tr>
!! <tr><td>grain_dens      </td><td>1.6    </td><td>real value                 </td><td>\copydoc interactions::grain_dens      </td></tr>
!! <tr><td>interactions_type</td><td>'none'</td><td>'none' \n 'aerodrag' \n 'aerodrag_ep' \n 'dragforce_dw' \n 'balsara'</td><td>\copydoc interactions::interactions_type</td></tr>
!! </table>
!! \n \n
!<
   subroutine init_interactions

      use bcast,      only: piernik_MPI_Bcast
      use constants,  only: PIERNIK_INIT_FLUIDS, cbuff_len
      use dataio_pub, only: die, code_progress, nh
      use fluidindex, only: flind
      use mpisetup,   only: master, slave, cbuff, lbuff, rbuff
      use units,      only: cm, gram

      implicit none

      namelist /INTERACTIONS/ collision_factor, cfl_interact, dragc_gas_dust, has_interactions, grain_size, grain_dens, interactions_type

      if (code_progress < PIERNIK_INIT_FLUIDS) call die("[interactions:init_interactions] fluids not initialized.")

      collision_factor  = 0.0
      cfl_interact      = 0.8
      dragc_gas_dust    = 0.0
      grain_size        = 10.0
      grain_dens        = 1.6

      has_interactions  = .false.

      interactions_type = 'none'

      if (master) then

         if (.not.nh%initialized) call nh%init()
         open(newunit=nh%lun, file=nh%tmp1, status="unknown")
         write(nh%lun,nml=INTERACTIONS)
         close(nh%lun)
         open(newunit=nh%lun, file=nh%par_file)
         nh%errstr=""
         read(unit=nh%lun, nml=INTERACTIONS, iostat=nh%ierrh, iomsg=nh%errstr)
         close(nh%lun)
         call nh%namelist_errh(nh%ierrh, "INTERACTIONS")
         read(nh%cmdl_nml,nml=INTERACTIONS, iostat=nh%ierrh)
         call nh%namelist_errh(nh%ierrh, "INTERACTIONS", .true.)
         open(newunit=nh%lun, file=nh%tmp2, status="unknown")
         write(nh%lun,nml=INTERACTIONS)
         close(nh%lun)
         call nh%compare_namelist()

         rbuff(1)  = collision_factor
         rbuff(2)  = cfl_interact
         rbuff(3)  = dragc_gas_dust
         rbuff(4)  = grain_size
         rbuff(5)  = grain_dens

         lbuff(1)  = has_interactions

         cbuff(1)  = interactions_type

      endif

      call piernik_MPI_Bcast(rbuff)
      call piernik_MPI_Bcast(lbuff)
      call piernik_MPI_Bcast(cbuff, cbuff_len)

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
      grain_dens_x_size = grain_dens * grain_size * gram / cm**2
      if (.not.allocated(epstein_factor)) allocate(epstein_factor(flind%fluids))

      epstein_factor = 0.0

!      if (grace_period_passed()) call interactions_grace_passed

   end subroutine init_interactions

   subroutine cleanup_interactions

      implicit none

      if (allocated(epstein_factor)) deallocate(epstein_factor)
      if (allocated(collfaq)) deallocate(collfaq)

   end subroutine cleanup_interactions

   subroutine interactions_grace_passed

      use constants,  only: DST
      use dataio_pub, only: warn, printinfo, die
#ifndef BALSARA
      use dataio_pub, only: msg
#endif /* !BALSARA */
      use fluidindex, only: flind
      use mpisetup,   only: master

      implicit none

      logical, save :: warned = .false.
      integer       :: i

      if (associated(flind%dst)) then
         do i = 1, flind%fluids
            if (flind%all_fluids(i)%fl%tag /= DST) then
               epstein_factor(flind%all_fluids(i)%fl%pos) = grain_dens_x_size / flind%all_fluids(i)%fl%cs !BEWARE iso assumed
               if (epstein_factor(flind%all_fluids(i)%fl%pos) <= 0.0) &
                  call warn("[interactions:interactions_grace_passed] epstein_factor <= 0.0, that's not good :/")
            else
               epstein_factor(flind%all_fluids(i)%fl%pos) = 0.0
            endif
         enddo
         ! has_interactions = .true.    !> \deprecated BEWARE: temporary hack,  switches on timestep_interactions, don't needed in implicit solver??
      else
         if (has_interactions .and. .not. warned .and. master) call warn("[interactions:interactions_grace_passed] Cannot initialize aerodynamical drag because dust does not exist.")
         warned = .true.
      endif

#ifndef BALSARA
      if (dragc_gas_dust > 0.0 .or. collision_factor > 0.0) then
         if (associated(flind%dst)) then
            if (master) call printinfo("[interactions:interactions_grace_passed] Initializing aerodynamical drag")
            allocate(collfaq(flind%fluids, flind%fluids))
            collfaq = collision_factor
            collfaq(:, flind%dst%pos) = dragc_gas_dust
            collfaq(flind%dst%pos, :) = dragc_gas_dust

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
            if (has_interactions .and. .not. warned .and. master) call warn("[interactions:interactions_grace_passed] Cannot initialize aerodynamical drag because dust does not exist.")
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

#ifndef BALSARA
   function fluid_interactions_aero_drag(dens, velx) result(acc)

      use fluidindex,       only: flind

      implicit none

      real, dimension(:,:), intent(in), pointer  :: dens
      real, dimension(:,:), intent(in), pointer  :: velx
      real, dimension(size(dens,1),size(dens,2)) :: acc

      acc(:, flind%dst%pos) = -dragc_gas_dust * (velx(:, flind%dst%pos) - velx(:, flind%neu%pos))
      acc(:, flind%neu%pos) = -acc(:, flind%dst%pos) * dens(:, flind%dst%pos) / dens(:, flind%neu%pos)
!      acc(:, flind%neu%pos) = -dragc_gas_dust * dens(:, flind%dst%pos) / dens(:, flind%neu%pos) * ( velx(:, flind%neu%pos) - velx(:, flind%dst%pos) )

   end function fluid_interactions_aero_drag

   function fluid_interactions_aero_drag_ep(dens, velx) result(acc)

      use fluidindex,       only: flind

      implicit none

      real, dimension(:,:), intent(in), pointer  :: dens
      real, dimension(:,:), intent(in), pointer  :: velx
      real, dimension(size(dens,1),size(dens,2)) :: acc

      acc(:, flind%dst%pos) = -dens(:, flind%neu%pos) / epstein_factor(flind%neu%pos) * (velx(:, flind%dst%pos) - velx(:, flind%neu%pos))
      acc(:, flind%neu%pos) = -acc(:, flind%dst%pos) * dens(:, flind%dst%pos) / dens(:, flind%neu%pos)

   end function fluid_interactions_aero_drag_ep
#endif /* !BALSARA */

!>
!! \details
!! Balsara Dinshaw S., Tilley David A., Rettig Terrence, Brittain Sean D. MNRAS (2009) 397: 24.
!! Tilley, David A., Balsara, Dinshaw S. MNRAS (2008) 389: 1058.
!<
   subroutine balsara_implicit_interactions(u1, vx, istep, sweep, i1, i2, cg)

      use constants,        only: half, one, zero, I_ONE, I_TWO, LO, HI, cs_i2_n, uh_n, RK2_1, RK2_2
      use dataio_pub,       only: msg, warn, die
      use fluidindex,       only: flind, iarr_all_dn, iarr_all_mx, iarr_all_swp
      use global,           only: dt
      use grid_cont,        only: grid_container
      use mpisetup,         only: master
      use named_array_list, only: qna, wna

      implicit none

      real, dimension(:,:), intent(inout)        :: u1
      real, dimension(:,:), intent(in)           :: vx
      integer, intent(in)                        :: istep
      integer(kind=4),               intent(in)  :: sweep              !< direction (x, y or z) we are doing calculations for
      integer,                       intent(in)  :: i1                 !< coordinate of sweep in the 1st remaining direction
      integer,                       intent(in)  :: i2                 !< coordinate of sweep in the 2nd remaining direction
      type(grid_container), pointer, intent(in)  :: cg                 !< current grid piece

      real, dimension(:), pointer                :: cs_iso2
      real, dimension(size(u1, 1), flind%fluids) :: vprim
      real, dimension(size(u1, 1), flind%all)    :: u0
      real, dimension(size(u1, 1))               :: delta, drag
      integer, dimension(2)                      :: tfl
      character(len=*), dimension(3), parameter  :: flns = ['ION', 'NEU', 'DST'] !< \todo move this to fluids modules
      logical, save                              :: initbalsara = .true.
      !>
      !! \deprecated BEWARE: this bit assumes that we have 2 fluids and \f$u_1 \equiv u_0 - \nabla F\f$
      !! \todo 2) half-time step should be \f$\le \frac{1}{2}\frac{c_s}{drag * \rho'_? |v'_d - v'_g|}\f$
      !! \todo 3) what if not isothermal?
      !! \todo 4) remove hardcoded integers - done
      !<
      if (epstein_factor(flind%neu%pos) <= zero) return

      tfl = [I_ONE, I_TWO] !> may become changeable in future

      if (initbalsara .and. master .and. (flind%fluids > I_TWO)) then
         write(msg,"(a,i2,4a)")"[interactions:balsara_implicit_interactions] ", flind%fluids, " fluids present, yet only two included to compute interactions: ", flns(flind%all_fluids(tfl(LO))%fl%tag), " & ", flns(flind%all_fluids(tfl(HI))%fl%tag)
         call warn(msg)
         initbalsara = .false.
      endif

      cs_iso2 => cg%q(qna%ind(cs_i2_n))%get_sweep(sweep,i1,i2)
      drag(:) = dt * half / grain_dens_x_size * sqrt( cs_iso2(:) + abs( vx(:, tfl(LO)) - vx(:, tfl(HI)) )**2)

      delta(:) = one + drag(:) * (u1(:, iarr_all_dn(tfl(LO))) + u1(:, iarr_all_dn(tfl(HI))))
      delta(:) = one/delta(:)

      ! here istep is tightly bound to canonical RK2
      select case (istep)
      case (RK2_1)
         vprim(:, tfl(LO)) =  delta(:)*( (1./u1(:, iarr_all_dn(tfl(LO))) + drag(:))*u1(:, iarr_all_mx(tfl(LO))) + drag(:)*u1(:, iarr_all_mx(tfl(HI))) )
         vprim(:, tfl(HI)) =  delta(:)*( (1./u1(:, iarr_all_dn(tfl(HI))) + drag(:))*u1(:, iarr_all_mx(tfl(HI))) + drag(:)*u1(:, iarr_all_mx(tfl(LO))) )
      case (RK2_2)
         u0(:, iarr_all_swp(sweep,:)) = transpose(cg%w(wna%ind(uh_n))%get_sweep(sweep,i1,i2))
         vprim(:, tfl(HI)) =  delta(:)*( drag(:)*( u0(:, iarr_all_dn(tfl(HI)))/u1(:, iarr_all_dn(tfl(HI)))*u0(:, iarr_all_mx(tfl(LO))) &
                                                 - u0(:, iarr_all_dn(tfl(LO)))/u1(:, iarr_all_dn(tfl(HI)))*u0(:, iarr_all_mx(tfl(HI))) &
                                                 + u1(:, iarr_all_mx(tfl(LO))) ) + u1(:, iarr_all_mx(tfl(HI))) * ( 1./u1(:, iarr_all_dn(tfl(HI))) + drag(:) )  )
         vprim(:, tfl(LO)) =  delta(:)*( drag(:)*( u0(:, iarr_all_dn(tfl(LO)))/u1(:, iarr_all_dn(tfl(LO)))*u0(:, iarr_all_mx(tfl(HI))) &
                                                 - u0(:, iarr_all_dn(tfl(HI)))/u1(:, iarr_all_dn(tfl(LO)))*u0(:, iarr_all_mx(tfl(LO))) &
                                                 + u1(:, iarr_all_mx(tfl(HI))) ) + u1(:, iarr_all_mx(tfl(LO))) * ( 1./u1(:, iarr_all_dn(tfl(LO))) + drag(:) )  )
      case default
         call die("[interactions:balsara_implicit_interactions] Unsupported substep")
      end select
      u1(:, iarr_all_mx) = u1(:, iarr_all_dn) * vprim(:,:)
      return
   end subroutine balsara_implicit_interactions

#if 0
! this function is not used anywhere

!>
!! \brief balsara_implicit_interactions rewritten to function
!! \details
!! Balsara Dinshaw S., Tilley David A., Rettig Terrence, Brittain Sean D. MNRAS (2009) 397: 24.
!! Tilley, David A., Balsara, Dinshaw S. MNRAS (2008) 389: 1058.
!<
   function balsara_based_interactions_function(u1, u0, vx, istep, sweep, i1, i2, cg) result(usrc)

      use constants,        only: half, one, zero, I_ONE, I_TWO, LO, HI, cs_i2_n, RK2_1, RK2_2
      use dataio_pub,       only: msg, warn, die
      use fluidindex,       only: iarr_all_dn, iarr_all_mx, flind
      use global,           only: dt
      use grid_cont,        only: grid_container
      use mpisetup,         only: master
      use named_array_list, only: qna

      implicit none

      real, dimension(:,:),          intent(in) :: u1
      real, dimension(:,:),          intent(in) :: u0
      real, dimension(:,:),          intent(in) :: vx
      integer,                       intent(in) :: istep
      integer(kind=4),               intent(in) :: sweep              !< direction (x, y or z) we are doing calculations for
      integer,                       intent(in) :: i1                 !< coordinate of sweep in the 1st remaining direction
      integer,                       intent(in) :: i2                 !< coordinate of sweep in the 2nd remaining direction
      type(grid_container), pointer, intent(in) :: cg                 !< current grid piece
      real, dimension(size(u1,1),size(u1,2))    :: usrc

      real, dimension(:), pointer                :: cs_iso2
      real, dimension(size(u1, 1))               :: delta, drag
      integer, dimension(2)                      :: tfl, inv
      integer                                    :: i
      character(len=*), dimension(3), parameter  :: flns = ['ION', 'NEU', 'DST'] !< \todo move this to fluids modules
      logical, save                              :: initbalsara = .true.
      !>
      !! \deprecated BEWARE: this bit assumes that we have 2 fluids and \f$u_1 \equiv u_0 - \nabla F\f$
      !! \todo 2) half-time step should be \f$\le \frac{1}{2}\frac{c_s}{drag * \rho'_? |v'_d - v'_g|}\f$
      !! \todo 3) what if not isothermal?
      !! \note BEWARE: it is function to replace balsara_implicit_interactions, yet it requires extensive testing
      !<
      usrc = 0.0
      if (epstein_factor(flind%neu%pos) <= zero) return

      tfl = [I_ONE, I_TWO] !> may become changeable in future
      inv = [tfl(HI), tfl(LO)]

      if (initbalsara .and. master .and. (flind%fluids > I_TWO)) then
         write(msg,"(a,i2,4a)")"[interactions:balsara_implicit_interactions] ", flind%fluids, " fluids present, yet only two included to compute interactions: ", flns(flind%all_fluids(tfl(LO))%fl%tag), " & ", flns(flind%all_fluids(tfl(HI))%fl%tag)
         call warn(msg)
         initbalsara = .false.
      endif

      cs_iso2 => cg%q(qna%ind(cs_i2_n))%get_sweep(sweep,i1,i2)
      drag(:) = half / grain_dens_x_size * sqrt( cs_iso2(:) + abs( vx(:, tfl(LO)) - vx(:, tfl(HI)) )**2)

      delta(:) = one + dt * drag(:) * (u1(:, iarr_all_dn(tfl(LO))) + u1(:, iarr_all_dn(tfl(HI))))
      delta(:) = one/delta(:)

      do i = I_ONE, I_TWO
         ! here istep is tightly bound to canonical RK2
         select case (istep)
         case (RK2_1)
            usrc(:,iarr_all_mx(tfl(i))) = drag(:) * ( u1(:, iarr_all_dn(tfl(i)))*u1(:, iarr_all_mx(inv(i))) - delta(:)*u1(:, iarr_all_dn(inv(i)))*u1(:, iarr_all_mx(tfl(i))) )
         case (RK2_2)
            usrc(:,iarr_all_mx(tfl(i))) = drag(:) * delta(:) * ( u0(:, iarr_all_dn(tfl(i)))*u0(:, iarr_all_mx(inv(i))) - u0(:, iarr_all_dn(inv(i)))*u0(:, iarr_all_mx(tfl(i))) &
                                                               + u1(:, iarr_all_dn(tfl(i)))*u1(:, iarr_all_mx(inv(i))) - u1(:, iarr_all_dn(inv(i)))*u1(:, iarr_all_mx(tfl(i))) )
         case default
            call die("[interactions:balsara_based_interactions_function] Unsupported substep")
         end select
      enddo
      ! usrc should be used to update u1 in internal_sources procedure: u1 = u1 + rk2factror*usrc*dt
      return
   end function balsara_based_interactions_function

#endif /* 0 */

   subroutine update_grain_size(new_size)
      use fluidindex, only: flind
      use constants,  only: DST
      use units,      only: gram, cm
      implicit none
      real, intent(in) :: new_size
      integer :: i

      grain_size = new_size
      grain_dens_x_size = grain_size * grain_dens * gram / cm**2

      epstein_factor(:) = 0.0
      do i = lbound(flind%all_fluids, 1), ubound(flind%all_fluids, 1)
         if (flind%all_fluids(i)%fl%tag /= DST) &
            epstein_factor(i) = grain_dens_x_size / flind%all_fluids(i)%fl%cs !BEWARE iso assumed
      enddo

   end subroutine update_grain_size

#ifndef BALSARA
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
            if (ifl /= jfl) acc(:, ifl) = acc(:, ifl) + collfaq(ifl,jfl) * dens(:, ifl) * dens(:, jfl) * (velx(:, jfl) - velx(:, ifl)) ! / dens(:, ifl) ?
         enddo
      enddo

   end function dragforce
#endif /* !BALSARA */

!>
!! \brief interface between sources and fluid_interactions
!! \todo do it better
!<
   function fluid_interactions_exec(nn,uu,velx) result(acc)

      use fluidindex, only: flind, iarr_all_dn

      implicit none

      integer(kind=4),                intent(in) :: nn                  !< array size
      real, dimension(nn, flind%all), intent(in) :: uu                  !< vector of conservative variables
      real, dimension(:,:), pointer,  intent(in) :: velx
      real, dimension(size(velx,1),size(velx,2)) :: acc
!locals
      real, dimension(nn, flind%fluids), target  :: density            !< gas density
      real, dimension(:,:),            pointer   :: dens

      dens => density
      density(:,:) = uu(:, iarr_all_dn)

      acc(:,:) = fluid_interactions(dens,velx)

   end function fluid_interactions_exec

end module interactions
