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
!! \brief Initialization of the neutral fluid
!!
!!
!!
!! In this module following namelist of parameters is specified:
!! \copydetails initneutral::init_neutral
!! \deprecated this module should not export any variables
!<

module initneutral
! pulled by ANY
   use fluidtypes, only: component_fluid
   implicit none

   public :: init_neutral, neutral_index, cleanup_neutral, neutral_fluid, &
      gamma_neu, cs_iso_neu, cs_iso_neu2, selfgrav_neu, idnn, imxn, imyn, imzn, ienn

   real                  :: gamma_neu             !< adiabatic index for the neutral gas component
   real                  :: cs_iso_neu            !< isothermal sound speed (p = cs_iso_neu<sup>2</sup>\f$\rho\f$), active only if neutral gas is \ref isothermal
   real                  :: cs_iso_neu2
   logical               :: selfgrav_neu          !< true if neutral gas is selfgravitating
   integer(kind=4)       :: idnn, imxn, imyn, imzn, ienn

   type, extends(component_fluid) :: neutral_fluid
      contains
         procedure, nopass :: get_tag
         procedure, pass :: get_cs => neu_cs
         procedure, nopass :: compute_flux => flux_neu
   end type neutral_fluid

contains

   real function neu_cs(this, cg, i, j, k)
      use grid_cont, only: grid_container
#ifndef ISO
      use func,      only: ekin
#endif /* !ISO */
      implicit none
      class(neutral_fluid), intent(in) :: this
      type(grid_container), pointer, intent(in) :: cg !< current grid container
      integer, intent(in) :: i, j, k

      real :: p
#ifdef ISO
      p  = cg%cs_iso2(i, j, k) * cg%u(this%idn, i, j, k)
      neu_cs = sqrt(cg%cs_iso2(i, j, k))
#else /* !ISO */
      p  = (cg%u(this%ien, i, j, k) - &
         &   ekin(cg%u(this%imx, i, j, k), cg%u(this%imy, i, j, k), cg%u(this%imz, i, j, k), cg%u(this%idn, i, j, k)) &
         & ) * this%gam_1
      neu_cs = sqrt(abs((this%gam * p) / cg%u(this%idn, i, j, k)))
#endif /* !ISO */
   end function neu_cs

   function get_tag() result(tag)
      use constants, only: idlen
      implicit none
      character(len=idlen)   :: tag

      tag = "NEU"
   end function get_tag

!>
!! \brief Routine to set parameters values from namelist FLUID_NEUTRAL
!!
!! \n \n
!! @b FLUID_NEUTRAL
!! \n \n
!! <table border="+1">
!! <tr><td width="150pt"><b>parameter</b></td><td width="135pt"><b>default value</b></td><td width="200pt"><b>possible values</b></td><td width="315pt"> <b>description</b></td></tr>
!! <tr><td>gamma_neu     </td><td>5./3.   </td><td>real value </td><td>\copydoc initneutral::gamma_neu    </td></tr>
!! <tr><td>cs_iso_neu    </td><td>1.0     </td><td>real value </td><td>\copydoc initneutral::cs_iso_neu   </td></tr>
!! <tr><td>selfgrav_neu  </td><td>.false. </td><td>logical    </td><td>\copydoc initneutral::selfgrav_neu </td></tr>
!! </table>
!! \n \n
!<
   subroutine init_neutral

      use dataio_pub,    only: par_file, ierrh, namelist_errh, compare_namelist, cmdl_nml, lun ! QA_WARN required for diff_nml
      use mpisetup,      only: rbuff, lbuff, comm, mpi_err, buffer_dim, master, slave, FIRST
      use mpi,           only: MPI_DOUBLE_PRECISION, MPI_LOGICAL

      implicit none

      namelist /FLUID_NEUTRAL/ gamma_neu, cs_iso_neu, selfgrav_neu

      gamma_neu    = 5./3.
      cs_iso_neu   = 1.0
      selfgrav_neu = .false.

      if (master) then

         diff_nml(FLUID_NEUTRAL)

         lbuff(1)   = selfgrav_neu

         rbuff(1)   = gamma_neu
         rbuff(2)   = cs_iso_neu

      endif

      call MPI_Bcast(rbuff,    buffer_dim, MPI_DOUBLE_PRECISION, FIRST, comm, mpi_err)
      call MPI_Bcast(lbuff,    buffer_dim, MPI_LOGICAL,          FIRST, comm, mpi_err)

      if (slave) then

         selfgrav_neu = lbuff(1)

         gamma_neu   = rbuff(1)
         cs_iso_neu  = rbuff(2)

      endif

      cs_iso_neu2      = cs_iso_neu**2

   end subroutine init_neutral

   subroutine neutral_index(flind)
      use constants,    only: NEU, xdim, ydim, zdim, ndims, I_ONE, I_TWO, I_THREE, I_FOUR
      use diagnostics,  only: ma1d, ma2d, my_allocate
      use fluidtypes,   only: var_numbers

      implicit none

      type(var_numbers), intent(inout) :: flind

      flind%neu%beg    = flind%all + I_ONE

      idnn = flind%all + I_ONE
      imxn = flind%all + I_TWO
      imyn = flind%all + I_THREE
      imzn = flind%all + I_FOUR

      flind%neu%idn = idnn
      flind%neu%imx = imxn
      flind%neu%imy = imyn
      flind%neu%imz = imzn

      flind%neu%all  = 4
      flind%all      = imzn
#ifndef ISO
      ienn          = imzn + I_ONE
      flind%all      = flind%all + I_ONE
      flind%neu%all  = flind%neu%all +I_ONE
      flind%neu%ien  = ienn
#endif /* !ISO */

      ma1d = [flind%neu%all]
      call my_allocate(flind%neu%iarr,       ma1d)
      ma2d = [ndims, flind%neu%all]
      call my_allocate(flind%neu%iarr_swp,   ma2d)

      !\deprecated repeated magic integers
      flind%neu%iarr(1:4)           = [idnn,imxn,imyn,imzn]
      flind%neu%iarr_swp(xdim, 1:4) = [idnn,imxn,imyn,imzn]
      flind%neu%iarr_swp(ydim, 1:4) = [idnn,imyn,imxn,imzn]
      flind%neu%iarr_swp(zdim, 1:4) = [idnn,imzn,imyn,imxn]

#ifndef ISO
      flind%neu%iarr(5)       = ienn
      flind%neu%iarr_swp(:,5) = ienn
      flind%neu%has_energy    = .true.

      flind%energ = flind%energ + I_ONE
#endif /* !ISO */

      flind%neu%end    = flind%all
      flind%components = flind%components + I_ONE
      flind%fluids     = flind%fluids + I_ONE
      flind%neu%pos    = flind%components
      if (selfgrav_neu)  flind%fluids_sg = flind%fluids_sg + I_ONE

      flind%neu%gam   = gamma_neu
      flind%neu%gam_1 = gamma_neu-1.0
      flind%neu%cs    = cs_iso_neu
      flind%neu%cs2   = cs_iso_neu**2
      flind%neu%tag   = NEU

      flind%neu%is_selfgrav   = selfgrav_neu
      flind%neu%is_magnetized = .false.
#ifndef ISO
      flind%neu%has_energy    = .true.
#endif /* !ISO */

   end subroutine neutral_index

   subroutine cleanup_neutral

      implicit none

   end subroutine cleanup_neutral

!==========================================================================================
!
! OPT: This routine may cost as much as 30% of rtvd. It seems that all the data fit well a 512kB L2 cache, but Ir:Dr:Dw is like 8:2:1
! OPT: \todo Try an explicit loop over RNG to check if we're better than the compiler
! OPT: similar treatment may be helpful for fluxionized.F90, fluxdust.F90 and fluxcosmicrays.F90
!
#define RNG 2:nm
!/*
!>
!! \brief Computation of %fluxes for the neutral fluid
!!
!!The flux functions for neutral fluid are given by
!!
!!\f[
!!  \vec{F}{(\vec{u})} =
!!  \left(\begin{array}{c}
!!    \rho v_x \\
!!    \rho v_x^2 + p \\
!!    \rho v_x v_y\\
!!    \rho v_x v_z\\
!!    (e + p)v_x
!!  \end{array}\right),
!!  \qquad
!!  \vec{G}{(\vec{u})} =
!!  \left(\begin{array}{c}
!!    \rho v_y \\
!!    \rho v_y v_x\\
!!    \rho v_y^2 + p\\
!!    \rho v_y v_z\\
!!    (e + p)v_y
!!  \end{array}\right),
!!\qquad
!!  \vec{H}{(\vec{u})} =
!!  \left(\begin{array}{c}
!!    \rho v_z \\
!!    \rho v_z v_x\\
!!    \rho v_z v_y \\
!!    \rho v_z^2 + p \\
!!    (e + p)v_z
!!  \end{array}\right),
!!\f]
!<
!*/
   subroutine flux_neu(flux, cfr, uu, n, vx, ps, bb, cs_iso2)

      use func,       only: ekin
      use fluidindex, only: idn, imx, imy, imz
#if defined(LOCAL_FR_SPEED) || !defined(ISO)
      use fluidindex, only: flind
#endif /* defined(LOCAL_FR_SPEED) || !defined(ISO) */
#ifdef LOCAL_FR_SPEED
      use constants,  only: small, half
      use global,     only: cfr_smooth
#endif /* LOCAL_FR_SPEED */
#ifdef GLOBAL_FR_SPEED
      use timestep,   only: c_all
#endif /* GLOBAL_FR_SPEED */
#ifndef ISO
      use fluidindex, only: ien
      use global,     only: smallp
#endif /* !ISO */

      implicit none

      integer(kind=4), intent(in)                  :: n         !< number of cells in the current sweep
      real, dimension(:,:), intent(inout), pointer :: flux      !< flux of neutral fluid
      real, dimension(:,:), intent(inout), pointer :: cfr       !< freezing speed for neutral fluid
      real, dimension(:,:), intent(in),    pointer :: uu        !< part of u for neutral fluid
      real, dimension(:),   intent(inout), pointer :: vx        !< velocity of neutral fluid for current sweep
      real, dimension(:),   intent(inout), pointer :: ps        !< pressure of neutral fluid for current sweep
      real, dimension(:,:), intent(in),    pointer :: bb        !< magnetic field x,y,z-components table
      real, dimension(:),   intent(in),    pointer :: cs_iso2   !< isothermal sound speed squared

      ! locals
      integer            :: nm
#ifdef LOCAL_FR_SPEED
      integer            :: i
      real               :: minvx     !<
      real               :: maxvx     !<
      real               :: amp       !<
#endif /* LOCAL_FR_SPEED */

      nm = n-1
      vx(RNG) = uu(imx,RNG)/uu(idn,RNG) ; vx(1) = vx(2); vx(n) = vx(nm)
#ifdef ISO
      ps(RNG)  = cs_iso2(RNG)*uu(idn,RNG) ; ps(1) = ps(2); ps(n) = ps(nm)
#else /* !ISO */
      ps(RNG)  = (uu(ien,RNG) - ekin(uu(imx,RNG),uu(imy,RNG),uu(imz,RNG),uu(idn,RNG)) )*(flind%neu%gam_1)
      ps(RNG)  = max(ps(RNG), smallp)
#endif /* !ISO */

      flux(idn,RNG)=uu(imx,RNG)
      flux(imx,RNG)=uu(imx,RNG)*vx(RNG)+ps(RNG)
      flux(imy,RNG)=uu(imy,RNG)*vx(RNG)
      flux(imz,RNG)=uu(imz,RNG)*vx(RNG)
#ifndef ISO
      flux(ien,RNG)=(uu(ien,RNG)+ps(RNG))*vx(RNG)
#endif /* !ISO */
      flux(:,1) = flux(:,2) ; flux(:,n) = flux(:,nm)

#ifdef LOCAL_FR_SPEED

      ! The freezing speed is now computed locally (in each cell)
      !  as in Trac & Pen (2003). This ensures much sharper shocks,
      !  but sometimes may lead to numerical instabilities
      minvx = minval(vx(RNG))
      maxvx = maxval(vx(RNG))
      amp   = half*(maxvx-minvx)
      !    c_fr  = 0.0
#ifdef ISO
      cfr(1,RNG) = sqrt(vx(RNG)**2+cfr_smooth*amp) + max(sqrt( abs(              ps(RNG))/uu(idn,RNG)),small)
#else /* !ISO */
      cfr(1,RNG) = sqrt(vx(RNG)**2+cfr_smooth*amp) + max(sqrt( abs(flind%neu%gam*ps(RNG))/uu(idn,RNG)),small)
#endif /* !ISO */
      !> \deprecated BEWARE: that is the cause of fast decreasing of timestep in galactic disk problem
      !>
      !! \todo find why is it so
      !! if such a treatment is OK then should be applied also in both cases of neutral and ionized gas
      !!    do i = 2,nm
      !!       cfr(1,i) = maxval( [c_fr(i-1), c_fr(i), c_fr(i+1)] )
      !!    enddo
      !<

      cfr(1,1) = cfr(1,2);  cfr(1,n) = cfr(1,nm)
      do i = 2, flind%neu%all
         cfr(i,:) = cfr(1,:)
      enddo
#endif /* LOCAL_FR_SPEED */

#ifdef GLOBAL_FR_SPEED
      ! The freezing speed is now computed globally
      !  (c=const for the whole domain) in subroutine 'timestep'

      !    cfr(:,:) = flind%neu%c   ! check which c_xxx is better
      cfr(:,:) = c_all
#endif /* GLOBAL_FR_SPEED */
      return
      if (.false.) write(0,*) bb, cs_iso2

   end subroutine flux_neu

end module initneutral
