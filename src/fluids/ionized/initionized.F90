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
!! \brief Initialization of the ionized fluid
!!
!!
!!
!! In this module following namelist of parameters is specified:
!! \copydetails initionized::init_ionized
!! \deprecated This module should not export any variables
!<

module initionized
! pulled by ANY
   use fluidtypes, only: component_fluid
   implicit none

   private
   public :: init_ionized, cleanup_ionized, ion_fluid

   real    :: gamma       !< adiabatic index for the ionized gas component
   real    :: cs_iso      !< isothermal sound speed (p = cs_iso<sup>2</sup>\f$\rho\f$), active only if ionized gas is \ref isothermal
   logical :: selfgrav    !< true if ionized gas is selfgravitating

   type, extends(component_fluid) :: ion_fluid
      contains
         procedure, nopass :: get_tag
         procedure, pass   :: get_cs => ion_cs
         procedure, pass   :: compute_flux => flux_ion
         procedure, pass   :: initialize_indices => initialize_ion_indices
   end type ion_fluid

contains

   subroutine initialize_ion_indices(this, flind)
      use constants,  only: ION
      use fluidtypes, only: var_numbers
      implicit none
      class(ion_fluid),     intent(inout) :: this
      type(var_numbers),    intent(inout) :: flind

      logical :: has_energy, is_magnetized
#ifdef ISO
      has_energy = .false.
#else /* !ISO */
      has_energy = .true.
#endif /* !ISO */
#ifdef MAGNETIZED
      is_magnetized = .true.
#else /* !MAGNETIZED */
      is_magnetized = .true.
#endif /* !MAGNETIZED */

      call this%set_fluid_index(flind, is_magnetized, selfgrav, has_energy, cs_iso, gamma, ION)

   end subroutine initialize_ion_indices

   real function ion_cs(this, i, j, k, u, b, cs_iso2)
      use constants, only: two
#ifndef ISO
      use func,      only: ekin
#endif /* !ISO */
#ifdef MAGNETIC
      use constants, only: xdim, ydim, zdim, half
      use domain,    only: dom
      use func,      only: emag
#else /* !MAGNETIC */
      use constants, only: zero
#endif /* !MAGNETIC */

      implicit none
      class(ion_fluid),                  intent(in) :: this
      integer,                           intent(in) :: i, j, k
      real, dimension(:,:,:,:), pointer, intent(in) :: u       !< pointer to array of fluid properties
      real, dimension(:,:,:,:), pointer, intent(in) :: b       !< pointer to array of magnetic fields (used for ionized fluid with MAGNETIC #defined)
      real, dimension(:,:,:),   pointer, intent(in) :: cs_iso2 !< pointer to array of isothermal sound speeds (used when ISO was #defined)

#ifdef MAGNETIC
      real :: bx, by, bz
#endif /* MAGNETIC */
      real :: pmag, p, ps

#ifdef MAGNETIC
      bx = half*(b(xdim,i,j,k) + b(xdim, i+dom%D_x, j,         k        ))
      by = half*(b(ydim,i,j,k) + b(ydim, i,         j+dom%D_y, k        ))
      bz = half*(b(zdim,i,j,k) + b(zdim, i,         j,         k+dom%D_z))

      pmag = emag(bx, by, bz)
#else /* !MAGNETIC */
      ! all_mag_boundaries has not been called so we cannot trust b(xdim, ie+dom%D_x:), b(ydim,:je+dom%D_y and b(zdim,:,:, ke+dom%D_z
      pmag = zero
#endif /* !MAGNETIC */

#ifdef ISO
      p  = cs_iso2(i, j, k) * u(this%idn, i, j, k)
      ps = p + pmag
      ion_cs = sqrt(abs((two * pmag + p) / u(this%idn, i, j, k)))
#else /* !ISO */
      ps = (u(this%ien, i, j, k) - &
         &   ekin(u(this%imx, i, j, k), u(this%imy, i, j, k), u(this%imz, i, j, k), u(this%idn, i, j, k)) &
         & ) * (this%gam_1) + (two - this%gam) * pmag
      p  = ps - pmag
      ion_cs = sqrt(abs((two * pmag + this%gam * p) / u(this%idn, i, j, k)))
#endif /* !ISO */
      if (.false.) print *, u(:, i, j, k), b(:, i, j, k), cs_iso2(i, j, k), this%cs
   end function ion_cs

   function get_tag() result(tag)
      use constants, only: idlen
      implicit none
      character(len=idlen) :: tag

      tag = "ION"
   end function get_tag

!>
!! \brief Routine to set parameters from namelist FLUID_IONIZED
!!
!! \n \n
!! @b FLUID_IONIZED
!! \n \n
!! <table border="+1">
!! <tr><td width="150pt"><b>parameter</b></td><td width="135pt"><b>default value</b></td><td width="200pt"><b>possible values</b></td><td width="315pt"> <b>description</b></td></tr>
!! <tr><td>gamma     </td><td>5./3.   </td><td>real value </td><td>\copydoc initionized::gamma    </td></tr>
!! <tr><td>cs_iso    </td><td>1.0     </td><td>real value </td><td>\copydoc initionized::cs_iso   </td></tr>
!! <tr><td>selfgrav  </td><td>.false. </td><td>logical    </td><td>\copydoc initionized::selfgrav </td></tr>
!! </table>
!! The list is active while \b "IONIZED" is defined.
!! \n \n
!<
   subroutine init_ionized

      use dataio_pub, only: nh ! QA_WARN required for diff_nml
      use mpisetup,   only: rbuff, lbuff, master, slave, piernik_MPI_Bcast

      implicit none

      namelist /FLUID_IONIZED/ gamma, cs_iso, selfgrav

      gamma    = 5./3.
      cs_iso   = 1.0
      selfgrav = .false.

      if (master) then

         if (.not.nh%initialized) call nh%init()
         open(newunit=nh%lun, file=nh%tmp1, status="unknown")
         write(nh%lun,nml=FLUID_IONIZED)
         close(nh%lun)
         open(newunit=nh%lun, file=nh%par_file)
         nh%errstr=""
         read(unit=nh%lun, nml=FLUID_IONIZED, iostat=nh%ierrh, iomsg=nh%errstr)
         close(nh%lun)
         call nh%namelist_errh(nh%ierrh, "FLUID_IONIZED")
         read(nh%cmdl_nml,nml=FLUID_IONIZED, iostat=nh%ierrh)
         call nh%namelist_errh(nh%ierrh, "FLUID_IONIZED", .true.)
         open(newunit=nh%lun, file=nh%tmp2, status="unknown")
         write(nh%lun,nml=FLUID_IONIZED)
         close(nh%lun)
         call nh%compare_namelist()

         lbuff(1)  = selfgrav

         rbuff(1)  = gamma
         rbuff(2)  = cs_iso

      endif

      call piernik_MPI_Bcast(rbuff)
      call piernik_MPI_Bcast(lbuff)

      if (slave) then

         selfgrav = lbuff(1)

         gamma    = rbuff(1)
         cs_iso   = rbuff(2)

      endif

   end subroutine init_ionized

   subroutine cleanup_ionized

      implicit none

   end subroutine cleanup_ionized

!/*
!>
!! \brief Computation of %fluxes for the ionized fluid
!!
!!The flux functions for ionized fluid are given by
!!
!!\f[
!!  \vec{F}{(\vec{u})} =
!!  \left(\begin{array}{c}
!!    \rho v_x \\
!!    \rho v_x^2 + p_* - B_x^2 \\
!!    \rho v_x v_y - B_x B_y\\
!!    \rho v_x v_z - B_x B_z\\
!!    (e + p)v_x - \vec{B} \cdot \vec{v} \; B_x
!!  \end{array}\right),
!!  \qquad
!!  \vec{G}{(\vec{u})} =
!!  \left(\begin{array}{c}
!!    \rho v_y \\
!!    \rho v_y v_x  - B_y B_x\\
!!    \rho v_y^2 + p_* - B_y^2 \\
!!    \rho v_y v_z  - B_y B_z\\
!!    (e + p)v_y - \vec{B} \cdot \vec{v} \; B_y
!!  \end{array}\right),
!!\qquad
!!  \vec{H}{(\vec{u})} =
!!  \left(\begin{array}{c}
!!    \rho v_z \\
!!    \rho v_z v_x  - B_z B_x\\
!!    \rho v_z v_y  - B_z B_y\\
!!    \rho v_z^2 + p_* - B_z^2 \\
!!    (e + p)v_z - \vec{B} \cdot \vec{v} \; B_z
!!  \end{array}\right),
!!\f]
!!where \f$p_* = p + B^2/2\f$, \f$e= e_{th} + \frac{1}{2} \rho v^2 + B^2/2\f$,  are the total pressure and total energy density,
!!while \f$e_{th}\f$ is thermal energy density and  \f$e_{mag} = B^2/2\f$ is the magnetic energy density.
!<
!*/
#define RNG 2:nm
   subroutine flux_ion(this, flux, cfr, uu, n, vx, ps, bb, cs_iso2, use_vx)

      use constants,    only: xdim, ydim, zdim, idn, imx, imy, imz
#ifndef ISO
      use constants,    only: ien
#endif /* !ISO */
      use dataio_pub,   only: die
      use func,         only: ekin, emag
#ifdef LOCAL_FR_SPEED
      use constants,    only: small, half
      use global,       only: cfr_smooth
#endif /* LOCAL_FR_SPEED */

      implicit none
      class(ion_fluid),     intent(in)             :: this
      integer(kind=4),      intent(in)             :: n         !< number of cells in the current sweep
      real, dimension(:,:), intent(inout), pointer :: flux      !< flux of ionized fluid
      real, dimension(:,:), intent(inout), pointer :: cfr       !< freezing speed for ionized fluid
      real, dimension(:,:), intent(in),    pointer :: uu        !< part of u for ionized fluid
      real, dimension(:),   intent(inout), pointer :: vx        !< velocity of ionized fluid for current sweep
      real, dimension(:),   intent(inout), pointer :: ps        !< pressure of ionized fluid for current sweep
      real, dimension(:,:), intent(in),    pointer :: bb        !< magnetic field x,y,z-components table
      real, dimension(:),   intent(in),    pointer :: cs_iso2   !< local isothermal sound speed squared (optional)
      logical,              intent(in)             :: use_vx    !< use provided vx instead of computing it

      ! locals
      real, dimension(n) :: p           !< thermal pressure of ionized fluid
      real, dimension(n) :: pmag        !< pressure of magnetic field
      integer            :: nm
#ifdef LOCAL_FR_SPEED
      integer            :: i
      real               :: minvx     !<
      real               :: maxvx     !<
      real               :: amp       !<
#endif /* LOCAL_FR_SPEED */

      nm = n-1
#ifdef MAGNETIC
      pmag(RNG)= emag(bb(RNG, xdim), bb(RNG, ydim), bb(RNG, zdim));  pmag(1) = pmag(2); pmag(n) = pmag(nm)
#else /* !MAGNETIC */
      pmag(:) = 0.0
#endif /* !MAGNETIC */
      if (.not. use_vx) then
         vx(RNG) = uu(RNG, imx) / uu(RNG, idn)
         vx(1) = vx(2)
         vx(n) = vx(nm)
      endif

#ifndef ISO
      if (associated(cs_iso2)) call die("[initionized:flux_ion] cs_iso2 should not be present")
#endif /* !ISO */

#ifdef ISO
      p(RNG)  = cs_iso2(RNG) * uu(RNG, idn)
      ps(RNG) = p(RNG) + pmag(RNG)
#else /* !ISO */
      ps(RNG) = (uu(RNG, ien) - ekin(uu(RNG, imx),uu(RNG, imy),uu(RNG, imz),uu(RNG, idn)) )*(this%gam_1) &
           & + (2.0 - this%gam)*pmag(RNG)
      p(RNG) = ps(RNG)- pmag(RNG);  p(1) = p(2); p(n) = p(nm)
#endif /* !ISO */
      ps(1) = ps(2); ps(n) = ps(nm)

      flux(RNG, idn)=uu(RNG, idn)*vx(RNG)
      flux(RNG, imx)=uu(RNG, imx)*vx(RNG)+ps(RNG) - bb(RNG, xdim)**2
      flux(RNG, imy)=uu(RNG, imy)*vx(RNG)-bb(RNG, ydim)*bb(RNG, xdim)
      flux(RNG, imz)=uu(RNG, imz)*vx(RNG)-bb(RNG, zdim)*bb(RNG, xdim)
#ifndef ISO
      flux(RNG, ien)=(uu(RNG, ien)+ps(RNG))*vx(RNG)-bb(RNG, xdim)*(bb(RNG, xdim)*uu(RNG, imx) &
                +bb(RNG, ydim)*uu(RNG, imy)+bb(RNG, zdim)*uu(RNG, imz))/uu(RNG, idn)
#endif /* !ISO */
      flux(1, :) = flux(2, :) ; flux(n, :) = flux(nm, :)

#ifdef LOCAL_FR_SPEED

      ! The freezing speed is now computed locally (in each cell)
      !  as in Trac & Pen (2003). This ensures much sharper shocks,
      !  but sometimes may lead to numerical instabilities
      minvx = minval(vx(RNG))
      maxvx = maxval(vx(RNG))
      amp   = half*(maxvx-minvx)
#ifdef ISO
      cfr(RNG, 1) = sqrt(vx(RNG)**2+cfr_smooth*amp) + max(sqrt( abs(2.0*pmag(RNG) +          p(RNG))/uu(RNG, idn)),small)
#else /* !ISO */
      cfr(RNG, 1) = sqrt(vx(RNG)**2+cfr_smooth*amp) + max(sqrt( abs(2.0*pmag(RNG) + this%gam*p(RNG))/uu(RNG, idn)),small)
#endif /* !ISO */
      !> \deprecated BEWARE: that is the cause of fast decreasing of timestep in galactic disk problem
      !>
      !! \todo find why is it so
      !! if such a treatment is OK then should be applied also in both cases of neutral and ionized gas
      !!    do i = 2,nm
      !!       cfr(1,i) = maxval( [c_fr(i-1), c_fr(i), c_fr(i+1)] )
      !!    enddo
      !<

      cfr(1,1) = cfr(2,1);  cfr(n, 1) = cfr(nm, 1)
      do i = 2, this%all
         cfr(:, i) = cfr(:, 1)
      enddo
#endif /* LOCAL_FR_SPEED */

#ifdef GLOBAL_FR_SPEED
      ! The freezing speed is now computed globally
      !  (c=const for the whole domain) in subroutine 'timestep'

      cfr(:,:) = this%c
#endif /* GLOBAL_FR_SPEED */

   end subroutine flux_ion

end module initionized
