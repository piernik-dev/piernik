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
      procedure, pass   :: get_mach => ion_mach
      procedure, pass   :: compute_flux => flux_ion
      procedure, pass   :: compute_pres => pres_ion
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
#ifdef MAGNETIC
      is_magnetized = .true.
#else /* !MAGNETIC */
      is_magnetized = .false.
#endif /* !MAGNETIC */

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
      use global,    only: cc_mag
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
      if (cc_mag) then
         pmag = emag(b(xdim,i,j,k), b(ydim,i,j,k), b(zdim,i,j,k))
      else
         bx = half*(b(xdim,i,j,k) + b(xdim, i+dom%D_x, j,         k        ))
         by = half*(b(ydim,i,j,k) + b(ydim, i,         j+dom%D_y, k        ))
         bz = half*(b(zdim,i,j,k) + b(zdim, i,         j,         k+dom%D_z))

         pmag = emag(bx, by, bz)
      endif
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

!>
!! \brief An estimate of (fast) Mach number based on upper estimate of fast magnetosonic speed
!! This is a bit complicated topic as fast magnetosonic speed depends locally on magnitude and orientation of magnetic field.
!! There are other characteristic speeds as well and each can be associated with respective Mach number.
!!
!! As this is the simplest approach, we use the same code for neutral fluid too.
!<

   real function ion_mach(this, i, j, k, u, b, cs_iso2)

      use func, only: sq_sum3

      implicit none

      class(ion_fluid),                  intent(in) :: this
      integer,                           intent(in) :: i, j, k
      real, dimension(:,:,:,:), pointer, intent(in) :: u       !< pointer to array of fluid properties
      real, dimension(:,:,:,:), pointer, intent(in) :: b       !< pointer to array of magnetic fields (used for ionized fluid with MAGNETIC #defined)
      real, dimension(:,:,:),   pointer, intent(in) :: cs_iso2 !< pointer to array of isothermal sound speeds (used when ISO was #defined)

      ion_mach = sqrt(sq_sum3(u(this%imx, i, j, k), u(this%imy, i, j, k), u(this%imz, i, j, k)))/u(this%idn, i, j, k) / this%get_cs(i, j, k, u, b, cs_iso2)

   end function ion_mach

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

      use bcast,      only: piernik_MPI_Bcast
      use dataio_pub, only: nh
      use mpisetup,   only: rbuff, lbuff, master, slave

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

         lbuff(1) = selfgrav

         rbuff(1) = gamma
         rbuff(2) = cs_iso

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
#define RNG2 2:nm
   subroutine flux_ion(this, flux, cfr, uu, n, vx, bb, cs_iso2)

      use constants, only: idn, imx, imy, imz
#ifndef ISO
      use constants, only: ien
#endif /* !ISO */
#ifdef MAGNETIC
      use constants, only: xdim, ydim, zdim
#endif /* MAGNETIC */
#ifdef LOCAL_FR_SPEED
      use constants, only: small, half
      use global,    only: cfr_smooth
#endif /* LOCAL_FR_SPEED */

      implicit none

      class(ion_fluid),     intent(in)             :: this
      integer(kind=4),      intent(in)             :: n         !< number of cells in the current sweep
      real, dimension(:,:), intent(inout), pointer :: flux      !< flux of ionized fluid
      real, dimension(:,:), intent(inout), pointer :: cfr       !< freezing speed for ionized fluid
      real, dimension(:,:), intent(in),    pointer :: uu        !< part of u for ionized fluid
      real, dimension(:),   intent(in),    pointer :: vx        !< velocity of ionized fluid for current sweep
      real, dimension(:,:), intent(in),    pointer :: bb        !< magnetic field x,y,z-components table
      real, dimension(:),   intent(in),    pointer :: cs_iso2   !< local isothermal sound speed squared (optional)

      ! locals
      real, dimension(n), target  :: ps        !< pressure of ionized fluid for current sweep
      real, dimension(:), pointer :: pps
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
      pps => ps
      call all_pres_ion(uu, n, bb, cs_iso2, this%gam, this%gam_1, pmag, p, pps)

      flux(RNG2, idn) = uu(RNG2, idn) * vx(RNG2)
#ifdef MAGNETIC
      flux(RNG2, imx) = uu(RNG2, imx) * vx(RNG2) + ps(RNG2) - bb(RNG2, xdim)**2
      flux(RNG2, imy) = uu(RNG2, imy) * vx(RNG2) - bb(RNG2, ydim) * bb(RNG2, xdim)
      flux(RNG2, imz) = uu(RNG2, imz) * vx(RNG2) - bb(RNG2, zdim) * bb(RNG2, xdim)
#ifndef ISO
      flux(RNG2, ien) = (uu(RNG2, ien) + ps(RNG2)) * vx(RNG2) - bb(RNG2, xdim) * (bb(RNG2, xdim) * uu(RNG2, imx) &
                                               + bb(RNG2, ydim) * uu(RNG2, imy) + bb(RNG2, zdim) * uu(RNG2, imz)) / uu(RNG2, idn)
#endif /* !ISO */
#else /* !MAGNETIC */
      flux(RNG2, imx) = uu(RNG2, imx) * vx(RNG2) + ps(RNG2)
      flux(RNG2, imy) = uu(RNG2, imy) * vx(RNG2)
      flux(RNG2, imz) = uu(RNG2, imz) * vx(RNG2)
#ifndef ISO
      flux(RNG2, ien) = (uu(RNG2, ien) + ps(RNG2)) * vx(RNG2)
#endif /* !ISO */
#endif /* !MAGNETIC */
      flux(1, :) = flux(2, :) ; flux(n, :) = flux(nm, :)

#ifdef LOCAL_FR_SPEED

      ! The freezing speed is now computed locally (in each cell)
      !  as in Trac & Pen (2003). This ensures much sharper shocks,
      !  but sometimes may lead to numerical instabilities
      minvx = minval(vx(RNG2))
      maxvx = maxval(vx(RNG2))
      amp   = half * (maxvx - minvx)
#ifdef ISO
      cfr(RNG2, 1) = sqrt(vx(RNG2)**2+cfr_smooth*amp) + max(sqrt( abs(2.0*pmag(RNG2) +          p(RNG2))/uu(RNG2, idn)),small)
#else /* !ISO */
      cfr(RNG2, 1) = sqrt(vx(RNG2)**2+cfr_smooth*amp) + max(sqrt( abs(2.0*pmag(RNG2) + this%gam*p(RNG2))/uu(RNG2, idn)),small)
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

   subroutine all_pres_ion(uu, n, bb, cs_iso2, gam, gam1, pmag, p, ps)

      use constants,  only: idn
#ifdef MAGNETIC
      use constants,  only: xdim, ydim, zdim
      use func,       only: emag
#endif /* MAGNETIC */
#ifndef ISO
      use constants,  only: imx, imy, imz, ien
      use dataio_pub, only: die
      use func,       only: ekin
#endif /* !ISO */

      implicit none

      integer(kind=4),      intent(in)             :: n         !< number of cells in the current sweep
      real, dimension(:,:), intent(in),    pointer :: uu        !< part of u for ionized fluid
      real, dimension(:,:), intent(in),    pointer :: bb        !< magnetic field x,y,z-components table
      real, dimension(:),   intent(in),    pointer :: cs_iso2   !< local isothermal sound speed squared (optional)
      real, dimension(:),   intent(inout), pointer :: ps        !< pressure of ionized fluid for current sweep
      real, dimension(:),   intent(inout)          :: p         !< thermal pressure of ionized fluid
      real, dimension(:),   intent(inout)          :: pmag      !< pressure of magnetic field
      real,                 intent(in)             :: gam, gam1

      ! locals
      integer :: nm

      nm = n - 1
#ifdef MAGNETIC
      pmag(RNG2) = emag(bb(RNG2, xdim), bb(RNG2, ydim), bb(RNG2, zdim));  pmag(1) = pmag(2); pmag(n) = pmag(nm)
#else /* !MAGNETIC */
      pmag(:) = 0.0
#endif /* !MAGNETIC */

#ifdef ISO
      p(RNG2)  = cs_iso2(RNG2) * uu(RNG2, idn)
      ps(RNG2) = p(RNG2) + pmag(RNG2)
#else /* !ISO */
      if (associated(cs_iso2)) call die("[initionized:all_pres_ion] cs_iso2 should not be associated")
      ps(RNG2) = (uu(RNG2, ien) - ekin(uu(RNG2, imx),uu(RNG2, imy),uu(RNG2, imz),uu(RNG2, idn)) )*gam1 + (2.0 - gam)*pmag(RNG2)
      p(RNG2) = ps(RNG2) - pmag(RNG2); p(1) = p(2); p(n) = p(nm)
#endif /* !ISO */
      ps(1) = ps(2); ps(n) = ps(nm)

#ifndef MAGNETIC
      return
      if (.false.) write(0,*) bb
#endif /* !MAGNETIC */
#ifdef ISO
      if (.false.) write(0,*) gam, gam1
#endif /* ISO */

   end subroutine all_pres_ion
#undef RNG2

   subroutine pres_ion(this, n, uu, bb, cs_iso2, ps)

      implicit none

      class(ion_fluid),     intent(in)             :: this
      integer(kind=4),      intent(in)             :: n        !< number of cells in the current sweep
      real, dimension(:,:), intent(in),    pointer :: uu       !< part of u for ionized fluid
      real, dimension(:,:), intent(in),    pointer :: bb       !< magnetic field x,y,z-components table
      real, dimension(:),   intent(in),    pointer :: cs_iso2  !< local isothermal sound speed squared (optional)
      real, dimension(:),   intent(inout), pointer :: ps       !< pressure of ionized fluid for current sweep

      ! locals
      real, dimension(n) :: p           !< thermal pressure of ionized fluid
      real, dimension(n) :: pmag        !< pressure of magnetic field

      call all_pres_ion(uu, n, bb, cs_iso2, this%gam, this%gam_1, pmag, p, ps)

   end subroutine pres_ion

end module initionized
