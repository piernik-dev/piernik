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

#ifndef COSM_RAYS
#error COSM_RAYS is required for 1cell
#endif /* COSM_RAYS */

#ifndef COSM_RAY_ELECTRONS
#error COSM_RAY_ELECTRONS is required for 1cell
#endif /* COSM_RAY_ELECTRONS */

module initproblem

   use constants,    only: fnamelen

   implicit none

   private
   public :: read_problem_par, problem_initial_conditions, problem_pointers

   integer(kind=4)    :: norm_step
   real               :: d0, p0, bx0, by0, bz0, x0, y0, z0, beta_cr, amp_cr1
   real               :: u_d0, u_b0, div_v0, ub_ampl, omega_d, omega_b, Btot, u_d_t
   character(len=fnamelen) :: outfile

   namelist /PROBLEM_CONTROL/ d0, p0, bx0, by0, bz0, x0, y0, z0, beta_cr, amp_cr1

   namelist /CRESP_ISOLATED/ u_d0, u_b0, div_v0, ub_ampl, omega_d, omega_b, Btot, u_d_t, outfile
contains

!-----------------------------------------------------------------------------

   subroutine problem_pointers

      implicit none

   end subroutine problem_pointers

!-----------------------------------------------------------------------------

   subroutine read_problem_par

      use constants,  only: cbuff_len
      use dataio_pub, only: nh
      use mpisetup,   only: cbuff, rbuff, master, slave, piernik_MPI_Bcast

      implicit none

      d0             = 1.0       !< density
      p0             = 1.0         !< pressure
      bx0            =   0.        !< Magnetic field component x
      by0            =   0.        !< Magnetic field component y
      bz0            =   0.        !< Magnetic field component z
      x0             = 0.0         !< x-position of the blob
      y0             = 0.0         !< y-position of the blob
      z0             = 0.0         !< z-position of the blob

      beta_cr        = 0.0         !< ambient level
      amp_cr1        = 1.0         !< amplitude of the blob

      outfile        = 'crs.dat'

      if (master) then

         if (.not.nh%initialized) call nh%init()
         open(newunit=nh%lun, file=nh%tmp1, status="unknown")
         write(nh%lun,nml=PROBLEM_CONTROL)
         close(nh%lun)
         open(newunit=nh%lun, file=nh%par_file)
         nh%errstr=""
         read(unit=nh%lun, nml=PROBLEM_CONTROL, iostat=nh%ierrh, iomsg=nh%errstr)
         close(nh%lun)
         call nh%namelist_errh(nh%ierrh, "PROBLEM_CONTROL")
         read(nh%cmdl_nml,nml=PROBLEM_CONTROL, iostat=nh%ierrh)
         call nh%namelist_errh(nh%ierrh, "PROBLEM_CONTROL", .true.)
         open(newunit=nh%lun, file=nh%tmp2, status="unknown")
         write(nh%lun,nml=PROBLEM_CONTROL)
         close(nh%lun)
         call nh%compare_namelist()

         rbuff(1)  = d0
         rbuff(2)  = p0
         rbuff(3)  = bx0
         rbuff(4)  = by0
         rbuff(5)  = bz0
         rbuff(6)  = x0
         rbuff(7)  = y0
         rbuff(8)  = z0
         rbuff(9)  = beta_cr
         rbuff(10) = amp_cr1

      endif

      call piernik_MPI_Bcast(rbuff)

      if (slave) then

         d0        = rbuff(1)
         p0        = rbuff(2)
         bx0       = rbuff(3)
         by0       = rbuff(4)
         bz0       = rbuff(5)
         x0        = rbuff(6)
         y0        = rbuff(7)
         z0        = rbuff(8)
         beta_cr   = rbuff(9)
         amp_cr1   = rbuff(10)

      endif

! Read the second namelist
      if (master) then

         if (.not.nh%initialized) call nh%init()
         open(newunit=nh%lun, file=nh%tmp1, status="unknown")
         write(nh%lun,nml=CRESP_ISOLATED)
         close(nh%lun)
         open(newunit=nh%lun, file=nh%par_file)
         nh%errstr=""
         read(unit=nh%lun, nml=CRESP_ISOLATED, iostat=nh%ierrh, iomsg=nh%errstr)
         close(nh%lun)
         call nh%namelist_errh(nh%ierrh, "CRESP_ISOLATED")
         read(nh%cmdl_nml,nml=CRESP_ISOLATED, iostat=nh%ierrh)
         call nh%namelist_errh(nh%ierrh, "CRESP_ISOLATED", .true.)
         open(newunit=nh%lun, file=nh%tmp2, status="unknown")
         write(nh%lun,nml=CRESP_ISOLATED)
         close(nh%lun)
         call nh%compare_namelist()

         rbuff(1)  = u_d0
         rbuff(2)  = u_b0
         rbuff(3)  = div_v0
         rbuff(4)  = ub_ampl
         rbuff(5)  = omega_d
         rbuff(6)  = omega_b
         rbuff(7)  = Btot
         rbuff(8)  = u_d_t

         cbuff(1)  = outfile

      endif

      call piernik_MPI_Bcast(rbuff)
      call piernik_MPI_Bcast(cbuff, cbuff_len)

      if (slave) then

         u_d0      = rbuff(1)
         u_b0      = rbuff(2)
         div_v0     = rbuff(3)
         ub_ampl   = rbuff(4)
         omega_d   = rbuff(5)
         omega_b   = rbuff(6)
         Btot      = rbuff(7)
         u_d_t     = rbuff(8)

         outfile   = trim(cbuff(1))

      endif

!  Rewind the outfile if it exists
      if (master) then
         open(newunit=nh%lun, file=outfile, status="replace", position="rewind")
         close(nh%lun)
      endif

   end subroutine read_problem_par

!-----------------------------------------------------------------------------

   subroutine problem_initial_conditions

      use cg_leaves,          only: leaves
      use cg_list,            only: cg_list_element
      use constants,          only: xdim, ydim, zdim, LO, HI, zero
      use domain,             only: dom
      use dataio_pub,         only: msg, printinfo
      use fluidindex,         only: flind
      use fluidtypes,         only: component_fluid
      use func,               only: ekin, emag
      use grid_cont,          only: grid_container
      use global,             only: skip_sweep, repeat_step
      use initcosmicrays,     only: iarr_crn, iarr_crs
      use user_hooks,         only: problem_customize_solution
#ifdef COSM_RAY_ELECTRONS
      use cresp_crspectrum,   only: cresp_get_scaled_init_spectrum
      use initcosmicrays,     only: iarr_cre_e, iarr_cre_n
      use initcrspectrum,     only: smallcree, cresp, cre_eff, use_cresp, adiab_active
#endif /* COSM_RAY_ELECTRONS */

      implicit none

      real                            :: cs_iso, decr
      type(cg_list_element),  pointer :: cgl
      type(grid_container),   pointer :: cg
#ifdef COSM_RAY_ELECTRONS
      real                            :: e_tot
#endif /* COSM_RAY_ELECTRONS */
      class(component_fluid), pointer :: fl
      integer                         :: i, j, k, icr, ipm, jpm, kpm

! Skip all the sweeps to simulate isolated case evolution -- allows to test CRESP algorithm without spatial transport
      skip_sweep  = [.true., .true., .true.]

      if (adiab_active) then
         repeat_step = .false.
         write (msg, *) "Parameter 'adiab_active' is TRUE, due to possible large variations of div(V) CFL warnings may show up: changing 'repeat_step' to FALSE (hard coded)"
         call printinfo(msg)
      endif

      problem_customize_solution => isolated_problem_customize_solution

      fl => flind%ion

      cs_iso = sqrt(p0/d0)
! NOTICE: Btot is preferred over bx0, by0 and bz0 if these are set as well and overwrites them
      if (Btot .gt. 0.) then  ! Distribute Btot evenly over all available directions
         if (dom%has_dir(xdim)) bx0 = sqrt(Btot**2) / (dom%D_x + dom%D_y + dom%D_z)
         if (dom%has_dir(ydim)) by0 = sqrt(Btot**2) / (dom%D_x + dom%D_y + dom%D_z)
         if (dom%has_dir(zdim)) bz0 = sqrt(Btot**2) / (dom%D_x + dom%D_y + dom%D_z)
      endif
! ignore B field in nonexistent direction
      if (.not.dom%has_dir(xdim)) bx0 = zero
      if (.not.dom%has_dir(ydim)) by0 = zero
      if (.not.dom%has_dir(zdim)) bz0 = zero

      cgl => leaves%first
      do while (associated(cgl))
         cg => cgl%cg

         call cg%set_constant_b_field([bx0, by0, bz0])  ! this acts only inside cg%ijkse box
         cg%u(fl%idn,:,:,:) = d0
         cg%u(fl%imx:fl%imz,:,:,:) = zero

#ifndef ISO
         do k = cg%ks, cg%ke
            do j = cg%js, cg%je
               do i = cg%is, cg%ie
                  cg%u(fl%ien,i,j,k) = p0/fl%gam_1 + &
                       &               ekin(cg%u(fl%imx,i,j,k), cg%u(fl%imy,i,j,k), cg%u(fl%imz,i,j,k), cg%u(fl%idn,i,j,k)) + &
                       &               emag(cg%b(xdim,i,j,k), cg%b(ydim,i,j,k), cg%b(zdim,i,j,k))
               enddo
            enddo
         enddo
#endif /* !ISO */

         cg%u(iarr_crs,:,:,:) = 0.0

! Spatial distribution: uniformly fill the whole domain
         do k = cg%ks, cg%ke
            do j = cg%js, cg%je
               do i = cg%is, cg%ie

                  decr = amp_cr1

                  cg%u(iarr_crn(1), i, j, k) = cg%u(iarr_crn(1), i, j, k) + amp_cr1 * decr
#ifdef COSM_RAY_ELECTRONS
                  if (use_cresp) then
                     e_tot = cre_eff * decr
                     if (e_tot > smallcree) then
                        cresp%n = 0.0 ;  cresp%e = 0.0
                        call cresp_get_scaled_init_spectrum(cresp%n, cresp%e, e_tot)
                        cg%u(iarr_cre_n,i,j,k) = cg%u(iarr_cre_n,i,j,k) + cresp%n
                        cg%u(iarr_cre_e,i,j,k) = cg%u(iarr_cre_e,i,j,k) + cresp%e
                     endif
                  endif
#endif /* COSM_RAY_ELECTRONS */
               enddo
            enddo
         enddo

         cgl => cgl%nxt
      enddo

      call append_cooling_terms(.true.)

   end subroutine problem_initial_conditions

   subroutine isolated_problem_customize_solution(forward)

      implicit none

      logical, intent(in)     :: forward
      logical                 :: set_tabs_for_cresp

      set_tabs_for_cresp = forward
      call append_cooling_terms(set_tabs_for_cresp)
      if (forward .eqv. .false.) then
         call printer
      endif

   end subroutine isolated_problem_customize_solution

   subroutine append_cooling_terms(set_tabs_for_cresp)

      use cg_leaves,       only: leaves
      use cg_list,         only: cg_list_element
      use constants,       only: zero
      use crhelpers,       only: divv_i, div_v
      use domain,          only: dom
      use dataio_pub,      only: msg, printinfo
      use fluidindex,      only: flind
      use fluidtypes,      only: component_fluid
      use grid_cont,       only: grid_container
      use global,          only: t
      use initcrspectrum,  only: adiab_active, synch_active

      implicit none

      type(cg_list_element),  pointer :: cgl
      type(grid_container),   pointer :: cg
      class(component_fluid), pointer :: fl
      real                            :: denom_dims, cos_omega_t
      integer                         :: i, j, k
      logical                         :: set_tabs_for_cresp

      fl => flind%ion
      cgl => leaves%first

      do while (associated(cgl))
         cg => cgl%cg

         if (synch_active) then
            call cg%set_constant_b_field([bx0, by0, bz0])  ! this acts only inside cg%ijkse box
         endif
         if (adiab_active) then
            denom_dims  = 1. / (dom%D_x + dom%D_y + dom%D_z)
            cos_omega_t = cos(omega_d * t)
            if (set_tabs_for_cresp .eqv. .true.) then

               do k = cg%ks, cg%ke
                  do j = cg%js, cg%je
                     do i = cg%is, cg%ie
                           cg%u(flind%ion%imx,i,j,k) = cg%u(flind%ion%idn,i,j,k) * cg%x(i) * (u_d0 + div_v0 * cos_omega_t) * denom_dims
                           cg%u(flind%ion%imy,i,j,k) = cg%u(flind%ion%idn,i,j,k) * cg%y(j) * (u_d0 + div_v0 * cos_omega_t) * denom_dims
                           cg%u(flind%ion%imz,i,j,k) = cg%u(flind%ion%idn,i,j,k) * cg%z(k) * (u_d0 + div_v0 * cos_omega_t) * denom_dims

                     enddo
                  enddo
               enddo

               call div_v(flind%ion%pos, cg)

            endif
         endif

#ifdef CRESP_VERBOSED
         write (msg, *) "Got values: u_d(numerical)", cg%q(divv_i)%point([i,j,k]), "u_d(t, set):", u_d0 + div_v0 * cos_omega_t
         if (adiab_active) call printinfo(msg)
#endif /* CRESP_VERBOSED */

         cgl => cgl%nxt
      enddo

   end subroutine append_cooling_terms

   subroutine printer

      use constants,       only: LO, HI
      use dataio_pub,      only: nh
      use global,          only: t, dt, repeat_step, cfl_violated
      use initcosmicrays,  only: ncre
      use initcrspectrum,  only: crel

      implicit none

      open(newunit=nh%lun, file=outfile, status="unknown", position="append")
      if (repeat_step .and. cfl_violated) then
         backspace(nh%lun)  ! rewind one line if step is repeated in order to keep consistent order of the data in crs file
      endif
      write(nh%lun, '(2e16.9, 3(1x,i8), 600(1x,ES18.9E3))') t, dt, ncre, crel%i_cut(LO), crel%i_cut(HI), crel%p, crel%f, crel%q
      close(nh%lun)

   end subroutine printer

end module initproblem
