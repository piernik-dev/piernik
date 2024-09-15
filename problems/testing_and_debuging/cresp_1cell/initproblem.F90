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

#ifndef CRESP
#error CRESP is required for 1cell
#endif /* CRESP */

module initproblem

   use constants, only: cbuff_len

   implicit none

   private
   public :: read_problem_par, problem_initial_conditions, problem_pointers

   real                     :: d0, p0, bx0, by0, bz0, amp_cr1, u_d0, u_b0, u_d_ampl, omega_d, Btot
   character(len=cbuff_len) :: outfile

   namelist /PROBLEM_CONTROL/ d0, p0, bx0, by0, bz0, amp_cr1, u_d0, u_b0, u_d_ampl, omega_d, Btot, outfile

contains

!-----------------------------------------------------------------------------

   subroutine problem_pointers

      use user_hooks, only: problem_customize_solution, user_reaction_to_redo_step

      implicit none

      problem_customize_solution => isolated_problem_customize_solution
      user_reaction_to_redo_step => printer_rewind_one_line

   end subroutine problem_pointers

!-----------------------------------------------------------------------------

   subroutine read_problem_par

      use bcast,      only: piernik_MPI_Bcast
      use dataio_pub, only: nh
      use mpisetup,   only: cbuff, rbuff, master, slave

      implicit none

      d0       = 1.0         !< density
      p0       = 1.0         !< pressure
      bx0      =   0.        !< Magnetic field component x
      by0      =   0.        !< Magnetic field component y
      bz0      =   0.        !< Magnetic field component z
      amp_cr1  = 1.0         !< amplitude of the blob
      u_b0     = 0.0         !< initial magnetic energy-density
      u_d0     = 0.0         !< initial magnitude of div_v
      u_d_ampl = -0.2        !< amplitude of variable u_d component (periodic), sums up with u_d0 to u_d
      omega_d  = 0.157079633 !< omega_d parameter for test with periodic adiabatic compression: u_d(t) = u_d0 + u_d_ampl * cos(omega_d * t)
      Btot     = 0.          !< Total amplitude of MF in micro Gauss

      outfile  = 'crs.dat'

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
         rbuff(6)  = amp_cr1
         rbuff(7)  = u_d0
         rbuff(8)  = u_b0
         rbuff(9)  = u_d_ampl
         rbuff(10) = omega_d
         rbuff(11) = Btot

         cbuff(1)  = outfile

      endif

      call piernik_MPI_Bcast(rbuff)
      call piernik_MPI_Bcast(cbuff, cbuff_len)

      if (slave) then

         d0       = rbuff(1)
         p0       = rbuff(2)
         bx0      = rbuff(3)
         by0      = rbuff(4)
         bz0      = rbuff(5)
         amp_cr1  = rbuff(6)
         u_d0     = rbuff(7)
         u_b0     = rbuff(8)
         u_d_ampl = rbuff(9)
         omega_d  = rbuff(10)
         Btot     = rbuff(11)

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

      use cg_leaves,        only: leaves
      use cg_list,          only: cg_list_element
      use constants,        only: xdim, ydim, zdim, zero, two
      use cresp_crspectrum, only: cresp_get_scaled_init_spectrum
      use domain,           only: dom
      use fluidindex,       only: flind
      use fluidtypes,       only: component_fluid
      use grid_cont,        only: grid_container
      use initcosmicrays,   only: iarr_crn, iarr_crs, iarr_cre_e, iarr_cre_n
      use initcrspectrum,   only: smallcree, cresp, cre_eff, use_cresp, f_synchIC, crel, total_init_cree
#ifndef ISO
      use func,             only: ekin, emag
#endif /* !ISO */

      implicit none

      class(component_fluid), pointer :: fl
      type(cg_list_element),  pointer :: cgl
      type(grid_container),   pointer :: cg
      real                            :: e_tot
      integer                         :: i, j, k

      fl => flind%ion

! NOTICE: Btot is preferred over bx0, by0 and bz0 if these are set as well and overwrites them, B field in nonexistent direction is ignored.
      if (Btot > zero .or. u_b0 > zero) then  ! Distribute Btot evenly over all available directions
         if (dom%has_dir(xdim)) bx0 = sqrt(Btot**2 / dom%eff_dim) + sqrt(two * u_b0 / f_synchIC / dom%eff_dim)
         if (dom%has_dir(ydim)) by0 = sqrt(Btot**2 / dom%eff_dim) + sqrt(two * u_b0 / f_synchIC / dom%eff_dim)
         if (dom%has_dir(zdim)) bz0 = sqrt(Btot**2 / dom%eff_dim) + sqrt(two * u_b0 / f_synchIC / dom%eff_dim)
      endif

      cgl => leaves%first
      do while (associated(cgl))
         cg => cgl%cg

         call cg%set_constant_b_field([bx0, by0, bz0])  ! this acts only inside cg%ijkse box
         cg%u(fl%idn,RNG) = d0
         cg%u(fl%imx:fl%imz,RNG) = zero

#ifndef ISO
         cg%u(fl%ien,RNG) = p0/fl%gam_1 + ekin(cg%u(fl%imx,RNG), cg%u(fl%imy,RNG), cg%u(fl%imz,RNG), cg%u(fl%idn,RNG)) + &
                          &               emag(cg%b(xdim,RNG), cg%b(ydim,RNG), cg%b(zdim,RNG))
#endif /* !ISO */

         cg%u(iarr_crs,RNG) = zero

! Spatial distribution: uniformly fill the whole domain
         cg%u(iarr_crn(1),RNG) = cg%u(iarr_crn(1),RNG) + amp_cr1**2

         if (use_cresp) then
            e_tot = cre_eff * amp_cr1
            if (e_tot > smallcree) then
               cresp%n = zero ;  cresp%e = zero
               call cresp_get_scaled_init_spectrum(cresp%n, cresp%e, e_tot)

               do k = cg%ks, cg%ke
                  do j = cg%js, cg%je
                     do i = cg%is, cg%ie
                        cg%u(iarr_cre_n,i,j,k) = cg%u(iarr_cre_n,i,j,k) + cresp%n
                        cg%u(iarr_cre_e,i,j,k) = cg%u(iarr_cre_e,i,j,k) + cresp%e
                     enddo
                  enddo
               enddo
            endif
         endif

         cgl => cgl%nxt
      enddo

      call append_cooling_terms

      crel%f = crel%f * e_tot / total_init_cree ! Scale f before saving initial spectrum
      call printer ! Save initial spectrum; CRESP is initialized at this point.

   end subroutine problem_initial_conditions

   subroutine isolated_problem_customize_solution(forward)

      implicit none

      logical, intent(in) :: forward

      if (forward) then
         call append_cooling_terms
      else
         call printer
      endif

   end subroutine isolated_problem_customize_solution

   subroutine append_cooling_terms

      use cg_leaves,      only: leaves
      use cg_list,        only: cg_list_element
      use constants,      only: LO, HI, three, xdim, ydim, zdim
      use crhelpers,      only: div_v
      use domain,         only: dom
      use fluidindex,     only: flind
      use fluidtypes,     only: component_fluid
      use grid_cont,      only: grid_container
      use global,         only: t
      use initcrspectrum, only: adiab_active, synch_active
#ifndef ISO
      use func,           only: ekin, emag
#endif /* !ISO */
#ifdef CRESP_VERBOSED
      use crhelpers,      only: divv_i
      use dataio_pub,     only: msg, printinfo
#endif /* CRESP_VERBOSED */

      implicit none

      type(cg_list_element),  pointer     :: cgl
      type(grid_container),   pointer     :: cg
      class(component_fluid), pointer     :: fl
      integer                             :: i, j, k
      real                                :: cos_f, cos_omega_t, denom_dims
#ifndef ISO
      real, dimension(:,:,:), allocatable :: int_ener
#endif /* !ISO */

      fl => flind%ion
      cgl => leaves%first

      do while (associated(cgl))
         cg => cgl%cg

#ifndef ISO
         allocate(int_ener(cg%lhn(xdim,LO):cg%lhn(xdim,HI), cg%lhn(ydim,LO):cg%lhn(ydim,HI), cg%lhn(zdim,LO):cg%lhn(zdim,HI)))
         int_ener = cg%u(fl%ien,:,:,:) - ekin(cg%u(fl%imx,:,:,:), cg%u(fl%imy,:,:,:), cg%u(fl%imz,:,:,:), cg%u(fl%idn,:,:,:)) - &
                  &                      emag(cg%b(xdim,:,:,:), cg%b(ydim,:,:,:), cg%b(zdim,:,:,:))
#endif /* !ISO */

         if (synch_active) call cg%set_constant_b_field([bx0, by0, bz0])  ! this acts only inside cg%ijkse box

         if (adiab_active) then
            denom_dims  = three / max(dom%eff_dim, 1)
            cos_omega_t = u_d0 + u_d_ampl * cos(omega_d * t)
            cos_f = cos_omega_t * denom_dims

            do i = cg%lhn(xdim, LO), cg%lhn(xdim, HI)
               cg%u(fl%imx,i,:,:) = cg%u(fl%idn,i,:,:) * cg%x(i) * cos_f
            enddo
            do j = cg%lhn(ydim, LO), cg%lhn(ydim, HI)
               cg%u(fl%imy,:,j,:) = cg%u(fl%idn,:,j,:) * cg%y(j) * cos_f
            enddo
            do k = cg%lhn(zdim, LO), cg%lhn(zdim, HI)
               cg%u(fl%imz,:,:,k) = cg%u(fl%idn,:,:,k) * cg%z(k) * cos_f
            enddo

            call div_v(fl%pos, cg)

#ifdef CRESP_VERBOSED
            write (msg, "(A,F10.7,A,F10.7)") "Adiabatic process: got u_d(0, 0, 0) values : u_d(numerical) = ", cg%q(divv_i)%point([0,0,0]) / three, " | u_d(t, set) = ", cos_omega_t
            call printinfo(msg)
#endif /* CRESP_VERBOSED */
         endif

#ifndef ISO
         cg%u(fl%ien,:,:,:) = int_ener + ekin(cg%u(fl%imx,:,:,:), cg%u(fl%imy,:,:,:), cg%u(fl%imz,:,:,:), cg%u(fl%idn,:,:,:)) + &
                  &                      emag(cg%b(xdim,:,:,:), cg%b(ydim,:,:,:), cg%b(zdim,:,:,:))
         deallocate(int_ener)
#endif /* !ISO */

         cgl => cgl%nxt
      enddo

   end subroutine append_cooling_terms

   subroutine printer

      use constants,      only: LO, HI
      use dataio_pub,     only: nh
      use global,         only: t, dt
      use initcosmicrays, only: nspc
      use initcrspectrum, only: crel

      implicit none

      open(newunit=nh%lun, file=outfile, status="unknown", position="append")
      write(nh%lun, '(2e16.9, 3(1x,i8), 600(1x,ES18.9E3))') t, dt, nspc, crel%i_cut(LO), crel%i_cut(HI), crel%p, crel%f, crel%q
      close(nh%lun)

   end subroutine printer

   subroutine printer_rewind_one_line

      use dataio_pub, only: nh

      implicit none

      open(newunit=nh%lun, file=outfile, status="unknown", position="append")
      backspace(nh%lun)  ! rewind one line if step is repeated in order to keep consistent order of the data in crs file
      close(nh%lun)

   end subroutine printer_rewind_one_line

end module initproblem
