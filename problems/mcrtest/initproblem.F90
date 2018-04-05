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

module initproblem

   implicit none

   private
   public :: read_problem_par, problem_initial_conditions, problem_pointers

   integer(kind=4)    :: norm_step
   real               :: t_sn
   real               :: d0, p0, bx0, by0, bz0, x0, y0, z0, r0, beta_cr, amp_cr, vxd0, vyd0, vzd0, expansion_cnst

   namelist /PROBLEM_CONTROL/ d0, p0, bx0, by0, bz0, x0, y0, z0, r0, vxd0, vyd0, vzd0, beta_cr, amp_cr, norm_step, expansion_cnst
contains

!-----------------------------------------------------------------------------

   subroutine problem_pointers

      implicit none

   end subroutine problem_pointers

!-----------------------------------------------------------------------------

   subroutine read_problem_par

      use constants,  only: I_TEN
      use dataio_pub, only: nh ! QA_WARN required for diff_nml
      use dataio_pub, only: die
      use domain,     only: dom
      use func,       only: operator(.equals.)
      use mpisetup,   only: ibuff, rbuff, master, slave, piernik_MPI_Bcast

      implicit none

      t_sn = 0.0

      d0           = 1.0e5     !< density
      p0           = 1.0       !< pressure
      bx0          =   0.      !< Magnetic field component x
      by0          =   0.      !< Magnetic field component y
      bz0          =   0.      !< Magnetic field component z
      x0           = 0.0       !< x-position of the blob
      y0           = 0.0       !< y-position of the blob
      z0           = 0.0       !< z-position of the blob
      r0           = 5.* minval(dom%L_(:)/dom%n_d(:), mask=dom%has_dir(:))  !< radius of the blob
      vxd0         = 0.0       !< initial velocity_x, refers to whole domain
      vyd0         = 0.0       !< initial velocity_y, refers to whole domain
      vzd0         = 0.0       !< initial velocity_z, refers to whole domain

      beta_cr      = 0.0       !< ambient level
      amp_cr       = 1.0       !< amplitude of the blob

      norm_step    = I_TEN     !< how often to compute the norm (in steps)

      expansion_cnst = 0.0

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
         rbuff(9)  = r0
         rbuff(10) = beta_cr
         rbuff(11) = amp_cr
         rbuff(12) = vxd0
         rbuff(13) = vyd0
         rbuff(14) = vzd0
         rbuff(15) = expansion_cnst

         ibuff(1)  = norm_step

      endif

      call piernik_MPI_Bcast(ibuff)
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
         r0        = rbuff(9)
         beta_cr   = rbuff(10)
         amp_cr    = rbuff(11)
         vxd0      = rbuff(12)
         vyd0      = rbuff(13)
         vzd0      = rbuff(14)
         expansion_cnst = rbuff(15)

         norm_step = int(ibuff(1), kind=4)

      endif

      if (r0 .equals. 0.) call die("[initproblem:read_problem_par] r0 == 0")

   end subroutine read_problem_par

!-----------------------------------------------------------------------------

   subroutine problem_initial_conditions

      use cg_leaves,      only: leaves
      use cg_list,        only: cg_list_element
      use constants,      only: xdim, ydim, zdim, LO, HI, pMAX
      use dataio_pub,     only: msg, warn, printinfo, die
      use domain,         only: dom, is_multicg
      use fluidindex,     only: flind
      use fluidtypes,     only: component_fluid
      use func,           only: ekin, emag, operator(.equals.), operator(.notequals.)
      use grid_cont,      only: grid_container
      use initcosmicrays, only: iarr_crn, iarr_crs, gamma_crn, K_crn_paral, K_crn_perp
      use mpisetup,       only: master, piernik_MPI_Allreduce
#ifdef COSM_RAYS_SOURCES
      use cr_data,        only: icr_H1, icr_C12, cr_table
#endif /* COSM_RAYS_SOURCES */
#ifdef COSM_RAY_ELECTRONS
     use initcosmicrays,   only: ncrn, iarr_cre_e, iarr_cre_n
     use initcrspectrum,   only: expan_order, taylor_coeff_2nd, taylor_coeff_3rd, e_small, cresp, cre_eff
     use cresp_crspectrum, only: e_tot_2_f_init_params, cresp_init_powl_spectrum, cresp_init_state, e_tot_2_en_powl_init_params
#endif /* COSM_RAY_ELECTRONS */

      implicit none

      class(component_fluid), pointer :: fl
      integer                         :: i, j, k, icr, ipm, jpm, kpm
      real                            :: cs_iso, xsn, ysn, zsn, r2, maxv, e_tot
      real                            :: sn_exp, sn_rdist2
      type(cg_list_element),  pointer :: cgl
      type(grid_container),   pointer :: cg
      logical        :: first_run = .true.
#ifndef COSM_RAYS_SOURCES
      integer, parameter              :: icr_H1 = 1, icr_C12 = 2
#endif /* !COSM_RAYS_SOURCES */

     if (first_run .eqv. .true.) then
      fl => flind%ion

      ! BEWARE: temporary fix
      xsn = x0
      ysn = y0
      zsn = z0

! Uniform equilibrium state

      cs_iso = sqrt(p0/d0)

      if (.not.dom%has_dir(xdim)) bx0 = 0. ! ignore B field in nonexistent direction
      if (.not.dom%has_dir(ydim)) by0 = 0.
      if (.not.dom%has_dir(zdim)) bz0 = 0.

      if ((bx0**2 + by0**2 + bz0**2 .equals. 0.) .and. (any(K_crn_paral(:) .notequals. 0.) .or. any(K_crn_perp(:) .notequals. 0.))) then
         call warn("[initproblem:problem_initial_conditions] No magnetic field is set, K_crn_* also have to be 0.")
         K_crn_paral(:) = 0.
         K_crn_perp(:)  = 0.
      endif


      cgl => leaves%first
      do while (associated(cgl))
         cg => cgl%cg

         call cg%set_constant_b_field([bx0, by0, bz0])
         cg%u(fl%idn, :, :, :) = d0
         cg%u(fl%imx:fl%imz, :, :, :) = 0.0

#ifndef ISO
         do k = cg%lhn(zdim,LO), cg%lhn(zdim,HI)
            do j = cg%lhn(ydim,LO), cg%lhn(ydim,HI)
               do i = cg%lhn(xdim,LO), cg%lhn(xdim,HI)
                  cg%u(fl%ien,i,j,k) = p0/fl%gam_1 + &
                       &               ekin(cg%u(fl%imx,i,j,k), cg%u(fl%imy,i,j,k), cg%u(fl%imz,i,j,k), cg%u(fl%idn,i,j,k)) + &
                       &               emag(cg%b(xdim,i,j,k), cg%b(ydim,i,j,k), cg%b(zdim,i,j,k))
               enddo
            enddo
         enddo
#endif /* !ISO */

#ifdef COSM_RAYS
         do icr = 1, flind%crs%all
            cg%u(iarr_crs(icr), :, :, :) =  beta_cr*fl%cs2 * cg%u(fl%idn, :, :, :)/(gamma_crn(icr)-1.0)
         enddo

! Explosions
         do icr = 1, flind%crn%all
            do k = cg%ks, cg%ke
               do j = cg%js, cg%je
                  do i = cg%is, cg%ie

                     do ipm=-1,1
                        do jpm=-1,1
                           do kpm=-1,1
                              r2 = (cg%x(i)-xsn+real(ipm)*dom%L_(xdim))**2+(cg%y(j)-ysn+real(jpm)*dom%L_(ydim))**2+(cg%z(k)-zsn+real(kpm)*dom%L_(zdim))**2
                              sn_rdist2 = r2/r0**2
                              sn_exp = 0.0
                              if(sn_rdist2 <= 10.0) then
                                 sn_exp = exp(-sn_rdist2)
                              endif
                              if (icr == cr_table(icr_H1)) then
                                 cg%u(iarr_crn(icr), i, j, k) = cg%u(iarr_crn(icr), i, j, k) + amp_cr*sn_exp
                              elseif (icr == cr_table(icr_C12)) then
                                 cg%u(iarr_crn(icr), i, j, k) = cg%u(iarr_crn(icr), i, j, k) + amp_cr*0.1*sn_exp ! BEWARE: magic number
                              else
                                 cg%u(iarr_crn(icr), i, j, k) = 0.0
                              endif
                           enddo
                        enddo
                     enddo

                  enddo
               enddo
            enddo
         enddo
#ifdef COSM_RAY_ELECTRONS
! Explosions @CRESP independent of cr nucleons
         do k = cg%ks, cg%ke
            do j = cg%js, cg%je
               do i = cg%is, cg%ie

                  e_tot = 0.0
                  do ipm=-1,1
                     do jpm=-1,1
                        do kpm=-1,1
                           r2 = (cg%x(i)-xsn+real(ipm)*dom%L_(xdim))**2+(cg%y(j)-ysn+real(jpm)*dom%L_(ydim))**2+(cg%z(k)-zsn+real(kpm)*dom%L_(zdim))**2
                           sn_rdist2 = r2/r0**2
                           sn_exp = 0.0
                           if(sn_rdist2 <= 10.0) then
                              sn_exp = exp(-sn_rdist2)
                           endif
                           e_tot = e_tot + amp_cr*sn_exp
                        enddo
                     enddo
                  enddo
                  cresp%n = cg%u(iarr_cre_n,i,j,k) ;  cresp%e = cg%u(iarr_cre_e,i,j,k) ; e_tot = e_tot * cre_eff
                  if (first_run) then
                     if (e_tot .gt. e_small) then     ! early phase - fill cells only when total passed energy is greater than e_small, amplitude computed from total explosion energy multiplied by factor cre_eff
                        call cresp_init_state(cresp%n,cresp%e, e_tot_2_f_init_params(e_tot))  ! initializes whole spectrum, accounts for "widening" due to e_small approximation
                     endif
                  else
                     call e_tot_2_en_powl_init_params(cresp%n, cresp%e, e_tot)
                  endif
                  cg%u(iarr_cre_n,i,j,k) = cg%u(iarr_cre_n,i,j,k) + cresp%n
                  cg%u(iarr_cre_e,i,j,k) = cg%u(iarr_cre_e,i,j,k) + cresp%e

               enddo
            enddo
         enddo
#endif /* COSM_RAY_ELECTRONS */
         cgl => cgl%nxt
      enddo
      cg => leaves%first%cg

      if (is_multicg) call die("[initproblem:problem_initial_conditions] multiple grid pieces per procesor not implemented yet") !nontrivial maxv

      do icr = 1, flind%crs%all
         maxv = maxval(cg%u(iarr_crs(icr),:,:,:))
         call piernik_MPI_Allreduce(maxv, pMAX)
         if (master) then
#ifdef COSM_RAY_ELECTRONS
            if (icr .lt. flind%cre%nbeg) then
               write(msg,*) '[initproblem:problem_initial_conditions] icr(nuc)  =',icr,' maxecr(nuc) =',maxv
            else if (icr .lt. flind%cre%ebeg .and. icr .ge. flind%cre%nbeg) then
               write(msg,*) '[initproblem:problem_initial_conditions] icr(cre_n)=',icr,' maxncr(cre) =',maxv
            else
               write(msg,*) '[initproblem:problem_initial_conditions] icr(cre_e)=',icr,' maxecr(cre) =',maxv
            endif
#else
            write(msg,*) '[initproblem:problem_initial_conditions] icr=',icr,' maxecr =',maxv
#endif /* COSM_RAY_ELECTRONS */
            call printinfo(msg)
         endif
      enddo

#endif /* COSM_RAYS */
#ifdef COSM_RAY_ELECTRONS
      write(msg,*) '[initproblem:problem_initial_conditions]: Taylor_exp._ord. (cresp)    = ', expan_order
      call printinfo(msg)
      write(msg,*) '[initproblem:problem_initial conditions]: Taylor_exp._coeff.(2nd,3rd) = ', taylor_coeff_2nd, taylor_coeff_3rd
      call printinfo(msg)
#endif /* COSM_RAY_ELECTRONS */

! Velocity field
      cgl => leaves%first
      do while (associated(cgl))
        cg => cgl%cg
        if (expansion_cnst .ne. 0.0 ) then ! adiabatic expansion / compression
         write(msg,*) '[initproblem:problem_initial_conditions] setting up expansion/compression, expansion_cnst=',expansion_cnst
#ifdef IONIZED
! Ionized
         do k = cg%lhn(zdim,LO), cg%lhn(zdim,HI)
            do j = cg%lhn(ydim,LO), cg%lhn(ydim,HI)
               do i = cg%lhn(xdim,LO), cg%lhn(xdim,HI)
                    cg%u(flind%ion%imx,i,j,k) = cg%u(flind%ion%idn,i,j,k) * (cg%x(i)-x0) * expansion_cnst  !< vxd0 * rho
                    cg%u(flind%ion%imy,i,j,k) = cg%u(flind%ion%idn,i,j,k) * (cg%y(j)-y0) * expansion_cnst  !< vyd0 * rho
                    cg%u(flind%ion%imz,i,j,k) = cg%u(flind%ion%idn,i,j,k) * (cg%z(k)-z0) * expansion_cnst  !< vzd0 * rho
                  enddo
                enddo
            enddo
#endif /* IONIZED */
        else

#ifdef IONIZED
! Ionized
            cg%u(flind%ion%imx,:,:,:) = vxd0 * cg%u(flind%ion%idn,:,:,:)
            cg%u(flind%ion%imy,:,:,:) = vyd0 * cg%u(flind%ion%idn,:,:,:)
            cg%u(flind%ion%imz,:,:,:) = vzd0 * cg%u(flind%ion%idn,:,:,:)
#endif /* IONIZED */

        endif
        cgl => cgl%nxt
      enddo

      first_run = .false.
    endif
   end subroutine problem_initial_conditions
end module initproblem
