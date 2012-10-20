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

module initproblem

! Initial condition for Keplerian disk
! Written by: M. Hanasz, March 2006

   use constants,    only: cbuff_len, dsetnamelen

   implicit none

   private
   public :: read_problem_par, init_prob, problem_pointers

   real                     :: d0, r_max, dout, alpha, r_in, r_out, f_in, f_out
   real                     :: dens_exp      !< exponent in profile density \f$\rho(R) = \rho_0 R^{-k}\f$
   real                     :: dens_amb      !< density of ambient medium (used for inner cutoff)
   real                     :: eps           !< dust to gas ratio
   real                     :: x_cut, a_cut  !< radius of inner disk cut-off
   real                     :: dens_max
   integer(kind=4)          :: cutoff_ncells !< width of cut-off profile
   real, save               :: T_inner = 0.0 !< Orbital period at the inner boundary, \todo save it to restart as an attribute
   !>
   !! \f$\tau\f$ in \f$\frac{Du}{Dt} = - \frac{u-u_0}{\tau}f(R)
   !! when initproblem::problem_customize_solution is used
   !<
   real                     :: dumping_coeff, drag_max, drag_min, amp_noise
   logical                  :: use_inner_orbital_period  !< use 1./T_inner as dumping_coeff
   integer(kind=4), parameter :: ngauss = 4
   real, dimension(ngauss)  :: gauss
   character(len=cbuff_len) :: mag_field_orient
   character(len=cbuff_len) :: densfile
   real, dimension(:), allocatable :: taus, tauf
   character(len=dsetnamelen), parameter :: inid_n = "u_0"

   namelist /PROBLEM_CONTROL/  alpha, d0, dout, r_max, mag_field_orient, r_in, r_out, f_in, f_out, &
      & dens_exp, eps, dens_amb, x_cut, cutoff_ncells, dumping_coeff, use_inner_orbital_period, &
      & drag_max, drag_min, densfile, amp_noise, gauss, a_cut, dens_max

contains
!-----------------------------------------------------------------------------
   subroutine problem_pointers

      !use dataio_user, only: user_vars_hdf5, user_reg_var_restart
      !use user_hooks,  only: finalize_problem, problem_customize_solution

      implicit none

      !user_reg_var_restart       => register_user_var

   end subroutine problem_pointers
!-----------------------------------------------------------------------------
   subroutine register_user_var

      use cg_list_global,   only: all_cg
      use constants,        only: AT_NO_B
      use named_array_list, only: wna

      implicit none

      integer(kind=4) :: dim4 !< BEWARE: workaround for gcc-4.6 bug

      dim4 = wna%lst(wna%fi)%dim4
      call all_cg%reg_var(inid_n, restart_mode = AT_NO_B, dim4 = dim4)

   end subroutine register_user_var
!-----------------------------------------------------------------------------
   subroutine read_problem_par

      use constants,             only: GEO_RPZ
      use dataio_pub,            only: nh      ! QA_WARN required for diff_nml
      use dataio_user,           only: user_vars_hdf5, user_reg_var_restart
      use domain,                only: dom
      use fluidboundaries_funcs, only: user_fluidbnd
      use gravity,               only: grav_pot_3d
      use mpisetup,              only: cbuff, rbuff, ibuff, lbuff, master, slave, piernik_MPI_Bcast
      use user_hooks,            only: problem_customize_solution, problem_grace_passed, problem_post_restart

      implicit none

      d0               = 1.0
      dout             = 1.0e-4
      r_max            = 1.0
      mag_field_orient = 'none'
      densfile         = ''
      alpha            = 1.0
      amp_noise        = 1.e-6
      gauss            = 0.0

      r_in             = 0.5
      f_in             = 10.0
      r_out            = 2.1
      f_out            = 0.0

      dens_exp         = 0.0
      dens_amb         = 1.e-3
      eps              = 1.0
      x_cut            = 2.5
      a_cut            = 15.0
      dens_max         = 500.0

      cutoff_ncells    = 8
      dumping_coeff    = 1.0

      use_inner_orbital_period = .false.

      if (master) then

         diff_nml(PROBLEM_CONTROL)

         ibuff(1) = cutoff_ncells

         lbuff(1) = use_inner_orbital_period

         cbuff(1) = mag_field_orient
         cbuff(2) = densfile

         rbuff(1) = d0
         rbuff(2) = dout
         rbuff(3) = r_max
         rbuff(4) = alpha
         rbuff(5) = r_in
         rbuff(6) = r_out
         rbuff(7) = f_in
         rbuff(8) = f_out
         rbuff(9)  = dens_exp
         rbuff(10) = eps
         rbuff(11) = dens_amb
         rbuff(12) = x_cut
         rbuff(13) = dumping_coeff
         rbuff(14) = drag_max
         rbuff(15) = drag_min
         rbuff(16) = amp_noise
         rbuff(17:20) = gauss
         rbuff(21) = a_cut
         rbuff(22) = dens_max

      endif

      call piernik_MPI_Bcast(cbuff, cbuff_len)
      call piernik_MPI_Bcast(rbuff)
      call piernik_MPI_Bcast(ibuff)
      call piernik_MPI_Bcast(lbuff)

      if (slave) then

         cutoff_ncells    = ibuff(1)

         use_inner_orbital_period = lbuff(1)

         mag_field_orient = cbuff(1)
         densfile         = cbuff(2)

         d0               = rbuff(1)
         dout             = rbuff(2)
         r_max            = rbuff(3)
         alpha            = rbuff(4)
         r_in             = rbuff(5)
         r_out            = rbuff(6)
         f_in             = rbuff(7)
         f_out            = rbuff(8)
         dens_exp         = rbuff(9)
         eps              = rbuff(10)
         dens_amb         = rbuff(11)
         x_cut            = rbuff(12)
         dumping_coeff    = rbuff(13)
         drag_max         = rbuff(14)
         drag_min         = rbuff(15)
         amp_noise        = rbuff(16)
         gauss            = rbuff(17:20)
         a_cut            = rbuff(21)
         dens_max         = rbuff(22)

      endif

      if (dom%geometry_type == GEO_RPZ) then ! BEWARE: cannot move this to problem_pointers because dom%geometry_type is set up in init_domain
         user_reg_var_restart       => register_user_var
         problem_customize_solution => problem_customize_solution_kepler
         user_fluidbnd => my_fbnd
         grav_pot_3d => my_grav_pot_3d
         user_vars_hdf5 => prob_vars_hdf5
         problem_grace_passed => add_random_noise
         problem_post_restart => kepler_problem_post_restart
!         problem_grace_passed => add_sine
      endif

   end subroutine read_problem_par
!-----------------------------------------------------------------------------
   subroutine add_sine

      use cg_list,     only: cg_list_element
      use cg_leaves,   only: leaves
      use constants,   only: dpi, xdim, zdim
      use dataio_pub,  only: printinfo, warn
      use domain,      only: dom
      use fluidindex,  only: flind
      use grid_cont,   only: grid_container
      use mpisetup,    only: master

      implicit none

      real :: kx, kz
      integer :: i, k
      type(cg_list_element), pointer :: cgl
      type(grid_container), pointer :: cg

      if (.not. associated(flind%dst)) then
         if (master) call warn("[initproblem:add_sine]: Cannot add sine wave perturbation to dust, because there is no dust.")
         return
      endif
      if (master) call printinfo("[initproblem:add_sine]: adding sine wave perturbation to dust")

      kx = 15.*dpi/dom%L_(xdim)
      kz =  5.*dpi/dom%L_(zdim)

      cgl => leaves%first
      do while (associated(cgl))
         cg => cgl%cg

         do i = cg%is, cg%ie
            do k = cg%ks, cg%ke
               cg%u(flind%dst%imx,i,:,k) = cg%u(flind%dst%imx,i,:,k) + amp_noise*sin(kx*cg%x(i) + kz*cg%z(k)) * cg%u(flind%dst%idn,i,:,k)
               cg%u(flind%dst%imy,i,:,k) = cg%u(flind%dst%imy,i,:,k) + amp_noise*sin(kx*cg%x(i) + kz*cg%z(k)) * cg%u(flind%dst%idn,i,:,k)
               cg%u(flind%dst%imz,i,:,k) = cg%u(flind%dst%imz,i,:,k) + amp_noise*sin(kx*cg%x(i) + kz*cg%z(k)) * cg%u(flind%dst%idn,i,:,k)
            enddo
         enddo
#ifdef DEBUG
         open(456, file="perturbation.dat", status="unknown")
         do k = cg%ks, cg%ke
            write(456,*) cg%z(k), sin(kz*cg%z(k))
         enddo
         close(456)
#endif /* DEBUG */

         cgl => cgl%nxt
      enddo

   end subroutine add_sine
!-----------------------------------------------------------------------------
   subroutine add_random_noise

      use cg_list,     only: cg_list_element
      use cg_leaves,   only: leaves
      use constants,   only: xdim, ydim, zdim
      use dataio_pub,  only: printinfo
      use grid_cont,   only: grid_container
      use fluidindex,  only: flind
      use mpisetup,    only: proc, master

      implicit none

      integer, dimension(:), allocatable    :: seed
      integer                               :: n, clock, i
      real, dimension(:,:,:,:), allocatable :: noise
      type(cg_list_element), pointer        :: cgl
      type(grid_container),  pointer        :: cg

      if (amp_noise <= 0.0) return
      if (master) call printinfo("[initproblem:add_random_noise]: adding random noise to dust")
      call random_seed(size=n)
      allocate(seed(n))
      call system_clock(count=clock)
      seed = clock*proc + 37 * [ (i-1, i = 1, n) ]
      call random_seed(put=seed)
      deallocate(seed)

      cgl => leaves%first
      do while (associated(cgl))
         cg => cgl%cg

         allocate(noise(3,cg%n_(xdim),cg%n_(ydim),cg%n_(zdim)))
         call random_number(noise)
         cg%u(flind%dst%imx,:,:,:) = cg%u(flind%dst%imx,:,:,:) +amp_noise*(1.0-2.0*noise(1,:,:,:)) * cg%u(flind%dst%idn,:,:,:)
         cg%u(flind%dst%imy,:,:,:) = cg%u(flind%dst%imy,:,:,:) +amp_noise*(1.0-2.0*noise(2,:,:,:)) * cg%u(flind%dst%idn,:,:,:)
         cg%u(flind%dst%imz,:,:,:) = cg%u(flind%dst%imz,:,:,:) +amp_noise*(1.0-2.0*noise(3,:,:,:)) * cg%u(flind%dst%idn,:,:,:)
         deallocate(noise)

         cgl => cgl%nxt
      enddo

   end subroutine add_random_noise
!-----------------------------------------------------------------------------
   subroutine init_prob

      use cg_list,          only: cg_list_element
      use cg_leaves,        only: leaves
      use cg_level_connected,    only: base_lev
      use constants,        only: dpi, xdim, ydim, zdim, GEO_XYZ, GEO_RPZ, DST, LO, HI
      use dataio_pub,       only: msg, printinfo, die
      use domain,           only: dom, is_multicg
      use fluidindex,       only: flind
      use fluidtypes,       only: component_fluid
      use gravity,          only: r_smooth, r_grav, n_gravr, ptmass, source_terms_grav, grav_pot2accel, grav_pot_3d
      use grid_cont,        only: grid_container
      use hydrostatic,      only: hydrostatic_zeq_densmid, set_default_hsparams, dprof
      use interactions,     only: epstein_factor
      use mpi,              only: MPI_DOUBLE_PRECISION
      use mpisetup,         only: master, comm, mpi_err, FIRST, proc
      use named_array_list, only: wna
      use units,            only: newtong, gram, cm, kboltz, mH

      implicit none

      integer                         :: i, j, k, kmid, p, middle_of_nx
      integer, dimension(1)           :: n_x_cut
      real                            :: xi, yj, zk, rc, vx, vy, vz, b0, sqr_gm, vr, vphi
      real                            :: csim2, gprim, H2

      real, dimension(:), allocatable :: grav, dens_prof, dens_cutoff, ln_dens_der, gdens
      class(component_fluid), pointer :: fl
      type(cg_list_element),  pointer :: cgl
      type(grid_container),   pointer :: cg

!   Secondary parameters
      call register_user_var

      cgl => leaves%first
      do while (associated(cgl))
         cg => cgl%cg

         if (is_multicg) call die("[initproblem:init_prob] multiple grid pieces per procesor not implemented yet") !nontrivial kmid, allocate

         sqr_gm = sqrt(newtong*ptmass)
         do k = 1, cg%n_(zdim)
            if (cg%z(k) < 0.0) kmid = k       ! the midplane is in between ksmid and ksmid+1
         enddo

         if (associated(flind%ion) .and. dom%geometry_type == GEO_XYZ) then
            fl => flind%ion
            csim2 = fl%cs2*(1.0+alpha)
            b0    = sqrt(2.*alpha*d0*fl%cs2)

            if (dom%has_dir(zdim)) call set_default_hsparams(cg)

            do j = 1, cg%n_(ydim)
               yj = cg%y(j)
               do i = 1, cg%n_(xdim)
                  xi = cg%x(i)
                  rc = sqrt(xi**2+yj**2)

                  if (dom%has_dir(zdim)) call hydrostatic_zeq_densmid(i, j, d0, csim2)

                  do k = 1, cg%n_(zdim)

                     vx = sqr_gm * (-yj)/(rc**2+r_smooth**2)**0.75
                     vy = sqr_gm * ( xi)/(rc**2+r_smooth**2)**0.75
                     vz = 0.0

                     cg%u(fl%idn,i,j,k) = min((rc/r_grav)**n_gravr,100.0)

                     if (dom%has_dir(zdim)) then
                        cg%u(fl%idn,i,j,k) = dprof(k)/cosh(cg%u(fl%idn,i,j,k))
                        cg%u(fl%idn,i,j,k) = max(cg%u(fl%idn,i,j,k), dout)
                     else
                        cg%u(fl%idn,i,j,k) = dout + (d0 - dout)/cosh(cg%u(fl%idn,i,j,k))
                     endif
                     cg%u(fl%idn,i,j,k) = cg%u(fl%idn,i,j,k)
                     cg%u(fl%imx,i,j,k) = vx*cg%u(fl%idn,i,j,k)
                     cg%u(fl%imy,i,j,k) = vy*cg%u(fl%idn,i,j,k)
                     cg%u(fl%imz,i,j,k) = vz*cg%u(fl%idn,i,j,k)
                     if (fl%ien > 0) then
                        cg%u(fl%ien,i,j,k) = fl%cs2/(fl%gam_1)*cg%u(fl%idn,i,j,k)
!                     cg%u(fl%ien,i,j,k) = max(cg%u(fl%ien,i,j,k), smallei)
                        cg%u(fl%ien,i,j,k) = cg%u(fl%ien,i,j,k) +0.5*(vx**2+vy**2+vz**2)*cg%u(fl%idn,i,j,k)
                     endif
                     if (trim(mag_field_orient) == 'toroidal') then
                        cg%b(xdim,i,j,k)   = -b0*sqrt(cg%u(fl%idn,i,j,k)/d0)*yj/rc
                        cg%b(ydim,i,j,k)   =  b0*sqrt(cg%u(fl%idn,i,j,k)/d0)*xi/rc
                        cg%b(zdim,i,j,k)   =  0.0
                     else if (trim(mag_field_orient) == 'vertical') then
                        cg%b(xdim,i,j,k)   =  0.0
                        cg%b(ydim,i,j,k)   =  0.0
                        cg%b(zdim,i,j,k)   =  b0
                     else if (trim(mag_field_orient) == 'none') then
                        cg%b(:,i,j,k)     =  0.0
                     endif

                     if (fl%ien > 0) cg%u(fl%ien,i,j,k)   = cg%u(fl%ien,i,j,k) + 0.5*sum(cg%b(:,i,j,k)**2,1)
                  enddo
               enddo
            enddo
         else if (dom%geometry_type == GEO_RPZ) then
            if (master) then
               call printinfo("------------------------------------------------------------------")
               call printinfo(" Assuming temperature profile for MMSN ")
               call printinfo(" T(R) = 150 ( R / 1 AU )^(-0.429) K")
               write(msg,'(A,F5.1,A)') " T(xmin) = ", mmsn_T(dom%edge(xdim, LO))," K"
               call printinfo(msg)
               write(msg,'(A,F5.1,A)') " T(xmax) = ", mmsn_T(dom%edge(xdim, HI))," K"
               call printinfo(msg)
               write(msg,'(A,F5.1,A)') " T_mean  = ", 0.5*(mmsn_T(dom%edge(xdim, LO))+mmsn_T(dom%edge(xdim, HI)))," K"
               call printinfo(msg)
               write(msg,'(A,F9.5)') " cs2(T_mean) = ", kboltz * 0.5*(mmsn_T(dom%edge(xdim, LO))+mmsn_T(dom%edge(xdim, HI))) / mH
               call printinfo(msg)
               if (associated(flind%neu)) then
                  write(msg,'(A,F12.3,A)') " T_real(cs2) = ", flind%neu%cs2*mH/kboltz, " K"
                  call printinfo(msg)
               endif
               call printinfo("------------------------------------------------------------------")
            endif
            call grav_pot_3d

            if (.not.allocated(grav)) allocate(grav(cg%n_(xdim)))
            if (.not.allocated(ln_dens_der)) allocate(ln_dens_der(cg%n_(xdim)))
            if (.not.allocated(dens_prof)) allocate(dens_prof(cg%n_(xdim)))
            if (.not.allocated(dens_cutoff)) allocate(dens_cutoff(cg%n_(xdim)))
            if (.not.allocated(tauf)) allocate(tauf(cg%n_(xdim)))
            if (.not.allocated(taus)) allocate(taus(cg%n_(xdim))) ! not deallocated

            call source_terms_grav
            call grav_pot2accel(xdim,1,1, cg%n_(xdim), grav, 1, cg)

            dens_prof(:) = d0 * cg%x(:)**(-dens_exp)  * gram / cm**2

            tauf(:) = epstein_factor(flind%neu%pos)/dens_prof(:)

            middle_of_nx = cg%n_(xdim)/2 + 1
            n_x_cut      = maxloc(cg%x, mask=cg%x<=x_cut)
#ifdef FGSL
            if (densfile /= "") then
               allocate(gdens(dom%n_d(xdim)+dom%nb*2))
               if (master) call read_dens_profile(densfile,gdens)
               call MPI_Bcast(gdens, size(gdens), MPI_DOUBLE_PRECISION, FIRST, comm, mpi_err)
               dens_prof(:) = gdens( base_lev%pse(proc)%c(cg%grid_id)%se(xdim, LO)+1:base_lev%pse(proc)%c(cg%grid_id)%se(xdim, HI)+1+dom%nb*2)
               deallocate(gdens)
            endif
#endif /* FGSL */
!              dens_prof = get_lcutoff2(cg%x(:), x_cut, a_cut)
!              dens_prof = dens_prof(:)*(1.0-get_lcutoff2(cg%x(:), x_cut, a_cut)) + dens_max*get_lcutoff2(cg%x(:), x_cut, a_cut)
!              dens_prof    = dens_prof * get_lcutoff(cutoff_ncells, int(middle_of_nx - n_x_cut(1), kind=4), cg%n_(xdim), 0.0, 1.0) + dens_amb

            !! \f$ v_\phi = \sqrt{R\left(c_s^2 \partial_R \ln\rho + \partial_R \Phi \right)} \f$
            ln_dens_der  = log(dens_prof)
            ln_dens_der(2:cg%n_(xdim))  = ( ln_dens_der(2:cg%n_(xdim)) - ln_dens_der(1:cg%n_(xdim)-1) ) / cg%dx
            ln_dens_der(1)        = ln_dens_der(2)
            T_inner               = dpi*cg%x(cg%is) / sqrt( abs(grav(cg%is)) * cg%x(cg%is) )
            write(msg,*) "T_inner = ", T_inner
            if (master) call printinfo(msg)
            write(msg,*) "III Kepler Law gives T = ", sqr_gm/dpi , " yr at 1 AU"
            if (master) call printinfo(msg)
#ifdef DEBUG
            open(143,file="dens_prof.dat",status="unknown")
            do p = 1, cg%n_(xdim)
               write(143,'(4(ES14.4,1X))') cg%x(p), dens_prof(p), sqrt( max(cg%x(p)*(flind%neu%cs2*ln_dens_der(p) + abs(grav(p))),0.0) ), &
                    sqrt( max(abs(grav(p)) * cg%x(p) - flind%neu%cs2*dens_exp,0.0))
            enddo
            close(143)
#endif /* DEBUG */

            do p = 1, flind%fluids
               fl => flind%all_fluids(p)%fl
               if (fl%tag /= DST .and. master) then
                  write(msg,'(A,F9.5)') "[initproblem:initprob] cs2 used = ", fl%cs2
                  call printinfo(msg)
               endif

               do j = 1, cg%n_(ydim)
                  yj = cg%y(j)
                  do i = 1, cg%n_(xdim)
                     xi = cg%x(i)
                     rc = xi + r_smooth

                     gprim = newtong*ptmass / xi**3
                     if (fl%cs > 0) then
                        H2 = 2.0*fl%cs2/gprim
                     else
                        H2 = 1.0
                     endif

                     vphi = 0.
                     do k = 1, cg%n_(zdim)
                        zk = cg%z(k)
!                     cg%u(fl%idn,i,j,k) = max(d0*(1./cosh((xi/r_max)**10)) * exp(-zk**2/H2),1.e-10))
                        cg%u(fl%idn,i,j,k) = dens_prof(i)
                        if (fl%tag == DST) cg%u(fl%idn,i,j,k) = eps * cg%u(fl%idn,i,j,k)

                        vr   = 0.0
                     ! that condition is not necessary since cs2 == 0.0 for dust
                        if (fl%tag /= DST) then
!                        vphi = sqrt( max(abs(grav(i)) * rc - fl%cs2*dens_exp,0.0))
                           vphi = sqrt( max(cg%x(i)*(fl%cs2*ln_dens_der(i) + abs(grav(i))),0.0) )
                        else
                           vphi = sqrt( max(abs(grav(i)) * rc, 0.0))
                        endif
                        vz   = 0.0

                        cg%u(fl%imx,i,j,k) = vr   * cg%u(fl%idn,i,j,k)
                        cg%u(fl%imy,i,j,k) = vphi * cg%u(fl%idn,i,j,k)
                        cg%u(fl%imz,i,j,k) = vz   * cg%u(fl%idn,i,j,k)
                        if (fl%ien > 0) then
                           cg%u(fl%ien,i,j,k) = fl%cs2/(fl%gam_1)*cg%u(fl%idn,i,j,k)
                           cg%u(fl%ien,i,j,k) = cg%u(fl%ien,i,j,k) + 0.5*(vr**2+vphi**2+vz**2)*cg%u(fl%idn,i,j,k)
                        endif
                     enddo
                     taus(i) = vphi/cg%x(i)*tauf(i) ! compiler complains that vphi may be used uninitialized here
                  enddo
               enddo

            enddo
            cg%w(wna%ind(inid_n))%arr(:,:,:,:) = cg%u(:,:,:,:)
            cg%b(:,:,:,:) = 0.0
            if (allocated(grav)) deallocate(grav)
            if (allocated(dens_prof)) deallocate(dens_prof)
#ifdef DEBUG
            open(123,file="tau.dat",status="unknown")
            do i = 1, cg%n_(xdim)
               write(123,*) cg%x(i), tauf(i), taus(i)
            enddo
            close(123)
#endif /* DEBUG */
         else
            call die("[initproblem:init_prob] I don't know what to do... :/")
         endif

         cgl => cgl%nxt
      enddo

   end subroutine init_prob
!-----------------------------------------------------------------------------
   real function mmsn_T(r)
      implicit none
      real, intent(in) :: r         ! [AU]
      real, parameter  :: T_0 = 150 ! [K]
      real, parameter  :: k   = 0.429

      mmsn_T = T_0 * r**(-k)
   end function mmsn_T
!-----------------------------------------------------------------------------
   subroutine kepler_problem_post_restart

      use cg_list,          only: cg_list_element
      use cg_leaves,        only: leaves
      use constants,        only: b0_n
      use fluidboundaries,  only: all_fluid_boundaries
      use named_array_list, only: wna
#ifdef TRACER
      use constants,        only: xdim, ydim, zdim
      use grid_cont,        only: grid_container
      use func,             only: resample_gauss
      use fluidindex,       only: flind
#endif /* TRACER */

      implicit none

      type(cg_list_element), pointer :: cgl
#ifdef TRACER
      type(grid_container), pointer :: cg
      integer :: i, j, k
#endif /* TRACER */

      call my_grav_pot_3d   ! reset gp, to get right values on the boundaries

      cgl => leaves%first
      do while (associated(cgl))
         cgl%cg%u  => cgl%cg%w(wna%ind(inid_n))%arr  ! BEWARE: Don't do things like that without parental supervision
         cgl%cg%b = 0.0
         cgl%cg%w(wna%ind(b0_n))%arr = 0.0
         cgl => cgl%nxt
      enddo

      call all_fluid_boundaries   ! all_fluid_boundaries properly set boundaries for %u pointer

      cgl => leaves%first
      do while (associated(cgl))
         cgl%cg%u  => cgl%cg%w(wna%fi)%arr ! Quick! Revert to sane state before anyone notices
         cgl => cgl%nxt
      enddo

#ifdef TRACER
      if (gauss(4) > 0.0) then
         cgl => leaves%first
         do while (associated(cgl))
            cg => cgl%cg

            do k = 1, cg%n_(zdim)
               do j = 1, cg%n_(ydim)
                  do i = 1, cg%n_(xdim)

                     cg%u(flind%trc%beg:flind%trc%end, i, j, k)   = &
                          resample_gauss( cg%x(i) - gauss(1), cg%y(j) - gauss(2), cg%z(k) - gauss(3), &
                                          cg%dl(xdim), cg%dl(ydim), cg%dl(zdim), &
                                          gauss(4), gauss(4), gauss(4), 10)

                  enddo
               enddo
            enddo
            cgl => cgl%nxt
         enddo
      endif

      call all_fluid_boundaries
#endif /* TRACER */

   end subroutine kepler_problem_post_restart
!-----------------------------------------------------------------------------
   subroutine problem_customize_solution_kepler(forward)

      use cg_list,          only: cg_list_element
      use cg_leaves,        only: leaves
      use constants,        only: xdim, ydim, zdim, I_ONE
      use dataio_pub,       only: die!, warn, msg
      use domain,           only: is_multicg
      use global,           only: dt, relax_time, smalld !, t, grace_period_passed
      use grid_cont,        only: grid_container
      use fluidboundaries,  only: all_fluid_boundaries
      use fluidindex,       only: flind!, iarr_all_mz, iarr_all_dn
      use mpisetup,         only: comm, mpi_err
      use mpi,              only: MPI_MAX, MPI_DOUBLE_PRECISION, MPI_IN_PLACE
      use named_array_list, only: wna
      ! use interactions,    only: dragc_gas_dust
#ifdef VERBOSE
!      use dataio_pub,       only: msg, printinfo
!      use mpisetup,         only: master
#endif /* VERBOSE */

      implicit none

      logical, intent(in)                     :: forward
      integer                                 :: i, j, k
      logical, save                           :: frun = .true.
      real, dimension(:,:), allocatable, save :: funcR
      logical, dimension(:,:,:), allocatable  :: adjust
      real, dimension(:,:,:), allocatable     :: vx_sign, vz_sign
      real, save                              :: x0, x1, y0, y1, a, b
      real                                    :: max_vx, mean_vy
      real, save                              :: max_vy = 100.0
      type(cg_list_element), pointer          :: cgl
      type(grid_container),  pointer          :: cg

!      if (grace_period_passed() .and. t <= x1) then
!         dragc_gas_dust = a*t+b
!#ifdef VERBOSE
!         if (master) write(msg,'(A,F6.1)') 'dragc_gas_dust = ', dragc_gas_dust
!         call printinfo(msg)
!#endif /* VERBOSE */
!      endif

      if (is_multicg) call die("[initproblem:problem_customize_solution_kepler] multiple grid pieces per procesor not implemented yet") !nontrivial

      cgl => leaves%first
      do while (associated(cgl))
         cg => cgl%cg
         if (.not.allocated(adjust)) allocate(adjust(cg%n_(xdim),cg%n_(ydim),cg%n_(zdim)))
         if (.not.allocated(vx_sign)) allocate(vx_sign(cg%n_(xdim),cg%n_(ydim),cg%n_(zdim)))
         if (.not.allocated(vz_sign)) allocate(vz_sign(cg%n_(xdim),cg%n_(ydim),cg%n_(zdim)))

         if (frun) then
            x0 = relax_time + 2.0
            x1 = x0 + 30.0
            y0 = drag_max
            y1 = drag_min
            a = (y0 - y1)/(x0 - x1)
            b = y0 - a*x0
            allocate(funcR(size(cg%u,dim=1), cg%n_(xdim)) )

            funcR(1,:) = -tanh((cg%x(:)-r_in+1.0)**f_in) + 1.0 + max( tanh((cg%x(:)-r_out+1.0)**f_out), 0.0)

            if (use_inner_orbital_period) then
               funcR(1,:) = funcR(1,:) / T_inner
            else
               funcR(1,:) = funcR(1,:) * dumping_coeff
            endif
#ifdef DEBUG
            open(212,file="funcR.dat",status="unknown")
            do j = 1, cg%n_(xdim)
               write(212,*) cg%x(j),funcR(1,j)
            enddo
            close(212)
#endif /* DEBUG */
            frun = .false.
            funcR(:,:) = spread(funcR(1,:),1,size(cg%u,dim=1))

            max_vy = maxval( abs(cg%u(flind%dst%imy,:,:,:))/cg%u(flind%dst%idn,:,:,:) )
            call MPI_Allreduce(MPI_IN_PLACE, max_vy, I_ONE, MPI_DOUBLE_PRECISION, MPI_MAX, comm, mpi_err)
         endif

         do j = 1, cg%n_(ydim)
            do k = 1, cg%n_(zdim)
               cg%u(:,:,j,k) = cg%u(:,:,j,k) - dt*(cg%u(:,:,j,k) - cg%w(wna%ind(inid_n))%arr(:,:,j,k))*funcR(:,:)
            enddo
         enddo

!         where ( cg%u(iarr_all_dn,:,:,:) < 2.0*smalld )
!            cg%u(iarr_all_mz,:,:,:) = cg%u(iarr_all_mz,:,:,:)*0.1
!         endwhere

         max_vx = maxval( abs(cg%u(flind%neu%imx,:,:,:))/cg%u(flind%neu%idn,:,:,:) )
         call MPI_Allreduce(MPI_IN_PLACE, max_vx, I_ONE, MPI_DOUBLE_PRECISION, MPI_MAX, comm, mpi_err)

         adjust = abs(cg%u(flind%dst%imx,:,:,:))/cg%u(flind%dst%idn,:,:,:) >= max_vx
         if ( any(adjust) ) then
!            write(msg,'(a,i6,a)') "[kepler_customize_solution]: ", count(adjust), &
!                " cells need adjustment"
!            call warn(msg)
            where (adjust)
               vx_sign = signum(cg%u(flind%dst%imx,:,:,:))
               vz_sign = signum(cg%u(flind%dst%imz,:,:,:))
            endwhere
            where (adjust)
               cg%u(flind%dst%idn,:,:,:) = max(cg%u(flind%dst%idn,:,:,:), 1.1*smalld)
               cg%u(flind%dst%imx,:,:,:) = vx_sign * cg%u(flind%neu%imx,:,:,:)/cg%u(flind%neu%idn,:,:,:) * cg%u(flind%dst%idn,:,:,:)
               cg%u(flind%dst%imz,:,:,:) = vz_sign * cg%u(flind%neu%imz,:,:,:)/cg%u(flind%neu%idn,:,:,:) * cg%u(flind%dst%idn,:,:,:)
            endwhere
         endif

         do i = 1, cg%n_(xdim)
            do j = 1, cg%n_(ydim)
               mean_vy = sum( cg%u(flind%dst%imy, i, j, :) /  cg%u(flind%dst%idn, i, j, :) ) / cg%n_(zdim)
               do k = 1, cg%n_(zdim)
                  if ( (abs(cg%u(flind%dst%imy, i, j, k) - mean_vy * cg%u(flind%dst%idn, i, j, k)) &
                     / (mean_vy * cg%u(flind%dst%idn, i, j, k)) >= 0.1) .and. &
                         cg%u(flind%dst%idn, i, j, k) < 10.0*smalld ) then
                     cg%u(flind%dst%idn, i ,j ,k) = max(cg%u(flind%dst%idn, i, j, k), 1.1*smalld)
                     cg%u(flind%dst%imy, i, j, k) = mean_vy * cg%u(flind%dst%idn, i, j, k)
                  endif
               enddo
            enddo
         enddo

!        where ( cg%u(flind%dst%imy,:,:,:) > max_vy*cg%u(flind%dst%idn,:,:,:) )
!           cg%u(flind%dst%idn,:,:,:) = max(cg%u(flind%dst%idn,:,:,:), 1.1*smalld)
!           cg%u(flind%dst%imy,:,:,:) = cg%u(flind%neu%imy,:,:,:)/cg%u(flind%neu%idn,:,:,:) * cg%u(flind%dst%idn,:,:,:)
!        endwhere

         if (allocated(adjust)) deallocate(adjust)
         if (allocated(vx_sign)) deallocate(vx_sign)
         if (allocated(vz_sign)) deallocate(vz_sign)

         cgl => cgl%nxt
      enddo

      call all_fluid_boundaries

      return
      if (forward) i = j ! suppress compiler warnings on unused arguments

   end subroutine problem_customize_solution_kepler
!-----------------------------------------------------------------------------
   subroutine my_grav_pot_3d

      use cg_list,     only: cg_list_element
      use cg_leaves,   only: leaves
      use constants,   only: xdim, zdim
      use gravity,     only: ptmass, sum_potential
      use grid_cont,   only: grid_container
      use units,       only: newtong

      implicit none

      logical, save                  :: frun = .true.
      real                           :: r2
      integer                        :: i, k
      type(cg_list_element), pointer :: cgl
      type(grid_container),  pointer :: cg

      if (frun) then
         cgl => leaves%first
         do while (associated(cgl))
            cg => cgl%cg

            do i = 1, cg%n_(xdim)
               do k = 1, cg%n_(zdim)
                  r2 = cg%x(i)**2! + cg%z(k)**2
                  cg%gp(i,:,k) = -newtong*ptmass / sqrt(r2)
               enddo
            enddo

            cgl => cgl%nxt
         enddo
      endif

      frun = .false.
      call sum_potential

   end subroutine my_grav_pot_3d
!-----------------------------------------------------------------------------
   subroutine my_fbnd(dir,side,cg)

      use constants,  only: xdim, LO
      use grid_cont,  only: grid_container

      implicit none

      integer(kind=4),               intent(in)    :: dir, side
      type(grid_container), pointer, intent(inout) :: cg

      if (dir == xdim) then
         if (side == LO) then
            call my_bnd_xl(cg)
         else
            call my_bnd_xr(cg)
         endif
      endif

   end subroutine my_fbnd
!-----------------------------------------------------------------------------
   subroutine my_bnd_xl(cg)

      use constants,  only: xdim, ydim, zdim
      use domain,     only: dom
      use gravity,    only: grav_pot2accel
      use grid_cont,  only: grid_container
      use fluidindex, only: iarr_all_dn, iarr_all_mx, iarr_all_my, iarr_all_mz, flind
#ifndef ISO
      use fluidindex, only: iarr_all_en
#endif /* ISO */

      implicit none

      type(grid_container), pointer, intent(inout) :: cg

      integer :: i
      real, dimension(cg%n_(xdim)) :: grav
      real, dimension(size(iarr_all_my), cg%n_(ydim), cg%n_(zdim)) :: vy,vym
      real, dimension(size(flind%all_fluids))    :: cs2_arr
      integer, dimension(size(flind%all_fluids)) :: ind_cs2

      do i = 1, size(flind%all_fluids)
         ind_cs2    = i
         cs2_arr(i) = flind%all_fluids(i)%fl%cs2
      enddo

      call grav_pot2accel(xdim,1,1, cg%n_(xdim), grav, 1, cg)

      do i = 1, dom%nb
         cg%u(iarr_all_dn,i,:,:) = cg%u(iarr_all_dn, cg%is,:,:)
         cg%u(iarr_all_mx,i,:,:) = min(0.0,cg%u(iarr_all_mx, cg%is,:,:))
!         do p = 1, size(flind%all_fluids)
!            cg%u(iarr_all_my(p),i,:,:) = sqrt( abs(grav(i)) * cg%x(i) - cs2_arr(p)*dens_exp) *  cg%u(iarr_all_dn(p),i,:,:)
!         enddo
         cg%u(iarr_all_my,i,:,:) = cg%u(iarr_all_my, cg%is,:,:)
         cg%u(iarr_all_mz,i,:,:) = cg%u(iarr_all_mz, cg%is,:,:)
#ifndef ISO
         cg%u(iarr_all_en,i,:,:) = cg%u(iarr_all_en, cg%is,:,:)
#endif /* !ISO */
      enddo

      do i = dom%nb,1,-1
         vym(:,:,:) = cg%u(iarr_all_my,i+2,:,:)/cg%u(iarr_all_dn,i+1,:,:)
         vy(:,:,:)  = cg%u(iarr_all_my,i+1,:,:)/cg%u(iarr_all_dn,i+1,:,:)
!         cg%u(iarr_all_my,i,:,:) = (vym(:,:,:) + (cg%x(i) - cg%x(i+2)) / (cg%x(i+1) - cg%x(i+2)) * (vy - vym))*cg%u(iarr_all_dn,i,:,:)
      enddo

   end subroutine my_bnd_xl
!-----------------------------------------------------------------------------
   subroutine my_bnd_xr(cg)

      use constants,        only: xdim
      use grid_cont,        only: grid_container
      use named_array_list, only: wna

      implicit none

      type(grid_container), pointer, intent(inout) :: cg

      cg%u(:, cg%ie+1:cg%n_(xdim),:,:) = cg%w(wna%ind(inid_n))%arr(:,cg%ie+1:cg%n_(xdim),:,:)
   end subroutine my_bnd_xr
!-----------------------------------------------------------------------------
   function get_lcutoff(width, dist, n, vmin, vmax) result(y)

      implicit none

      integer(kind=4), intent(in) :: width  !< width of tanh profile [cells]
      integer(kind=4), intent(in) :: dist   !< distance between the expected position of the profile and the middle cell of the domain [cells]
      integer(kind=4), intent(in) :: n      !< length of the profile array
      real,            intent(in) :: vmin   !< minimal value of the profile
      real,            intent(in) :: vmax   !< maximum value of the profile
      real, dimension(n)  :: y, x

      real, parameter     :: kstep = 1.0 !< iteration step for tanh fit
      real                :: dv, k
      integer             :: nn, i

      x = [(real(i), i=0,n-1)] / real(n) * 10.0 - 5.0

      dv = vmax - vmin
      nn = huge(1) ; k = 0.0

      do while (width < nn)
         k = k + kstep
         nn = get_ncells(x,k)
      enddo

      y = 0.5*dv*(tanh(x*k) + 1.0) + vmin

      if (dist < 0) then
         y = eoshift(y,dim=1,shift=dist-width/2,boundary=vmin)
      else
         y = eoshift(y,dim=1,shift=dist+width/2,boundary=vmax)
      endif
   end function get_lcutoff

   function get_lcutoff2(x, x0, a) result (y)
      use constants,   only: pi
      implicit none
      real, intent(in), dimension(:) :: x
      real, intent(in) :: x0, a
      real, dimension(size(x)) :: y

      y = atan(-(x-x0)*a)/pi + 0.5

   end function get_lcutoff2
!-----------------------------------------------------------------------------
   integer function get_ncells(x,k)
      implicit none
      real, intent(in) :: k
      real, intent(in), dimension(:) :: x
      real, dimension(size(x))       :: y
      y = tanh(x*k)
      get_ncells = count(y > -0.99 .and. y < 0.99)
   end function get_ncells
!-----------------------------------------------------------------------------
   subroutine prob_vars_hdf5(var,tab, ierrh, cg)

      use grid_cont,  only: grid_container
      use interactions, only: epstein_factor
      use fluidindex,   only: flind

      implicit none

      character(len=*), intent(in)                    :: var
      real(kind=4), dimension(:,:,:), intent(inout)   :: tab
      integer, intent(inout)                          :: ierrh
      type(grid_container), pointer, intent(in)       :: cg

      ierrh = 0
      select case (trim(var))
         case ("tauf")
            tab(:,:,:) = real(epstein_factor(flind%neu%pos) / cg%u(flind%neu%idn,cg%is:cg%ie,cg%js:cg%je,cg%ks:cg%ke), 4)
         case default
            ierrh = -1
      end select

   end subroutine prob_vars_hdf5
!-----------------------------------------------------------------------------
#ifdef FGSL
   subroutine read_dens_profile(densfile,gdens)

      use dataio_pub, only: printinfo, msg, warn
      use domain,     only: dom
      use fgsl,       only: fgsl_size_t, fgsl_interp_accel, fgsl_interp, fgsl_int, fgsl_char, fgsl_strmax, fgsl_interp_cspline, &
           &                fgsl_interp_accel_alloc, fgsl_interp_alloc, fgsl_interp_name, fgsl_interp_init, fgsl_interp_eval, fgsl_interp_free, fgsl_interp_accel_free
      !, fgsl_spline
      use domain,     only: dom

      implicit none

      character(len=*), intent(in)                   :: densfile
      real, dimension(:), intent(out)                :: gdens

      type(fgsl_interp_accel) :: acc
      type(fgsl_interp) :: a_interp
      integer(fgsl_int) :: fstatus
      character(kind=fgsl_char,len=fgsl_strmax) :: fname
      real, dimension(:), allocatable :: x,y
      real                            :: xi
      integer(fgsl_size_t)            :: n, nmax
      integer(kind=4) :: i, nxd

      write(msg,*) "[initproblem:read_dens_profile] Reading ", trim(densfile)
      open(1,file=densfile, status="old", form='unformatted')
         read(1) i
         n = i
         allocate(y(n), x(n))
         read(1) x
         read(1) y
      close(1)

      nmax = n
      nxd = int(size(gdens) - 2*dom%nb, kind=4)

      if (nxd == n) then
         call printinfo("[initproblem:read_dens_profile] Saved profile has required dimension \o/")
         gdens(dom%nb+1:dom%nb+nxd) = y(:)
      else
         call warn("[initproblem:read_dens_profile] Saved profile has different dimension :/")
         write(msg,'(A,I5,A,I5,A)') "[initproblem:read_dens_profile] Performing spline interpolation from",n," to ",nxd," cells."
         call printinfo(msg)
         acc = fgsl_interp_accel_alloc()
         a_interp = fgsl_interp_alloc(fgsl_interp_cspline,nmax)
         fname = fgsl_interp_name(a_interp)
         fstatus = fgsl_interp_init(a_interp,x,y,nmax)
         do i = 1, nxd
            xi = (1-i /real(nxd))*x(1) + (i/real(nxd))*x(n)
            gdens(dom%nb+i) = fgsl_interp_eval(a_interp, x, y, xi, acc)
         enddo

         call fgsl_interp_free(a_interp)
         call fgsl_interp_accel_free(acc)
      endif

      do i = 1, dom%nb
         gdens(i)           = gdens(dom%nb+1)
         gdens(nxd+dom%nb+i) = gdens(nxd+dom%nb)
      enddo
      deallocate(x,y)
      return
   end subroutine read_dens_profile
#endif /* FGSL */
!-----------------------------------------------------------------------------
   elemental function signum(a) result (b)
      implicit none
      real, intent(in) :: a
      real :: b
      if (a > 0.0) then
         b = 1.0
      else
         b = -1.0
      endif
   end function signum
end module initproblem
! vim: set tw=120:
