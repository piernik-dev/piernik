! $Id: initproblem.F90 4545 2011-06-21 13:55:05Z gawrysz $
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

   real                     :: rho0, cs0, R0,  r_max, dout, alpha, r_in, r_out, f_in, f_out
   real                     :: growth_time         !< time after which particles will reach initproblem::final_grain_size
   real                     :: initial_grain_size   !< particle size after global::relax_time
   real                     :: final_grain_size    !< particle size after initproblem::growth_time
   real                     :: eps                 !< dust to gas ratio
   !>
   !! \f$\tau\f$ in \f$\frac{Du}{Dt} = - \frac{u-u_0}{\tau}f(R)
   !! when initproblem::problem_customize_solution is used
   !<
   real                     :: dumping_coeff, amp_noise
   logical                  :: use_inner_orbital_period  !< use 1./T_inner as dumping_coeff
   character(len=cbuff_len) :: densfile
   character(len=dsetnamelen), parameter :: inid_n = "u_0"

   namelist /PROBLEM_CONTROL/  rho0, cs0, R0, r_in, r_out, f_in, f_out, eps, dumping_coeff, &
      &  use_inner_orbital_period, amp_noise, growth_time, initial_grain_size, final_grain_size

contains
!-----------------------------------------------------------------------------
   subroutine problem_pointers

      use dataio_user,         only: user_reg_var_restart, user_vars_hdf5
      use gravity,             only: grav_pot_3d, grav_type
!     use fluidboundaries_funcs, only: user_fluidbnd
      use user_hooks,          only: problem_customize_solution, problem_grace_passed

      implicit none

      user_reg_var_restart       => register_user_var
      problem_customize_solution => problem_customize_solution_kepler
!     user_fluidbnd => my_fbnd
      grav_pot_3d => my_grav_pot_3d
      user_vars_hdf5 => prob_vars_hdf5
      problem_grace_passed => add_random_noise
      grav_type => my_grav_ptmass_pure
!     problem_grace_passed => add_sine

   end subroutine problem_pointers
!-----------------------------------------------------------------------------
   subroutine read_problem_par

      use dataio_pub,          only: ierrh, par_file, namelist_errh, compare_namelist, cmdl_nml, lun      ! QA_WARN required for diff_nml
      use mpisetup,            only: rbuff, lbuff, master, slave, piernik_MPI_Bcast

      implicit none

      rho0             = 1.0
      R0               = 1.0e-4
      cs0              = 1.0
      amp_noise        = 1.e-6

      r_in             = 0.5
      f_in             = 10.0
      r_out            = 2.1
      f_out            = 0.0

      eps              = 1.0
      dumping_coeff    = 0.0

      use_inner_orbital_period = .true.

      if (master) then

         diff_nml(PROBLEM_CONTROL)

         lbuff(1) = use_inner_orbital_period

         rbuff(1)  = rho0
         rbuff(2)  = R0
         rbuff(3)  = cs0
         rbuff(4)  = r_in
         rbuff(5)  = r_out
         rbuff(6)  = f_in
         rbuff(7)  = f_out
         rbuff(8)  = eps
         rbuff(9)  = amp_noise
         rbuff(10) = dumping_coeff
         rbuff(11) = growth_time
         rbuff(12) = initial_grain_size
         rbuff(13) = final_grain_size

      endif

      call piernik_MPI_Bcast(rbuff)
      call piernik_MPI_Bcast(lbuff)

      if (slave) then

         use_inner_orbital_period = lbuff(1)

         rho0               = rbuff(1)
         R0                 = rbuff(2)
         cs0                = rbuff(3)
         r_in               = rbuff(4)
         r_out              = rbuff(5)
         f_in               = rbuff(6)
         f_out              = rbuff(7)
         eps                = rbuff(8)
         amp_noise          = rbuff(9)
         dumping_coeff      = rbuff(10)
         growth_time        = rbuff(11)
         initial_grain_size = rbuff(12)
         final_grain_size   = rbuff(13)

      endif

   end subroutine read_problem_par
!-----------------------------------------------------------------------------
   subroutine register_user_var

      use cg_list_global, only: all_cg
      use constants,      only: AT_NO_B
      use named_array_list, only: wna

      implicit none

      integer(kind=4) :: dim4

      dim4 = wna%lst(wna%fi)%dim4
      call all_cg%reg_var(inid_n, restart_mode = AT_NO_B, dim4 = dim4)

   end subroutine register_user_var
!-----------------------------------------------------------------------------
   subroutine add_sine

      use constants,   only: dpi, xdim, zdim
      use dataio_pub,  only: printinfo, warn
      use domain,      only: dom
      use fluidindex,  only: flind
      use cg_list_bnd, only: leaves
      use cg_list,     only: cg_list_element
      use grid_cont,   only: grid_container
      use mpisetup,    only: master

      implicit none

      real, parameter :: amp = 1.e-6
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
               cg%u(flind%dst%imx,i,:,k) = cg%u(flind%dst%imx,i,:,k) + amp*sin(kx*cg%x(i) + kz*cg%z(k)) * cg%u(flind%dst%idn,i,:,k)
               cg%u(flind%dst%imy,i,:,k) = cg%u(flind%dst%imy,i,:,k) + amp*sin(kx*cg%x(i) + kz*cg%z(k)) * cg%u(flind%dst%idn,i,:,k)
               cg%u(flind%dst%imz,i,:,k) = cg%u(flind%dst%imz,i,:,k) + amp*sin(kx*cg%x(i) + kz*cg%z(k)) * cg%u(flind%dst%idn,i,:,k)
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

      use constants,   only: xdim, ydim, zdim
      use dataio_pub,  only: printinfo
      use cg_list_bnd, only: leaves
      use cg_list,     only: cg_list_element
      use grid_cont,   only: grid_container
      use fluidindex,  only: flind
      use mpisetup,    only: proc, master

      implicit none

      integer, dimension(:), allocatable  :: seed
      integer :: n, clock, i
      real, dimension(:,:,:,:), allocatable :: noise
      type(cg_list_element), pointer :: cgl
      type(grid_container), pointer :: cg

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
         cg%u(flind%dst%imx,:,:,:) = cg%u(flind%dst%imx,:,:,:) +amp_noise*noise(1,:,:,:) * cg%u(flind%dst%idn,:,:,:)
         cg%u(flind%dst%imy,:,:,:) = cg%u(flind%dst%imy,:,:,:) +amp_noise*noise(2,:,:,:) * cg%u(flind%dst%idn,:,:,:)
         cg%u(flind%dst%imz,:,:,:) = cg%u(flind%dst%imz,:,:,:) +amp_noise*noise(3,:,:,:) * cg%u(flind%dst%idn,:,:,:)
         deallocate(noise)

         cgl => cgl%nxt
      enddo

   end subroutine add_random_noise
!-----------------------------------------------------------------------------
   subroutine init_prob

      use cg_list,       only: cg_list_element
      use cg_list_bnd,   only: leaves
      use constants,     only: DST, GEO_RPZ, xdim, ydim, zdim
      use dataio_pub,    only: msg, printinfo, die
      use cart_comm,     only: cdd
      use domain,        only: dom, is_multicg
      use fluidindex,    only: flind
      use fluidtypes,    only: component_fluid
      use func,          only: ekin
      use global,        only: smalld
      use gravity,       only: ptmass, grav_pot2accel
      use grid_cont,     only: grid_container
      use hydrostatic  , only: hydrostatic_zeq_densmid, set_default_hsparams, dprof
      use mpi,           only: MPI_COMM_NULL
      use mpisetup,      only: master
      use named_array  , only: wna
      use units,         only: newtong

      implicit none

      integer :: i, j, k, p
      real    :: zk, R, vz, sqr_gm, vr
      real    :: rho, cs2
      real, dimension(:), allocatable :: vphi, grav, ln_dens_der

      class(component_fluid), pointer  :: fl
      type(cg_list_element), pointer :: cgl
      type(grid_container), pointer :: cg

!   Secondary parameters
      if (cdd%comm3d == MPI_COMM_NULL) call die("[initproblem:init_prob] comm3d == MPI_COMM_NULL not implemented") !pcoords
      if (dom%geometry_type /= GEO_RPZ) call die("[initproblem:init_prob] only cylindrical geometry supported")

      sqr_gm = sqrt(newtong*ptmass)

      call register_user_var

      cgl => leaves%first
      do while (associated(cgl))
         cg => cgl%cg

         if (is_multicg) call die("[initproblem:init_prob] multiple grid pieces per procesor not implemented yet") !nontrivial kmid, allocate

         do p = 1, flind%fluids
            fl => flind%all_fluids(p)%fl
            if (fl%tag /= DST .and. master) then
               write(msg,'(A,F9.5)') "[init_problem:initprob] cs2 used = ", fl%cs2
               call printinfo(msg)
            endif

            call set_default_hsparams(cg)

            do i = 1, cg%n_(xdim)
               R = cg%x(i)

               rho = max(rho0*(R0/R)**(1.75),smalld)   ! bylo 1.5
               cs2 = (cs0*sqrt(R0/R))**2

               cg%cs_iso2(i,:,:) = cs2

               do j = 1, cg%n_(ydim)
                  call hydrostatic_zeq_densmid(i, j, rho, cs2)
                  do k = 1, cg%n_(zdim)
                     zk = cg%z(k)
                     cg%u(fl%idn,i,j,k) = max(dprof(k)/(1.0+eps), smalld)
                     if (fl%tag == DST) cg%u(fl%idn,i,j,k) = max(eps * cg%u(flind%neu%idn,i,j,k), smalld)

                  enddo
               enddo
            enddo

            allocate(ln_dens_der(cg%n_(xdim)), grav(cg%n_(xdim)), vphi(cg%n_(xdim)))

            do j = 1, cg%n_(ydim)
               do k = 1, cg%n_(zdim)


                  ln_dens_der  = log(cg%u(fl%idn,:,j,k))
                  if (dom%has_dir(xdim)) then
                     call grav_pot2accel(xdim, j, k, cg%n_(xdim), grav, 1, cg)
                     ln_dens_der(2:cg%n_(xdim))  = ( ln_dens_der(2:cg%n_(xdim)) - ln_dens_der(1:cg%n_(xdim)-1) ) / cg%dx
                     ln_dens_der(1)        = ln_dens_der(2)
                     if (fl%tag /= DST) then
                        vphi = sqrt( max(cg%x(:)*(cg%cs_iso2(:,j,k)*ln_dens_der(:) + abs(grav(:))),0.0) )
                     else
                        vphi = sqrt( max(abs(grav(:)) * cg%x(:), 0.0))
                     endif
                  else
                     vphi(:) = sqrt( newtong*ptmass / cg%x(:) )
                  endif

                  vr = 0.0
                  vz = 0.0

                  cg%u(fl%imx,:,j,k) = vr     *cg%u(fl%idn,:,j,k)
                  cg%u(fl%imy,:,j,k) = vphi(:)*cg%u(fl%idn,:,j,k)
                  cg%u(fl%imz,:,j,k) = vz     *cg%u(fl%idn,:,j,k)

               enddo
            enddo
            deallocate(ln_dens_der, grav, vphi)
            if (fl%has_energy) then
               cg%u(fl%ien,:,:,:) = fl%cs2/(fl%gam_1)*cg%u(fl%idn,:,:,:)
               cg%u(fl%ien,:,:,:) = cg%u(fl%ien,:,:,:) + ekin( cg%u(fl%imx,:,:,:), cg%u(fl%imy,:,:,:), &
                  & cg%u(fl%imz,:,:,:), cg%u(fl%idn,:,:,:))
            endif
         enddo

         cg%w(wna%ind(inid_n))%arr(:,:,:,:) = cg%u(:,:,:,:)
         cg%b(:,:,:,:) = 0.0
         cgl => cgl%nxt
      enddo
      open(12,file="vel_profile.dat",status="unknown")
         do i = 1, cg%n_(xdim)
            write(12,'(3(E12.5,1X))') cg%x(i), cg%u(flind%all_fluids(1:2)%fl%imy,i,max(cg%n_(ydim)/2,1),max(cg%n_(zdim)/2,1)) / &
                &  cg%u(flind%all_fluids(1:2)%fl%idn,i,max(cg%n_(ydim)/2,1),max(cg%n_(zdim)/2,1))
         enddo
      close(12)

   end subroutine init_prob
!-----------------------------------------------------------------------------
   subroutine problem_customize_solution_kepler

      use constants,       only: dpi, xdim, ydim, zdim
      use dataio_pub,      only: die
      use domain,          only: is_multicg, dom
      use fluidboundaries, only: all_fluid_boundaries
      use fluidindex,      only: iarr_all_dn, iarr_all_mz
      use cg_list,         only: cg_list_element
      use global,          only: t, grace_period_passed, relax_time, smalld !, dt
      use gravity,         only: ptmass
      use cg_list_bnd,     only: leaves
      use grid_cont,       only: grid_container
      use interactions,    only: update_grain_size
      use named_array_list, only: wna
      use units,           only: newtong
#ifdef VERBOSE
!      use dataio_pub,      only: msg, printinfo
!      use mpisetup,        only: master
#endif /* VERBOSE */

      implicit none

      integer                               :: j, k
      logical, save                         :: frun = .true.
      real, dimension(:,:), allocatable, save :: funcR
      real, save :: x0, x1, y0, y1, a, b
      type(cg_list_element), pointer :: cgl
      type(grid_container), pointer :: cg

      if (is_multicg) call die("[initproblem:problem_customize_solution_kepler] multiple grid pieces per procesor not implemented yet") !nontrivial

      cgl => leaves%first
      do while (associated(cgl))
         cg => cgl%cg

         if (frun) then
            x0 = relax_time
            x1 = x0 + growth_time
            y0 = initial_grain_size
            y1 = final_grain_size
            a = (y0 - y1)/(x0 - x1)
            b = y0 - a*x0

            allocate(funcR(size(cg%u,dim=1), cg%n_(xdim)) )     ! BEWARE: will fail with multiple cgs

            funcR(1,:) = 0.0
            funcR(1,:) = funcR(1,:) - tanh((cg%x(:)-r_in+1.0)**f_in) + 1.0
            where (cg%x >= r_out)
               funcR(1,:) = funcR(1,:) + max( tanh((cg%x(:)-r_out+1.0)**f_out), 0.0)
            endwhere

            if (use_inner_orbital_period) then
               funcR(1,:) = funcR(1,:) * sqrt( newtong*ptmass/cg%x(dom%nb)**3 ) / dpi
            else
               funcR(1,:) = funcR(1,:) * dumping_coeff
            endif
            open(212,file="funcR.dat",status="unknown")
            write(212,*) "# ",sqrt( newtong*ptmass/cg%x(dom%nb)**3 ) / dpi
            do j = 1, cg%n_(xdim)
               write(212,*) cg%x(j),funcR(1,j)
            enddo
            close(212)
            frun = .false.
            funcR(:,:) = spread(funcR(1,:),1,size(cg%u,dim=1))
         endif

         if (grace_period_passed()) call update_grain_size(a*t+b)
!         do j = 1, cg%n_(ydim)
!            do k = 1, cg%n_(zdim)
!               cg%u(iarr_all_dn,:,j,k) = cg%u(iarr_all_dn,:,j,k) - dt*(cg%u(iarr_all_dn,:,j,k) - cg%w(wna%ind(inid_n))%arr(:,:,j,k))*funcR(:,:)
!            enddo
!         enddo
         do j = 1, cg%n_(ydim)
            do k = 1, cg%n_(zdim)
               cg%u(:,:,j,k) = (1.-funcR(:,:))*cg%u(iarr_all_dn,:,j,k) + cg%w(wna%ind(inid_n))%arr(:,:,j,k)*funcR(:,:)
            enddo
         enddo

         where ( cg%u(iarr_all_dn,:,:,:) < 1.05*smalld )
            cg%u(iarr_all_mz,:,:,:) = cg%u(iarr_all_mz,:,:,:)*0.1
         endwhere

         cgl => cgl%nxt
      enddo

      call all_fluid_boundaries

   end subroutine problem_customize_solution_kepler
!-----------------------------------------------------------------------------
   subroutine my_grav_pot_3d

      use gravity,     only: sum_potential
      use cg_list_bnd, only: leaves
      use cg_list,     only: cg_list_element
      use grid_cont,   only: grid_container
      use types,       only: axes

      implicit none

      logical, save :: frun = .true.
      type(cg_list_element), pointer :: cgl
      type(grid_container), pointer  :: cg
      type(axes)                     :: ax

      if (frun) then
         cgl => leaves%first
         do while (associated(cgl))
            cg => cgl%cg
            if (.not.allocated(ax%x)) allocate(ax%x(size(cg%x)))
            if (.not.allocated(ax%y)) allocate(ax%y(size(cg%y)))
            if (.not.allocated(ax%z)) allocate(ax%z(size(cg%z)))
            ax%x = cg%x
            ax%y = cg%y
            ax%z = cg%z

            call my_grav_ptmass_pure(cg%gp,ax,.false.)

            if (allocated(ax%x)) deallocate(ax%x)
            if (allocated(ax%y)) deallocate(ax%y)
            if (allocated(ax%z)) deallocate(ax%z)

            cgl => cgl%nxt
         enddo
      endif

      frun = .false.
      call sum_potential

   end subroutine my_grav_pot_3d
!-----------------------------------------------------------------------------
   subroutine my_grav_ptmass_pure(gp, ax, flatten)

      use units,   only: newtong
      use gravity, only: ptmass, ptm_x, ptm_z
      use types,   only: axes

      implicit none

      real, dimension(:,:,:), pointer       :: gp
      type(axes), intent(in)                :: ax
      logical, intent(in), optional         :: flatten

      integer             :: i, j
      real                :: GM, R2

      GM        = newtong*ptmass


      do i = lbound(gp,1), ubound(gp,1)
         R2 = (ax%x(i) - ptm_x)**2
         do j = lbound(gp,2), ubound(gp,2)
            gp(i,j,:) = -GM / sqrt( (ax%z(:) - ptm_z)**2 + R2 )
         enddo
      enddo
      if (.false.) write(6,*) flatten

   end subroutine my_grav_ptmass_pure
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

      use domain,     only: dom
      use grid_cont,  only: grid_container
      use gravity,    only: grav_pot2accel
      use fluidindex, only: iarr_all_dn, iarr_all_mx, iarr_all_my, iarr_all_mz, flind
#ifndef ISO
      use fluidindex, only: iarr_all_en
#endif /* ISO */
      use constants,  only: xdim, ydim, zdim

      implicit none

      type(grid_container), pointer, intent(in) :: cg

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

      use constants,   only: xdim
      use grid_cont,   only: grid_container
      use named_array_list, only: wna

      implicit none

      type(grid_container), pointer, intent(in) :: cg

      cg%u(:, cg%ie+1:cg%n_(xdim),:,:) = cg%w(wna%ind(inid_n))%arr(:,cg%ie+1:cg%n_(xdim),:,:)

   end subroutine my_bnd_xr
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
end module initproblem
