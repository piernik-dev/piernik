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

   use constants,    only: cbuff_len

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
   real, target, allocatable, dimension(:,:,:,:) :: den0, mtx0, mty0, mtz0, ene0
   real, target, allocatable, dimension(:,:,:)   :: harr
   integer, parameter       :: dname_len = 10

   namelist /PROBLEM_CONTROL/  rho0, cs0, R0, r_in, r_out, f_in, f_out, eps, dumping_coeff, &
      &  use_inner_orbital_period, amp_noise, growth_time, initial_grain_size, final_grain_size

contains
!-----------------------------------------------------------------------------
   subroutine problem_pointers
      use dataio_user,         only: problem_write_restart, problem_read_restart, user_vars_hdf5
      use gravity,             only: grav_pot_3d, grav_type
!     use fluidboundaries_pub, only: user_bnd_xl, user_bnd_xr
      use types,               only: problem_customize_solution, problem_grace_passed
      implicit none
      problem_write_restart => write_initial_fld_to_restart
      problem_read_restart  => read_initial_fld_from_restart
      problem_customize_solution => problem_customize_solution_kepler
!     user_bnd_xl => my_bnd_xl
!     user_bnd_xr => my_bnd_xr
      grav_pot_3d => my_grav_pot_3d
      user_vars_hdf5 => prob_vars_hdf5
      problem_grace_passed => add_random_noise
      grav_type => my_grav_ptmass_pure
!     problem_grace_passed => add_sine

   end subroutine problem_pointers
!-----------------------------------------------------------------------------
   subroutine read_problem_par

      use dataio_pub,          only: ierrh, par_file, namelist_errh, compare_namelist, cmdl_nml      ! QA_WARN required for diff_nml
      use mpi,                 only: MPI_DOUBLE_PRECISION, MPI_LOGICAL
      use mpisetup,            only: rbuff, lbuff, buffer_dim, master, slave, comm, ierr

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

      call MPI_Bcast(rbuff,           buffer_dim, MPI_DOUBLE_PRECISION, 0, comm, ierr)
      call MPI_Bcast(lbuff,           buffer_dim, MPI_LOGICAL,          0, comm, ierr)

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
   subroutine add_sine

      use constants,  only: dpi, xdim, zdim
      use dataio_pub, only: printinfo, warn
      use domain,     only: dom
      use fluidindex, only: flind
      use grid,       only: cga
      use grid_cont,  only: cg_list_element, grid_container
      use mpisetup,   only: master

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

      call cga%get_root(cgl)
      do while (associated(cgl))
         cg => cgl%cg

         do i = cg%is, cg%ie
            do k = cg%ks, cg%ke
               cg%u%arr(flind%dst%imx,i,:,k) = cg%u%arr(flind%dst%imx,i,:,k) + amp*sin(kx*cg%x(i) + kz*cg%z(k)) * cg%u%arr(flind%dst%idn,i,:,k)
               cg%u%arr(flind%dst%imy,i,:,k) = cg%u%arr(flind%dst%imy,i,:,k) + amp*sin(kx*cg%x(i) + kz*cg%z(k)) * cg%u%arr(flind%dst%idn,i,:,k)
               cg%u%arr(flind%dst%imz,i,:,k) = cg%u%arr(flind%dst%imz,i,:,k) + amp*sin(kx*cg%x(i) + kz*cg%z(k)) * cg%u%arr(flind%dst%idn,i,:,k)
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

      use dataio_pub, only: printinfo
      use grid,       only: cga
      use grid_cont,  only: cg_list_element, grid_container
      use fluidindex, only: flind
      use mpisetup,   only: proc, master

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

      call cga%get_root(cgl)
      do while (associated(cgl))
         cg => cgl%cg

         allocate(noise(3,cg%nx,cg%ny,cg%nz))
         call random_number(noise)
         cg%u%arr(flind%dst%imx,:,:,:) = cg%u%arr(flind%dst%imx,:,:,:) +amp_noise -2.0*amp_noise*noise(1,:,:,:) * cg%u%arr(flind%dst%idn,:,:,:)
         cg%u%arr(flind%dst%imy,:,:,:) = cg%u%arr(flind%dst%imy,:,:,:) +amp_noise -2.0*amp_noise*noise(2,:,:,:) * cg%u%arr(flind%dst%idn,:,:,:)
         cg%u%arr(flind%dst%imz,:,:,:) = cg%u%arr(flind%dst%imz,:,:,:) +amp_noise -2.0*amp_noise*noise(3,:,:,:) * cg%u%arr(flind%dst%idn,:,:,:)
         deallocate(noise)

         cgl => cgl%nxt
      enddo

   end subroutine add_random_noise
!-----------------------------------------------------------------------------
   subroutine init_prob

      use constants,    only: DST, GEO_RPZ, xdim
      use global,       only: smalld
      use dataio_pub,   only: msg, printinfo, die
      use domain,       only: geometry_type, cdd
      use fluidindex,   only: flind
      use fluidtypes,   only: component_fluid
      use gravity,      only: ptmass, grav_pot2accel
      use grid,         only: cga
      use grid_cont,    only: cg_list_element, grid_container
      use mpi,          only: MPI_COMM_NULL
      use mpisetup,     only: master, comm
      use units,        only: newtong
      use hydrostatic,  only: hydrostatic_zeq_densmid
      use func,         only: ekin

      implicit none

      integer :: i, j, k, p
      real    :: zk, R, vz, sqr_gm, vr
      real    :: rho, cs2
      real, dimension(:), allocatable :: vphi, grav, ln_dens_der

      type(component_fluid), pointer  :: fl
      type(cg_list_element), pointer :: cgl
      type(grid_container), pointer :: cg

!   Secondary parameters
      if (cdd%comm3d == MPI_COMM_NULL) call die("[initproblem:init_prob] comm3d == MPI_COMM_NULL not implemented") !pcoords
      if (geometry_type /= GEO_RPZ) call die("[initproblem:init_prob] only cylindrical geometry supported")

      sqr_gm = sqrt(newtong*ptmass)

      call cga%get_root(cgl)

      do while (associated(cgl))
         cg => cgl%cg

         if (ubound(cga%cg_all(:), dim=1) > 1) call die("[initproblem:init_prob] multiple grid pieces per procesor not implemented yet") !nontrivial kmid, allocate

         if (.not.allocated(den0)) allocate(den0(flind%fluids, cg%nx, cg%ny, cg%nz))
         if (.not.allocated(mtx0)) allocate(mtx0(flind%fluids, cg%nx, cg%ny, cg%nz))
         if (.not.allocated(mty0)) allocate(mty0(flind%fluids, cg%nx, cg%ny, cg%nz))
         if (.not.allocated(mtz0)) allocate(mtz0(flind%fluids, cg%nx, cg%ny, cg%nz))
         if (.not.allocated(ene0)) allocate(ene0(flind%fluids, cg%nx, cg%ny, cg%nz))

         do p = 1, flind%fluids
            fl => flind%all_fluids(p)
            if (fl%tag /= DST .and. master) then
               write(msg,'(A,F9.5)') "[init_problem:initprob] cs2 used = ", fl%cs2
               call printinfo(msg)
            endif


            do i = 1, cg%nx
               R = cg%x(i)

               rho = max(rho0*(R0/R)**(1.5),smalld)
               cs2 = (cs0*sqrt(R0/R))**2

               cg%cs_iso2%arr(i,:,:) = cs2

               do j = 1, cg%ny
                  call hydrostatic_zeq_densmid(i, j, rho, cs2, cg=cg)
                  do k = 1, cg%nz
                     zk = cg%z(k)
                     cg%u%arr(fl%idn,i,j,k) = max(cg%dprof(k)/(1.0+eps), smalld)
                     if (fl%tag == DST) cg%u%arr(fl%idn,i,j,k) = eps * cg%u%arr(flind%neu%idn,i,j,k)

                  enddo
               enddo
            enddo

            allocate(ln_dens_der(cg%nx), grav(cg%nx), vphi(cg%nx))

            do j = 1, cg%ny
               do k = 1, cg%nz

                  call grav_pot2accel(xdim, j, k, cg%nx, grav, 1, cg)

                  ln_dens_der  = log(cg%u%arr(fl%idn,:,j,k))
                  ln_dens_der(2:cg%nx)  = ( ln_dens_der(2:cg%nx) - ln_dens_der(1:cg%nx-1) ) / cg%dx
                  ln_dens_der(1)        = ln_dens_der(2)

                  vr = 0.0
                  !vphi = sqrt( newtong*ptmass / R )
                  vz = 0.0

                  if (fl%tag /= DST) then
                     vphi = sqrt( max(cg%x(:)*(cg%cs_iso2%arr(:,j,k)*ln_dens_der(:) + abs(grav(:))),0.0) )
                  else
                     vphi = sqrt( max(abs(grav(:)) * cg%x(:), 0.0))
                  endif


                  cg%u%arr(fl%imx,:,j,k) = vr     *cg%u%arr(fl%idn,:,j,k)
                  cg%u%arr(fl%imy,:,j,k) = vphi(:)*cg%u%arr(fl%idn,:,j,k)
                  cg%u%arr(fl%imz,:,j,k) = vz     *cg%u%arr(fl%idn,:,j,k)

               enddo
            enddo
            deallocate(ln_dens_der, grav, vphi)
            if (fl%has_energy) then
               cg%u%arr(fl%ien,:,:,:) = fl%cs2/(fl%gam_1)*cg%u%arr(fl%idn,:,:,:)
               cg%u%arr(fl%ien,:,:,:) = cg%u%arr(fl%ien,:,:,:) + ekin( cg%u%arr(fl%imx,:,:,:), cg%u%arr(fl%imy,:,:,:), &
                  & cg%u%arr(fl%imz,:,:,:), cg%u%arr(fl%idn,:,:,:))
               ene0(p,:,:,:)   = cg%u%arr(fl%ien,:,:,:)
            else
               ene0(p,:,:,:)   = 0.0
            endif

            den0(p,:,:,:) = cg%u%arr(fl%idn,:,:,:)
            mtx0(p,:,:,:) = cg%u%arr(fl%imx,:,:,:)
            mty0(p,:,:,:) = cg%u%arr(fl%imy,:,:,:)
            mtz0(p,:,:,:) = cg%u%arr(fl%imz,:,:,:)
         enddo
         cg%b%arr(:,:,:,:) = 0.0
         cgl => cgl%nxt
      enddo
#ifdef DEBUG
      open(12,file="vel_profile.dat",status="unknown")
         do i = 1, cg%nx
            write(12,'(3(E12.5,1X))') cg%x(i), cg%u%arr(flind%all_fluids(1:2)%imy,i,max(cg%ny/2,1),max(cg%nz/2,1)) / &
                &  cg%u%arr(flind%all_fluids(1:2)%idn,i,max(cg%ny/2,1),max(cg%nz/2,1))
         enddo
      close(12)
#endif /* DEBUG */

   end subroutine init_prob
!-----------------------------------------------------------------------------
   subroutine write_initial_fld_to_restart(file_id, cg)

      use constants,   only: AT_NO_B
      use grid_cont,   only: grid_container
      use hdf5,        only: HID_T
      use dataio_hdf5, only: write_arr_to_restart
      use fluidindex,  only: flind

      implicit none

      integer(HID_T),intent(in)  :: file_id
      type(grid_container), pointer, intent(in) :: cg

      integer :: i
      character(len=dname_len) :: dname
      real, dimension(:,:,:), pointer :: p3d

      if (associated(p3d)) nullify(p3d)
      do i = lbound(den0,1), ubound(den0,1)
         if (allocated(den0)) then
            write(dname,'(2a)') flind%all_fluids(i)%get_tag(), '_den0'
            p3d => den0(i, :, :, :)
            call write_arr_to_restart(file_id, p3d, AT_NO_B, dname, cg)
            nullify(p3d)
         endif
         if (allocated(mtx0)) then
            write(dname,'(2a)') flind%all_fluids(i)%get_tag(), '_mtx0'
            p3d => mtx0(i, :, :, :)
            call write_arr_to_restart(file_id, p3d, AT_NO_B, dname, cg)
            nullify(p3d)
         endif
         if (allocated(mty0)) then
            write(dname,'(2a)') flind%all_fluids(i)%get_tag(), '_mty0'
            p3d => mty0(i, :, :, :)
            call write_arr_to_restart(file_id, p3d, AT_NO_B, dname, cg)
            nullify(p3d)
         endif
         if (allocated(mtz0)) then
            write(dname,'(2a)') flind%all_fluids(i)%get_tag(), '_mtz0'
            p3d => mtz0(i, :, :, :)
            call write_arr_to_restart(file_id, p3d, AT_NO_B, dname, cg)
            nullify(p3d)
         endif
         if (allocated(ene0)) then
            write(dname,'(2a)') flind%all_fluids(i)%get_tag(), '_ene0'
            p3d => ene0(i, :, :, :)
            call write_arr_to_restart(file_id, p3d, AT_NO_B, dname, cg)
            nullify(p3d)
         endif
      enddo

   end subroutine write_initial_fld_to_restart
!-----------------------------------------------------------------------------
   subroutine read_initial_fld_from_restart(file_id, cg)

      use constants,   only: AT_NO_B
      use hdf5,        only: HID_T
      use grid_cont,   only: grid_container
      use fluidindex,  only: flind
      use dataio_hdf5, only: read_arr_from_restart

      implicit none

      integer(HID_T),intent(in) :: file_id
      type(grid_container), pointer, intent(in) :: cg

      character(len=dname_len) :: dname
      real, dimension(:,:,:), pointer :: p3d
      integer :: i

      ! /todo First query for existence of den0, vlx0 and vly0, then allocate
      if (.not.allocated(den0)) allocate(den0(flind%fluids, cg%nx, cg%ny, cg%nz))
      if (.not.allocated(mtx0)) allocate(mtx0(flind%fluids, cg%nx, cg%ny, cg%nz))
      if (.not.allocated(mty0)) allocate(mty0(flind%fluids, cg%nx, cg%ny, cg%nz))
      if (.not.allocated(mtz0)) allocate(mtz0(flind%fluids, cg%nx, cg%ny, cg%nz))
      if (.not.allocated(ene0)) allocate(ene0(flind%fluids, cg%nx, cg%ny, cg%nz))
      if (.not.allocated(harr)) allocate(harr(cg%nx, cg%ny, cg%nz))

      if (.not.associated(p3d)) p3d => harr(:,:,:)
      do i=1, flind%fluids
         write(dname,'(2a)') flind%all_fluids(i)%get_tag(), '_den0'
         call read_arr_from_restart(file_id, p3d, AT_NO_B, dname, cg)
         den0(i,:,:,:) = harr(:,:,:)

         write(dname,'(2a)') flind%all_fluids(i)%get_tag(), '_mtx0'
         call read_arr_from_restart(file_id, p3d, AT_NO_B, dname, cg)
         mtx0(i,:,:,:) = harr(:,:,:)

         write(dname,'(2a)') flind%all_fluids(i)%get_tag(), '_mty0'
         call read_arr_from_restart(file_id, p3d, AT_NO_B, dname, cg)
         mty0(i,:,:,:) = harr(:,:,:)

         write(dname,'(2a)') flind%all_fluids(i)%get_tag(), '_mtz0'
         call read_arr_from_restart(file_id, p3d, AT_NO_B, dname, cg)
         mtz0(i,:,:,:) = harr(:,:,:)

         write(dname,'(2a)') flind%all_fluids(i)%get_tag(), '_ene0'
         call read_arr_from_restart(file_id, p3d, AT_NO_B, dname, cg)
         ene0(i,:,:,:) = harr(:,:,:)
      enddo
      if (associated(p3d)) nullify(p3d)
      if (allocated(harr)) deallocate(harr)

   end subroutine read_initial_fld_from_restart
!-----------------------------------------------------------------------------
   subroutine problem_customize_solution_kepler

      use dataio_pub,      only: die
      use grid,            only: cga
      use grid_cont,       only: cg_list_element, grid_container
      use fluidboundaries, only: all_fluid_boundaries
      use fluidindex,      only: iarr_all_dn, iarr_all_mx, iarr_all_my, iarr_all_mz
      use global,          only: dt, t, grace_period_passed, relax_time
      use constants,       only: dpi
      use units,           only: newtong
      use gravity,         only: ptmass
#ifndef ISO
      use fluidindex,      only: iarr_all_en
#endif /* ISO */
      use interactions,    only: update_grain_size
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

      if (ubound(cga%cg_all(:), dim=1) > 1) call die("[initproblem:problem_customize_solution_kepler] multiple grid pieces per procesor not implemented yet") !nontrivial

      call cga%get_root(cgl)
      do while (associated(cgl))
         cg => cgl%cg

         if (frun) then
            x0 = relax_time
            x1 = x0 + growth_time
            y0 = initial_grain_size
            y1 = final_grain_size
            a = (y0 - y1)/(x0 - x1)
            b = y0 - a*x0

            allocate(funcR(size(iarr_all_dn), cg%nx) )

            funcR(1,:) = 0.0
            funcR(1,:) = funcR(1,:) - tanh((cg%x(:)-r_in+1.0)**f_in) + 1.0
            where (cg%x >= r_out)
               funcR(1,:) = funcR(1,:) + max( tanh((cg%x(:)-r_out+1.0)**f_out), 0.0)
            endwhere

            if (use_inner_orbital_period) then
               funcR(1,:) = funcR(1,:) * sqrt( newtong*ptmass/cg%x(:)**3 ) / dpi
            else
               funcR(1,:) = funcR(1,:) * dumping_coeff
            endif
#ifdef DEBUG
            open(212,file="funcR.dat",status="unknown")
            do j = 1, cg%nx
               write(212,*) cg%x(j),funcR(1,j)
            enddo
            close(212)
#endif /* DEBUG */
            frun = .false.
            funcR(:,:) = spread(funcR(1,:),1,size(iarr_all_dn))
         endif
         funcR = 0.0

         if (grace_period_passed()) call update_grain_size(a*t+b)
         do j = 1, cg%ny
            do k = 1, cg%nz
               cg%u%arr(iarr_all_dn,:,j,k) = cg%u%arr(iarr_all_dn,:,j,k) - dt*(cg%u%arr(iarr_all_dn,:,j,k) - den0(:,:,j,k))*funcR(:,:)
               cg%u%arr(iarr_all_mx,:,j,k) = cg%u%arr(iarr_all_mx,:,j,k) - dt*(cg%u%arr(iarr_all_mx,:,j,k) - mtx0(:,:,j,k))*funcR(:,:)
               cg%u%arr(iarr_all_my,:,j,k) = cg%u%arr(iarr_all_my,:,j,k) - dt*(cg%u%arr(iarr_all_my,:,j,k) - mty0(:,:,j,k))*funcR(:,:)
               cg%u%arr(iarr_all_mz,:,j,k) = cg%u%arr(iarr_all_mz,:,j,k) - dt*(cg%u%arr(iarr_all_mz,:,j,k) - mtz0(:,:,j,k))*funcR(:,:)
#ifndef ISO
               cg%u%arr(iarr_all_en,:,j,k) = cg%u%arr(iarr_all_en,:,j,k) - dt*(cg%u%arr(iarr_all_en,:,j,k) - ene0(:,:,j,k)*funcR(:,:)
#endif /* !ISO */
            enddo
         enddo
!         where ( cg%u%arr(iarr_all_dn,:,:,:) < 2.*smalld )
!            cg%u%arr(iarr_all_mx,:,:,:) = cg%u%arr(iarr_all_mx,:,:,:)*0.1
!            cg%u%arr(iarr_all_mz,:,:,:) = cg%u%arr(iarr_all_mz,:,:,:)*0.1
!         endwhere

         cgl => cgl%nxt
      enddo

      call all_fluid_boundaries

   end subroutine problem_customize_solution_kepler
!-----------------------------------------------------------------------------
   subroutine my_grav_pot_3d

      use gravity,   only: sum_potential, grav_pot_3d_bnd
      use grid,      only: cga
      use grid_cont, only: cg_list_element, grid_container
      use types,     only: axes

      implicit none

      logical, save :: frun = .true.
      type(cg_list_element), pointer :: cgl
      type(grid_container), pointer  :: cg
      type(axes)                     :: ax

      if (frun) then
         call cga%get_root(cgl)
         do while (associated(cgl))
            cg => cgl%cg
            if (.not.allocated(ax%x)) allocate(ax%x(size(cg%x)))
            if (.not.allocated(ax%y)) allocate(ax%y(size(cg%y)))
            if (.not.allocated(ax%z)) allocate(ax%z(size(cg%z)))
            ax%x = cg%x
            ax%y = cg%y
            ax%z = cg%z

            call my_grav_ptmass_pure(cg%gp%arr,ax,.false.)

            if (allocated(ax%x)) deallocate(ax%x)
            if (allocated(ax%y)) deallocate(ax%y)
            if (allocated(ax%z)) deallocate(ax%z)

            cgl => cgl%nxt
         enddo
      endif

      frun = .false.
      call grav_pot_3d_bnd
      call sum_potential

   end subroutine my_grav_pot_3d

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


      do i = 1, ubound(gp,1)
         R2 = (ax%x(i) - ptm_x)**2
         do j = 1, ubound(gp,2)
            gp(i,j,:) = -GM / sqrt( (ax%z(:) - ptm_z)**2 + R2 )
         enddo
      enddo
      if (.false.) write(6,*) flatten
   end subroutine my_grav_ptmass_pure
!-----------------------------------------------------------------------------
   subroutine my_bnd_xl(cg)

      use grid_cont,  only: grid_container
      use gravity,    only: grav_pot2accel
      use fluidindex, only: iarr_all_dn, iarr_all_mx, iarr_all_my, iarr_all_mz, flind
#ifndef ISO
      use fluidindex, only: iarr_all_en
#endif /* ISO */
      use constants,  only: xdim

      implicit none

      type(grid_container), pointer, intent(in) :: cg

      integer :: i
      real, dimension(cg%nx) :: grav
      real, dimension(size(iarr_all_my), cg%ny, cg%nz) :: vy,vym
      real, dimension(size(flind%all_fluids))    :: cs2_arr
      integer, dimension(size(flind%all_fluids)) :: ind_cs2

      do i = 1, size(flind%all_fluids)
         ind_cs2    = i
         cs2_arr(i) = flind%all_fluids(i)%cs2
      enddo

      call grav_pot2accel(xdim,1,1, cg%nx, grav, 1, cg)

      do i = 1, cg%nb
         cg%u%arr(iarr_all_dn,i,:,:) = cg%u%arr(iarr_all_dn, cg%is,:,:)
         cg%u%arr(iarr_all_mx,i,:,:) = min(0.0,cg%u%arr(iarr_all_mx, cg%is,:,:))
!         do p = 1, size(flind%all_fluids)
!            cg%u%arr(iarr_all_my(p),i,:,:) = sqrt( abs(grav(i)) * cg%x(i) - cs2_arr(p)*dens_exp) *  cg%u%arr(iarr_all_dn(p),i,:,:)
!         enddo
         cg%u%arr(iarr_all_my,i,:,:) = cg%u%arr(iarr_all_my, cg%is,:,:)
         cg%u%arr(iarr_all_mz,i,:,:) = cg%u%arr(iarr_all_mz, cg%is,:,:)
#ifndef ISO
         cg%u%arr(iarr_all_en,i,:,:) = cg%u%arr(iarr_all_en, cg%is,:,:)
#endif /* !ISO */
      enddo

      do i = cg%nb,1,-1
         vym(:,:,:) = cg%u%arr(iarr_all_my,i+2,:,:)/cg%u%arr(iarr_all_dn,i+1,:,:)
         vy(:,:,:)  = cg%u%arr(iarr_all_my,i+1,:,:)/cg%u%arr(iarr_all_dn,i+1,:,:)
!         cg%u%arr(iarr_all_my,i,:,:) = (vym(:,:,:) + (cg%x(i) - cg%x(i+2)) / (cg%x(i+1) - cg%x(i+2)) * (vy - vym))*cg%u%arr(iarr_all_dn,i,:,:)
      enddo

   end subroutine my_bnd_xl
!-----------------------------------------------------------------------------
   subroutine my_bnd_xr(cg)

      use grid_cont,  only: grid_container
      use fluidindex, only: iarr_all_dn, iarr_all_mx, iarr_all_my, iarr_all_mz
#ifndef ISO
      use fluidindex, only: iarr_all_en
#endif /* ISO */

      implicit none

      type(grid_container), pointer, intent(in) :: cg

      cg%u%arr(iarr_all_dn, cg%ie+1:cg%nx,:,:) = den0(:, cg%ie+1:cg%nx,:,:)
      cg%u%arr(iarr_all_mx, cg%ie+1:cg%nx,:,:) = mtx0(:, cg%ie+1:cg%nx,:,:)
      cg%u%arr(iarr_all_my, cg%ie+1:cg%nx,:,:) = mty0(:, cg%ie+1:cg%nx,:,:)
      cg%u%arr(iarr_all_mz, cg%ie+1:cg%nx,:,:) = mtz0(:, cg%ie+1:cg%nx,:,:)
#ifndef ISO
      cg%u%arr(iarr_all_en, cg%ie+1:cg%nx,:,:) = ene0(:, cg%ie+1:cg%nx,:,:)
#endif /* !ISO */
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
            tab(:,:,:) = real(epstein_factor(flind%neu%pos) / cg%u%arr(flind%neu%idn,cg%is:cg%ie,cg%js:cg%je,cg%ks:cg%ke), 4)
         case default
            ierrh = -1
      end select

   end subroutine prob_vars_hdf5
!-----------------------------------------------------------------------------
end module initproblem
