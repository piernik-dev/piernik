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
!! \brief Definition of the Wengen 4 test problem.
!!
!! \details Look here: http://users.camk.edu.pl/gawrysz/test4/ for detailed description of initial conditions and a catalogue of results.
!<

module initproblem

   use constants,    only: dsetnamelen, cbuff_len, ndims

   implicit none

   private
   public :: read_problem_par, problem_initial_conditions, problem_pointers

   !namelist parameters
   character(len=cbuff_len) :: input_file                     !< File with initial conditions
   real                     :: gamma_loc                      !< gamma used for calculating initial T distribution
   real                     :: mass_mul                       !< density scaling factor, default 1.0, try also 1.02
   real                     :: ambient_density                !< modify velocities below this density
   real                     :: cs_mul                         !< temperature scaling factor, implemented for debugging only
   real                     :: damp_factor                    !< Set 1. to clear velocities in the ambient medium, 0. does nothing
   integer(kind=4)          :: divine_intervention_type       !< select type of every-step solution alteration
   real                     :: mincs2, maxcs2                 !< extreme soundspeed values found in the IC file
   real                     :: r_in                           !< inner radius of d_i_t = 3, for r < r_in we enforce den0, vlx0, vly0
   real                     :: r_out                          !< outer radius of d_i_t = 3, for r > r_out we enforce den0, vlx0, vly0
   ! BEWARE: small value of f_{in,out{ (<10.) may smear d_i_t over many cells around r_{in,out}, large value (~100) is equal to imposing step function at r_{in,out}
   real                     :: f_in                           !< smoothing factor for cutoff at r_in
   real                     :: f_out                          !< smoothing factor for cutoff at r_out
   real                     :: alfasupp                       !< scaling factor for d_i_t = 3, in most cases should be = 1
   logical                  :: fake_ic                        !< Skip reading the IC file (useful only for debugging, or running under valgrind)
   real                     :: T_disk                         !< temperature of the disk
   real                     :: mean_mol_weight                !< mean molecular weight

   namelist /PROBLEM_CONTROL/  input_file, gamma_loc, mass_mul, ambient_density, cs_mul, damp_factor, divine_intervention_type, mincs2, maxcs2, &
      &                        r_in, r_out, f_in, f_out, alfasupp, fake_ic, T_disk, mean_mol_weight

   ! private data
   integer, parameter :: ic_nx = 512, ic_ny = 512, ic_nz = 52 !< initial conditions size
   real, parameter    :: ic_xysize = 8.                       !< X- and Y- size of the domain covered by the IC
   real, parameter    :: ic_zsize = (ic_xysize*ic_nz)/ic_nx   !< Z-size of the domain covered by the IC
   real, parameter    :: ic_dx = ic_xysize/ic_nx              !< dx=dy=dz in the IC
   real, dimension(ndims) :: starpos, starvel                 !< the primary star initial position and velocity
   real, allocatable, dimension(:, :, :, :) :: ic_data        !< Storage for local part of the IC file
   enum, bind(C)
      enumerator :: D0, VX0, VY0
   end enum
   character(len=dsetnamelen), dimension(D0:VY0), parameter :: q_n = [ "den0", "vlx0", "vly0" ] !< Names of initial condition (t=0.) arrays used for divine_intervention_type = 3

contains

!> \brief Pointers to user-provided routines

   subroutine problem_pointers

      use dataio_user, only: user_attrs_wr
      use user_hooks,  only: problem_customize_solution, cleanup_problem

      implicit none

      problem_customize_solution => problem_customize_solution_wt4
      user_attrs_wr              => problem_initial_conditions_attrs
      cleanup_problem            => cleanup_wt4

   end subroutine problem_pointers

!> \brief Read problem parameters

   subroutine read_problem_par

      use cg_list_global, only: all_cg
      use constants,      only: AT_NO_B
      use dataio_pub,     only: nh      ! QA_WARN required for diff_nml
      use mpisetup,       only: rbuff, cbuff, ibuff, lbuff, master, slave, piernik_MPI_Bcast

      implicit none

      integer :: i
!      integer, parameter :: maxsub = 10  !< upper limit for subsampling

      ! namelist default parameter values
      input_file      = './test4-512.alt'

      gamma_loc       = 1.4
      mass_mul        = 1.0
      cs_mul          = 1.0
      ambient_density = 1e-5
      damp_factor     = 0.9
      mincs2          = 8.725322e-4
      maxcs2          = 5.8168972e-3
      r_in            = 0.5
      r_out           = 3.3
      f_in            = 10.0
      f_out           = 50.0
      alfasupp        = 1.0
      T_disk          = 20.0
      mean_mol_weight = 2.0

      divine_intervention_type = 2

      fake_ic = .false.

      starpos(:)      = 0.0 ! for test4-512 [ -0.00190265, 0.0379506,   0.00083884  ]
      starvel(:)      = 0.0 ! for test4-512 [  0.00172268, 0.00178423, -6.20918e-05 ]

      if (master) then

         diff_nml(PROBLEM_CONTROL)

         cbuff(1) =  input_file

         rbuff(1) = gamma_loc
         rbuff(2) = mass_mul
         rbuff(3) = ambient_density
         rbuff(4) = cs_mul
         rbuff(5) = damp_factor
         rbuff(6) = mincs2
         rbuff(7) = maxcs2
         rbuff(8) = r_in
         rbuff(9) = r_out
         rbuff(10) = f_in
         rbuff(11) = f_out
         rbuff(12) = alfasupp
         rbuff(13) = T_disk
         rbuff(14) = mean_mol_weight

         ibuff(1) = divine_intervention_type

         lbuff(1) = fake_ic

      endif

      call piernik_MPI_Bcast(cbuff, cbuff_len)
      call piernik_MPI_Bcast(ibuff)
      call piernik_MPI_Bcast(rbuff)
      call piernik_MPI_Bcast(lbuff)

      if (slave) then

         input_file   = trim(cbuff(1))

         gamma_loc       = rbuff(1)
         mass_mul        = rbuff(2)
         ambient_density = rbuff(3)
         cs_mul          = rbuff(4)
         damp_factor     = rbuff(5)
         mincs2          = rbuff(6)
         maxcs2          = rbuff(7)
         r_in            = rbuff(8)
         r_out           = rbuff(9)
         f_in            = rbuff(10)
         f_out           = rbuff(11)
         alfasupp        = rbuff(12)
         T_disk          = rbuff(13)
         mean_mol_weight = rbuff(14)

         divine_intervention_type = ibuff(1)

         fake_ic = lbuff(1)

      endif

      if (mass_mul < 0.) mass_mul = 1.

      do i = D0, VY0
         call all_cg%reg_var(q_n(i), restart_mode = AT_NO_B)
      enddo

   end subroutine read_problem_par

!>
!! \brief Read the file with initial conditions.
!!
!! \details Read the initial conditions on the master process and send it to all other processes.
!! The data is quite large in size and should be kept in memory because it is not possible to determine how many times the problem_initial_conditions will be called
!! (especially in setups employing AMR).
!<

   subroutine read_IC_file

      use cg_leaves,   only: leaves
      use constants,   only: ndims
      use dataio_pub,  only: msg, die
      use grid_cont,   only: grid_container
      use mpi,         only: MPI_DOUBLE_PRECISION
      use mpisetup,    only: proc, master, FIRST, LAST, comm, status, mpi_err

      implicit none

      integer                             :: i, j, k, v, pe, ostat
      type(grid_container), pointer       :: cg
      enum, bind(C)
         enumerator :: DEN0 = 1, VELX0, VELZ0 = VELX0+ndims-1, ENER0
      end enum

      cg => leaves%first%cg

      if (allocated(ic_data)) call die("[initproblem:read_IC_file] ic_data already allocated")
      allocate(ic_data(ic_nx, ic_ny, ic_nz, DEN0:ENER0))

      if (master) then
         open(1, file=input_file, status='old', iostat=ostat)
         if (ostat /= 0) then
            write(msg,'(3a,i4)')"[initproblem:read_IC_file] cannot read ic_data from file '",input_file,"' at PE#",proc
            call die(msg)
         endif
      endif

      do v = DEN0, ENER0
         if (master) then ! read the quantities, then send to everyone interested
            do k = 1, ic_nz
               do j = 1, ic_ny
                  do i = 1, ic_nx
                     read(1,*) ic_data(i, j, k, v)
                  enddo
               enddo
            enddo
            do pe = FIRST+1, LAST
               call MPI_Send(ic_data, size(ic_data), MPI_DOUBLE_PRECISION, pe, pe, comm, mpi_err)
            enddo
         else
            call MPI_Recv(ic_data, size(ic_data), MPI_DOUBLE_PRECISION, FIRST, proc, comm, status, mpi_err)
         endif
      enddo

      if (master) close(1)

      if (mass_mul /= 1.0) ic_data(:, :, :, DEN0) = ic_data(:, :, :, DEN0) * mass_mul

      do v = VELX0, VELZ0 ! convert velocity to momentum
         ic_data(:, :, :, v) = ic_data(:, :, :, v) * ic_data(:, :, :, DEN0)
      enddo

      ! U = ( kB * T ) / (mean_mol_weight * (gamma - 1))
      ! cs2 = (gamma) * kB * T / mean_mol_weight.
      !   => cs2 = U * (gamma) * (gamma - 1)
      ic_data(:, :, :, ENER0) = ic_data(:, :, :, ENER0) * (gamma_loc - 1.0) * cs_mul ! * gamma_loc

      ! BEWARE: Until we have decent hdf5 restart that would be problem
      ! dependent, i.e. >=gcc-4.5 / >=ifort 10.1, things like
      ! that must be present in problem.par
      !mincs2 = minval(ic_data(:, :, :, ENER0))
      !maxcs2 = maxval(ic_data(:, :, :, ENER0))

   end subroutine read_IC_file

!> \brief deallocate working arrays

   subroutine cleanup_wt4

      implicit none

      if (allocated(ic_data)) deallocate(ic_data)

   end subroutine cleanup_wt4

!> \bried initialize fluids with the initial conditions data

   subroutine problem_initial_conditions

      use cg_list,          only: cg_list_element
      use cg_leaves,        only: leaves
      use constants,        only: small, GEO_XYZ, GEO_RPZ
      use dataio_pub,       only: warn, printinfo, msg, die
      use domain,           only: dom
      use global,           only: smalld
      use grid_cont,        only: grid_container
      use fluidindex,       only: flind
      use fluidtypes,       only: component_fluid
      use mpisetup,         only: master
      use named_array_list, only: qna
      use units,            only: kboltz, mH

      implicit none

      real, parameter                 :: beat_dx = 1e-5
      integer                         :: i, j, k, iic, jic, kic
      type(cg_list_element),  pointer :: cgl
      type(grid_container),   pointer :: cg
      class(component_fluid), pointer :: fl
      real, dimension(:,:,:), pointer :: q0
      real                            :: vr, vphi

      if (.not. allocated(ic_data) .and. .not. fake_ic) call read_IC_file

      fl => flind%neu
      cgl => leaves%first
      do while (associated(cgl))
         cg => cgl%cg

         if (master) then
            if (maxval(cg%dl(:), mask=dom%has_dir(:)) > ic_dx) then
               write(msg,'(a)')     "[initproblem:problem_initial_conditions] Too low resolution"
               call warn(msg)
            endif
            if (abs(ic_dx/cg%dx-anint(ic_dx/cg%dx)) > beat_dx) then
               write(msg,'(a,f8.4)')"[initproblem:problem_initial_conditions] X-direction requires interpolation ic_dx/dx= ", ic_dx/cg%dx
               call warn(msg)
            endif
            if (abs(ic_dx/cg%dy-anint(ic_dx/cg%dy)) > beat_dx) then
               write(msg,'(a,f8.4)')"[initproblem:problem_initial_conditions] Y-direction requires interpolation ic_dx/dy= ", ic_dx/cg%dy
               call warn(msg)
            endif
            if (abs(ic_dx/cg%dz-anint(ic_dx/cg%dz)) > beat_dx) then
               write(msg,'(a,f8.4)')"[initproblem:problem_initial_conditions] Z-direction requires interpolation ic_dx/dz= ", ic_dx/cg%dz
               call warn(msg)
            endif
         endif

         if (fake_ic) then
            cg%u(fl%idn, cg%is:cg%ie, cg%js:cg%je, cg%ks:cg%ke) = 1.
            cg%u(fl%imx, cg%is:cg%ie, cg%js:cg%je, cg%ks:cg%ke) = 0.
            cg%u(fl%imy, cg%is:cg%ie, cg%js:cg%je, cg%ks:cg%ke) = 0.
            cg%u(fl%imz, cg%is:cg%ie, cg%js:cg%je, cg%ks:cg%ke) = 0.
            cg%cs_iso2(cg%is:cg%ie, cg%js:cg%je, cg%ks:cg%ke) = 1e-2
         else
            do k = cg%ks, cg%ke
               kic = nint((cg%z(k) + ic_zsize/2.)/ic_dx)
               select case (dom%geometry_type)
                  case (GEO_XYZ)
                     do j = cg%js, cg%je
                        jic = nint((cg%y(j) + ic_xysize/2.)/ic_dx)
                        do i = cg%is, cg%ie
                           iic = nint((cg%x(i) + ic_xysize/2.)/ic_dx)
                           call set_point(i, j, k, iic, jic, kic)
                        enddo
                     enddo
                  case (GEO_RPZ)
                     do j = cg%js, cg%je
                        do i = cg%is, cg%ie
                           iic = nint((cg%x(i)*cos(cg%y(j)) + ic_xysize/2.)/ic_dx)
                           jic = nint((cg%x(i)*sin(cg%y(j)) + ic_xysize/2.)/ic_dx)
                           call set_point(i, j, k, iic, jic, kic)
                           vr   =  cg%u(fl%imx, i, j, k)*cos(cg%y(j)) + cg%u(fl%imy, i, j, k)*sin(cg%y(j))
                           vphi = -cg%u(fl%imx, i, j, k)*sin(cg%y(j)) + cg%u(fl%imy, i, j, k)*cos(cg%y(j))
                           cg%u(fl%imx, i, j, k) = vr
                           cg%u(fl%imy, i, j, k) = vphi
                        enddo
                     enddo
                  case default
                     call die("[initproblem:problem_initial_conditions] geometry not supported.")
               end select
            enddo
         endif

         do i = 1, dom%nb
            cg%u(:, cg%is-i,:,:)    = cg%u(:, cg%is,:,:)
            cg%u(:, cg%ie+i,:,:)    = cg%u(:, cg%ie,:,:)
            cg%cs_iso2(cg%is-i,:,:) = cg%cs_iso2(cg%is,:,:)
            cg%cs_iso2(cg%ie+i,:,:) = cg%cs_iso2(cg%ie,:,:)

            cg%u(:,:,cg%js-i,:)     = cg%u(:,:, cg%js,:)
            cg%u(:,:,cg%je+i,:)     = cg%u(:,:, cg%je,:)
            cg%cs_iso2(:,cg%js-i,:) = cg%cs_iso2(:, cg%js,:)
            cg%cs_iso2(:,cg%je+i,:) = cg%cs_iso2(:, cg%je,:)

            cg%u(:,:,:, cg%ks-i)    = cg%u(:,:,:, cg%ks)
            cg%u(:,:,:, cg%ke+i)    = cg%u(:,:,:, cg%ke)
            cg%cs_iso2(:,:,cg%ks-i) = cg%cs_iso2(:,:, cg%ks)
            cg%cs_iso2(:,:,cg%ke+i) = cg%cs_iso2(:,:, cg%ke)
         enddo
         if (master ) then
            write(msg,'(2(a,g15.7))') '[initproblem:problem_initial_conditionslem]: minval(dens)    = ', minval(cg%u(fl%idn,:,:,:)),      ' maxval(dens)    = ', maxval(cg%u(fl%idn,:,:,:))
            call printinfo(msg, .true.)
            write(msg,'(2(a,g15.7))') '[initproblem:problem_initial_conditionslem]: minval(cs_iso2) = ', minval(cg%cs_iso2(:,:,:)), ' maxval(cs_iso2) = ', maxval(cg%cs_iso2(:,:,:))
            call printinfo(msg, .true.)
         endif

         cg%b(:, :, :, :) = 0.0

         do i = D0, VY0
            q0 => cg%q(qna%ind(q_n(i)))%arr
            select case (i)
               case (D0)
                  q0 = cg%u(fl%idn,:,:,:)
               case (VX0)
                  q0 = cg%u(fl%imx,:,:,:) / cg%u(fl%idn,:,:,:)
               case (VY0)
                  q0 = cg%u(fl%imy,:,:,:) / cg%u(fl%idn,:,:,:)
               case default
                  call die("[initproblem:problem_initial_conditionslem] Illegal quantity")
            end select
         enddo
         ! It would be cool to dump a restart file here but this would make a cyclic dependency

         cgl => cgl%nxt
      enddo

#ifndef UMUSCL
      if (master ) call warn("[initproblem:problem_initial_conditionslem]: Without UMUSCL you'll likely get Monet-like density maps.")
#endif /* !UMUSCL */

   contains

      subroutine set_point(i, j, k, iic, jic, kic)

         implicit none

         integer, intent(in) :: i, j, k, iic, jic, kic

         if ( iic >= lbound(ic_data, dim=1) .and. iic <= ubound(ic_data, dim=1) .and. &
              jic >= lbound(ic_data, dim=2) .and. jic <= ubound(ic_data, dim=2) .and. &
              kic >= lbound(ic_data, dim=3) .and. kic <= ubound(ic_data, dim=3) ) then
            cg%u(fl%idn, i, j, k) = ic_data(iic, jic, kic, 1) ! simple injection
            cg%u(fl%imx, i, j, k) = ic_data(iic, jic, kic, 2)
            cg%u(fl%imy, i, j, k) = ic_data(iic, jic, kic, 3)
            cg%u(fl%imz, i, j, k) = ic_data(iic, jic, kic, 4)
            ! cs_iso2_arr(i, j, k) = ic_data(iic, jic, kic, 5)
            cg%cs_iso2(  i, j, k) = (gamma_loc) * kboltz * T_disk / mean_mol_weight / mH
         else
            cg%u(fl%idn, i, j, k)     = smalld
            cg%u(fl%imx, i, j, k)     = small
            cg%u(fl%imy, i, j, k)     = small
            cg%u(fl%imz, i, j, k)     = small
            cg%cs_iso2(  i, j, k)     = mincs2
         endif

      end subroutine set_point

   end subroutine problem_initial_conditions

!> \brief Add some attributes to the datafiles

   subroutine problem_initial_conditions_attrs(file_id)

      use hdf5,  only: HID_T, SIZE_T
      use h5lt,  only: h5ltset_attribute_double_f
      use units, only: fpiG

      implicit none

      integer(HID_T),intent(in) :: file_id
      integer(SIZE_T)           :: bufsize
      integer(kind=4)           :: error

      bufsize = 1
      call h5ltset_attribute_double_f(file_id, "/", "fpiG", [fpiG], bufsize, error)

   end subroutine problem_initial_conditions_attrs

!> \brief modify the density and velocity fields to provide kind of boundary conditions enforced far from domain boundaries

   subroutine problem_customize_solution_wt4(forward)

      use cg_list,          only: cg_list_element
      use cg_leaves,        only: leaves
      use constants,        only: xdim, ydim, zdim, LO, HI, GEO_XYZ, GEO_RPZ
      use dataio_pub,       only: warn, die
      use domain,           only: dom
      use fluidindex,       only: flind
      use fluidtypes,       only: component_fluid
      use grid_cont,        only: grid_container
      use named_array_list, only: qna

      implicit none

      logical, intent(in)               :: forward
      integer                           :: i, j, k
      real, allocatable, dimension(:)   :: mod_str
      real, parameter                   :: max_ambient = 100. ! do not modify solution if density is above max_ambient * ambient_density
      real, allocatable, dimension(:,:) :: alf
      real                              :: rc, ambient_density_min
      type(cg_list_element),  pointer   :: cgl
      type(grid_container),   pointer   :: cg
      class(component_fluid), pointer   :: fl
      real, dimension(:,:,:), pointer   :: den0, vlx0, vly0

      fl => flind%neu
      cgl => leaves%first
      do while (associated(cgl))
         cg => cgl%cg

         allocate(mod_str(cg%is:cg%ie))

         select case (divine_intervention_type)
            case (1)                                                                                ! crude
               if (dom%geometry_type /= GEO_XYZ) call die("[initproblem:problem_customize_solution_wt4] Non-cartesian geometry not supported (divine_intervention_type=1).")! remapping required
               where (cg%u(fl%idn, cg%is:cg%ie, cg%js:cg%je, cg%ks:cg%ke) < ambient_density)
                  cg%u(fl%imx, cg%is:cg%ie, cg%js:cg%je, cg%ks:cg%ke) = (1. - damp_factor) * cg%u(fl%imx, cg%is:cg%ie, cg%js:cg%je, cg%ks:cg%ke)
                  cg%u(fl%imy, cg%is:cg%ie, cg%js:cg%je, cg%ks:cg%ke) = (1. - damp_factor) * cg%u(fl%imy, cg%is:cg%ie, cg%js:cg%je, cg%ks:cg%ke)
                  cg%u(fl%imz, cg%is:cg%ie, cg%js:cg%je, cg%ks:cg%ke) = (1. - damp_factor) * cg%u(fl%imz, cg%is:cg%ie, cg%js:cg%je, cg%ks:cg%ke)
                  cg%cs_iso2(cg%is:cg%ie, cg%js:cg%je, cg%ks:cg%ke) = mincs2
               elsewhere
                  cg%cs_iso2(cg%is:cg%ie, cg%js:cg%je, cg%ks:cg%ke) = maxcs2
               endwhere
            case (2)                                                                                ! smooth
               if (dom%geometry_type /= GEO_XYZ) call die("[initproblem:problem_customize_solution_wt4] Non-cartesian geometry not supported (divine_intervention_type=2).")
               ambient_density_min = ambient_density / max_ambient
               do k = cg%ks, cg%ke
                  do j = cg%js, cg%je
                     mod_str(cg%is:cg%ie) = max(0., (1. + 1./max_ambient) * ambient_density_min / (max(0., cg%u(fl%idn, cg%is:cg%ie, j, k)) + ambient_density_min) - 1./max_ambient)
                     ! ifort can have memory leaks on WHERE - let's provide explicit loop for this crappy compiler
#ifdef __INTEL_COMPILER
                     do i = cg%is, cg%ie
                        if (mod_str(i) > max_ambient**(-2)) then
                           cg%u(fl%idn,     i, j, k) = cg%u(fl%idn, i, j, k) + ambient_density_min * mod_str(i)
                           cg%u(fl%imx,     i, j, k) = cg%u(fl%imx, i, j, k) * (1. - damp_factor   * mod_str(i))
                           cg%u(fl%imy,     i, j, k) = cg%u(fl%imy, i, j, k) * (1. - damp_factor   * mod_str(i))
                           cg%u(fl%imz,     i, j, k) = cg%u(fl%imz, i, j, k) * (1. - damp_factor   * mod_str(i))
                           cg%cs_iso2(i, j, k) = maxcs2           -  (maxcs2-mincs2)    * mod_str(i)
                        endif
                     enddo
#else /* !__INTEL_COMPILER */
                     where (mod_str(cg%is:cg%ie) > max_ambient**(-2))
                        cg%u(fl%idn,     cg%is:cg%ie, j, k) = cg%u(fl%idn, cg%is:cg%ie, j, k) + ambient_density_min * mod_str(cg%is:cg%ie)
                        cg%u(fl%imx,     cg%is:cg%ie, j, k) = cg%u(fl%imx, cg%is:cg%ie, j, k) * (1. - damp_factor   * mod_str(cg%is:cg%ie))
                        cg%u(fl%imy,     cg%is:cg%ie, j, k) = cg%u(fl%imy, cg%is:cg%ie, j, k) * (1. - damp_factor   * mod_str(cg%is:cg%ie))
                        cg%u(fl%imz,     cg%is:cg%ie, j, k) = cg%u(fl%imz, cg%is:cg%ie, j, k) * (1. - damp_factor   * mod_str(cg%is:cg%ie))
                        cg%cs_iso2(cg%is:cg%ie, j, k) = maxcs2               -  (maxcs2-mincs2)    * mod_str(cg%is:cg%ie)
                     endwhere
#endif /* !__INTEL_COMPILER */
                  enddo
               enddo
            case (3)
               den0 => cg%q(qna%ind(q_n(D0)))%arr
               vlx0 => cg%q(qna%ind(q_n(VX0)))%arr
               vly0 => cg%q(qna%ind(q_n(VY0)))%arr
               allocate(alf(cg%lhn(xdim, LO):cg%lhn(xdim, HI), cg%lhn(ydim, LO):cg%lhn(ydim, HI)))
               do j = cg%lhn(ydim, LO), cg%lhn(ydim, HI)
                  do i = cg%lhn(xdim, LO), cg%lhn(xdim, HI)
                     select case (dom%geometry_type)
                        case (GEO_XYZ)
                           rc = sqrt(cg%x(i)**2 + cg%y(j)**2)
                        case (GEO_RPZ)
                           rc = cg%x(i)
                        case default
                           rc = 0 ! suppress compiler warning
                           call die("[initproblem:problem_customize_solution_wt4] geometry not supported (divine_intervention_type=3).")
                     end select
                     alf(i,j) = -alfasupp*0.5*(tanh((rc-r_in)/r_in*f_in)-1.)
                     alf(i,j) = alf(i,j) + alfasupp*0.5*(tanh((rc-r_out)/r_out*f_out) + 1.)
                  enddo
               enddo
               do k =  cg%lhn(zdim, LO), cg%lhn(zdim, HI)
                  cg%u(fl%idn, :, :, k) = (1. - alf(:,:))*cg%u(fl%idn, :, :, k) + alf*den0(:, :, k)
                  cg%u(fl%imx, :, :, k) = (1. - alf(:,:))*cg%u(fl%imx, :, :, k) + alf*den0(:, :, k) * vlx0(:, :, k)
                  cg%u(fl%imy, :, :, k) = (1. - alf(:,:))*cg%u(fl%imy, :, :, k) + alf*den0(:, :, k) * vly0(:, :, k)
               enddo
               do k = cg%ks, cg%ke
                  do j = cg%js, cg%je
                     mod_str(cg%is:cg%ie) = max(0., (1. + 1./max_ambient) * ambient_density / (max(0., cg%u(fl%idn, cg%is:cg%ie, j, k)) + ambient_density) - 1./max_ambient)
                     where (mod_str(cg%is:cg%ie) > max_ambient**(-2))
                        cg%u(fl%imz,     cg%is:cg%ie, j, k) = cg%u(fl%imz, cg%is:cg%ie, j, k) * (1. - damp_factor * mod_str(cg%is:cg%ie))
                     endwhere
                  enddo
               enddo
               deallocate(alf)
            case default
               call warn("[initproblem:problem_customize_solution_wt4] Unknown divine_intervention_type")
         end select

         deallocate(mod_str)

         cgl => cgl%nxt
      enddo

      return
      if (.false. .and. forward) i = j ! suppress compiler warnings on unused arguments

   end subroutine problem_customize_solution_wt4

end module initproblem
