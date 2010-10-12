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
!    Initial implemetation of PIERNIK code was based on TVD split MHD code by
!    Ue-Li Pen
!        see: Pen, Arras & Wong (2003) for algorithm and
!             http://www.cita.utoronto.ca/~pen/MHD
!             for original source code "mhd.f90"
!
!    For full list of developers see $PIERNIK_HOME/license/pdt.txt
!
#include "piernik.def"
#include "macros.h"

module initproblem

   use mpisetup,    only: cbuff_len
   use problem_pub, only: problem_name, run_id

   integer, parameter :: ic_nx = 512, ic_ny = 512, ic_nz = 52 !< initial conditions size
   integer, parameter :: ic_vars = 5                          !< number of quantities in the IC
   real, parameter    :: ic_xysize = 8.                       !< X- and Y- size of the domain covered by the IC
   real, parameter    :: ic_zsize = (ic_xysize*ic_nz)/ic_nx   !< Z-size of the domain covered by the IC
   real, parameter    :: ic_dx = ic_xysize/ic_nx              !< dx=dy=dz in the IC
   real               :: mass_mul                             !< density scaling factor, default 1.0, try also 1.02
   real               :: cs_mul                               !< temperature scaling factor, implemented for debugging only
   real               :: gamma_loc                            !< gamma used for calculating initial T distribution
   real, dimension(3) :: starpos, starvel                     !< the primary star initial position and velocity
   real               :: mincs2, maxcs2                       !< extreme soundspeed values found in the IC file
   real               :: ambient_density                      !< modify velocities below this density
   real               :: damp_factor                          !< Set 1. to clear velocities in the ambient medium, 0. does nothing
   integer            :: divine_intervention_type             !< select type of every-step solution alteration
   logical            :: fake_ic                              !< Skip reading the IC file (useful only for debugging, or running under valgrind)
   character(len=cbuff_len) :: input_file                     !< File with initial conditions

   real, allocatable, dimension(:, :, :, :) :: ic_data        !< Storage for local part of the IC file
   integer :: ic_is, ic_ie, ic_js, ic_je, ic_ks, ic_ke        !< range  of IC file covering local domain
   ! BEWARE: following arrays can be allocated in various places, yet deallocated nowhere
   !  they are used only for divine_intervention_type = 3
   real, allocatable, dimension(:, :, :), target  :: den0     !< Density distribution at t = 0
   real, allocatable, dimension(:, :, :), target  :: vlx0     !< Distribution of X component of velocity at t = 0
   real, allocatable, dimension(:, :, :), target  :: vly0     !< Distribution of Y component of velocity at t = 0
   real                 :: r_in                               !< inner radius of d_i_t = 3, for r < r_in we enforce den0, vlx0, vly0
   real                 :: r_out                              !< outer radius of d_i_t = 3, for r > r_out we enforce den0, vlx0, vly0
   ! BEWARE: small value of f_{in,out{ (<10.) may smear d_i_t over many cells around r_{in,out}, large value (~100) is equal to imposing step function at r_{in,out}
   real                 :: f_in                               !< smoothing factor for cutoff at r_in
   real                 :: f_out                              !< smoothing factor for cutoff at r_out
   real                 :: alfasupp                           !< scaling factor for d_i_t = 3, in most cases should be = 1
   real                 :: T_disk                             !< temperature of the disk
   real                 :: mean_mol_weight                    !< mean molecular weight

   namelist /PROBLEM_CONTROL/  problem_name, run_id, input_file, gamma_loc, mass_mul, ambient_density, cs_mul, damp_factor, divine_intervention_type, mincs2, maxcs2, &
      &                        r_in, r_out, f_in, f_out, alfasupp, fake_ic, T_disk, mean_mol_weight

contains

!-----------------------------------------------------------------------------

   subroutine read_problem_par

      use grid,          only: xmin, xmax, ymin, ymax, zmin, zmax
      use mpisetup,      only: ierr, rbuff, cbuff, ibuff, lbuff, proc, &
           &                    MPI_CHARACTER, MPI_DOUBLE_PRECISION, MPI_INTEGER, MPI_LOGICAL, &
           &                    buffer_dim, comm, smalld
      use constants,     only: pi
      use dataio_public, only: ierrh, msg, par_file, die, namelist_errh, compare_namelist

      implicit none

!      integer, parameter :: maxsub = 10  !< upper limit for subsampling

      ! namelist default parameter values
      problem_name    = 'wengen4'
      input_file      = './test4-512.alt'
      run_id          = 'tst'

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

      if (proc == 0) then

         diff_nml(PROBLEM_CONTROL)

         cbuff(1) =  problem_name
         cbuff(2) =  run_id
         cbuff(3) =  input_file

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

      call MPI_Bcast(cbuff, cbuff_len*buffer_dim, MPI_CHARACTER,        0, comm, ierr)
      call MPI_Bcast(ibuff,           buffer_dim, MPI_INTEGER,          0, comm, ierr)
      call MPI_Bcast(rbuff,           buffer_dim, MPI_DOUBLE_PRECISION, 0, comm, ierr)
      call MPI_Bcast(lbuff,           buffer_dim, MPI_LOGICAL,          0, comm, ierr)

      if (proc /= 0) then

         problem_name = trim(cbuff(1))
         run_id       = cbuff(2)(1:3)
         input_file   = trim(cbuff(3))

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

   end subroutine read_problem_par

!-----------------------------------------------------------------------------

   subroutine read_IC_file

      use grid,          only: xminb, xmaxb, yminb, ymaxb, zminb, zmaxb
      use mpisetup,      only: proc, nproc, comm3d, status, ierr, MPI_INTEGER, MPI_DOUBLE_PRECISION
      use dataio_public, only: msg, die

      implicit none

      integer, parameter                  :: margin = 1
      integer                             :: i, j, k, v, pe, ostat
      real, allocatable, dimension(:,:,:) :: ic_v
      integer, parameter                  :: NDIM = 3
      integer, dimension(2*NDIM)          :: ic_rng

      ! calculate index ranges for the subset of IC file covering local domain with a safety margin for interpolation
      ic_is = min(ic_nx, max(1,     1+floor((xminb + ic_xysize/2.)/ic_dx) - margin) )
      ic_ie = max(1,     min(ic_nx, ceiling((xmaxb + ic_xysize/2.)/ic_dx) + margin) )
      ic_js = min(ic_ny, max(1,     1+floor((yminb + ic_xysize/2.)/ic_dx) - margin) )
      ic_je = max(1,     min(ic_ny, ceiling((ymaxb + ic_xysize/2.)/ic_dx) + margin) )
      ic_ks = min(ic_nz, max(1,     1+floor((zminb + ic_zsize/2. )/ic_dx) - margin) )
      ic_ke = max(1,     min(ic_nz, ceiling((zmaxb + ic_zsize/2. )/ic_dx) + margin) )

      if (allocated(ic_data)) call die("[initproblem:read_IC_file] ic_data already allocated")
      allocate(ic_data(ic_is:ic_ie, ic_js:ic_je, ic_ks:ic_ke, ic_vars))

      if (proc == 0) then
         open(1, file=input_file, status='old', iostat=ostat)
         if (ostat /= 0) then
            write(msg,'(3a,i4)')"[initproblem:read_IC_file] cannot read ic_data from file '",input_file,"' at PE#",proc
            call die(msg)
         endif
         if (allocated(ic_v)) deallocate(ic_v)
         allocate(ic_v(ic_nx, ic_ny, ic_nz))
      endif

      do v = 1, ic_vars
         if (proc == 0) then ! read the quantities, then send to everyone interested
            do k = 1, ic_nz
               do j = 1, ic_ny
                  do i = 1, ic_nx
                     read(1,*) ic_v(i,j,k)
                  enddo
               enddo
            enddo
            ic_data(ic_is:ic_ie, ic_js:ic_je, ic_ks:ic_ke, v) = ic_v(ic_is:ic_ie, ic_js:ic_je, ic_ks:ic_ke)
            do pe = 1, nproc-1
               call MPI_Recv( ic_rng, 2*NDIM, MPI_INTEGER, pe, pe, comm3d, status, ierr)
               call MPI_Send(      ic_v(ic_rng(1):ic_rng(2), ic_rng(3):ic_rng(4), ic_rng(5):ic_rng(6)), &
                    &         size(ic_v(ic_rng(1):ic_rng(2), ic_rng(3):ic_rng(4), ic_rng(5):ic_rng(6))), &
                    &         MPI_DOUBLE_PRECISION, pe, pe, comm3d, ierr)
            enddo
         else
            call MPI_Send( [ ic_is, ic_ie, ic_js, ic_je, ic_ks, ic_ke ], 2*NDIM, MPI_INTEGER, 0, proc, comm3d, ierr)
            call MPI_Recv(      ic_data(ic_is:ic_ie, ic_js:ic_je, ic_ks:ic_ke, v), &
                 &         size(ic_data(ic_is:ic_ie, ic_js:ic_je, ic_ks:ic_ke, v)), &
                 &         MPI_DOUBLE_PRECISION, 0, proc, comm3d, status, ierr)
         endif
      enddo

      if (allocated(ic_v)) deallocate(ic_v)
      if (proc == 0) close(1)

      if (mass_mul /= 1.0) ic_data(:, :, :, 1) = ic_data(:, :, :, 1) * mass_mul

      do v = 2, 4 ! convert velocity to momentum
         ic_data(:, :, :, v) = ic_data(:, :, :, v) * ic_data(:, :, :, 1)
      enddo

      ! U = ( kB * T ) / (mean_mol_weight * (gamma - 1))
      ! cs2 = (gamma) * kB * T / mean_mol_weight.
      !   => cs2 = U * (gamma) * (gamma - 1)
      ic_data(:, :, :, 5) = ic_data(:, :, :, 5) * (gamma_loc - 1.0) * cs_mul ! * gamma_loc

      ! BEWARE: Until we have decent hdf5 restart that would be problem
      ! dependent, i.e. >=gcc-4.5 / >=ifort 10.1, things like
      ! that must be present in problem.par
      !mincs2 = minval(ic_data(:, :, :, 5))
      !maxcs2 = maxval(ic_data(:, :, :, 5))

   end subroutine read_IC_file

!-----------------------------------------------------------------------------

   subroutine init_prob

      use mpisetup,      only: proc, smalld
      use arrays,        only: u, b, cs_iso2_arr
      use grid,          only: is, ie, js, je, ks, ke, nx, ny, nz, nb, x, y, z, dx, dy, dz
      use initionized,   only: idni, imxi, imyi, imzi
      use list_hdf5,     only: additional_attrs, problem_write_restart, problem_read_restart
      use constants,     only: small, kboltz, mH
      use dataio_public, only: die, warn, printinfo, msg
      use types,         only: problem_customize_solution

      implicit none

      real, parameter :: beat_dx = 1e-5
      integer :: i, j, k, iic, jic, kic

      if (proc == 0) then
         if (max(dx, dy, dz) > ic_dx) then
            write(msg,'(a)')     "[initproblem:init_prob] Too low resolution" ! call die
            call warn(msg)
         endif
         if (abs(ic_dx/dx-anint(ic_dx/dx)) > beat_dx) then
            write(msg,'(a,f8.4)')"[initproblem:init_prob] X-direction requires interpolation ic_dx/dx= ", ic_dx/dx
            call warn(msg)
         endif
         if (abs(ic_dx/dy-anint(ic_dx/dy)) > beat_dx) then
            write(msg,'(a,f8.4)')"[initproblem:init_prob] Y-direction requires interpolation ic_dx/dy= ", ic_dx/dy
            call warn(msg)
         endif
         if (abs(ic_dx/dz-anint(ic_dx/dz)) > beat_dx) then
            write(msg,'(a,f8.4)')"[initproblem:init_prob] Z-direction requires interpolation ic_dx/dz= ", ic_dx/dz
            call warn(msg)
         endif
      endif

      if (fake_ic) then
         u(idni, is:ie, js:je, ks:ke) = 1.
         u(imxi, is:ie, js:je, ks:ke) = 0.
         u(imyi, is:ie, js:je, ks:ke) = 0.
         u(imzi, is:ie, js:je, ks:ke) = 0.
         cs_iso2_arr(is:ie, js:je, ks:ke) = 1e-2
      else
         call read_IC_file
         do k = ks, ke
            kic = nint((z(k) + ic_zsize/2.)/ic_dx)
            do j = js, je
               jic = nint((y(j) + ic_xysize/2.)/ic_dx)
               do i = is, ie
                  iic = nint((x(i) + ic_xysize/2.)/ic_dx)
                  if (iic >= ic_is .and. iic <= ic_ie .and. jic >= ic_js .and. jic <= ic_je .and. kic >= ic_ks .and. kic <= ic_ke) then
                     u(idni, i, j, k)     = ic_data(iic, jic, kic, 1) ! simple injection
                     u(imxi, i, j, k)     = ic_data(iic, jic, kic, 2)
                     u(imyi, i, j, k)     = ic_data(iic, jic, kic, 3)
                     u(imzi, i, j, k)     = ic_data(iic, jic, kic, 4)
      !               cs_iso2_arr(i, j, k) = ic_data(iic, jic, kic, 5)
                     cs_iso2_arr(i, j, k) = (gamma_loc) * kboltz * T_disk / mean_mol_weight / mH
                  else
                     u(idni, i, j, k)     = smalld
                     u(imxi, i, j, k)     = small
                     u(imyi, i, j, k)     = small
                     u(imzi, i, j, k)     = small
                     cs_iso2_arr(i, j, k) = mincs2
                  endif
               enddo
            enddo
         enddo
         if (allocated(ic_data)) deallocate(ic_data)
      endif

      do i = 1,nb
         u(:,i,:,:)               = u(:,nb+1,:,:)
         u(:,nx-nb+i,:,:)         = u(:,nx-nb,:,:)
         cs_iso2_arr(i,:,:)       = cs_iso2_arr(nb+1,:,:)
         cs_iso2_arr(nx-nb+i,:,:) = cs_iso2_arr(nx-nb,:,:)

         u(:,:,i,:)               = u(:,:,nb+1,:)
         u(:,:,ny-nb+i,:)         = u(:,:,ny-nb,:)
         cs_iso2_arr(:,i,:)       = cs_iso2_arr(:,nb+1,:)
         cs_iso2_arr(:,ny-nb+i,:) = cs_iso2_arr(:,ny-nb,:)

         u(:,:,:,i)               = u(:,:,:,nb+1)
         u(:,:,:,nz-nb+i)         = u(:,:,:,nz-nb)
         cs_iso2_arr(:,:,i)       = cs_iso2_arr(:,:,nb+1)
         cs_iso2_arr(:,:,nz-nb+i) = cs_iso2_arr(:,:,nz-nb)
      enddo
      if (proc == 0 ) then
         write(msg,'(2(a,g15.7))') '[initproblem:init_problem]: minval(dens)    = ', minval(u(idni,:,:,:)),      ' maxval(dens)    = ', maxval(u(idni,:,:,:))
         call printinfo(msg, .true.)
         write(msg,'(2(a,g15.7))') '[initproblem:init_problem]: minval(cs_iso2) = ', minval(cs_iso2_arr(:,:,:)), ' maxval(cs_iso2) = ', maxval(cs_iso2_arr(:,:,:))
         call printinfo(msg, .true.)
      endif

      b(:, 1:nx, 1:ny, 1:nz) = 0.0
      additional_attrs      => init_prob_attrs
      problem_write_restart => write_initial_fld_to_restart
      problem_read_restart  => read_initial_fld_from_restart

      ! BEWARE: den0, vlx0 and vly0 are used only with divine_intervention_type = 3
      if (.not.allocated(den0)) allocate(den0(nx,ny,nz))
      if (.not.allocated(vlx0)) allocate(vlx0(nx,ny,nz))
      if (.not.allocated(vly0)) allocate(vly0(nx,ny,nz))

      den0 = u(idni,:,:,:)
      vlx0 = u(imxi,:,:,:) / den0
      vly0 = u(imyi,:,:,:) / den0

      ! It would be cool to dump a restart file here but this would make a cyclic dependency

      problem_customize_solution => problem_customize_solution_wt4

#ifndef UMUSCL
      if (proc == 0 ) call warn("[initproblem:init_problem]: Without UMUSCL you'll likely get Monet-like density maps.")
#endif /* !UMUSCL */

      return
   end subroutine init_prob

!-----------------------------------------------------------------------------

   subroutine init_prob_attrs(file_id)

      use hdf5, only: HID_T, SIZE_T
      use h5lt, only: h5ltset_attribute_double_f
      use constants, only: fpiG

      implicit none

      integer(HID_T),intent(in)  :: file_id
      integer(SIZE_T) :: bufsize = 1
      integer :: error

      call h5ltset_attribute_double_f(file_id, "/", "fpiG", [fpiG], bufsize,error)

   end subroutine init_prob_attrs

!-----------------------------------------------------------------------------

   subroutine write_initial_fld_to_restart(file_id)

      use hdf5,        only: HID_T
      use grid,        only: nx, ny, nz
      use dataio_hdf5, only: write_3darr_to_restart

      implicit none

      integer(HID_T),intent(in)  :: file_id

      if ( divine_intervention_type == 3) then
        if (allocated(den0)) call write_3darr_to_restart(den0(:,:,:), file_id, "den0", nx, ny, nz)
        if (allocated(vlx0)) call write_3darr_to_restart(vlx0(:,:,:), file_id, "vlx0", nx, ny, nz)
        if (allocated(vly0)) call write_3darr_to_restart(vly0(:,:,:), file_id, "vly0", nx, ny, nz)
      endif

   end subroutine write_initial_fld_to_restart

!-----------------------------------------------------------------------------

   subroutine read_initial_fld_from_restart(file_id)

      use hdf5,        only: HID_T
      use grid,        only: nx, ny, nz
      use dataio_hdf5, only: read_3darr_from_restart

      implicit none

      integer(HID_T),intent(in) :: file_id

      real, dimension(:,:,:), pointer :: p3d

      ! /todo First query for existence of den0, vlx0 and vly0, then allocate
      if (divine_intervention_type == 3) then
         if (.not.allocated(den0)) allocate(den0(nx,ny,nz))
         if (.not.allocated(vlx0)) allocate(vlx0(nx,ny,nz))
         if (.not.allocated(vly0)) allocate(vly0(nx,ny,nz))

         if (.not.associated(p3d)) p3d => den0(:,:,:)
         call read_3darr_from_restart(file_id,"den0",p3d,nx,ny,nz)
         if (associated(p3d)) nullify(p3d)
         if (.not.associated(p3d)) p3d => vlx0(:,:,:)
         call read_3darr_from_restart(file_id,"vlx0",p3d,nx,ny,nz)
         if (associated(p3d)) nullify(p3d)
         if (.not.associated(p3d)) p3d => vly0(:,:,:)
         call read_3darr_from_restart(file_id,"vly0",p3d,nx,ny,nz)
         if (associated(p3d)) nullify(p3d)
      endif

   end subroutine read_initial_fld_from_restart

!-----------------------------------------------------------------------------

   subroutine problem_customize_solution_wt4

     use mpisetup,    only: proc
     use arrays,      only: u, cs_iso2_arr
     use grid,        only: is, ie, js, je, ks, ke, nx, ny, nz, x, y
     use initionized, only: idni, imxi, imyi, imzi

     implicit none

     integer :: i, j, k
     real, allocatable, dimension(:) :: mod_str
     real, parameter :: max_ambient = 100. ! do not modify solution if density is above max_ambient * ambient_density
     real, allocatable, dimension(:,:) :: alf
     real            :: rc, ambient_density_min

     if (.not. allocated(mod_str)) allocate(mod_str(is:ie)) !BEWARE not deallocated anywhere yet

     select case (divine_intervention_type)
        case (1)                                                                                ! crude
           where (u(idni, is:ie, js:je, ks:ke) < ambient_density)
              u(imxi, is:ie, js:je, ks:ke) = (1. - damp_factor) * u(imxi, is:ie, js:je, ks:ke)
              u(imyi, is:ie, js:je, ks:ke) = (1. - damp_factor) * u(imyi, is:ie, js:je, ks:ke)
              u(imzi, is:ie, js:je, ks:ke) = (1. - damp_factor) * u(imzi, is:ie, js:je, ks:ke)
              cs_iso2_arr(is:ie, js:je, ks:ke) = mincs2
           elsewhere
              cs_iso2_arr(is:ie, js:je, ks:ke) = maxcs2
           endwhere
        case (2)                                                                                ! smooth
           ambient_density_min = ambient_density / max_ambient
           do k = ks, ke
              do j = js, je
                 mod_str(is:ie) = max(0., (1. + 1./max_ambient) * ambient_density_min / (max(0., u(idni, is:ie, j, k)) + ambient_density_min) - 1./max_ambient)
                 ! ifort can have memory leaks on WHERE - let's provide explicit loop for this crappy compiler
                 ! The __IFORT__ macro has to be defined manually, e.g. in appropriate compiler.in file
#ifndef __IFORT__
                 where (mod_str(is:ie) > max_ambient**(-2))
                    u(idni,     is:ie, j, k) = u(idni, is:ie, j, k) + ambient_density_min * mod_str(is:ie)
                    u(imxi,     is:ie, j, k) = u(imxi, is:ie, j, k) * (1. - damp_factor   * mod_str(is:ie))
                    u(imyi,     is:ie, j, k) = u(imyi, is:ie, j, k) * (1. - damp_factor   * mod_str(is:ie))
                    u(imzi,     is:ie, j, k) = u(imzi, is:ie, j, k) * (1. - damp_factor   * mod_str(is:ie))
                    cs_iso2_arr(is:ie, j, k) = maxcs2               -  (maxcs2-mincs2)    * mod_str(is:ie)
                 endwhere
#else /* !__IFORT__ */
                 do i = is, ie
                    if (mod_str(i) > max_ambient**(-2)) then
                       u(idni,     i, j, k) = u(idni, i, j, k) + ambient_density_min * mod_str(i)
                       u(imxi,     i, j, k) = u(imxi, i, j, k) * (1. - damp_factor   * mod_str(i))
                       u(imyi,     i, j, k) = u(imyi, i, j, k) * (1. - damp_factor   * mod_str(i))
                       u(imzi,     i, j, k) = u(imzi, i, j, k) * (1. - damp_factor   * mod_str(i))
                       cs_iso2_arr(i, j, k) = maxcs2           -  (maxcs2-mincs2)    * mod_str(i)
                    endif
                 enddo
#endif /* !__IFORT__ */
              enddo
           enddo
        case (3)
          if (.not. allocated(alf)) then
             allocate(alf(nx,ny))
             do i = 1,nx
                do j = 1,ny
                   rc = sqrt(x(i)**2 + y(j)**2)
                   alf(i,j) = -alfasupp*0.5*(tanh((rc-r_in)/r_in*f_in)-1.)
                   alf(i,j) = alf(i,j) + alfasupp*0.5*(tanh((rc-r_out)/r_out*f_out) + 1.)
                enddo
             enddo
          endif
          do k = 1,nz
             u(idni, :, :, k) = (1. - alf(:,:))*u(idni, :, :, k) + alf*den0(:, :, k)
             u(imxi, :, :, k) = (1. - alf(:,:))*u(imxi, :, :, k) + alf*den0(:, :, k) * vlx0(:, :, k)
             u(imyi, :, :, k) = (1. - alf(:,:))*u(imyi, :, :, k) + alf*den0(:, :, k) * vly0(:, :, k)
          enddo
          do k = ks, ke
             do j = js, je
                mod_str(is:ie) = max(0., (1. + 1./max_ambient) * ambient_density / (max(0., u(idni, is:ie, j, k)) + ambient_density) - 1./max_ambient)
                where (mod_str(is:ie) > max_ambient**(-2))
                   u(imzi,     is:ie, j, k) = u(imzi, is:ie, j, k) * (1. - damp_factor * mod_str(is:ie))
                endwhere
             enddo
          enddo
     end select

   end subroutine problem_customize_solution_wt4
end module initproblem
