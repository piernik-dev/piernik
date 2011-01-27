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

   use mpisetup,    only: cbuff_len
   implicit none

   private
   public :: read_problem_par, init_prob

   real                     :: d0, r_max, dout, alpha, r_in, r_out, f_in, f_out
   real                     :: dens_exp      !< exponent in profile density \f$\rho(R) = \rho_0 R^{-k}\f$
   real                     :: dens_amb      !< density of ambient medium (used for inner cutoff)
   real                     :: eps           !< dust to gas ratio
   real                     :: x_cut         !< radius of inner disk cut-off
   integer                  :: cutoff_ncells !< width of cut-off profile
   real, save               :: T_inner = 0.0 !< Orbital period at the inner boundary, \todo save it to restart as an attribute
   !>
   !! \f$\tau\f$ in \f$\frac{Du}{Dt} = - \frac{u-u_0}{\tau}f(R)
   !! when initproblem::problem_customize_solution is used
   !<
   real                     :: dumping_coeff
   logical                  :: use_inner_orbital_period  !< use 1./T_inner as dumping_coeff
   character(len=cbuff_len) :: mag_field_orient
   real, target, allocatable, dimension(:,:,:,:) :: den0, mtx0, mty0, mtz0, ene0
   integer, parameter       :: dname_len = 10

   namelist /PROBLEM_CONTROL/  alpha, d0, dout, r_max, mag_field_orient, r_in, r_out, f_in, f_out, &
      & dens_exp, eps, dens_amb, x_cut, cutoff_ncells, dumping_coeff, use_inner_orbital_period

contains
!-----------------------------------------------------------------------------
   subroutine read_problem_par
      use dataio_pub,          only: ierrh, par_file, namelist_errh, compare_namelist, cmdl_nml      ! QA_WARN required for diff_nml
      use mpisetup,            only: cbuff, rbuff, ibuff, lbuff, buffer_dim, master, slave, comm, ierr
      use mpi,                 only: MPI_CHARACTER, MPI_DOUBLE_PRECISION, MPI_INTEGER, MPI_LOGICAL
      use gravity,             only: grav_pot_3d
      use grid,                only: geometry
      use types,               only: problem_customize_solution, problem_grace_passed
      use list_hdf5,           only: problem_write_restart, problem_read_restart
      use fluidboundaries_pub, only: user_bnd_xl, user_bnd_xr
      implicit none

      d0               = 1.0
      dout             = 1.0e-4
      r_max            = 1.0
      mag_field_orient = 'none'
      alpha            = 1.0

      r_in             = 0.5
      f_in             = 10.0
      r_out            = 2.1
      f_out            = 0.0

      dens_exp         = 0.0
      dens_amb         = 1.e-3
      eps              = 1.0
      x_cut            = 0.6

      cutoff_ncells    = 8.0
      dumping_coeff    = 1.0

      use_inner_orbital_period = .false.

      if (master) then

         diff_nml(PROBLEM_CONTROL)

         ibuff(1) = cutoff_ncells

         lbuff(1) = use_inner_orbital_period

         cbuff(1) = mag_field_orient

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

      endif

      call MPI_Bcast(cbuff, cbuff_len*buffer_dim, MPI_CHARACTER,        0, comm, ierr)
      call MPI_Bcast(rbuff,           buffer_dim, MPI_DOUBLE_PRECISION, 0, comm, ierr)
      call MPI_Bcast(ibuff,           buffer_dim, MPI_INTEGER         , 0, comm, ierr)
      call MPI_Bcast(lbuff,           buffer_dim, MPI_LOGICAL         , 0, comm, ierr)

      if (slave) then

         cutoff_ncells    = ibuff(1)

         use_inner_orbital_period = lbuff(1)

         mag_field_orient = cbuff(1)

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

      endif

      if (geometry=="cylindrical") then
         problem_write_restart => write_initial_fld_to_restart
         problem_read_restart  => read_initial_fld_from_restart
         problem_customize_solution => problem_customize_solution_kepler
         user_bnd_xl => my_bnd_xl
         user_bnd_xr => my_bnd_xr
         grav_pot_3d => my_grav_pot_3d
!         problem_grace_passed => add_random_noise
         problem_grace_passed => add_sine
      endif

   end subroutine read_problem_par
!-----------------------------------------------------------------------------
   subroutine add_sine
      use mpisetup,     only: master, dom
      use dataio_pub,   only: printinfo
      use fluidindex,   only: flind
      use grid,         only: cg
      use arrays,       only: u
      use fluidindex,   only: flind
      use constants,    only: dpi
      implicit none
      real, parameter :: amp = 1.e-6
      real :: kx, kz
      integer :: i, k

      if (master) call printinfo("[initproblem:add_sine]: adding sine wave perturbation to dust")

      kx = 15.*dpi/dom%Lx
      kz =  5.*dpi/dom%Lz

      do i = cg%is, cg%ie
         do k = cg%ks, cg%ke
            u(flind%dst%imx,i,:,k) = u(flind%dst%imx,i,:,k) + amp*sin(kx*cg%x(i) + kz*cg%z(k)) * u(flind%dst%idn,i,:,k)
            u(flind%dst%imy,i,:,k) = u(flind%dst%imy,i,:,k) + amp*sin(kx*cg%x(i) + kz*cg%z(k)) * u(flind%dst%idn,i,:,k)
            u(flind%dst%imz,i,:,k) = u(flind%dst%imz,i,:,k) + amp*sin(kx*cg%x(i) + kz*cg%z(k)) * u(flind%dst%idn,i,:,k)
         enddo
      enddo
#ifdef DEBUG
      open(456, file="perturbation.dat", status="unknown")
      do k = cg%ks, cg%ke
         write(456,*) cg%z(k), sin(kz*cg%z(k))
      enddo
      close(456)
#endif /* DEBUG */

   end subroutine add_sine

   subroutine add_random_noise
      use grid,         only: cg
      use arrays,       only: u
      use fluidindex,   only: flind
      use mpisetup,     only: proc, master
      use dataio_pub,   only: printinfo

      implicit none
      integer, dimension(:), allocatable  :: seed
      integer :: n, clock, i
      real, dimension(:,:,:,:), allocatable :: noise
      real, parameter :: amp = 1.e-6

      if (master) call printinfo("[initproblem:add_random_noise]: adding random noise to dust")
      call random_seed(size=n)
      allocate(seed(n))
      call system_clock(count=clock)
      seed = clock*proc + 37 * (/ (i-1, i = 1, n) /)
      call random_seed(put=seed)
      deallocate(seed)

      allocate(noise(3,cg%nx,cg%ny,cg%nz))
      call random_number(noise)
      u(flind%dst%imx,:,:,:) = u(flind%dst%imx,:,:,:) +amp -2.0*amp*noise(1,:,:,:) * u(flind%dst%idn,:,:,:)
      u(flind%dst%imy,:,:,:) = u(flind%dst%imy,:,:,:) +amp -2.0*amp*noise(2,:,:,:) * u(flind%dst%idn,:,:,:)
      u(flind%dst%imz,:,:,:) = u(flind%dst%imz,:,:,:) +amp -2.0*amp*noise(3,:,:,:) * u(flind%dst%idn,:,:,:)
      deallocate(noise)
   end subroutine add_random_noise
!-----------------------------------------------------------------------------
   subroutine init_prob

      use dataio_pub,          only: msg, printinfo
      use types,               only: component_fluid
      use arrays,              only: u, b, dprof
      use constants,           only: newtong, gram, cm, kboltz, mH, dpi
      use fluidindex,          only: ibx, iby, ibz, flind
      use gravity,             only: r_smooth, r_grav, n_gravr, ptmass, source_terms_grav, grav_pot2accel, grav_pot_3d
      use grid,                only: cg, geometry
      use hydrostatic,         only: hydrostatic_zeq_densmid
      use mpisetup,            only: zdim, has_dir, dom, master
      use types,               only: component_fluid
      use dataio_pub,          only: die

      implicit none

      integer :: i, j, k, kmid, p, middle_of_nx
      integer, dimension(1) :: n_x_cut
      real    :: xi, yj, zk, rc, vx, vy, vz, b0, sqr_gm, vr, vphi
      real    :: csim2, gprim, H2

      real, dimension(:), allocatable :: grav, dens_prof, dens_cutoff, ln_dens_der
      type(component_fluid), pointer  :: fl

!   Secondary parameters

      sqr_gm = sqrt(newtong*ptmass)
      do k = 1, cg%nz
         if (cg%z(k) < 0.0) kmid = k       ! the midplane is in between ksmid and ksmid+1
      enddo

      if (associated(flind%ion) .and. geometry=='cartesian') then
         fl => flind%ion
         csim2 = fl%cs2*(1.0+alpha)
         b0    = sqrt(2.*alpha*d0*fl%cs2)

         do j = 1, cg%ny
            yj = cg%y(j)
            do i = 1, cg%nx
               xi = cg%x(i)
               rc = sqrt(xi**2+yj**2)

               if (has_dir(zdim)) call hydrostatic_zeq_densmid(i, j, d0, csim2)

               do k = 1, cg%nz

                  vx = sqr_gm * (-yj)/(rc**2+r_smooth**2)**0.75
                  vy = sqr_gm * ( xi)/(rc**2+r_smooth**2)**0.75
                  vz = 0.0

                  u(fl%idn,i,j,k) = min((rc/r_grav)**n_gravr,100.0)

                  if (has_dir(zdim)) then
                     u(fl%idn,i,j,k) = dprof(k)/cosh(u(fl%idn,i,j,k))
                     u(fl%idn,i,j,k) = max(u(fl%idn,i,j,k), dout)
                  else
                     u(fl%idn,i,j,k) = dout + (d0 - dout)/cosh(u(fl%idn,i,j,k))
                  endif
                  u(fl%idn,i,j,k) = u(fl%idn,i,j,k)
                  u(fl%imx,i,j,k) = vx*u(fl%idn,i,j,k)
                  u(fl%imy,i,j,k) = vy*u(fl%idn,i,j,k)
                  u(fl%imz,i,j,k) = vz*u(fl%idn,i,j,k)
                  if (fl%ien > 0) then
                     u(fl%ien,i,j,k) = fl%cs2/(fl%gam_1)*u(fl%idn,i,j,k)
!                     u(fl%ien,i,j,k) = max(u(fl%ien,i,j,k), smallei)
                     u(fl%ien,i,j,k) = u(fl%ien,i,j,k) +0.5*(vx**2+vy**2+vz**2)*u(fl%idn,i,j,k)
                  endif
                  if (trim(mag_field_orient) .eq. 'toroidal') then
                     b(ibx,i,j,k)   = -b0*sqrt(u(fl%idn,i,j,k)/d0)*yj/rc
                     b(iby,i,j,k)   =  b0*sqrt(u(fl%idn,i,j,k)/d0)*xi/rc
                     b(ibz,i,j,k)   =  0.0
                  else if (trim(mag_field_orient) .eq. 'vertical') then
                     b(ibx,i,j,k)   =  0.0
                     b(iby,i,j,k)   =  0.0
                     b(ibz,i,j,k)   =  b0
                  else if (trim(mag_field_orient) .eq. 'none') then
                     b(:,i,j,k)     =  0.0
                  endif

                  if (fl%ien > 0) u(fl%ien,i,j,k)   = u(fl%ien,i,j,k) + 0.5*sum(b(:,i,j,k)**2,1)
               enddo
            enddo
         enddo
      else if (geometry=='cylindrical') then
         if (master) then
            call printinfo("------------------------------------------------------------------")
            call printinfo(" Assuming temperature profile for MMSN ")
            call printinfo(" T(R) = 150 ( R / 1 AU )^(-0.429) K")
            write(msg,'(A,F5.1,A)') " T(xmin) = ", mmsn_T(dom%xmin)," K"
            call printinfo(msg)
            write(msg,'(A,F5.1,A)') " T(xmax) = ", mmsn_T(dom%xmax)," K"
            call printinfo(msg)
            write(msg,'(A,F5.1,A)') " T_mean  = ", 0.5*(mmsn_T(dom%xmin)+mmsn_T(dom%xmax))," K"
            call printinfo(msg)
            write(msg,'(A,F9.5)') " cs2(T_mean) = ", kboltz * 0.5*(mmsn_T(dom%xmin)+mmsn_T(dom%xmax)) / mH
            call printinfo(msg)
            write(msg,'(A,ES12.3,A)') " T_real(cs2) = ", flind%neu%cs2*mH/kboltz, " K"
            call printinfo(msg)
            call printinfo("------------------------------------------------------------------")
         endif
         call grav_pot_3d

         if (.not.allocated(den0)) allocate(den0(flind%fluids, cg%nx, cg%ny, cg%nz))
         if (.not.allocated(mtx0)) allocate(mtx0(flind%fluids, cg%nx, cg%ny, cg%nz))
         if (.not.allocated(mty0)) allocate(mty0(flind%fluids, cg%nx, cg%ny, cg%nz))
         if (.not.allocated(mtz0)) allocate(mtz0(flind%fluids, cg%nx, cg%ny, cg%nz))
         if (.not.allocated(ene0)) allocate(ene0(flind%fluids, cg%nx, cg%ny, cg%nz))

         if (.not.allocated(grav)) allocate(grav(cg%nx))
         if (.not.allocated(ln_dens_der)) allocate(ln_dens_der(cg%nx))
         if (.not.allocated(dens_prof)) allocate(dens_prof(cg%nx))
         if (.not.allocated(dens_cutoff)) allocate(dens_cutoff(cg%nx))

         call source_terms_grav
         call grav_pot2accel('xsweep',1,1, cg%nx, grav, 1)

         dens_prof(:) = d0 * cg%x(:)**(-dens_exp)  * gram / cm**2

         middle_of_nx = cg%nx/2 + 1
         n_x_cut      = maxloc(cg%x, mask=cg%x<=x_cut)
!         dens_prof    = min(dens_prof, get_lcutoff(cutoff_ncells, (middle_of_nx - n_x_cut(1)), cg%nx, dens_amb, dens_prof(n_x_cut(1))) )
         dens_prof    = dens_prof * get_lcutoff(cutoff_ncells, (middle_of_nx - n_x_cut(1)), cg%nx, 0.0, 1.0) + dens_amb

         !! \f$ v_\phi = \sqrt{R\left(c_s^2 \partial_R \ln\rho + \partial_R \Phi \right)} \f$
         ln_dens_der  = log(dens_prof)
         ln_dens_der(2:cg%nx)  = ( ln_dens_der(2:cg%nx) - ln_dens_der(1:cg%nx-1) ) / cg%dx
         ln_dens_der(1)        = ln_dens_der(2)
         T_inner               = dpi*cg%x(cg%is) / sqrt( abs(grav(cg%is)) * cg%x(cg%is) )
         write(msg,*) "III Kepler Law gives T = ", sqr_gm/dpi , " yr at 1 AU"
         if (master) call printinfo(msg)
#ifdef DEBUG
         open(143,file="dens_prof.dat",status="unknown")
         do p = 1, cg%nx
            write(143,'(4(ES14.4,1X))') cg%x(p), dens_prof(p), sqrt( max(cg%x(p)*(flind%neu%cs2*ln_dens_der(p) + abs(grav(p))),0.0) ), &
               sqrt( max(abs(grav(p)) * cg%x(p) - flind%neu%cs2*dens_exp,0.0))
         enddo
         close(143)
#endif /* DEBUG */

         do p = 1, flind%fluids
            fl => flind%all_fluids(p)
            if (fl%tag /= "DST" .and. master) then
               write(msg,'(A,F9.5)') "[init_problem:initprob] cs2 used = ", fl%cs2
               call printinfo(msg)
            endif

            do j = 1, cg%ny
               yj = cg%y(j)
               do i = 1, cg%nx
                  xi = cg%x(i)
                  rc = xi + r_smooth

                  gprim = newtong*ptmass / xi**3
                  if (fl%cs > 0) then
                     H2 = 2.0*fl%cs2/gprim
                  else
                     H2 = 1.0
                  endif

                  do k = 1, cg%nz
                     zk = cg%z(k)
!                     u(fl%idn,i,j,k) = max(d0*(1./cosh((xi/r_max)**10)) * exp(-zk**2/H2),1.e-10))
                     u(fl%idn,i,j,k) = dens_prof(i)
                     if (fl%tag == "DST") u(fl%idn,i,j,k) = eps * u(fl%idn,i,j,k)

                     vr   = 0.0
                     ! that condition is not necessary since cs2 == 0.0 for dust
                     if (fl%tag /= "DST") then
!                        vphi = sqrt( max(abs(grav(i)) * rc - fl%cs2*dens_exp,0.0))
                         vphi = sqrt( max(cg%x(i)*(fl%cs2*ln_dens_der(i) + abs(grav(i))),0.0) )
                     else
                        vphi = sqrt( max(abs(grav(i)) * rc, 0.0))
                     endif
                     vz   = 0.0

                     u(fl%imx,i,j,k) = vr   * u(fl%idn,i,j,k)
                     u(fl%imy,i,j,k) = vphi * u(fl%idn,i,j,k)
                     u(fl%imz,i,j,k) = vz   * u(fl%idn,i,j,k)
                     if (fl%ien > 0) then
                        u(fl%ien,i,j,k) = fl%cs2/(fl%gam_1)*u(fl%idn,i,j,k)
                        u(fl%ien,i,j,k) = u(fl%ien,i,j,k) + 0.5*(vr**2+vphi**2+vz**2)*u(fl%idn,i,j,k)
                        ene0(p,i,j,k)   = u(fl%ien,i,j,k)
                     else
                        ene0(p,i,j,k)   = 0.0
                     endif
                  enddo
               enddo
            enddo

            den0(p,:,:,:) = u(fl%idn,:,:,:)
            mtx0(p,:,:,:) = u(fl%imx,:,:,:)
            mty0(p,:,:,:) = u(fl%imy,:,:,:)
            mtz0(p,:,:,:) = u(fl%imz,:,:,:)
         enddo
         b = 0.0
         if (allocated(grav)) deallocate(grav)
         if (allocated(dens_prof)) deallocate(dens_prof)
      else
         call die("[initproblem:init_prob] I don't know what to do... :/")
      endif

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
   subroutine write_initial_fld_to_restart(file_id)

      use hdf5,        only: HID_T
      use grid,        only: cg
      use dataio_hdf5, only: write_3darr_to_restart
      use fluidindex,  only: flind

      implicit none

      integer(HID_T),intent(in)  :: file_id
      integer :: i
      character(len=dname_len) :: dname

      do i = LBOUND(den0,1), UBOUND(den0,1)
         write(dname,'(2a)') flind%all_fluids(i)%tag, '_den0'
         if (allocated(den0)) call write_3darr_to_restart(den0(i,:,:,:), file_id, dname, cg%nx, cg%ny, cg%nz)
         write(dname,'(2a)') flind%all_fluids(i)%tag, '_mtx0'
         if (allocated(mtx0)) call write_3darr_to_restart(mtx0(i,:,:,:), file_id, dname, cg%nx, cg%ny, cg%nz)
         write(dname,'(2a)') flind%all_fluids(i)%tag, '_mty0'
         if (allocated(mty0)) call write_3darr_to_restart(mty0(i,:,:,:), file_id, dname, cg%nx, cg%ny, cg%nz)
         write(dname,'(2a)') flind%all_fluids(i)%tag, '_mtz0'
         if (allocated(mtz0)) call write_3darr_to_restart(mtz0(i,:,:,:), file_id, dname, cg%nx, cg%ny, cg%nz)
         write(dname,'(2a)') flind%all_fluids(i)%tag, '_ene0'
         if (allocated(ene0)) call write_3darr_to_restart(ene0(i,:,:,:), file_id, dname, cg%nx, cg%ny, cg%nz)
      enddo

   end subroutine write_initial_fld_to_restart
!-----------------------------------------------------------------------------
   subroutine read_initial_fld_from_restart(file_id)

      use hdf5,        only: HID_T
      use grid,        only: cg
      use fluidindex,  only: flind
      use dataio_hdf5, only: read_3darr_from_restart
      use fluidindex,  only: flind

      implicit none

      integer(HID_T),intent(in) :: file_id

      character(len=dname_len) :: dname
      real, dimension(:,:,:), pointer :: p3d
      integer :: i

      ! /todo First query for existence of den0, vlx0 and vly0, then allocate
      if (.not.allocated(den0)) allocate(den0(flind%fluids, cg%nx, cg%ny, cg%nz))
      if (.not.allocated(mtx0)) allocate(mtx0(flind%fluids, cg%nx, cg%ny, cg%nz))
      if (.not.allocated(mty0)) allocate(mty0(flind%fluids, cg%nx, cg%ny, cg%nz))
      if (.not.allocated(mtz0)) allocate(mtz0(flind%fluids, cg%nx, cg%ny, cg%nz))
      if (.not.allocated(ene0)) allocate(ene0(flind%fluids, cg%nx, cg%ny, cg%nz))

      do i=1, flind%fluids
         write(dname,'(2a)') flind%all_fluids(i)%tag, '_den0'
         if (.not.associated(p3d)) p3d => den0(i,:,:,:)
         call read_3darr_from_restart(file_id,dname,p3d, cg%nx, cg%ny, cg%nz)
         if (associated(p3d)) nullify(p3d)

         write(dname,'(2a)') flind%all_fluids(i)%tag, '_mtx0'
         if (.not.associated(p3d)) p3d => mtx0(i,:,:,:)
         call read_3darr_from_restart(file_id,dname,p3d, cg%nx, cg%ny, cg%nz)
         if (associated(p3d)) nullify(p3d)

         write(dname,'(2a)') flind%all_fluids(i)%tag, '_mty0'
         if (.not.associated(p3d)) p3d => mty0(i,:,:,:)
         call read_3darr_from_restart(file_id,dname,p3d, cg%nx, cg%ny, cg%nz)
         if (associated(p3d)) nullify(p3d)

         write(dname,'(2a)') flind%all_fluids(i)%tag, '_mtz0'
         if (.not.associated(p3d)) p3d => mtz0(i,:,:,:)
         call read_3darr_from_restart(file_id,dname,p3d, cg%nx, cg%ny, cg%nz)
         if (associated(p3d)) nullify(p3d)

         write(dname,'(2a)') flind%all_fluids(i)%tag, '_ene0'
         if (.not.associated(p3d)) p3d => ene0(i,:,:,:)
         call read_3darr_from_restart(file_id,dname,p3d, cg%nx, cg%ny, cg%nz)
         if (associated(p3d)) nullify(p3d)
      enddo

   end subroutine read_initial_fld_from_restart
!-----------------------------------------------------------------------------
   subroutine problem_customize_solution_kepler
      use mpisetup,        only: dt, t, grace_period_passed, relax_time
      use arrays,          only: u
      use grid,            only: cg
      use fluidboundaries, only: all_fluid_boundaries
      use fluidindex,      only: iarr_all_dn, iarr_all_mx, iarr_all_my, iarr_all_mz
#ifndef ISO
      use fluidindex,      only: iarr_all_en
#endif /* ISO */
      use interactions,    only: dragc_gas_dust
#ifdef VERBOSE
      use dataio_pub,      only: msg, printinfo
      use mpisetup,        only: master
#endif /* VERBOSE */
      implicit none
      integer                               :: j, k
      logical, save                         :: frun = .true.
      real, dimension(:,:), allocatable, save :: funcR
      real, save :: x0, x1, y0, y1, a, b

      if (grace_period_passed() .and. t <= x1) then
         dragc_gas_dust = a*t+b
#ifdef VERBOSE
         if (master) write(msg,'(A,F6.1)') 'dragc_gas_dust = ', dragc_gas_dust
         call printinfo(msg)
#endif /* VERBOSE */
      endif

      if (frun) then
         x0 = relax_time + 2.0
         x1 = x0 + 30.0
         y0 = 100.0
         y1 = 1.0
         a = (y0 - y1)/(x0 - x1)
         b = y0 - a*x0
         allocate(funcR(size(iarr_all_dn), cg%nx) )

         funcR(1,:) = -tanh((cg%x(:)-r_in+1.0)**f_in) + 1.0 + max( tanh((cg%x(:)-r_out+1.0)**f_out), 0.0)

         if (use_inner_orbital_period) then
            funcR(1,:) = funcR(1,:) / T_inner
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

      do j = 1, cg%ny
         do k = 1, cg%nz
            u(iarr_all_dn,:,j,k) = u(iarr_all_dn,:,j,k) - dt*(u(iarr_all_dn,:,j,k) - den0(:,:,j,k))*funcR(:,:)
            u(iarr_all_mx,:,j,k) = u(iarr_all_mx,:,j,k) - dt*(u(iarr_all_mx,:,j,k) - mtx0(:,:,j,k))*funcR(:,:)
            u(iarr_all_my,:,j,k) = u(iarr_all_my,:,j,k) - dt*(u(iarr_all_my,:,j,k) - mty0(:,:,j,k))*funcR(:,:)
            u(iarr_all_mz,:,j,k) = u(iarr_all_mz,:,j,k) - dt*(u(iarr_all_mz,:,j,k) - mtz0(:,:,j,k))*funcR(:,:)
#ifndef ISO
            u(iarr_all_en,:,j,k) = u(iarr_all_en,:,j,k) - dt*(u(iarr_all_en,:,j,k) - ene0(:,:,j,k)*funcR(:,:)
#endif
         enddo
      enddo
      call all_fluid_boundaries
   end subroutine problem_customize_solution_kepler
!-----------------------------------------------------------------------------
   subroutine my_grav_pot_3d
      use constants, only: newtong
      use gravity,   only: ptmass, sum_potential
      use arrays,    only: gp
      use grid,      only: cg
      implicit none
      logical, save :: frun = .true.
      real          :: r2
      integer       :: i, k

      if (frun) then
         do i = 1, cg%nx
            do k = 1, cg%nz
               r2 = cg%x(i)**2! + cg%z(k)**2
               gp(i,:,k) = -newtong*ptmass / sqrt(r2)
            enddo
         enddo
      endif

      frun = .false.
      call sum_potential

   end subroutine my_grav_pot_3d
!-----------------------------------------------------------------------------
   subroutine my_bnd_xl
      use grid,         only: cg
      use arrays,       only: u
      use gravity,      only: grav_pot2accel
      use fluidindex,   only: iarr_all_dn, iarr_all_mx, iarr_all_my, iarr_all_mz, flind
#ifndef ISO
      use fluidindex,   only: iarr_all_en
#endif /* ISO */
      implicit none
      integer :: i, p
      real, dimension(cg%nx) :: grav
      real, dimension(size(iarr_all_my), cg%ny, cg%nz) :: vy,vym
      real, dimension(size(flind%all_fluids))    :: cs2_arr
      integer, dimension(size(flind%all_fluids)) :: ind_cs2

      do i = 1, size(flind%all_fluids)
         ind_cs2    = i
         cs2_arr(i) = flind%all_fluids(i)%cs2
      enddo

      call grav_pot2accel('xsweep',1,1, cg%nx, grav, 1)

      do i = 1, cg%nb
         u(iarr_all_dn,i,:,:) = u(iarr_all_dn, cg%is,:,:)
         u(iarr_all_mx,i,:,:) = max(0.0,u(iarr_all_mx, cg%is,:,:))
         do p = 1, size(flind%all_fluids)
            u(iarr_all_my(p),i,:,:) = sqrt( abs(grav(i)) * cg%x(i) - cs2_arr(p)*dens_exp) *  u(iarr_all_dn(p),i,:,:)
         enddo
         u(iarr_all_mz,i,:,:) = u(iarr_all_mz, cg%is,:,:)
#ifndef ISO
         u(iarr_all_en,i,:,:) = u(iarr_all_en, cg%is,:,:)
#endif
      enddo

      do i = cg%nb,1,-1
         vym(:,:,:) = u(iarr_all_my,i+2,:,:)/u(iarr_all_dn,i+1,:,:)
         vy(:,:,:)  = u(iarr_all_my,i+1,:,:)/u(iarr_all_dn,i+1,:,:)
         u(iarr_all_my,i,:,:) = (vym(:,:,:) + (cg%x(i) - cg%x(i+2)) / (cg%x(i+1) - cg%x(i+2)) * (vy - vym))*u(iarr_all_dn,i,:,:)
      enddo

   end subroutine my_bnd_xl
!-----------------------------------------------------------------------------
   subroutine my_bnd_xr
      use grid,   only: cg
      use arrays, only: u
      use fluidindex,  only: iarr_all_dn, iarr_all_mx, iarr_all_my, iarr_all_mz
#ifndef ISO
      use fluidindex,  only: iarr_all_en
#endif /* ISO */
      implicit none

      u(iarr_all_dn, cg%ie+1:cg%nx,:,:) = den0(:, cg%ie+1:cg%nx,:,:)
      u(iarr_all_mx, cg%ie+1:cg%nx,:,:) = mtx0(:, cg%ie+1:cg%nx,:,:)
      u(iarr_all_my, cg%ie+1:cg%nx,:,:) = mty0(:, cg%ie+1:cg%nx,:,:)
      u(iarr_all_mz, cg%ie+1:cg%nx,:,:) = mtz0(:, cg%ie+1:cg%nx,:,:)
#ifndef ISO
      u(iarr_all_en, cg%ie+1:cg%nx,:,:) = ene0(:, cg%ie+1:cg%nx,:,:)
#endif
   end subroutine my_bnd_xr
!-----------------------------------------------------------------------------
   function get_lcutoff(width,dist,n,vmin,vmax) result(y)
      implicit none
      integer, intent(in) :: width  !< width of tanh profile [cells]
      integer, intent(in) :: dist   !< distance between the expected position of the profile and the middle cell of the domain [cells]
      integer, intent(in) :: n      !< length of the profile array
      real, intent(in)    :: vmin   !< minimal value of the profile
      real, intent(in)    :: vmax   !< maximum value of the profile
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
end module initproblem
