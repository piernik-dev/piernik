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
   character(len=cbuff_len) :: mag_field_orient
   real, target, allocatable, dimension(:,:,:,:) :: den0, mtx0, mty0, mtz0, ene0
   integer, parameter       :: dname_len = 10

   namelist /PROBLEM_CONTROL/  alpha, d0, dout, r_max, mag_field_orient, r_in, r_out, f_in, f_out

contains
!-----------------------------------------------------------------------------
   subroutine read_problem_par
      use dataio_pub,          only: ierrh, par_file, namelist_errh, compare_namelist      ! QA_WARN required for diff_nml
      use mpisetup,            only: cbuff, rbuff, buffer_dim, master, slave, comm, ierr
      use mpi,                 only: MPI_CHARACTER, MPI_DOUBLE_PRECISION
      use gravity,             only: grav_pot_3d
      use grid,                only: geometry
      use types,               only: problem_customize_solution
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

      if (master) then

         diff_nml(PROBLEM_CONTROL)


         cbuff(3) = mag_field_orient

         rbuff(1) = d0
         rbuff(2) = dout
         rbuff(3) = r_max
         rbuff(4) = alpha
         rbuff(5) = r_in
         rbuff(6) = r_out
         rbuff(7) = f_in
         rbuff(8) = f_out

      endif

      call MPI_Bcast(cbuff, cbuff_len*buffer_dim, MPI_CHARACTER,        0, comm, ierr)
      call MPI_Bcast(rbuff,           buffer_dim, MPI_DOUBLE_PRECISION, 0, comm, ierr)

      if (slave) then

         mag_field_orient = cbuff(3)

         d0               = rbuff(1)
         dout             = rbuff(2)
         r_max            = rbuff(3)
         alpha            = rbuff(4)
         r_in             = rbuff(5)
         r_out            = rbuff(6)
         f_in             = rbuff(7)
         f_out            = rbuff(8)

      endif

      if (geometry=="cylindrical") then
         problem_write_restart => write_initial_fld_to_restart
         problem_read_restart  => read_initial_fld_from_restart
         problem_customize_solution => problem_customize_solution_kepler
         user_bnd_xl => my_bnd_xl
         user_bnd_xr => my_bnd_xr
         grav_pot_3d => my_grav_pot_3d
      endif
   end subroutine read_problem_par
!-----------------------------------------------------------------------------
   subroutine init_prob
      use types,               only: component_fluid
      use arrays,              only: u, b, dprof
      use constants,           only: newtong
      use fluidindex,          only: ibx, iby, ibz, nvar
      use gravity,             only: r_smooth, r_grav, n_gravr, ptmass, source_terms_grav, grav_pot2accel, grav_pot_3d
      use grid,                only: x, y, z, nx, ny, nz, zdim, has_dir, geometry
      use hydrostatic,         only: hydrostatic_zeq_densmid
      use mpisetup,            only: smalld, smallei
      use types,               only: component_fluid
      use dataio_pub,          only: die

      implicit none

      integer :: i, j, k, kmid, p
      real    :: xi, yj, zk, rc, vx, vy, vz, b0, sqr_gm, vr, vphi
      real    :: csim2, gprim, H2

      real, dimension(:), allocatable :: grav
      type(component_fluid), pointer  :: fl


!   Secondary parameters

      sqr_gm = sqrt(newtong*ptmass)
      do k = 1, nz
         if (z(k) < 0.0) kmid = k       ! the midplane is in between ksmid and ksmid+1
      enddo

      if (associated(nvar%ion) .and. geometry=='cartesian') then
         fl => nvar%ion
         csim2 = fl%cs2*(1.0+alpha)
         b0    = sqrt(2.*alpha*d0*fl%cs2)

         do j = 1,ny
            yj = y(j)
            do i = 1,nx
               xi = x(i)
               rc = sqrt(xi**2+yj**2)

               if (has_dir(zdim)) call hydrostatic_zeq_densmid(i, j, d0, csim2)

               do k = 1,nz

                  vx = sqr_gm * (-yj)/(rc**2+r_smooth**2)**0.75
                  vy = sqr_gm * ( xi)/(rc**2+r_smooth**2)**0.75
                  vz = 0.0

                  u(fl%idn,i,j,k) = min((rc/r_grav)**n_gravr,100.0)

                  if (has_dir(zdim)) then
                     u(fl%idn,i,j,k) = dout + (dprof(k)-dout)/cosh(u(fl%idn,i,j,k))
                  else
                     u(fl%idn,i,j,k) = dout + (d0 - dout)/cosh(u(fl%idn,i,j,k))
                  endif
                  u(fl%idn,i,j,k) = max(u(fl%idn,i,j,k), smalld)
                  u(fl%imx,i,j,k) = vx*u(fl%idn,i,j,k)
                  u(fl%imy,i,j,k) = vy*u(fl%idn,i,j,k)
                  u(fl%imz,i,j,k) = vz*u(fl%idn,i,j,k)
                  if (fl%ien > 0) then
                     u(fl%ien,i,j,k) = fl%cs2/(fl%gam_1)*u(fl%idn,i,j,k)
                     u(fl%ien,i,j,k) = max(u(fl%ien,i,j,k), smallei)
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

         call grav_pot_3d

         if (.not.allocated(den0)) allocate(den0(nvar%fluids,nx,ny,nz))
         if (.not.allocated(mtx0)) allocate(mtx0(nvar%fluids,nx,ny,nz))
         if (.not.allocated(mty0)) allocate(mty0(nvar%fluids,nx,ny,nz))
         if (.not.allocated(mtz0)) allocate(mtz0(nvar%fluids,nx,ny,nz))
         if (.not.allocated(ene0)) allocate(ene0(nvar%fluids,nx,ny,nz))

         if (.not.allocated(grav)) allocate(grav(nx))

         do p = 1, nvar%fluids
            fl => nvar%all_fluids(p)
            call source_terms_grav
            call grav_pot2accel('xsweep',1,1, nx, grav, 1)


            do j = 1,ny
               yj = y(j)
               do i = 1,nx
                  xi = x(i)
                  rc = xi + r_smooth

                  gprim = newtong*ptmass / xi**3
                  if (fl%cs > 0) then
                     H2 = 2.0*fl%cs2/gprim
                  else
                     H2 = 1.0
                  endif

                  do k = 1,nz
                     zk = z(k)
                     u(fl%idn,i,j,k) = max(d0*(1./cosh((xi/r_max)**10)) * exp(-zk**2/H2),smalld)

                     vr   = 0.0
!                     vphi = sqrt( max(abs(grav(i)) * rc - fl%cs2,0.0))
                     vphi = sqrt( max(abs(grav(i)) * rc, 0.0))
                     vz   = 0.0

                     u(fl%imx,i,j,k) = vr   * u(fl%idn,i,j,k)
                     u(fl%imy,i,j,k) = vphi * u(fl%idn,i,j,k)
                     u(fl%imz,i,j,k) = vz   * u(fl%idn,i,j,k)
                     if (fl%ien > 0) then
                        u(fl%ien,i,j,k) = fl%cs2/(fl%gam_1)*u(fl%idn,i,j,k)
                        u(fl%ien,i,j,k) = max(u(fl%ien,i,j,k), smallei)
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
      else
         call die("[initproblem:init_prob] I don't know what to do... :/")
      endif
      return
   end subroutine init_prob
!-----------------------------------------------------------------------------
   subroutine write_initial_fld_to_restart(file_id)

      use hdf5,        only: HID_T
      use grid,        only: nx, ny, nz
      use dataio_hdf5, only: write_3darr_to_restart
      use fluidindex,  only: nvar

      implicit none

      integer(HID_T),intent(in)  :: file_id
      integer :: i
      character(len=dname_len) :: dname

      do i = LBOUND(den0,1), UBOUND(den0,1)
         write(dname,'(2a)') nvar%all_fluids(i)%tag, '_den0'
         if (allocated(den0)) call write_3darr_to_restart(den0(i,:,:,:), file_id, dname, nx, ny, nz)
         write(dname,'(2a)') nvar%all_fluids(i)%tag, '_mtx0'
         if (allocated(mtx0)) call write_3darr_to_restart(mtx0(i,:,:,:), file_id, dname, nx, ny, nz)
         write(dname,'(2a)') nvar%all_fluids(i)%tag, '_mty0'
         if (allocated(mty0)) call write_3darr_to_restart(mty0(i,:,:,:), file_id, dname, nx, ny, nz)
         write(dname,'(2a)') nvar%all_fluids(i)%tag, '_mtz0'
         if (allocated(mtz0)) call write_3darr_to_restart(mtz0(i,:,:,:), file_id, dname, nx, ny, nz)
         write(dname,'(2a)') nvar%all_fluids(i)%tag, '_ene0'
         if (allocated(ene0)) call write_3darr_to_restart(ene0(i,:,:,:), file_id, dname, nx, ny, nz)
      enddo

   end subroutine write_initial_fld_to_restart
!-----------------------------------------------------------------------------
   subroutine read_initial_fld_from_restart(file_id)

      use hdf5,        only: HID_T
      use grid,        only: nx, ny, nz
      use fluidindex,  only: nvar
      use dataio_hdf5, only: read_3darr_from_restart
      use fluidindex,  only: nvar

      implicit none

      integer(HID_T),intent(in) :: file_id

      character(len=dname_len) :: dname
      real, dimension(:,:,:), pointer :: p3d
      integer :: i

      ! /todo First query for existence of den0, vlx0 and vly0, then allocate
      if (.not.allocated(den0)) allocate(den0(nvar%fluids,nx,ny,nz))
      if (.not.allocated(mtx0)) allocate(mtx0(nvar%fluids,nx,ny,nz))
      if (.not.allocated(mty0)) allocate(mty0(nvar%fluids,nx,ny,nz))
      if (.not.allocated(mtz0)) allocate(mtz0(nvar%fluids,nx,ny,nz))
      if (.not.allocated(ene0)) allocate(ene0(nvar%fluids,nx,ny,nz))

      do i=1, nvar%fluids
         write(dname,'(2a)') nvar%all_fluids(i)%tag, '_den0'
         if (.not.associated(p3d)) p3d => den0(i,:,:,:)
         call read_3darr_from_restart(file_id,dname,p3d,nx,ny,nz)
         if (associated(p3d)) nullify(p3d)

         write(dname,'(2a)') nvar%all_fluids(i)%tag, '_mtx0'
         if (.not.associated(p3d)) p3d => mtx0(i,:,:,:)
         call read_3darr_from_restart(file_id,dname,p3d,nx,ny,nz)
         if (associated(p3d)) nullify(p3d)

         write(dname,'(2a)') nvar%all_fluids(i)%tag, '_mty0'
         if (.not.associated(p3d)) p3d => mty0(i,:,:,:)
         call read_3darr_from_restart(file_id,dname,p3d,nx,ny,nz)
         if (associated(p3d)) nullify(p3d)

         write(dname,'(2a)') nvar%all_fluids(i)%tag, '_mtz0'
         if (.not.associated(p3d)) p3d => mtz0(i,:,:,:)
         call read_3darr_from_restart(file_id,dname,p3d,nx,ny,nz)
         if (associated(p3d)) nullify(p3d)

         write(dname,'(2a)') nvar%all_fluids(i)%tag, '_ene0'
         if (.not.associated(p3d)) p3d => ene0(i,:,:,:)
         call read_3darr_from_restart(file_id,dname,p3d,nx,ny,nz)
         if (associated(p3d)) nullify(p3d)
      enddo

   end subroutine read_initial_fld_from_restart
!-----------------------------------------------------------------------------
   subroutine problem_customize_solution_kepler
      use mpisetup,        only: dt
      use arrays,          only: u
      use grid,            only: x, nx, ny, nz
      use fluidboundaries, only: all_fluid_boundaries
      use fluidindex,      only: iarr_all_dn, iarr_all_mx, iarr_all_my, iarr_all_mz
#ifndef ISO
      use fluidindex,      only: iarr_all_en
#endif /* ISO */
      implicit none
      integer                               :: i, j, k
      logical, save                         :: frun = .true.
      real, dimension(:,:), allocatable, save :: funcR

      if (frun) then
         allocate(funcR(size(iarr_all_dn),nx) )

         funcR(1,:) = -tanh((x(:)-r_in+1.0)**f_in) + 1.0
         funcR(1,:) = alpha*funcR(1,:)
         open(212,file="funcR.dat",status="unknown")
         do i = 1, nx
            write(212,*) x(i),funcR(1,i)
         enddo
         close(212)
         frun = .false.
         funcR(:,:) = spread(funcR(1,:),1,size(iarr_all_dn))
      endif

      do j = 1, ny
         do k = 1, nz
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
      use grid,      only: x, z, nx, nz
      implicit none
      logical, save :: frun = .true.
      real          :: r2
      integer       :: i, k

      if (frun) then
         do i = 1, nx
            do k = 1, nz
               r2 = x(i)**2 + z(k)**2
               gp(i,:,k) = -newtong*ptmass / sqrt(r2)
            enddo
         enddo
      endif

      frun = .false.
      call sum_potential

   end subroutine my_grav_pot_3d
!-----------------------------------------------------------------------------
   subroutine my_bnd_xl
      use grid,         only: nb,nx,ny,nz,x
      use arrays,       only: u
      use gravity,      only: grav_pot2accel
      use fluidindex,   only: iarr_all_dn, iarr_all_mx, iarr_all_my, iarr_all_mz, nvar
#ifndef ISO
      use fluidindex,   only: iarr_all_en
#endif /* ISO */
      implicit none
      integer :: i, p
      real, dimension(nx) :: grav
      real, dimension(size(iarr_all_my),ny,nz) :: vy,vym
      real, dimension(size(nvar%all_fluids))    :: cs2_arr
      integer, dimension(size(nvar%all_fluids)) :: ind_cs2

      do i = 1, size(nvar%all_fluids)
         ind_cs2    = i
         cs2_arr(i) = nvar%all_fluids(i)%cs2
      enddo

      call grav_pot2accel('xsweep',1,1, nx, grav, 1)

      do i = 1,nb
         u(iarr_all_dn,i,:,:) = u(iarr_all_dn,nb+1,:,:)
         u(iarr_all_mx,i,:,:) = max(0.0,u(iarr_all_mx,nb+1,:,:))
         do p = 1, size(nvar%all_fluids)
            u(iarr_all_my(p),i,:,:) = sqrt( abs(grav(i)) * x(i) - cs2_arr(p)) *  u(iarr_all_dn(p),i,:,:)
         enddo
         u(iarr_all_mz,i,:,:) = u(iarr_all_mz,nb+1,:,:)
#ifndef ISO
         u(iarr_all_en,i,:,:) = u(iarr_all_en,nb+1,:,:)
#endif
      enddo

      do i = nb,1,-1
         vym(:,:,:) = u(iarr_all_my,i+2,:,:)/u(iarr_all_dn,i+1,:,:)
         vy(:,:,:)  = u(iarr_all_my,i+1,:,:)/u(iarr_all_dn,i+1,:,:)
!         u(iarr_all_my,i,:,:) = (vym(:,:,:) + (x(i) - x(i+2)) / (x(i+1) - x(i+2)) * (vy - vym))*u(iarr_all_dn,i,:,:)
      enddo

   end subroutine my_bnd_xl
!-----------------------------------------------------------------------------
   subroutine my_bnd_xr
      use grid,   only: nb, nxb
      use arrays, only: u
      use fluidindex,  only: iarr_all_dn, iarr_all_mx, iarr_all_my, iarr_all_mz
#ifndef ISO
      use fluidindex,  only: iarr_all_en
#endif /* ISO */
      implicit none

      u(iarr_all_dn,nxb+nb+1:nxb+2*nb,:,:) = den0(:,nxb+nb+1:nxb+2*nb,:,:)
      u(iarr_all_mx,nxb+nb+1:nxb+2*nb,:,:) = mtx0(:,nxb+nb+1:nxb+2*nb,:,:)
      u(iarr_all_my,nxb+nb+1:nxb+2*nb,:,:) = mty0(:,nxb+nb+1:nxb+2*nb,:,:)
      u(iarr_all_mz,nxb+nb+1:nxb+2*nb,:,:) = mtz0(:,nxb+nb+1:nxb+2*nb,:,:)
#ifndef ISO
      u(iarr_all_en,nxb+nb+1:nxb+2*nb,:,:) = ene0(:,nxb+nb+1:nxb+2*nb,:,:)
#endif
   end subroutine my_bnd_xr
!-----------------------------------------------------------------------------
end module initproblem
