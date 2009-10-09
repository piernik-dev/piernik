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

!>
!! \brief Module containing all main subroutines and functions that govern gravity force in the code
!<
module gravity

   use constants

   character(LEN=9) :: gp_status
   real    :: g_z, g_y, dg_dz, r_gc
   real    :: ptmass, ptm_x, ptm_y, ptm_z, r_smooth
   integer :: nsub
   real    :: h_grav,  r_grav
   integer :: n_gravh, n_gravr, n_gravr2
   real    :: tune_zeq, tune_zeq_bnd


   contains

!>
!! \brief Routine that sets the initial values of gravity parameters from namelist @c GRAVITY
!<
   subroutine init_grav
      use errh, only : namelist_errh
      use mpisetup
      implicit none
      integer :: ierrh
      character(LEN=100) :: par_file, tmp_log_file

      namelist /GRAVITY/ g_z,g_y, dg_dz, r_gc,  &
                     ptmass,ptm_x,ptm_y,ptm_z,r_smooth, &
                     nsub, tune_zeq, tune_zeq_bnd,      &
                     h_grav, r_grav, n_gravr, n_gravr2, n_gravh

      par_file = trim(cwd)//'/problem.par'
      tmp_log_file = trim(cwd)//'/tmp.log'

      g_z     = 0.0
      g_y     = 0.0
      dg_dz   = 0.0
      r_gc    = 8500
      ptmass  = 0.0
      ptm_x   = 0.0
      ptm_y   = 0.0
      ptm_z   = 0.0
      r_smooth= 0.0
      nsub    = 10
      tune_zeq     = 1.0
      tune_zeq_bnd = 1.0
      h_grav = 1.e6
      r_grav = 1.e6
      n_gravr = 0
      n_gravr2= 0
      n_gravh = 0

      if(proc .eq. 0) then
         open(1,file=par_file)
            read(unit=1,nml=GRAVITY,iostat=ierrh)
            call namelist_errh(ierrh,'GRAVITY')
         close(1)
         open(3, file=tmp_log_file, position='append')
            write(unit=3,nml=GRAVITY)
         close(3)

         ibuff(1)  = nsub
         ibuff(2)  = n_gravr
         ibuff(3)  = n_gravr2
         ibuff(4)  = n_gravh

         rbuff(1)  = g_z
         rbuff(2)  = g_y
         rbuff(3)  = dg_dz
         rbuff(4)  = r_gc
         rbuff(5)  = ptmass
         rbuff(6)  = ptm_x
         rbuff(7)  = ptm_y
         rbuff(8)  = ptm_z
         rbuff(9)  = tune_zeq
         rbuff(10) = tune_zeq_bnd
         rbuff(11) = r_smooth
         rbuff(12) = h_grav
         rbuff(13) = r_grav

         call MPI_BCAST(ibuff,    buffer_dim, MPI_INTEGER,          0, comm, ierr)
         call MPI_BCAST(rbuff,    buffer_dim, MPI_DOUBLE_PRECISION, 0, comm, ierr)

     else

         call MPI_BCAST(ibuff,    buffer_dim, MPI_INTEGER,          0, comm, ierr)
         call MPI_BCAST(rbuff,    buffer_dim, MPI_DOUBLE_PRECISION, 0, comm, ierr)

         nsub                = ibuff(1)
         n_gravr             = ibuff(2)
         n_gravr2            = ibuff(3)
         n_gravh             = ibuff(4)

         g_z                 = rbuff(1)
         g_y                 = rbuff(2)
         dg_dz               = rbuff(3)
         r_gc                = rbuff(4)
         ptmass              = rbuff(5)
         ptm_x               = rbuff(6)
         ptm_y               = rbuff(7)
         ptm_z               = rbuff(8)
         tune_zeq            = rbuff(9)
         tune_zeq_bnd        = rbuff(10)
         r_smooth            = rbuff(11)
         h_grav              = rbuff(12)
         r_grav              = rbuff(13)

      endif


   end subroutine init_grav

!--------------------------------------------------------------------------
!>
!! \brief Routine that compute values of gravitational potential filling in gp array and setting gp_status character string
!<
   subroutine grav_pot_3d
      use arrays, only     : gp
      use grid, only       : nx,ny,nz,x,y,z
      use mpisetup, only   : smalld
      use initfluids, only : cs_iso2
#ifdef GRAV_USER
      use gravity_user, only : grav_pot_user
#endif /* GRAV_USER */

      implicit none
      integer :: i, j, k
      real    :: r2, rc, fr

      gp_status = ''

#ifdef GRAV_NULL
      gp(:,:,:) = 0.0      ! do nothing

#elif defined (GRAV_UNIFORM)
      do i = 1, nz
         gp(:,:,i) = -g_z*z(i)
      enddo

#elif defined (GRAV_LINEAR)
      do i = 1, nz
         gp(:,:,i) = -0.5 * dg_dz * z(i)**2
      enddo

#elif defined (GRAV_PTMASS)
       do i = 1, nx
          do j = 1, ny
             do k = 1, nz
               rc = dsqrt(x(i)**2+y(j)**2)
               r2 = x(i)**2 + y(j)**2 + z(k)**2
               fr = min( (rc/r_grav)**n_gravr , 100.0)
               fr = max(1./cosh(fr),smalld/100.)
               gp(i,j,k) = -newtong*ptmass / dsqrt(r2 + r_smooth**2)
               gp(i,j,k) = gp(i,j,k) - cs_iso2 * dlog(fr) ! *d0
             enddo
          enddo
       enddo

#elif defined (GRAV_PTMASSPURE)
       do i = 1, nx
          do j = 1, ny
             do k = 1, nz
               rc = dsqrt(x(i)**2+y(j)**2)
               r2 = (x(i) - ptm_x)**2 + (y(j) - ptm_y)**2 + (z(k) - ptm_z)**2
               gp(i,j,k) = -newtong*ptmass / dsqrt(r2 + r_smooth**2)
             enddo
          enddo
       enddo

#elif defined (GRAV_PTFLAT)
       do i = 1, nx
          do j = 1, ny
             rc = dsqrt(x(i)**2+y(j)**2)
             fr = min( (rc/r_grav)**n_gravr , 100.0)
             fr = max(1./cosh(fr),smalld/100.)
             gp(i,j,:) = -newtong*ptmass / dsqrt(rc**2 + r_smooth**2)
             gp(i,j,:) = gp(i,j,:) - cs_iso2 * dlog(fr) ! *d0
          enddo
       enddo

#elif defined (GRAV_USER)
       call grav_pot_user()

#else /* GRAV_(SPECIFIED) */
       gp_status = 'undefined'
!#warning 'GRAV declared, but gravity model undefined in grav_pot'
! niektore modele grawitacji realizowane sa za pomoca przyspieszenia
! (np. 'galactic') z ktorego liczony jest potencjal
#endif /* GRAV_(SPECIFIED) */
!-----------------------

      if(gp_status .eq. 'undefined') then
         call grav_accel2pot
      endif

   end subroutine grav_pot_3d

!--------------------------------------------------------------------------
!>
!! \brief Routine that compute values of gravitational acceleration
!! \param sweep string of characters that points out the current sweep direction
!! |param i1 integer, number of column in the first direction after one pointed out by sweep
!! \param i2 integer, number of column in the second direction after one pointed out by sweep
!! \param xsw 1D position array in the direction pointed out by sweep
!! \param n number of elements of xsw array
!! \return grav 1D array of gravitational acceleration values computed for positions from xsw
!<

   subroutine grav_accel(sweep, i1,i2, xsw, n, grav)
      use grid, only : nb
      use arrays, only : gp
      use grid, only :   x,y,z,dl,xdim,ydim,zdim,nx,ny,nz
      use grid, only :   is,ie,js,je,ks,ke, xr,yr,zr
#if defined GRAV_PTMASS || defined GRAV_PTFLAT
      use mpisetup, only : smalld
#endif /* GRAV_PTMASS || GRAV_PTFLAT */
#if defined GRAV_ACC_USER
      use gravity_user, only : grav_accel_user
#endif /* GRAV_ACC_USER  */
      implicit none
      character, intent(in) :: sweep*6
      integer, intent(in)   :: i1, i2
      integer, intent(in)   :: n
      real, dimension(:)    :: xsw
      real, dimension(n),intent(out) :: grav
      real                  :: x1, x2
#if defined GRAV_PTMASS || defined GRAV_PTFLAT
      real, dimension(n)    :: r, r32
#endif /* GRAV_PTMASS || GRAV_PTFLAT */

      select case (sweep)
         case('xsweep')
            x1  = y(i1)
            x2  = z(i2)
         case('ysweep')
            x1  = z(i1)
            x2  = x(i2)
         case('zsweep')
            x1  = x(i1)
            x2  = y(i2)
      end select

! simplified, z component only of Galactic gravitational acceleration from Ferriere'98
      if (sweep == 'zsweep') then
         grav = 3.23e8 * (  &
           (-4.4e-9 * exp(-(r_gc-r_gc_sun)/(4.9*kpc)) * xsw/sqrt(xsw**2+(0.2*kpc)**2)) &
           -( 1.7e-9 * (r_gc_sun**2 + (2.2*kpc)**2)/(r_gc**2 + (2.2*kpc)**2)*xsw/kpc) )
!          -Om*(Om+G) * Z * (kpc ?) ! in the transition region between rigid
!                                   ! and flat rotation F'98: eq.(36)

      else
         grav=0.0
      endif

      if (n_gravh .ne. 0) then
         grav(:)=grav(:)/cosh((xsw(:)/h_grav)**n_gravh )
      endif

   end subroutine grav_accel

!>
!! \brief Routine that compute values of gravitational acceleration using gravitational potential array gp
!! \param sweep string of characters that points out the current sweep direction
!! |param i1 integer, number of column in the first direction after one pointed out by sweep
!! \param i2 integer, number of column in the second direction after one pointed out by sweep
!! \param n number of elements of returned array grav
!! \return grav 1D array of gravitational acceleration values computed for positions from xsw
!<
   subroutine grav_pot2accel(sweep, i1,i2, n, grav)
      use arrays, only : gp
      use grid, only : dl, xdim, ydim, zdim
      implicit none
      character, intent(in)          :: sweep*6
      integer, intent(in)            :: i1, i2
      integer, intent(in)            :: n
      real, dimension(n),intent(out) :: grav

! Gravitational accelleration is computed on right cell boundaries

      select case(sweep)
         case('xsweep')
            grav(2:n-1) = 0.5*(gp(1:n-2,i1,i2) - gp(3:n,i1,i2))/dl(xdim)
         case('ysweep')
            grav(2:n-1) = 0.5*(gp(i2,1:n-2,i1) - gp(i2,3:n,i1))/dl(ydim)
         case('zsweep')
            grav(2:n-1) = 0.5*(gp(i1,i2,1:n-2) - gp(i1,i2,3:n))/dl(zdim)
      end select
      grav(1) = grav(2); grav(n) = grav(n-1)

   end subroutine grav_pot2accel

!--------------------------------------------------------------------------
!>
!! \brief Routine that uses gravity acceleration given in grav_accel to compute values of gravitational potential filling in gp array
!<
   subroutine grav_accel2pot
      use mpisetup
      use arrays, only : gp
      use grid, only  : dl,xdim,ydim,zdim,is,ie,js,je,ks,ke
      use grid, only  : nb,nx,ny,nz,zr,yr,xr

      implicit none
      integer               :: i, j, k, ip, pgpmax
      real, allocatable     :: gpwork(:,:,:)
      real gravrx(nx), gravry(ny), gravrz(nz)
      real                  :: gp_max
      integer, dimension(3) :: loc_gp_max
      integer               :: proc_gp_max
      integer               :: px, py, pz, pc(3)
      real dgpx_proc, dgpx_all(0:nproc-1), &
           dgpy_proc, dgpy_all(0:nproc-1), &
           dgpz_proc, dgpz_all(0:nproc-1), &
           dgpx(0:pxsize-1,0:pysize-1,0:pzsize-1), &
           dgpy(0:pxsize-1,0:pysize-1,0:pzsize-1), &
           dgpz(0:pxsize-1,0:pysize-1,0:pzsize-1), &
           ddgp(0:pxsize-1,0:pysize-1,0:pzsize-1)


! Gravitational potential gp(i,j,k) is defined in cell centers
! First we find gp(:,:,:) in each block separately assuming:
! Instead of gp, gpwork is used to not change already computed gp when
! disk, bulge and halo parts of galactic potential are computed.
! gravpart is used to distinguish which part has to be computed.

      allocate(gpwork(nx,ny,nz))
      gpwork(1,1,1) = 0.0

      call grav_accel('xsweep',1,1, xr(:), nx, gravrx)
      do i = 1, nx-1
         gpwork(i+1,1,1) = gpwork(i,1,1) - gravrx(i)*dl(xdim)
      enddo

      do i=1,nx
         call grav_accel('ysweep',1, i, yr(:), ny, gravry)
         do j = 1, ny-1
            gpwork(i,j+1,1) = gpwork(i,j,1) - gravry(j)*dl(ydim)
         enddo
      enddo

      do i=1,nx
         do j=1,ny
            call grav_accel('zsweep', i, j, zr(:), nz, gravrz)
            do k = 1, nz-1
               gpwork(i,j,k+1) = gpwork(i,j,k) - gravrz(k)*dl(zdim)
            enddo
         enddo
      enddo

      dgpx_proc = gpwork(is,1,1)-gpwork(1,1,1)
      dgpy_proc = gpwork(1,js,1)-gpwork(1,1,1)
      dgpz_proc = gpwork(1,1,ks)-gpwork(1,1,1)

      call MPI_GATHER ( dgpx_proc, 1, MPI_DOUBLE_PRECISION, &
                        dgpx_all,  1, MPI_DOUBLE_PRECISION, &
                        0, comm3d,err )
      call MPI_GATHER ( dgpy_proc, 1, MPI_DOUBLE_PRECISION, &
                        dgpy_all,  1, MPI_DOUBLE_PRECISION, &
                       0, comm3d,err )
      call MPI_GATHER ( dgpz_proc, 1, MPI_DOUBLE_PRECISION, &
                        dgpz_all,  1, MPI_DOUBLE_PRECISION, &
                        0, comm3d,err )


      if(proc .eq. 0) then

         do ip = 0, nproc-1
            call MPI_CART_COORDS(comm3d, ip, ndims, pc, ierr)

            px = pc(1)
            py = pc(2)
            pz = pc(3)

            dgpx(px,py,pz) = dgpx_all(ip)
            dgpy(px,py,pz) = dgpy_all(ip)
            dgpz(px,py,pz) = dgpz_all(ip)

         enddo

         ddgp(0,0,0) = 0.0
         do i = 1, pxsize-1
            ddgp(i,0,0) = ddgp(i-1,0,0) + dgpx(i-1,0,0)
         enddo

         do i=0, pxsize-1
            do j = 1, pysize-1
               ddgp(i,j,0) = ddgp(i,j-1,0) + dgpy(i,j-1,0)
            enddo
         enddo

         do i=0, pxsize-1
            do j=0, pysize-1
               do k = 1, pzsize-1
                  ddgp(i,j,k)= ddgp(i,j,k-1) + dgpz(i,j,k-1)
               enddo
            enddo
         enddo

      endif

      call MPI_BCAST(ddgp, nproc, MPI_DOUBLE_PRECISION, 0, comm, ierr)

      px = pcoords(1)
      py = pcoords(2)
      pz = pcoords(3)

      gpwork = gpwork + ddgp(px,py,pz)

      gp_max      = maxval(gpwork(is:ie,js:je,ks:ke))
      loc_gp_max  = maxloc(gpwork(is:ie,js:je,ks:ke)) &
                  + (/nb,nb,nb/)

      call mpifind(gp_max, 'max', loc_gp_max, proc_gp_max)
      pgpmax = proc_gp_max

      call MPI_BCAST(gp_max, 1, MPI_DOUBLE_PRECISION, pgpmax, comm, ierr)
      gpwork = gpwork - gp_max

      gp=gpwork
      deallocate(gpwork)

   end subroutine grav_accel2pot


end module gravity
