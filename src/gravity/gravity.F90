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

   subroutine init_grav
      use errh, only : namelist_errh
      use mpisetup
      implicit none
      integer :: ierrh
      character(LEN=100) :: par_file, tmp_log_file

      namelist /GRAVITY/ g_z,g_y, dg_dz, r_gc, &
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
   subroutine grav_pot_3d
      use arrays, only : gp
      use grid, only : nx,ny,nz,z
      implicit none
      integer i, j
      real, allocatable :: gpot(:)
      allocate(gpot(nz))

      do j = 1,ny
         do i = 1,nx
            call grav_pot('zsweep', i,j, z(:), nz, gpot, gp_status)
            if(gp_status .eq. 'undefined') exit
            gp(i,j,:) = gpot
         enddo
      enddo

      if(gp_status .eq. 'undefined') then
!         gravpart = 'default'
         call grav_accel2pot
      endif

      deallocate(gpot)

   end subroutine grav_pot_3d

!--------------------------------------------------------------------------
   subroutine grav_pot(sweep, i1,i2, xsw, n, gpot,status,temp_log)
#if defined GRAV_PTMASS || defined GRAV_PTFLAT || defined GRAV_PTMASSPURE
      use start,  only : csim2,smalld
      use grid, only : x,y,z
#endif /* GRAV_PTMASS || GRAV_PTFLAT || GRAV_PTMASSPURE */
#ifdef GRAV_GALACTIC
      use grid, only : z
#endif /* GRAV_GALACTIC */
#ifdef GRAV_USER
      use gravity_user, only : grav_pot_user
#endif /* GRAV_USER */

      implicit none
      character, intent(in) :: sweep*6
      integer, intent(in)   :: i1, i2, n
      logical, optional     :: temp_log
      real, dimension(n)    :: xsw
      real, dimension(n),intent(out) :: gpot
      character, intent(inout) :: status*9
#if defined GRAV_PTMASS || defined GRAV_PTFLAT || defined GRAV_PTMASSPURE
      real                  :: x1, x2
      real, dimension(n)    :: x3, rc, fr
#endif /* GRAV_PTMASS || GRAV_PTFLAT || GRAV_PTMASSPURE */
#ifdef GRAV_GALACTIC
      real                  :: x1
#endif /* GRAV_GALACTIC */
      status = ''

#ifdef GRAV_NULL
      gpot = 0.0
      ! do nothing
#elif defined (GRAV_UNIFORM)
      select case (sweep)
         case('xsweep')
            gpot = -g_z*z(i2)
         case('ysweep')
            gpot = -g_z*z(i1)
         case('zsweep')
            gpot = -g_z*xsw
      end select
#elif defined (GRAV_LINEAR)
      select case (sweep)
         case('xsweep')
            gpot = -0.5*dg_dz*z(i2)**2
         case('ysweep')
            gpot = -0.5*dg_dz*z(i1)**2
         case('zsweep')
            gpot = -0.5*dg_dz*xsw**2
      end select

#elif defined (GRAV_PTMASS)
      select case (sweep)
         case('xsweep')
            x1  = y(i1)-ptm_y
            x2  = z(i2)-ptm_z
            x3  = xsw - ptm_x
         case('ysweep')
            x1  = z(i1)-ptm_z
            x2  = x(i2)-ptm_x
            x3  = xsw - ptm_y
         case('zsweep')
            x1  = x(i1)-ptm_x
            x2  = y(i2)-ptm_y
            x3  = xsw - ptm_z
         end select
         rc = dsqrt(x1**2+x2**2)
         fr = min( (rc/r_grav)**n_gravr , 100.0)
         fr = max(1./cosh(fr),smalld/100.)
         gpot = -newtong*ptmass/dsqrt(x1**2+x2**2+x3**2+r_smooth**2)
         gpot = gpot - csim2*dlog(fr) ! *d0
#elif defined (GRAV_PTMASSPURE)
         select case (sweep)
            case('xsweep')
               x1  = y(i1)-ptm_y
               x2  = z(i2)-ptm_z
               x3  = xsw - ptm_x
            case('ysweep')
               x1  = z(i1)-ptm_z
               x2  = x(i2)-ptm_x
               x3  = xsw - ptm_y
            case('zsweep')
               x1  = x(i1)-ptm_x
               x2  = y(i2)-ptm_y
               x3  = xsw - ptm_z
         end select
         rc = dsqrt(x1**2+x2**2)
         gpot = -newtong*ptmass/dsqrt(x1**2+x2**2+x3**2+r_smooth**2)
#elif defined (GRAV_PTFLAT)
         select case (sweep)
            case('xsweep')
               x1  = y(i1)-ptm_y
               x2  = 0.0
               x3  = xsw - ptm_x
            case('ysweep')
               x1  = 0.0
               x2  = x(i2)-ptm_x
               x3  = xsw - ptm_y
            case('zsweep')
               x1  = x(i1)-ptm_x
               x2  = y(i2)-ptm_y
               x3  = 0.0
         end select
         rc = dsqrt(x1**2+x2**2)
         fr(:) = min((rc(:)/r_grav)**n_gravr,100.0)
         fr = 1./cosh(fr)+smalld
         gpot = -newtong*ptmass/dsqrt(x1**2+x2**2+x3**2+r_smooth**2)
         gpot = gpot - csim2*dlog(fr) ! *d0
#elif defined (GRAV_USER)
         call grav_pot_user(gpot,sweep,i1,i2,xsw,n,status,temp_log)

#else /* GRAV_(SPECIFIED) */
            status = 'undefined'
!#warning 'GRAV declared, but gravity model undefined in grav_pot'
! niektore modele grawitacji realizowane sa za pomoca przyspieszenia
! (np. 'galactic') z ktorego liczony jest potencjal
#endif /* GRAV_(SPECIFIED) */

   end subroutine grav_pot

!--------------------------------------------------------------------------

   subroutine grav_accel(sweep, i1,i2, xsw, n, grav)
      use grid, only : nb
      use arrays, only : gp
      use grid, only :   x,y,z,dl,xdim,ydim,zdim,nx,ny,nz
      use grid, only :   is,ie,js,je,ks,ke, xr,yr,zr
#if defined GRAV_PTMASS || defined GRAV_PTFLAT
      use start, only : csim2,smalld
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


#ifdef GRAV_UNIFORM
!      case ('uniform')
      if (sweep == 'zsweep') then
         grav = g_z
      else
         grav = 0.0
      endif
#elif defined (GRAV_LINEAR)
      if (sweep == 'zsweep') then
         grav = dg_dz*xsw
      else
         grav = 0.0
      endif
#elif defined (GRAV_PTMASS)
      if (sweep == 'xsweep') then
         r  = sqrt((xsw-ptm_x)*(xsw-ptm_x) + (x1-ptm_y)*(x1-ptm_y)  &
               + (x2-ptm_z)*(x2-ptm_z) + r_smooth*r_smooth)
      elseif(sweep == 'ysweep') then
         r  = sqrt((xsw-ptm_y)*(xsw-ptm_y) + (x1-ptm_z)*(x1-ptm_z)  &
               + (x2-ptm_x)*(x2-ptm_x) + r_smooth*r_smooth)
      else
         r  = sqrt((xsw-ptm_z)*(xsw-ptm_z) + (x1-ptm_x)*(x1-ptm_x)  &
               + (x2-ptm_y)*(x2-ptm_y) + r_smooth*r_smooth)
      endif
      r32  = r*r*r
      grav = -newtong*ptmass*xsw/r32
#elif defined (GRAV_PTFLAT)
      if (sweep == 'xsweep') then
         r  = sqrt((xsw-ptm_x)*(xsw-ptm_x) + (x1-ptm_y)*(x1-ptm_y)  &
               + r_smooth*r_smooth)
         r32  = r*r*r
         grav = -newtong*ptmass*xsw/r32
      elseif(sweep == 'ysweep') then
         r  = sqrt((xsw-ptm_y)*(xsw-ptm_y) + (x2-ptm_x)*(x2-ptm_x)  &
               + r_smooth*r_smooth)
         r32  = r*r*r
         grav = -newtong*ptmass*xsw/r32
      else
         grav = 0.0
      endif
#elif defined (GRAV_GALACTIC)
! simplified, z component only of Galactic gravitational acceleration from Ferriere'98
      if (sweep == 'zsweep') then
         grav = cmps2 * (  &
           (-4.4e-9 * exp(-(r_gc-r_gc_sun)/(4.9*kpc)) * xsw/sqrt(xsw**2+(0.2*kpc)**2)) &
           -( 1.7e-9 * (r_gc_sun**2 + (2.2*kpc)**2)/(r_gc**2 + (2.2*kpc)**2)*xsw/kpc) )
!          -Om*(Om+G) * Z * (kpc ?) ! in the transition region between rigid
!                                   ! and flat rotation F'98: eq.(36)

!           write(*,*) cmps2, kpc, r_gc, r_gc_sun


!          write(*,*) grav
!	  stop
      else
         grav=0.0
      endif

#elif defined (GRAV_ACC_USER)
      call grav_accel_user(grav,sweep,i1,i2,xsw,n)

#elif defined (GRAV_NULL)
      grav=0.0
#else /* GRAV_(SPECIFIED) */
!#error 'GRAV declared, but gravity model undefined'
#endif /* GRAV_(SPECIFIED) */

      if (n_gravh .ne. 0) then
         grav(:)=grav(:)/cosh((xsw(:)/h_grav)**n_gravh )
      endif

   end subroutine grav_accel

!--------------------------------------------------------------------------

!  subroutine grav_profile
!
!    implicit none
!
!        if(n_gravh .ne. 0) then
!          gravity_profile = 1./cosh((z/(h_grav))**n_gravh)
!        else
!         gravity_profile(:) = 1.
!       endif
!
!       write(*,*) gravity_profile
!
!  end subroutine  grav_profile
!
!*****************************************************************************


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
            grav(1:n-1) = (gp(1:n-1,i1,i2) - gp(2:n,i1,i2))/dl(xdim)
         case('ysweep')
            grav(1:n-1) = (gp(i2,1:n-1,i1) - gp(i2,2:n,i1))/dl(ydim)
         case('zsweep')
            grav(1:n-1) = (gp(i1,i2,1:n-1) - gp(i1,i2,2:n))/dl(zdim)
      end select

   end subroutine grav_pot2accel

!--------------------------------------------------------------------------

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
