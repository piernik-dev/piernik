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

!>
!! \brief [DW] Module containing a subroutine that arranges %hydrostatic equilibrium in the vertical (z) direction
!<
module hydrostatic
! pulled by GRAV
   implicit none

   private
#ifdef GRAV
   public :: hydrostatic_zeq_coldens, hydrostatic_zeq_densmid, gprofs, nstot, zs, dzs, outh_bnd
#endif /* GRAV */

   real, allocatable, dimension(:), save :: zs, gprofs
   real,    save :: dzs
   integer, save :: nstot
#ifndef NEW_HYDROSTATIC
   real,    save :: dmid
#endif /* !NEW_HYDROSTATIC */

contains
#ifdef GRAV
!>
!! \brief Routine that arranges %hydrostatic equilibrium in the vertical (z) direction
!<
   subroutine hydrostatic_main

      use arrays,     only: dprof
      use dataio_pub, only: die
      use gravity,    only: nsub !, gp_status
      use grid,       only: nz, zl, zr

      implicit none

      real, allocatable, dimension(:) :: dprofs
      integer :: ksub, ksmid, k
#ifndef NEW_HYDROSTATIC
      real    :: factor
#endif /* !NEW_HYDROSTATIC */

      allocate(dprofs(nstot))

      ksmid = 0
      do ksub=1, nstot
         if (zs(ksub) .lt. 0.0) ksmid = ksub      ! the midplane is in between
      enddo                                  ! ksmid and ksmid+1
      if (ksmid == 0) call die("[hydrostatic:hydrostatic_main] ksmid not set")

#ifdef NEW_HYDROSTATIC
      if (ksmid .lt. nstot) then
         dprofs(ksmid+1) = 1.0
         do ksub=ksmid+1, nstot-1
            dprofs(ksub+1) = dprofs(ksub)*(1.0 + 0.5*(gprofs(ksub)+gprofs(ksub+1))*dzs)
         enddo
      endif

      if (ksmid .gt. 1) then
         dprofs(ksmid) = 1.0
         do ksub=ksmid, 2, -1
            dprofs(ksub-1) = dprofs(ksub)*(1.0 - 0.5*(gprofs(ksub)+gprofs(ksub-1))*dzs)
         enddo
      endif

      dprof(:) =0.0
      do k=1,nz
         do ksub=1, nstot
            if (zs(ksub) .gt. zl(k) .and. zs(ksub) .lt. zr(k)) then
               dprof(k) = dprof(k) + dprofs(ksub)/real(nsub)
            endif
         enddo
      enddo
#else /* !NEW_HYDROSTATIC */
      if (ksmid .lt. nstot) then
         dprofs(ksmid+1) = dmid
         do ksub=ksmid+1, nstot-1
            factor = (2.0 + dzs*gprofs(ksub))/(2.0 - dzs*gprofs(ksub))
            dprofs(ksub+1) = factor * dprofs(ksub)
         enddo
      endif

      if (ksmid .gt. 1) then
         dprofs(ksmid) = dmid
         do ksub=ksmid, 2, -1
            factor = (2.0 - dzs*gprofs(ksub))/(2.0 + dzs*gprofs(ksub))
            dprofs(ksub-1) = factor * dprofs(ksub)
         enddo
      endif

      dprof(:) =0.0
      do k=1,nz
         do ksub=1, nstot
            if (zs(ksub) .gt. zl(k) .and. zs(ksub) .lt. zr(k)) then
               dprof(k) = dprof(k) + dprofs(ksub)/real(nsub)
            endif
         enddo
      enddo
#endif /* !NEW_HYDROSTATIC */

      if (allocated(dprofs)) deallocate(dprofs)

      return
   end subroutine hydrostatic_main

   subroutine get_gprofs_accel(iia,jja)
      use gravity, only: tune_zeq, grav_accel
      use grid,    only: nx, ny
      implicit none
      integer, intent(in) :: iia, jja
      integer :: ia, ja


      ia = min(nx,max(1, iia))
      ja = min(ny,max(1, jja))
      call grav_accel('zsweep',ia, ja, zs, nstot, gprofs)
      gprofs = tune_zeq*gprofs

   end subroutine get_gprofs_accel

   subroutine get_gprofs_gparray(iia,jja)
      use arrays,  only: gp
      use gravity, only: tune_zeq
      use grid,    only: nz, z
      implicit none
      integer, intent(in) :: iia, jja
      integer :: ksub, k
      real, allocatable, dimension(:)  :: gpots

      allocate(gpots(nstot))
      k = 1; gpots(:) = 0.0
      do ksub=1, nstot
         if (zs(ksub) >= z(min(k+1,nz))) k = k + 1
         if (zs(ksub) >= z(min(k,nz-1)) .and. zs(ksub) < z(min(k+1,nz))) then
            gpots(ksub) = gp(iia,jja,k) + (zs(ksub) - z(k)) * &
             (gp(iia,jja,min(k+1,nz)) - gp(iia,jja,k)) / (z(min(k+1,nz)) - z(min(k,nz-1)))
         endif
      enddo
!         call grav_pot('zsweep', ia,ja, zs, nstot, gpots,gp_status,.true.)
      gprofs(1:nstot-1) = (gpots(1:nstot-1) - gpots(2:nstot))/dzs
      gprofs = tune_zeq*gprofs
      if (allocated(gpots)) deallocate(gpots)

   end subroutine get_gprofs_gparray

   subroutine hydrostatic_zeq_coldens(iia,jja,coldens,csim2)
      use arrays,  only: dprof
      implicit none
      integer, intent(in) :: iia, jja
      real,    intent(in) :: coldens, csim2
      real :: sdprof

      sdprof = 1.0
      call hydrostatic_zeq_densmid(iia,jja,sdprof,csim2)
      sdprof = sum(dprof)
      dprof = dprof * coldens / sdprof

   end subroutine hydrostatic_zeq_coldens

   subroutine hydrostatic_zeq_densmid(iia,jja,d0,csim2)
      use arrays,     only: dprof
      use constants,  only: small
      use dataio_pub, only: die
      implicit none
      integer, intent(in) :: iia, jja
      real,    intent(in) :: d0, csim2

      if (d0 .le. small) then
         call die("[hydrostatic:hydrostatic_zeq_densmid] d0 must be /= 0")
      endif
#ifndef NEW_HYDROSTATIC
      dmid = d0
#endif /* !NEW_HYDROSTATIC */

      call start_hydrostatic(iia,jja,csim2)
#ifdef NEW_HYDROSTATIC
      dprof = dprof * d0
#endif /* NEW_HYDROSTATIC */
      call finish_hydrostatic

   end subroutine hydrostatic_zeq_densmid

   subroutine start_hydrostatic(iia,jja,csim2)
      use dataio_pub, only: die
      use gravity,    only: get_gprofs, gprofs_target, nsub
      use grid,       only: zmin, zmax, zdim, nzt, dl, nb
      implicit none
      integer, intent(in) :: iia, jja
      real,    intent(in) :: csim2
      integer :: ksub
      if (.not.associated(get_gprofs)) then
         select case (gprofs_target)
            case ('accel')
               get_gprofs => get_gprofs_accel
            case ('gparr')
               get_gprofs => get_gprofs_gparray
            case default
               call die("[hydrostatic:start_hydrostatic] get_gprofs'' target has not been specified")
         end select
      endif
      nstot = nsub * nzt
      dzs = (zmax-zmin)/real(nstot-2*nb*nsub)
      allocate(zs(nstot), gprofs(nstot))
      do ksub=1, nstot
         zs(ksub) = zmin-nb*dl(zdim) + dzs/2 + (ksub-1)*dzs
      enddo
      call get_gprofs(iia,jja)
      gprofs = gprofs / csim2
      call hydrostatic_main
   end subroutine start_hydrostatic

   subroutine finish_hydrostatic
      implicit none
      if (allocated(zs))     deallocate(zs)
      if (allocated(gprofs)) deallocate(gprofs)
   end subroutine finish_hydrostatic

   subroutine outh_bnd(kb,kk,minmax)
      use gravity,             only: grav_accel, nsub, tune_zeq_bnd
      use arrays,              only: u
      use grid,                only: nx, ny, z
      use mpisetup,            only: smalld
      use initfluids,          only: gamma, cs_iso2
      use fluidindex,          only: nvar, iarr_all_dn, iarr_all_mx, iarr_all_my, iarr_all_mz
#ifndef ISO
      use fluidindex,          only: iarr_all_en
      use mpisetup,            only: smallei
#endif /* !ISO */
#ifdef COSM_RAYS
      use fluidindex,          only: iarr_all_crs
      use initcosmicrays,      only: smallecr
#endif /* COSM_RAYS */

      implicit none
      integer, intent(in)           :: kb, kk
      character(len=*), intent(in)  :: minmax

      integer                             :: ksub, ifluid, i, j
      real, dimension(nvar%fluids,nx,ny)  :: db, csi2b
#ifndef ISO
      real, dimension(nvar%fluids,nx,ny)  :: ekb, eib
#endif /* !ISO */
      real, dimension(nsub+1)             :: zs, gprofs
      real, dimension(nvar%fluids,nsub+1) :: dprofs
      real, dimension(nvar%fluids)        :: factor
      real                                :: dzs,z1,z2


      db = u(iarr_all_dn,:,:,kb)
      db = max(db,smalld)
#ifdef ISO
      csi2b = cs_iso2
#else /* !ISO */
      ekb = 0.5*(u(iarr_all_mx,:,:,kb)**2+u(iarr_all_my,:,:,kb)**2+u(iarr_all_mz,:,:,kb)**2)/db
      eib = u(iarr_all_en,:,:,kb) - ekb
      eib = max(eib,smallei)
      do ifluid=1,nvar%fluids
         csi2b(ifluid,:,:) = (gamma(ifluid)-1.0)*eib(ifluid,:,:)/db(ifluid,:,:)
      enddo
#endif /* !ISO */
      z1 = z(kb)
      z2 = z(kk)
      dzs = (z2-z1)/real(nsub)

      do ksub=1, nsub+1
         zs(ksub) = z1 + dzs/2 + (ksub-1)*dzs
      enddo

      do j=1,ny
         do i=1,nx

            call grav_accel('zsweep',i,j, zs, nsub, gprofs)
            gprofs=tune_zeq_bnd * gprofs

            dprofs(:,1) = db(:,i,j)
            do ksub=1, nsub
               factor = (1.0 + 0.5*dzs*gprofs(ksub)/csi2b(:,i,j))  &
                        /(1.0 - 0.5*dzs*gprofs(ksub)/csi2b(:,i,j))
               dprofs(:,ksub+1) = factor * dprofs(:,ksub)
            enddo

            db(:,i,j)  = dprofs(:,nsub+1)
            db(:,i,j)  = max(db(:,i,j), smalld)

            u(iarr_all_dn,i,j,kk)      =     db(:,i,j)
            u(iarr_all_mx,i,j,kk)      =     u(iarr_all_mx,i,j,kb)
            u(iarr_all_my,i,j,kk)      =     u(iarr_all_my,i,j,kb)
            u(iarr_all_mz,i,j,kk)      =     u(iarr_all_mz,i,j,kb)
! zakomentowac nastepna linie jesli warunek diodowy nie ma byc stosowany razem z hydrostatycznym
!           if (minmax == 'max') then
!              u(iarr_all_mz,i,j,kk)          =     max(u(iarr_all_mz,i,j,kk),0.0)
!           else
!              u(iarr_all_mz,i,j,kk)          =     min(u(iarr_all_mz,i,j,kk),0.0)
!           endif
            if (.false.) print *, minmax
#ifndef ISO
            eib(:,i,j) = csi2b(:,i,j)*db(:,i,j)/(gamma-1)
            eib(:,i,j) = max(eib(:,i,j), smallei)
            u(iarr_all_en,i,j,kk)      =     ekb(:,i,j) + eib(:,i,j)
#endif /* !ISO */
#ifdef COSM_RAYS
            u(iarr_all_crs,i,j,kk)     =     smallecr
#endif /* COSM_RAYS */
         enddo ! i
      enddo ! j
   end subroutine outh_bnd

#endif /* GRAV */
end module hydrostatic
