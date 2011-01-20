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
!! \brief (DW) Module containing a subroutine that arranges %hydrostatic equilibrium in the vertical (z) direction
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
   subroutine hydrostatic_main(sd)

      use arrays,     only: dprof
      use dataio_pub, only: die
      use gravity,    only: nsub !, gp_status
      use grid,       only: cg

      implicit none

      real, intent(out), optional     :: sd
      real, allocatable, dimension(:) :: dprofs
      integer :: ksub, ksmid, k
#ifndef NEW_HYDROSTATIC
      real    :: factor
#endif /* !NEW_HYDROSTATIC */

      allocate(dprofs(nstot))

      ksmid = 0
      do ksub=1, nstot
         if (zs(ksub) < 0.0) ksmid = ksub      ! the midplane is in between
      enddo                                  ! ksmid and ksmid+1
      if (ksmid == 0) call die("[hydrostatic:hydrostatic_main] ksmid not set")

#ifdef NEW_HYDROSTATIC
      if (ksmid < nstot) then
         dprofs(ksmid+1) = 1.0
         do ksub=ksmid+1, nstot-1
            dprofs(ksub+1) = dprofs(ksub)*(1.0 + 0.5*(gprofs(ksub)+gprofs(ksub+1))*dzs)
         enddo
      endif

      if (ksmid > 1) then
         dprofs(ksmid) = 1.0
         do ksub=ksmid, 2, -1
            dprofs(ksub-1) = dprofs(ksub)*(1.0 - 0.5*(gprofs(ksub)+gprofs(ksub-1))*dzs)
         enddo
      endif
#else /* !NEW_HYDROSTATIC */
      if (ksmid < nstot) then
         dprofs(ksmid+1) = dmid
         do ksub=ksmid+1, nstot-1
            factor = (2.0 + dzs*gprofs(ksub))/(2.0 - dzs*gprofs(ksub))
            dprofs(ksub+1) = factor * dprofs(ksub)
         enddo
      endif

      if (ksmid > 1) then
         dprofs(ksmid) = dmid
         do ksub=ksmid, 2, -1
            factor = (2.0 - dzs*gprofs(ksub))/(2.0 + dzs*gprofs(ksub))
            dprofs(ksub-1) = factor * dprofs(ksub)
         enddo
      endif
#endif /* !NEW_HYDROSTATIC */

      dprof(:) =0.0
      if (present(sd)) sd = 0.0
      do k=1, cg%nz
         do ksub=1, nstot
            if (zs(ksub) > cg%zl(k) .and. zs(ksub) < cg%zr(k)) then
               dprof(k) = dprof(k) + dprofs(ksub)/real(nsub)
            endif
            if (present(sd)) then
               if (zs(ksub) .gt. cg%zmin .and. zs(ksub) .lt. cg%zmax) then
                  sd = sd + dprofs(ksub)/real(nsub)
               endif
            endif
         enddo
      enddo

      if (allocated(dprofs)) deallocate(dprofs)

      return
   end subroutine hydrostatic_main

   subroutine get_gprofs_accel(iia,jja)
      use gravity, only: tune_zeq, grav_accel
      use grid,    only: cg
      implicit none
      integer, intent(in) :: iia, jja
      integer :: ia, ja

      ia = min(cg%nx,max(1, iia))
      ja = min(cg%ny,max(1, jja))
      call grav_accel('zsweep',ia, ja, zs, nstot, gprofs)
      gprofs = tune_zeq*gprofs

   end subroutine get_gprofs_accel

!>
!! \brief Routine that has to offer a z-sweep of gravity potential with extended z-grid
!! \deprecated probably now the routine should have different name than gparray which got from the commented part of code
!<
   subroutine get_gprofs_gparray(iia,jja)

      use arrays,  only: gp
      use gravity, only: tune_zeq, grav_type
      use grid,    only: cg
      use types,   only: axes

      implicit none

      integer, intent(in)                  :: iia, jja
      integer                              :: ksub, k
      real, allocatable, dimension(:,:,:)  :: gpots
      type(axes)                           :: ax

      allocate(gpots(1,1,nstot))
! BEWARE: good approximation only for pzsize == 1 (eventuall for pzsize == 2 still works well, but never for higher psize)
!      k = 1; gpots(:) = 0.0
!      do ksub=1, nstot
!         if (zs(ksub) >= cg%z(min(k+1, cg%nz))) k = k + 1
!         if (zs(ksub) >= cg%z(min(k, cg%nz-1)) .and. zs(ksub) < cg%z(min(k+1, cg%nz))) then
!            gpots(ksub) = gp(iia,jja,k) + (zs(ksub) - cg%z(k)) * &
!             (gp(iia,jja,min(k+1, cg%nz)) - gp(iia,jja,k)) / (cg%z(min(k+1, cg%nz)) - cg%z(min(k, cg%nz-1)))
!         endif
!      enddo
!         call grav_pot('zsweep', ia,ja, zs, nstot, gpots,gp_status,.true.)
      if (.not.allocated(ax%x)) allocate(ax%x(1))
      if (.not.allocated(ax%y)) allocate(ax%y(1))
      if (.not.allocated(ax%z)) allocate(ax%z(nstot))
      ax%x = cg%x(iia)
      ax%y = cg%y(jja)
      ax%z = zs
      call grav_type(gpots,ax)
      gprofs(1:nstot-1) = (gpots(1,1,1:nstot-1) - gpots(1,1,2:nstot))/dzs
      gprofs(nstot) = 0. ! or maybe gprofs(nstot-1) ?
      gprofs = tune_zeq*gprofs
      if (allocated(gpots)) deallocate(gpots)
      if (allocated(ax%x))  deallocate(ax%x)
      if (allocated(ax%y))  deallocate(ax%y)
      if (allocated(ax%z))  deallocate(ax%z)

   end subroutine get_gprofs_gparray

   subroutine hydrostatic_zeq_coldens(iia,jja,coldens,csim2)

      use arrays,   only: dprof
      use grid,     only: cg

      implicit none

      integer, intent(in)   :: iia, jja
      real,    intent(in)   :: coldens, csim2
      real                  :: sdprof, sd
      integer               :: comm1d
      logical, dimension(3) :: remain

      sdprof = 1.0
      call hydrostatic_zeq_densmid(iia,jja,sdprof,csim2,sd)
      sdprof = sd * cg%dz
      dprof = dprof * coldens / sdprof

   end subroutine hydrostatic_zeq_coldens

   subroutine hydrostatic_zeq_densmid(iia,jja,d0,csim2,sd)

      use arrays,     only: dprof
      use constants,  only: small
      use dataio_pub, only: die

      implicit none

      integer, intent(in) :: iia, jja
      real,    intent(in) :: d0, csim2
      real,    intent(inout), optional :: sd

      if (d0 <= small) then
         call die("[hydrostatic:hydrostatic_zeq_densmid] d0 must be /= 0")
      endif
#ifndef NEW_HYDROSTATIC
      dmid = d0
#endif /* !NEW_HYDROSTATIC */

      call start_hydrostatic(iia,jja,csim2,sd)
#ifdef NEW_HYDROSTATIC
      dprof = dprof * d0
#endif /* NEW_HYDROSTATIC */
      call finish_hydrostatic

   end subroutine hydrostatic_zeq_densmid

   subroutine start_hydrostatic(iia,jja,csim2,sd)

      use dataio_pub, only: die
      use gravity,    only: get_gprofs, gprofs_target, nsub
      use grid,       only: cg
      use mpisetup,   only: zdim

      implicit none

      integer, intent(in) :: iia, jja
      real,    intent(in) :: csim2
      real,    intent(inout), optional :: sd
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
      nstot = nsub * cg%nzt
      dzs = (cg%zmax-cg%zmin)/real(nstot-2*cg%nb*nsub)
      allocate(zs(nstot), gprofs(nstot))
      do ksub=1, nstot
         zs(ksub) = cg%zmin-cg%nb*cg%dl(zdim) + dzs/2 + (ksub-1)*dzs
      enddo
      call get_gprofs(iia,jja)
      gprofs = gprofs / csim2
      call hydrostatic_main(sd)

   end subroutine start_hydrostatic

   subroutine finish_hydrostatic

      implicit none

      if (allocated(zs))     deallocate(zs)
      if (allocated(gprofs)) deallocate(gprofs)

   end subroutine finish_hydrostatic

   subroutine outh_bnd(kb,kk,minmax)

      use dataio_pub,          only: die
      use gravity,             only: grav_accel, nsub, tune_zeq_bnd
      use arrays,              only: u
      use grid,                only: cg
      use mpisetup,            only: smalld
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

      integer                             :: ksub, i, j
      real, dimension(nvar%fluids, cg%nx, cg%ny) :: db, csi2b
#ifndef ISO
      integer                             :: ifluid
      real, dimension(nvar%fluids, cg%nx, cg%ny) :: ekb, eib
#endif /* !ISO */
      real, dimension(nsub+1)             :: zs, gprofs
      real, dimension(nvar%fluids,nsub+1) :: dprofs
      real, dimension(nvar%fluids)        :: factor
      real                                :: dzs,z1,z2

      if (.not.associated(grav_accel)) call die("[hydrostatic:outh_bnd] grav_accel not associated")

      db = u(iarr_all_dn,:,:,kb)
      db = max(db,smalld)
#ifdef ISO
      csi2b = maxval(nvar%all_fluids(:)%cs2)   !> \deprecated BEWARE should be fluid dependent
#else /* !ISO */
      ekb = 0.5*(u(iarr_all_mx,:,:,kb)**2+u(iarr_all_my,:,:,kb)**2+u(iarr_all_mz,:,:,kb)**2)/db
      eib = u(iarr_all_en,:,:,kb) - ekb
      eib = max(eib,smallei)
      do ifluid=1,nvar%fluids
         csi2b(ifluid,:,:) = (nvar%all_fluids(ifluid)%gam_1)*eib(ifluid,:,:)/db(ifluid,:,:)
      enddo
#endif /* !ISO */
      z1 = cg%z(kb)
      z2 = cg%z(kk)
      dzs = (z2-z1)/real(nsub)

      do ksub=1, nsub+1
         zs(ksub) = z1 + dzs/2 + (ksub-1)*dzs
      enddo

      do j=1, cg%ny
         do i=1, cg%nx

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
            !> \deprecated to use outh together with outd user should manually interfere in the code of outh_bnd routine
! zakomentowac nastepna linie jesli warunek diodowy nie ma byc stosowany razem z hydrostatycznym
!           if (minmax == 'max') then
!              u(iarr_all_mz,i,j,kk)          =     max(u(iarr_all_mz,i,j,kk),0.0)
!           else
!              u(iarr_all_mz,i,j,kk)          =     min(u(iarr_all_mz,i,j,kk),0.0)
!           endif
            if (.false.) print *, minmax
#ifndef ISO
            eib(:,i,j) = csi2b(:,i,j)*db(:,i,j)/(nvar%all_fluids(:)%gam_1)
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
