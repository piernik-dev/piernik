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
!! \details There are two routines to call to set hydrostatic equilibrium:
!! @n hydrostatic_zeq_coldens that fixes column density,
!! @n hydrostatic_zeq_densmid that fixes density value in the midplane.
!! @n Additionally there is also outh_bnd routine to keep hydrostatic equilibrium on the boundaries.
!<
module hydrostatic
! pulled by GRAV
   implicit none

   private
#ifdef GRAV
   public :: hydrostatic_zeq_coldens, hydrostatic_zeq_densmid, gprofs, nstot, zs, dzs, outh_bnd
#endif /* GRAV */

   real, allocatable, dimension(:), save :: zs        !< array of z-positions of subgrid cells centers
   real, allocatable, dimension(:), save :: gprofs    !< array of gravitational acceleration in a column of subgrid
   real,                            save :: dzs       !< length of the subgrid cell in z-direction
   integer,                         save :: nstot     !< total number of subgrid cells in a column through all z-blocks
   real,                            save :: dmid      !< density value in a midplane (fixed for hydrostatic_zeq_densmid, overwritten by hydrostatic_zeq_coldens)

   interface
      subroutine hzeqscheme(ksub, up, factor)
         implicit none
         integer, intent(in)  :: ksub
         real,    intent(in)  :: up
         real,    intent(out) :: factor
      end subroutine hzeqscheme
   end interface

   procedure(hzeqscheme), pointer :: hzeq_scheme => NULL()

contains
#ifdef GRAV
!>
!! \brief Routine that establishes hydrostatic equilibrium for fixed column density
!! \details Routine calls the routine of the case of fixed plane density value and use the correction for column density.
!! To properly use this routine it is important to make sure that get_gprofs pointer has been associated. See details of start_hydrostatic routine.
!! \param iia x index of z-column
!! \param jja y index of z-column
!! \param coldens column density value for given x and y coordinates
!! \param csim2 sqare of sound velocity
!<
   subroutine hydrostatic_zeq_coldens(iia,jja,coldens,csim2)

      use arrays,   only: dprof  !! \todo convert to intent(inout)

      implicit none

      integer, intent(in)   :: iia, jja
      real,    intent(in)   :: coldens, csim2
      real                  :: sdprof, sd

      sdprof = 1.0
      call hydrostatic_zeq_densmid(iia,jja,sdprof,csim2,sd)
      dprof = dprof * coldens / sd

   end subroutine hydrostatic_zeq_coldens

!>
!! \brief Routne that establishes hydrostatic equilibrium for fixed plane density value
!! \details To properly use this routine it is important to make sure that get_gprofs pointer has been associated. See details of start_hydrostatic routine.
!! \param iia x-coordinate of z-column
!! \param jja y-coordinate of z-column
!! \param d0 plane density value for given x and y coordinates
!! \param csim2 sqare of sound velocity
!! \param sd optional variable to give a sum of dprofs array from hydrostatic_main routine
!<
   subroutine hydrostatic_zeq_densmid(iia,jja,d0,csim2,sd)

      use constants,  only: small
      use dataio_pub, only: die

      implicit none

      integer, intent(in) :: iia, jja
      real,    intent(in) :: d0, csim2
      real,    intent(inout), optional :: sd

      if (d0 <= small) then
         call die("[hydrostatic:hydrostatic_zeq_densmid] d0 must be /= 0")
      endif
      dmid = d0

      call start_hydrostatic(iia,jja,csim2,sd)
      call finish_hydrostatic

   end subroutine hydrostatic_zeq_densmid

!>
!! \brief Routine that arranges %hydrostatic equilibrium in the vertical (z) direction
!<
   subroutine hydrostatic_main(sd)

      use arrays,     only: dprof
      use dataio_pub, only: die
      use gravity,    only: nsub
      use grid,       only: cg
      use mpisetup,   only: dom

      implicit none

      real, intent(out), optional     :: sd
      real, allocatable, dimension(:) :: dprofs
      integer :: ksub, ksmid, k
      real    :: factor

      allocate(dprofs(nstot))

      ksmid = 0
#ifdef NEW_HYDROSTATIC
      ksmid = minloc(abs(gprofs),1)          ! generally the midplane is where gravity is 0,  practically we want the least gravity absolute value
      hzeq_scheme => hzeq_scheme_v2
#else /* !NEW_HYDROSTATIC */
      ksmid = maxloc(zs,1,mask=(zs < 0.0))   ! the midplane is in between ksmid and ksmid+1
      hzeq_scheme => hzeq_scheme_v1
#endif /* !NEW_HYDROSTATIC */
      if (ksmid == 0) call die("[hydrostatic:hydrostatic_main] ksmid not set")

      if (ksmid < nstot) then
         dprofs(ksmid+1) = dmid
         do ksub=ksmid+1, nstot-1
            call hzeq_scheme(ksub,  1.0, factor)
            dprofs(ksub+1) = factor * dprofs(ksub)
         enddo
      endif

      if (ksmid > 1) then
         dprofs(ksmid) = dmid
         do ksub=ksmid, 2, -1
            call hzeq_scheme(ksub, -1.0, factor)
            dprofs(ksub-1) = factor * dprofs(ksub)
         enddo
      endif

      dprof(:) =0.0
      do k=1, cg%nz
         do ksub=1, nstot
            if (zs(ksub) > cg%zl(k) .and. zs(ksub) < cg%zr(k)) then
               dprof(k) = dprof(k) + dprofs(ksub)/real(nsub)
            endif
         enddo
      enddo
      if (present(sd)) then
         sd = 0.0
         do ksub=1, nstot
            if (zs(ksub) > dom%zmin .and. zs(ksub) < dom%zmax) then
               sd = sd + dprofs(ksub)*dzs
            endif
         enddo
      endif

      if (allocated(dprofs)) deallocate(dprofs)

      return
   end subroutine hydrostatic_main

   subroutine hzeq_scheme_v1(ksub, up, factor)
      implicit none
      integer, intent(in)  :: ksub
      real,    intent(in)  :: up
      real,    intent(out) :: factor
      factor = (2.0 + up*gprofs(ksub))/(2.0 - up*gprofs(ksub))
   end subroutine hzeq_scheme_v1

   subroutine hzeq_scheme_v2(ksub, up, factor)
      implicit none
      integer, intent(in)  :: ksub
      real,    intent(in)  :: up
      real,    intent(out) :: factor
      factor = gprofs(ksub)+gprofs(ksub+nint(up))
      factor = (4.0 + up*factor)/(4.0 - up*factor)
   end subroutine hzeq_scheme_v2

   subroutine get_gprofs_accel(iia,jja)

      use gravity,   only: tune_zeq, grav_accel
      use grid,      only: cg
      use constants, only: zdim

      implicit none

      integer, intent(in) :: iia, jja
      integer :: ia, ja

      ia = min(cg%nx,max(1, iia))
      ja = min(cg%ny,max(1, jja))
      call grav_accel(zdim, ia, ja, zs, nstot, gprofs)
      gprofs = tune_zeq*gprofs

   end subroutine get_gprofs_accel

!>
!! \brief Routine that has to offer a z-sweep of gravity potential with extended z-grid
!! \deprecated probably now the routine should have different name than gparray which got from the commented part of code
!! \warning in case of moving 'use types, only: axes'' behind use gravity there could be gcc(4.5) internal compiler error: in fold_convert_loc, at fold-const.c:2792 (solved in >=gcc-4.6)
!<
   subroutine get_gprofs_gparray(iia,jja)

      use types,   only: axes
      use arrays,  only: gp
      use gravity, only: tune_zeq, grav_type
      use grid,    only: cg

      implicit none

      integer, intent(in)                  :: iia, jja
      real, allocatable, dimension(:,:,:)  :: gpots
      type(axes)                           :: ax
      integer                              :: nstot1

      nstot1 = nstot + 1
      allocate(gpots(1,1,nstot1))
      if (.not.allocated(ax%x)) allocate(ax%x(1))
      if (.not.allocated(ax%y)) allocate(ax%y(1))
      if (.not.allocated(ax%z)) allocate(ax%z(nstot1))
      ax%x          = cg%x(iia)
      ax%y          = cg%y(jja)
      ax%z(1:nstot) = zs - 0.5*dzs
      ax%z(nstot1)  = ax%z(nstot) + dzs
      call grav_type(gpots,ax)
      gprofs(1:nstot) = (gpots(1,1,1:nstot) - gpots(1,1,2:nstot1))/dzs
      gprofs = tune_zeq*gprofs
      if (allocated(gpots)) deallocate(gpots)
      if (allocated(ax%x))  deallocate(ax%x)
      if (allocated(ax%y))  deallocate(ax%y)
      if (allocated(ax%z))  deallocate(ax%z)
   end subroutine get_gprofs_gparray

!>
!! \brief Routine that prepares data for constructing hydrostatic equilibrium by hydrostatic_main routine
!! \details It is important to have get_gprofs pointer associated to a proper routine that gives back the column of nsub*nzt elements of gravitational acceleration in z direction. In the most common cases the gprofs_target parameter from GRAVITY namelist may be used. When it is set to 'accel' or 'gparr' the pointer is associated to get_gprofs_accel or get_gprofs_gparray routines, respectively.
!! \note After calling this routine gprofs is multiplied by dzs/csim2 which are assumed to be constant. This is done for optimizing the hydrostatic_main routine.
!! \param iia x-coordinate of z-column
!! \param jja y-coordinate of z-column
!! \param csim2 sqare of sound velocity
!! \param sd optional variable to give a sum of dprofs array from hydrostatic_main routine
!<
   subroutine start_hydrostatic(iia,jja,csim2,sd)

      use dataio_pub, only: die
      use gravity,    only: get_gprofs, gprofs_target, nsub
      use grid,       only: cg
      use constants,  only: zdim
      use mpisetup,   only: dom

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
      nstot = nsub * dom%nzt
      dzs = (dom%zmax-dom%zmin)/real(nstot-2*cg%nb*nsub)
      allocate(zs(nstot), gprofs(nstot))
      do ksub=1, nstot
         zs(ksub) = dom%zmin-cg%nb*cg%dl(zdim) + (real(ksub)-0.5)*dzs
      enddo
      call get_gprofs(iia,jja)
      gprofs = gprofs / csim2 *dzs
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
      use grid,                only: cg
      use mpisetup,            only: smalld
      use constants,           only: zdim
      use fluidindex,          only: flind, iarr_all_dn, iarr_all_mx, iarr_all_my, iarr_all_mz
#ifndef ISO
      use fluidindex,          only: iarr_all_en
      use mpisetup,            only: smallei
#endif /* !ISO */
#ifdef COSM_RAYS
      use fluidindex,          only: iarr_all_crs
      use initcosmicrays,      only: smallecr
#endif /* COSM_RAYS */

      implicit none

      integer,          intent(in)                :: kb, kk
      character(len=*), intent(in)                :: minmax

      integer                                     :: ksub, i, j
      real, dimension(flind%fluids, cg%nx, cg%ny) :: db, csi2b
#ifndef ISO
      integer                                     :: ifluid
      real, dimension(flind%fluids, cg%nx, cg%ny) :: ekb, eib
#endif /* !ISO */
      real, dimension(nsub+1)                     :: zs, gprofs
      real, dimension(flind%fluids,nsub+1)        :: dprofs
      real, dimension(flind%fluids)               :: factor
      real                                        :: dzs,z1,z2

      if (.not.associated(grav_accel)) call die("[hydrostatic:outh_bnd] grav_accel not associated")

      db = cg%u%arr(iarr_all_dn,:,:,kb)
      db = max(db,smalld)
#ifdef ISO
      csi2b = maxval(flind%all_fluids(:)%cs2)   !> \deprecated BEWARE should be fluid dependent
#else /* !ISO */
      ekb = 0.5*(cg%u%arr(iarr_all_mx,:,:,kb)**2+cg%u%arr(iarr_all_my,:,:,kb)**2+cg%u%arr(iarr_all_mz,:,:,kb)**2)/db
      eib = cg%u%arr(iarr_all_en,:,:,kb) - ekb
      eib = max(eib,smallei)
      do ifluid=1,flind%fluids
         csi2b(ifluid,:,:) = (flind%all_fluids(ifluid)%gam_1)*eib(ifluid,:,:)/db(ifluid,:,:)
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

            call grav_accel(zdim, i, j, zs, nsub, gprofs)
            gprofs=tune_zeq_bnd * gprofs

            dprofs(:,1) = db(:,i,j)
            do ksub=1, nsub
               factor = (1.0 + 0.5*dzs*gprofs(ksub)/csi2b(:,i,j))  &
                        /(1.0 - 0.5*dzs*gprofs(ksub)/csi2b(:,i,j))
               dprofs(:,ksub+1) = factor * dprofs(:,ksub)
            enddo

            db(:,i,j)  = dprofs(:,nsub+1)
            db(:,i,j)  = max(db(:,i,j), smalld)

            cg%u%arr(iarr_all_dn,i,j,kk)      =     db(:,i,j)
            cg%u%arr(iarr_all_mx,i,j,kk)      =     cg%u%arr(iarr_all_mx,i,j,kb)
            cg%u%arr(iarr_all_my,i,j,kk)      =     cg%u%arr(iarr_all_my,i,j,kb)
            cg%u%arr(iarr_all_mz,i,j,kk)      =     cg%u%arr(iarr_all_mz,i,j,kb)
            !> \deprecated to use outh together with outd user should manually interfere in the code of outh_bnd routine
! zakomentowac nastepna linie jesli warunek diodowy nie ma byc stosowany razem z hydrostatycznym
!           if (minmax == 'max') then
!              cg%u%arr(iarr_all_mz,i,j,kk)          =     max(cg%u%arr(iarr_all_mz,i,j,kk),0.0)
!           else
!              cg%u%arr(iarr_all_mz,i,j,kk)          =     min(cg%u%arr(iarr_all_mz,i,j,kk),0.0)
!           endif
            if (.false.) print *, minmax
#ifndef ISO
            eib(:,i,j) = csi2b(:,i,j)*db(:,i,j)/(flind%all_fluids(:)%gam_1)
            eib(:,i,j) = max(eib(:,i,j), smallei)
            cg%u%arr(iarr_all_en,i,j,kk)      =     ekb(:,i,j) + eib(:,i,j)
#endif /* !ISO */
#ifdef COSM_RAYS
            cg%u%arr(iarr_all_crs,i,j,kk)     =     smallecr
#endif /* COSM_RAYS */
         enddo ! i
      enddo ! j
   end subroutine outh_bnd

#endif /* GRAV */
end module hydrostatic
