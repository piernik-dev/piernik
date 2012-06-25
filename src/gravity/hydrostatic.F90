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
   use grid_cont, only: grid_container

   implicit none

   private
#ifdef GRAV
   public :: set_default_hsparams, hydrostatic_zeq_coldens, hydrostatic_zeq_densmid, cleanup_hydrostatic, outh_bnd
   public :: dprof, gprofs, nstot, zs, dzs, hsmin, hsbn, hstn, hsl, sdlim, hscg
#endif /* GRAV */

   real, allocatable, dimension(:), save :: zs        !< array of z-positions of subgrid cells centers
   real, allocatable, dimension(:), save :: gprofs    !< array of gravitational acceleration in a column of subgrid
   real, allocatable, dimension(:), save :: dprof     !< Array used for storing density during calculation of hydrostatic equilibrium
   real,                            save :: dzs       !< length of the subgrid cell in z-direction
   integer(kind=4),                 save :: nstot     !< total number of subgrid cells in a column through all z-blocks
   real,                            save :: dmid      !< density value in a midplane (fixed for hydrostatic_zeq_densmid, overwritten by hydrostatic_zeq_coldens)
   real,                            save :: hsmin     !< lower position limit
   integer(kind=4),                 save :: hsbn      !< number of cells in proceeded block
   integer,                         save :: hstn      !< number of cells in proceeded domain
   real, allocatable, dimension(:), save :: hsl       !< lower borders of cells of proceeded block
   real, dimension(2),              save :: sdlim     !< edges for sd sum
   type(grid_container), pointer,   save :: hscg

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
   subroutine hydrostatic_zeq_coldens(iia, jja, coldens, csim2)

      implicit none

      integer, intent(in)    :: iia, jja
      real,    intent(in)    :: coldens, csim2
      real                   :: sdprof, sd

      sdprof = 1.0
      call hydrostatic_zeq_densmid(iia, jja, sdprof, csim2, sd)
      dprof = dprof * coldens / sd

   end subroutine hydrostatic_zeq_coldens

!>
!! \brief Routine that establishes hydrostatic equilibrium for fixed plane density value
!! \details To properly use this routine it is important to make sure that get_gprofs pointer has been associated. See details of start_hydrostatic routine.
!! \param iia x-coordinate of z-column
!! \param jja y-coordinate of z-column
!! \param d0 plane density value for given x and y coordinates
!! \param csim2 sqare of sound velocity
!! \param sd optional variable to give a sum of dprofs array from hydrostatic_main routine
!<
   subroutine hydrostatic_zeq_densmid(iia, jja, d0, csim2, sd)

      use constants,  only: small
      use dataio_pub, only: die

      implicit none

      integer,        intent(in)    :: iia, jja
      real,           intent(in)    :: d0, csim2
      real, optional, intent(inout) :: sd

      if (d0 <= small) call die("[hydrostatic:hydrostatic_zeq_densmid] d0 must be /= 0")
      dmid = d0

      call start_hydrostatic(iia, jja, csim2, sd)
      call finish_hydrostatic

   end subroutine hydrostatic_zeq_densmid

!>
!! \brief Routine to set up sizes of arrays used in hydrostatic module. Settings depend on cg structure.
!! \details Routine has to be called before the firs usage of hydrostatic_zeq_coldens/densmid if there is no other equivalent user settings.
!<
   subroutine set_default_hsparams(cg)

      use constants,   only: zdim, LO, HI, I_ONE
      use diagnostics, only: my_allocate, my_deallocate
      use domain,      only: dom
      use gravity,     only: nsub
      use grid_cont,   only: grid_container

      implicit none

      type(grid_container), pointer, intent(in) :: cg

      hscg => cg

      nstot = nsub * dom%n_t(zdim)
      dzs   = (dom%edge(zdim, HI)-dom%edge(zdim, LO))/real(nstot-2*dom%nb*nsub)
      hsmin = dom%edge(zdim, LO)-dom%nb*cg%dl(zdim)
      hsbn  = cg%n_(zdim)
      hstn  = dom%n_t(zdim)
      sdlim = dom%edge(zdim,:)
      if (allocated(dprof)) call my_deallocate(dprof)
      call my_allocate(dprof, [cg%n_(zdim)], "dprof")
      if (allocated(hsl)) call my_deallocate(hsl)
      call my_allocate(hsl, [hsbn+I_ONE], "hsl")
      hsl(1:hsbn) = cg%zl(1:hsbn)
      hsl(hsbn+1) = cg%zr(hsbn)

   end subroutine set_default_hsparams

!>
!! \brief Routine to clean up after the last usage of hydrostatic routines
!<
   subroutine cleanup_hydrostatic

      use diagnostics, only: my_deallocate

      implicit none

      if (allocated(dprof)) call my_deallocate(dprof)
      if (allocated(hsl))   call my_deallocate(hsl)
      if (associated(hscg)) nullify(hscg)

   end subroutine cleanup_hydrostatic

!>
!! \brief Routine that arranges %hydrostatic equilibrium in the vertical (z) direction
!<
   subroutine hydrostatic_main(sd)

      use constants,  only: LO, HI
      use dataio_pub, only: die
      use gravity,    only: nsub

      implicit none

      real, optional, intent(out)     :: sd
      real, allocatable, dimension(:) :: dprofs
      integer                         :: ksub, ksmid, k
      real                            :: factor

      allocate(dprofs(nstot))

      ksmid = 0
#ifdef HYDROSTATIC_V2
      ksmid = minloc(abs(gprofs),1)          ! generally the midplane is where gravity is 0,  practically we want the least gravity absolute value
      hzeq_scheme => hzeq_scheme_v2
#else /* !HYDROSTATIC_V2 */
      ksmid = maxloc(zs,1,mask=(zs < 0.0))   ! the midplane is in between ksmid and ksmid+1
      hzeq_scheme => hzeq_scheme_v1
#endif /* !HYDROSTATIC_V2 */
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

      dprof(:) = 0.0
      do k=1, hsbn
         do ksub=1, nstot
            if (zs(ksub) > hsl(k) .and. zs(ksub) < hsl(k+1)) then
               dprof(k) = dprof(k) + dprofs(ksub)/real(nsub)
            endif
         enddo
      enddo

      if (present(sd)) then
         sd = 0.0
         do ksub=1, nstot
            if (zs(ksub) > sdlim(LO) .and. zs(ksub) < sdlim(HI)) sd = sd + dprofs(ksub)*dzs
         enddo
      endif

      if (allocated(dprofs)) deallocate(dprofs)

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

   subroutine get_gprofs_accel(iia, jja)

      use constants, only: xdim, ydim, zdim
      use gravity,   only: tune_zeq, grav_accel

      implicit none

      integer, intent(in) :: iia, jja
      integer             :: ia, ja

      ia = min(hscg%n_(xdim), int(max(1, iia), kind=4))
      ja = min(hscg%n_(ydim), int(max(1, jja), kind=4))
      call grav_accel(zdim, ia, ja, zs, nstot, gprofs)
      gprofs = tune_zeq*gprofs

   end subroutine get_gprofs_accel

!>
!! \brief Routine that has to offer a z-sweep of external gravity potential with extended z-grid
!! \warning in case of moving 'use types, only: axes'' behind use gravity there could be gcc(4.5) internal compiler error: in fold_convert_loc, at fold-const.c:2792 (solved in >=gcc-4.6)
!<
   subroutine get_gprofs_extgp(iia, jja)

      use constants, only: half
      use gravity,   only: tune_zeq, grav_type
      use types,     only: axes

      implicit none

      integer, intent(in)             :: iia, jja
      real, dimension(:,:,:), pointer :: gpots
      type(axes)                      :: ax
      integer                         :: nstot1

      nstot1 = nstot + 1
      allocate(gpots(1,1,nstot1))
      if (.not.allocated(ax%x)) allocate(ax%x(1))
      if (.not.allocated(ax%y)) allocate(ax%y(1))
      if (.not.allocated(ax%z)) allocate(ax%z(nstot1))
      ax%x          = hscg%x(iia)
      ax%y          = hscg%y(jja)
      ax%z(1:nstot) = zs - half*dzs
      ax%z(nstot1)  = ax%z(nstot) + dzs
      call grav_type(gpots,ax)
      gprofs(1:nstot) = (gpots(1,1,1:nstot) - gpots(1,1,2:nstot1))/dzs
      gprofs = tune_zeq*gprofs
      if (associated(gpots)) deallocate(gpots)
      if (allocated(ax%x))   deallocate(ax%x)
      if (allocated(ax%y))   deallocate(ax%y)
      if (allocated(ax%z))   deallocate(ax%z)
   end subroutine get_gprofs_extgp

!>
!! \brief Routine that prepares data for constructing hydrostatic equilibrium by hydrostatic_main routine
!! \details It is important to have get_gprofs pointer associated to a proper routine that gives back the column of nsub*nzt elements of gravitational acceleration in z direction. In the most common cases the gprofs_target parameter from GRAVITY namelist may be used. When it is set to 'accel' or 'extgp' the pointer is associated to get_gprofs_accel or get_gprofs_extgp routines, respectively.
!! \note After calling this routine gprofs is multiplied by dzs/csim2 which are assumed to be constant. This is done for optimizing the hydrostatic_main routine.
!! \param iia x-coordinate of z-column
!! \param jja y-coordinate of z-column
!! \param csim2 sqare of sound velocity
!! \param sd optional variable to give a sum of dprofs array from hydrostatic_main routine
!<
   subroutine start_hydrostatic(iia, jja, csim2, sd)

      use constants,  only: half
      use dataio_pub, only: die
      use gravity,    only: get_gprofs, gprofs_target

      implicit none

      integer,        intent(in)    :: iia, jja
      real,           intent(in)    :: csim2
      real, optional, intent(inout) :: sd
      integer                       :: ksub

      if (.not.associated(get_gprofs)) then
         select case (gprofs_target)
            case ('accel')
               get_gprofs => get_gprofs_accel
            case ('extgp')
               get_gprofs => get_gprofs_extgp
            case default
               call die("[hydrostatic:start_hydrostatic] get_gprofs'' target has not been specified")
         end select
      endif
      allocate(zs(nstot), gprofs(nstot))
      do ksub=1, nstot
         zs(ksub) = hsmin + (real(ksub)-half)*dzs
      enddo
      call get_gprofs(iia, jja)
      gprofs = gprofs / csim2 *dzs
      call hydrostatic_main(sd)

   end subroutine start_hydrostatic

   subroutine finish_hydrostatic

      implicit none

      if (allocated(zs))     deallocate(zs)
      if (allocated(gprofs)) deallocate(gprofs)

   end subroutine finish_hydrostatic

   !>
   !! \todo this procedure is incompatible with cg%cs_iso2
   !<

   subroutine outh_bnd(side, cg, diode)

      use constants,      only: xdim, ydim, zdim, half, HI, INT4
      use dataio_pub,     only: die
      use domain,         only: dom
      use fluidindex,     only: flind, iarr_all_dn, iarr_all_mx, iarr_all_my, iarr_all_mz
      use func,           only: ekin
      use global,         only: smalld
      use gravity,        only: nsub, get_gprofs, tune_zeq_bnd
      use grid_cont,      only: grid_container
#ifndef ISO
      use fluidindex,     only: iarr_all_en
      use global,         only: smallei
#endif /* !ISO */
#ifdef COSM_RAYS
      use fluidindex,     only: iarr_all_crs
      use initcosmicrays, only: smallecr
#endif /* COSM_RAYS */

      implicit none

      integer(kind=4),               intent(in)    :: side
      type(grid_container), pointer, intent(inout) :: cg
      logical,                       intent(in)    :: diode

      integer(kind=4)                              :: ib, ssign, kb, kk
      integer                                      :: ksub, i, j, lksub
      real, dimension(:,:,:), allocatable          :: db, csi2b, dbr
      real, dimension(:,:),   allocatable          :: dprofs
      real, dimension(flind%fluids)                :: factor
#ifndef ISO
      real, dimension(:,:,:), allocatable          :: ekb, eib
      integer                                      :: ifluid
#endif /* !ISO */

      if (.not.associated(get_gprofs)) call die("[hydrostatic:outh_bnd] get_gprofs not associated")

      hscg => cg
      nstot = int(3*nsub/2+1,kind=4)
      allocate(zs(nstot), gprofs(nstot), dprofs(flind%fluids,nstot))

      if (any([allocated(db), allocated(csi2b), allocated(dbr)])) call die("[hydrostatic:outh_bnd] db, dbr or csi2b already allocated")
      allocate(db(flind%fluids, cg%n_(xdim), cg%n_(ydim)), csi2b(flind%fluids, cg%n_(xdim), cg%n_(ydim)), dbr(flind%fluids, cg%n_(xdim), cg%n_(ydim)))
#ifndef ISO
      if (any([allocated(ekb), allocated(eib)])) call die("[hydrostatic:outh_bnd] ekb or eib already allocated")
      allocate(ekb(flind%fluids, cg%n_(xdim), cg%n_(ydim)), eib(flind%fluids, cg%n_(xdim), cg%n_(ydim)))
#endif /* !ISO */

      ssign = 2_INT4*side - 3_INT4
      dzs = (cg%z(cg%ijkse(zdim,side)+ssign)-cg%z(cg%ijkse(zdim,side)))/real(nsub)
      dbr = 1.0
      do ib=0_INT4, dom%nb
         kb = cg%ijkse(zdim,side)+ssign*(ib-1_INT4)
         kk = kb + ssign
         zs(:) = cg%z(kb) + dzs*(real([(ksub,ksub=1,nstot)])+real(nsub-3)*half)

         db = cg%u(iarr_all_dn,:,:,kb)
         db = max(db,smalld)
#ifdef ISO
         csi2b = maxval(flind%all_fluids(:)%cs2)   !> \deprecated BEWARE should be fluid dependent
#else /* !ISO */
         ekb = ekin(cg%u(iarr_all_mx,:,:,kb),cg%u(iarr_all_my,:,:,kb),cg%u(iarr_all_mz,:,:,kb),db)
         eib = cg%u(iarr_all_en,:,:,kb) - ekb
         eib = max(eib,smallei)
         do ifluid=1,flind%fluids
            csi2b(ifluid,:,:) = (flind%all_fluids(ifluid)%gam_1)*eib(ifluid,:,:)/db(ifluid,:,:)
         enddo
#endif /* !ISO */

         do j=1, cg%n_(ydim)
            do i=1, cg%n_(xdim)

               call get_gprofs(i,j)
               gprofs = tune_zeq_bnd * gprofs
               dprofs(:,1) = dbr(:,i,j)
               do ksub=1, nstot-1
                  factor = (2.0 + dzs*gprofs(ksub)/csi2b(:,i,j)) / (2.0 - dzs*gprofs(ksub)/csi2b(:,i,j))     !> \todo use hzeq_scheme here
                  dprofs(:,ksub+1) = factor * dprofs(:,ksub)
               enddo

               db(:,i,j) = 0.0
               lksub = 0
               do ksub=1, nstot
                  if (zs(ksub) > cg%zl(kk) .and. zs(ksub) < cg%zr(kk)) then
                     db(:,i,j) = db(:,i,j) + dprofs(:,ksub)/real(nsub)
                     lksub = ksub
                  endif
               enddo
               if (ib == 0_INT4) dprofs(:,lksub) = dprofs(:,lksub) * cg%u(iarr_all_dn,i,j,kk) / db(:,i,j)
               dbr(:,i,j) = dprofs(:,lksub)

               db(:,i,j)  = max(db(:,i,j), smalld)
#ifndef ISO
               eib(:,i,j) = csi2b(:,i,j)*db(:,i,j)/(flind%all_fluids(:)%gam_1)
               eib(:,i,j) = max(eib(:,i,j), smallei)
#endif /* !ISO */
            enddo
         enddo

         if (ib /= 0_INT4) then
            cg%u(iarr_all_dn,:,:,kk) = db(:,:,:)
            cg%u(iarr_all_mx,:,:,kk) = cg%u(iarr_all_mx,:,:,kb)
            cg%u(iarr_all_my,:,:,kk) = cg%u(iarr_all_my,:,:,kb)
            cg%u(iarr_all_mz,:,:,kk) = cg%u(iarr_all_mz,:,:,kb)
            if (diode) then
               if (side == HI) then
                  cg%u(iarr_all_mz,:,:,kk) = max(cg%u(iarr_all_mz,:,:,kk),0.0)
               else
                  cg%u(iarr_all_mz,:,:,kk) = min(cg%u(iarr_all_mz,:,:,kk),0.0)
               endif
            endif
#ifndef ISO
            ekb = ekin(cg%u(iarr_all_mx,:,:,kk),cg%u(iarr_all_my,:,:,kk),cg%u(iarr_all_mz,:,:,kk),db)
            cg%u(iarr_all_en,:,:,kk) = ekb + eib
#endif /* !ISO */
#ifdef COSM_RAYS
            cg%u(iarr_all_crs,:,:,kk) = smallecr
#endif /* COSM_RAYS */
         endif
      enddo

      deallocate(db, dbr,csi2b,zs,gprofs)
#ifndef ISO
      deallocate(ekb,eib)
#endif /* !ISO */

   end subroutine outh_bnd

#endif /* GRAV */
end module hydrostatic
