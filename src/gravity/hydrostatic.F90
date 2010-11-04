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
#include "piernik.def"

!>
!! \brief [DW] Module containing a subroutine that arranges %hydrostatic equilibrium in the vertical (z) direction
!<
module hydrostatic

   implicit none

   private
#ifdef GRAV
   public :: hydrostatic_zeq_coldens, hydrostatic_zeq_densmid, gprofs, nstot, zs, dzs
#endif /* GRAV */

   real, allocatable, dimension(:), save :: zs, gprofs
   real,    save :: dzs
   integer, save :: nstot
   logical, save :: hstarted = .false.
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

      if (.not.hstarted) call die("[hydrostatic:hydrostatic_main] procedure used before initializing with start_hydrostatic")
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
            dprofs(ksub+1) = dprofs(ksub)*(0.5*(gprofs(ksub)+gprofs(ksub+1))*dzs + 1.0)
         enddo
      endif

      if (ksmid .gt. 1) then
         dprofs(ksmid) = 1.0
         do ksub=ksmid, 2, -1
            dprofs(ksub-1) = dprofs(ksub)*(0.5*(gprofs(ksub)+gprofs(ksub-1))*dzs + 1.0)
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
      use gravity, only: get_gprofs
      implicit none
      integer, intent(in) :: iia, jja
      real,    intent(in) :: coldens, csim2
      real :: sdprof

      call start_hydrostatic
      call get_gprofs(iia,jja)
      gprofs = gprofs / csim2
      call hydrostatic_main
      sdprof = sum(dprof)
      dprof = dprof * coldens / sdprof
      call finish_hydrostatic

   end subroutine hydrostatic_zeq_coldens

   subroutine hydrostatic_zeq_densmid(iia,jja,d0,csim2)
      use arrays,     only: dprof
      use constants,  only: small
      use dataio_pub, only: die
      use gravity,    only: get_gprofs
      implicit none
      integer, intent(in) :: iia, jja
      real,    intent(in) :: d0, csim2

      if (d0 .le. small) then
         call die("[hydrostatic:hydrostatic_zeq_densmid] d0 must be /= 0")
      endif
#ifndef NEW_HYDROSTATIC
      dmid = d0
#endif /* !NEW_HYDROSTATIC */

      call start_hydrostatic
      call get_gprofs(iia,jja)
      gprofs = gprofs / csim2
      call hydrostatic_main
      dprof = dprof * d0
      call finish_hydrostatic

   end subroutine hydrostatic_zeq_densmid

   subroutine start_hydrostatic
      use dataio_pub, only: die
      use gravity,    only: get_gprofs, gprofs_target, nsub
      use grid,       only: zmin, zmax, zdim, nzt, dl, nb
      implicit none
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
      hstarted = .true.
   end subroutine start_hydrostatic

   subroutine finish_hydrostatic
      implicit none
      if (allocated(zs))     deallocate(zs)
      if (allocated(gprofs)) deallocate(gprofs)
      hstarted = .false.
   end subroutine finish_hydrostatic
#endif /* GRAV */
end module hydrostatic
