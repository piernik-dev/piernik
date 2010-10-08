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
!! \brief [DW] Module containing a subroutine that arranges %hydrostatic equilibrium in the vertical (z) direction
!<
module hydrostatic
#ifdef GRAV
  contains
!>
!! \brief Routine that arranges %hydrostatic equilibrium in the vertical (z) direction
!! \param iia integer, number of a column in the x direction
!! \param jja integer, number of a column in the y direction
!! \param d0 real, density value in the midplane
!! \param csim2 real, square value of sound speed
!! \param dprof array of reals, computed density distribution in vertical direction returned by the routine
!<
    subroutine hydrostatic_zeq(iia,jja, d0, csim2, dprof)

      use errh, only: die
      use constants, only: small
      use mpisetup, only: proc
      use grid, only: nx,ny,nz,dl,zdim,z,zl,zr,nzt,nb,zmin,zmax
      use gravity, only: grav_accel,gp_status,nsub,tune_zeq
      use arrays, only: gp
#ifndef ISO
      use arrays, only: eprof
#endif /* ISO */
      implicit none
      real, intent(inout)              :: d0
      integer, intent(in)              :: iia, jja
      real, dimension(nz), intent(out) :: dprof
      real, intent(in)                 :: csim2

      real, allocatable, dimension(:)  :: zs, dprofs, gprofs, gpots
      integer :: ksub, ksmid, k, ia, ja, nstot, iter, itermx
      real    :: dzs, factor, dmid

!      csim2 = cs_iso2*(1.0+alpha)

      ia = min(nx,max(1, iia))
      ja = min(ny,max(1, jja))

      nstot=nsub*nzt

      allocate(zs(nstot), dprofs(nstot), gprofs(nstot), gpots(nstot))
      itermx = 20
      if(d0 .gt. small) then
         dmid = d0
         iter = 0
      else
         call die("[hydrostatic:hydrostatic_zeq] d0 must be /= 0")
         dmid = 0. ! just for suppressing compiler warning
      endif

      dzs = (zmax-zmin)/real(nstot-2*nb*nsub)

      ksmid = 0
      do ksub=1, nstot
        zs(ksub) = zmin-nb*dl(zdim) + dzs/2 + (ksub-1)*dzs  !
        if(zs(ksub) .lt. 0.0) ksmid = ksub      ! the midplane is in between
      enddo                                  ! ksmid and ksmid+1
      if (ksmid == 0) call die("[hydrostatic:hydrostatic_zeq] ksmid not set")

      if(gp_status .eq. 'undefined') then
        call grav_accel('zsweep',ia, ja, zs, nstot, gprofs)
      else
        k = 1; gpots(:) = 0.0
        do ksub=1, nstot
           if(zs(ksub) >= z(min(k+1,nz))) k = k + 1
           if(zs(ksub) >= z(min(k,nz-1)) .and. zs(ksub) < z(min(k+1,nz))) then
               gpots(ksub) = gp(iia,jja,k) + (zs(ksub) - z(k)) * &
                (gp(iia,jja,min(k+1,nz)) - gp(iia,jja,k)) / (z(min(k+1,nz)) - z(min(k,nz-1)))
           endif
        enddo
!        call grav_pot('zsweep', ia,ja, zs, nstot, gpots,gp_status,.true.)
        gprofs(1:nstot-1) = (gpots(1:nstot-1) - gpots(2:nstot))/dzs
      endif
      gprofs = tune_zeq*gprofs

      if(ksmid .lt. nstot) then
        dprofs(ksmid+1) = dmid
        do ksub=ksmid+1, nstot-1
          factor = (1.0 + 0.5*dzs*gprofs(ksub)/csim2)  &
                  /(1.0 - 0.5*dzs*gprofs(ksub)/csim2)
          dprofs(ksub+1) = factor * dprofs(ksub)
        enddo
      endif

      if(ksmid .gt. 1) then
        dprofs(ksmid) = dmid
        do ksub=ksmid, 2, -1
          factor = (1.0 - 0.5*dzs*gprofs(ksub)/csim2)  &
                  /(1.0 + 0.5*dzs*gprofs(ksub)/csim2)
          dprofs(ksub-1) = factor * dprofs(ksub)
        enddo
      endif

      dprof(:) =0.0
      do k=1,nz
        do ksub=1, nstot
          if(zs(ksub) .gt. zl(k) .and. zs(ksub) .lt. zr(k)) then
            dprof(k) = dprof(k) + dprofs(ksub)/real(nsub)
          endif
        enddo
      enddo

    if(allocated(zs)) deallocate(zs)
    if(allocated(dprofs)) deallocate(dprofs)
    if(allocated(gprofs)) deallocate(gprofs)
    if(allocated(gpots)) deallocate(gpots)

    return


      end subroutine hydrostatic_zeq

#endif /* GRAV */
end module hydrostatic
