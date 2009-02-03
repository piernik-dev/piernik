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

module hydrostatic
#ifdef GRAV
  contains

    subroutine hydrostatic_zeq(iia,jja, d0, dprof)
      use constants
      use start, only : csim2
      use mpisetup, only : proc
      use grid, only : nx,ny,nz,dl,zdim,z,zl,zr,nzt,nb,zmin,zmax

      use gravity, only  : grav_accel,grav_pot,gp_status,nsub,tune_zeq
#ifndef ISO
      use arrays, only   : eprof
#endif /* ISO */
      implicit none
      real, intent(inout)              :: d0
      integer, intent(in)              :: iia, jja
      real, dimension(nz), intent(out) :: dprof

      integer nstot
      real, allocatable ::  zs(:), dprofs(:), gprofs(:), gpots(:)
      integer ksub, ksmid, k, ia, ja
      real dzs, factor

      real dmid, ddmid
      integer iter, itermx

      ia = min(nx,max(1, iia))
      ja = min(ny,max(1, jja))

      nstot=nsub*nzt
   
      allocate(zs(nstot), dprofs(nstot), gprofs(nstot), gpots(nstot))
      itermx = 20
      if(d0 .gt. small) then
        dmid = d0
        iter = 0
      else
        if(proc .eq.0)  write(*,*) '"d0" must be /= 0'
         stop
      endif


      dzs = (zmax-zmin)/real(nstot-2*nb*nsub)

      do ksub=1, nstot
        zs(ksub) = zmin-nb*dl(zdim) + dzs/2 + (ksub-1)*dzs  !
        if(zs(ksub) .lt. 0.0) ksmid = ksub      ! the midplane is in between
      enddo                                  ! ksmid and ksmid+1

      if(gp_status .eq. 'undefined') then
        call grav_accel('zsweep',ia, ja, zs, nstot, gprofs)
      else
        gp_status = 'hydrozeq'
        call grav_pot('zsweep', ia,ja, zs, nstot, gpots,gp_status,.true.)
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


    deallocate(zs,dprofs,gprofs,gpots)


    return


      end subroutine hydrostatic_zeq

#endif /* GRAV */
end module hydrostatic
