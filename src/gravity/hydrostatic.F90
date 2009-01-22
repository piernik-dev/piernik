! $Id$
#include "piernik.def"

module hydrostatic

! Written by M. Hanasz March-May 2006


#ifdef GRAV
  contains

    subroutine hydrostatic_zeq(iia,jja, d0, c_si, alpha, dprof)
      use constants
      use mpi_setup, only : proc
      use grid, only : nx,ny,nz,dl,zdim,z,zl,zr,nzt,nb,zmin,zmax

      use gravity, only  : grav_accel,grav_pot,gp_status,nsub,tune_zeq
#ifndef ISO
      use arrays, only   : eprof
#endif /* ISO */
      implicit none
      real, intent(inout)              :: d0
      real, intent(in)                 :: c_si, alpha
      integer, intent(in)              :: iia, jja
      real, dimension(nz), intent(out) :: dprof

      integer nstot
      real, allocatable ::  zs(:), dprofs(:), gprofs(:), gpots(:)
      integer ksub, ksmid, k, ia, ja
      real dzs, factor

      real dmid, ddmid, csim2
      integer iter, itermx

      ia = min(nx,max(1, iia))
      ja = min(ny,max(1, jja))

      nstot=nsub*nzt
   
      csim2 = c_si**2 * (1.0 + alpha)

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
