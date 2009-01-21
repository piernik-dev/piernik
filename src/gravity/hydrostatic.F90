! $Id$
#include "piernik.def"

module hydrostatic

! Written by M. Hanasz March-May 2006


#ifdef GRAV
  contains

    subroutine hydrostatic_zeq(iia,jja, d0, dprof)
      use start, only    : nsub,tune_zeq,csim2, col_dens, proc
      use constants, only : k_B, hydro_mass, small, pc
      use grid, only   : nx,ny,nz,dl,zdim,z,zl,zr, nzt
      use grid, only : nb,zmin,zmax

      use gravity, only  : grav_accel,grav_pot,gp_status
#ifndef ISO
      use arrays, only   : eprof
!      use start, only    : c_si
#endif /* ISO */
      implicit none
      real, intent(inout)              :: d0
      integer, intent(in)              :: iia, jja
      real, dimension(nz), intent(out) :: dprof

      integer nstot
      real, allocatable ::  zs(:), dprofs(:), gprofs(:), gpots(:)
      integer ksub, ksmid, k, ia, ja
      real dzs, factor

#if defined GALAXY || defined GALACTIC_DISK
      real cd, cdold,  dold, a, b
#endif /* GALAXY || GALACTIC_DISK */
      real dmid, ddmid, dcol_dens
      integer iter, itermx

      ia = min(nx,max(1, iia))
      ja = min(ny,max(1, jja))

      nstot=nsub*nzt

      allocate(zs(nstot), dprofs(nstot), gprofs(nstot), gpots(nstot))
#ifdef GALACTIC_DISK
      col_dens = d0
#endif /* GALACTIC_DISK */
      itermx = 20
      if(col_dens .gt. small) then
        dmid = 1.0
        ddmid = 1.0
        dcol_dens = 0.001*col_dens
        iter = 1
      elseif(d0 .gt. small) then
        dmid = d0
        iter = 0
      else
        if(proc .eq.0)  write(*,*) 'One of "d0" or "col_dens" must be .ne. 0'
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

#if defined GALAXY || defined GALACTIC_DISK
100   continue
#endif /* GALAXY || GALACTIC_DISK */

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


#if defined GALAXY || defined GALACTIC_DISK

      cd = sum(dprofs(:)) * dzs * pc

      if(col_dens .ne. 0.0) then
        dmid = d0
      else
        col_dens = cd
      endif

      if(abs(cd - col_dens) .gt. dcol_dens) then
        if(iter .eq. 1) then
          dold = dmid
          cdold = cd
          dmid = dmid+ddmid
        else
          a = (cd - cdold)/(dmid - dold)
          b = cd - a*dmid
          dold = dmid
          cdold = cd
          dmid = (col_dens - b)/a
        endif
        d0 = dmid
!       if(proc .eq.0)  write(*,888) d0, cd ,iter
        iter = iter+1
        if (iter .gt. itermx) stop
        goto 100
      endif

!      if(proc .eq.0)  write(*,888) d0, cd ,iter

#ifndef ISO
!      eprof(:) = c_si**2/(gamma-1.0) * dprof(:)	! fluid number has to be given to use gamma here
#endif /* ISO */

!888  format('Midplane density =',f10.4,2x,'Column density =', e10.4,2x,'iter=',i4 )
#endif /* GALAXY || GALACTIC_DISK */

    deallocate(zs,dprofs,gprofs,gpots)


    return


      end subroutine hydrostatic_zeq

!--------------------------------------------------------------------------
#ifdef NOT_WORKING
    subroutine hydro_thermal_zeq(ia, ja, hfl, T0, dprof,eprof,tprof,bprof)
      use constants, only : k_B, hydro_mass, small
      use grid, only : nb,zmin,zmax
      use start, only    : nsub,tune_zeq,proc,col_dens, &
         G_sup1, G_uv1, gamma, alpha
      use arrays, only   : nx,ny,nz,dl,zdim,z,zl,zr, nzt
      use thermal, only  : d_temp_dz, cool
     use gravity, only  : grav_accel

      implicit none
      integer                          :: ia, ja, hfl
      real, intent(in)                 :: T0
      real, dimension(nz), intent(out) :: dprof,eprof,tprof,bprof
      real gprof(nz)

      integer nstot
      real, allocatable ::  zs(:), gprofs(:), dprofs(:), eprofs(:), tprofs(:), cprofs(:), bprofs(:), cfuncs(:)
      integer ksub, ksmid, k
      real dzs, gravz, tmp, dtmpdz,zz(1),gg(1)

      ia = max(1, ia)
      ia = min(nx,ia)
      ja = max(1, ja)
      ja = min(ny,ja)

      nstot=nsub*nzt
      allocate(zs(nstot), gprofs(nstot), dprofs(nstot), eprofs(nstot), tprofs(nstot), cprofs(nstot), bprofs(nstot),cfuncs(nstot))

          call grav_accel('zsweep',ia, ja, z, nz, gprof)
!         write(*,*) gprof

      dzs = (zmax-zmin)/real(nstot-2*nb*nsub)

! Search for the index ksmid of the cell, whose center lies just below the midplane

      do ksub=1, nstot
        zs(ksub) = zmin-nb*dl(zdim) + dzs/2 + (ksub-1)*dzs
        if(zs(ksub) .lt. 0.0) ksmid = ksub       ! the midplane is in between
      enddo                                  ! ksmid and ksmid+1

! Integration of the hydro-thermal equilibrium up to the top
! of the subdivided z-grid

      if(ksmid .lt. nstot) then
        tprofs(ksmid+1) = T0
        do ksub=ksmid+1, nstot-1

          zz(1) = zs(ksub)
          call grav_accel('zsweep',ia, ja, zz, 1, gg)
          gravz = tune_zeq*gg(1)

          tmp =  tprofs(ksub)
          call d_temp_dz(gravz,tmp,dtmpdz)
          tmp = tmp + dtmpdz*dzs/2.

          zz(1) =zs(ksub) + dzs/2.
          call grav_accel('zsweep',ia, ja, zz, 1, gg)
          gravz = tune_zeq*gg(1)

          call d_temp_dz(gravz,tmp,dtmpdz)
          tprofs(ksub+1) = tprofs(ksub) + dtmpdz*dzs

!         write(*,"(i8,5(1x,e10.4))") ksub, zs(ksub), gravz,dtmpdz, tprofs(ksub)
        enddo
      endif

! Integration of the hydro-thermal equilibrium down to the bottom of the z-subgrid

      if(ksmid .gt. 1) then
        tprofs(ksmid) = T0
        do ksub=ksmid, 2, -1

          zz(1) = zs(ksub)
          call grav_accel('zsweep',ia, ja, zz, 1, gg)
          gravz = tune_zeq*gg(1)

          tmp =  tprofs(ksub)
          call d_temp_dz(gravz,tmp,dtmpdz)
          tmp = tmp - dtmpdz*dzs/2.

          zz(1) =zs(ksub) - dzs/2.
          call grav_accel('zsweep',ia, ja, zz, 1, gg)
          gravz = tune_zeq*gg(1)

          call d_temp_dz(gravz,tmp,dtmpdz)
          tprofs(ksub-1) = tprofs(ksub) - dtmpdz*dzs

!         write(*,"(i8,5(1x,e10.4))") ksub, zs(ksub), gravz,dtmpdz, tprofs(ksub)
        enddo
      endif

! Computation of all the relavant quantities on the subgrid

      call cool(nstot, tprofs, cfuncs)
      dprofs(:) = G_sup1/(cfuncs(:) - G_uv1)
      cprofs(:) = sqrt(tprofs(:)*k_B / hydro_mass)
      eprofs(:) = cprofs(:)**2/(gamma(hfl)-1.0) * dprofs(:)
      bprofs(:) = sqrt(2.*alpha*dprofs(:)*cprofs(:)**2)

! Remapping of the subgrid quantities onto a basic z-grid

      dprof(:) =0.0
      eprof(:) =0.0
      tprof(:) =0.0
!     cprof(:) =0.0
      bprof(:) =0.0
      do k=1,nz
        do ksub=1, nstot
          if(zs(ksub) .gt. zl(k) .and. zs(ksub) .lt. zr(k)) then
            dprof(k) = dprof(k) + dprofs(ksub)/real(nsub)
            eprof(k) = eprof(k) + eprofs(ksub)/real(nsub)
            tprof(k) = tprof(k) + tprofs(ksub)/real(nsub)
!           cprof(k) = cprof(k) + cprofs(ksub)/real(nsub)
            bprof(k) = bprof(k) + bprofs(ksub)/real(nsub)

          endif
        enddo
      enddo

!      do k=1, nz
!       write(*,"(i8,10(1x,e10.4))") k, z(k), tprof(k),  dprof(k), eprof(k), bprof(k)
!      enddo
!      stop

    deallocate(zs,gprofs,tprofs,dprofs,cfuncs)

    return


      end subroutine hydro_thermal_zeq
#endif /* NOT_WORKING */
#endif /* GRAV */
end module hydrostatic
