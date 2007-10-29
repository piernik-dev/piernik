module hydrostatic

! Written by M. Hanasz March-May 2006

  use start
  use arrays
  use grid
#ifdef GRAV
  use gravity
#endif
  use thermal

#ifdef GRAV
  contains

    subroutine hydrostatic_zeq(iia,jja, d0, dprof)
      
      implicit none
      real, intent(in)                 :: d0
      integer, intent(in)              :: iia, jja
      real, dimension(nz), intent(out) :: dprof
      
      integer nstot
      real, allocatable ::  zs(:), dprofs(:), gprofs(:), gpots(:)
      integer ks, ksmid, k, ia, ja
      real z1,z2, dzs, factor

      ia = min(nx,max(1, iia))
      ja = min(ny,max(1, jja))
            
      nstot=nsub*nzt
      
      allocate(zs(nstot), dprofs(nstot), gprofs(nstot), gpots(nstot))
          
      dzs = (zmax-zmin)/real(nstot-2*nb*nsub)
      
      do ks=1, nstot
        zs(ks) = zmin-nb*dl(zdim) + dzs/2 + (ks-1)*dzs  !
        if(zs(ks) .lt. 0.0) ksmid = ks       ! the midplane is in between 
      enddo                                  ! ksmid and ksmid+1
      
      if(gp_status .eq. 'undefined') then
        call grav_accel('zsweep',ia, ja, zs, nstot, gprofs)
      else
        call grav_pot('zsweep', ia,ja, zs, nstot, gpots,gp_status)
        gprofs(1:nstot-1) = (gpots(1:nstot-1) - gpots(2:nstot))/dzs
      endif

      gprofs = tune_zeq*gprofs
      if(ksmid .lt. nstot) then 
        dprofs(ksmid+1) = d0
        do ks=ksmid+1, nstot-1
          factor = (1.0 + 0.5*dzs*gprofs(ks)/csim2)  &
                  /(1.0 - 0.5*dzs*gprofs(ks)/csim2)     
          dprofs(ks+1) = factor * dprofs(ks)  
        enddo
      endif
        
      if(ksmid .gt. 1) then
        dprofs(ksmid) = d0
        do ks=ksmid, 2, -1
          factor = (1.0 - 0.5*dzs*gprofs(ks)/csim2)  &
                  /(1.0 + 0.5*dzs*gprofs(ks)/csim2)     
          dprofs(ks-1) = factor * dprofs(ks)  
        enddo
      endif
            
      dprof(:) =0.0      
      do k=1,nz
        do ks=1, nstot
          if(zs(ks) .gt. zl(k) .and. zs(ks) .lt. zr(k)) then
            dprof(k) = dprof(k) + dprofs(ks)/real(nsub)
          endif
        enddo
      enddo

!!      if(ia.eq.5 .and.ja.eq.5) then
!        write(*,*)
!        write(*,*) 'proc=',proc
!        do k=1,nz
!         write(*,999) k, z(k), dprof(k)        
!        enddo
!       stop
!!      endif


    deallocate(zs,dprofs,gprofs,gpots)

    return


999 format((1x,i4),10(1x,e10.4))

      end subroutine hydrostatic_zeq

!--------------------------------------------------------------------------

    subroutine hydro_thermal_zeq(ia, ja, T0, dprof,eprof,tprof,bprof)
      
      implicit none
      integer                          :: ia, ja
      real, intent(in)                 :: T0
      real, dimension(nz), intent(out) :: dprof,eprof,tprof,bprof
      real gprof(nz)
      
      integer nstot
      real, allocatable ::  zs(:), gprofs(:), dprofs(:), eprofs(:), tprofs(:), cprofs(:), bprofs(:), cfuncs(:)
      integer ks, ksmid, k
      real z1,z2, dzs, factor, gravz, tmp, dtmpdz,zz(1),gg(1)
      real csim2  
      real tune

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

      do ks=1, nstot
        zs(ks) = zmin-nb*dl(zdim) + dzs/2 + (ks-1)*dzs  
        if(zs(ks) .lt. 0.0) ksmid = ks       ! the midplane is in between 
      enddo                                  ! ksmid and ksmid+1
      
! Integration of the hydro-thermal equilibrium up to the top 
! of the subdivided z-grid
      
      if(ksmid .lt. nstot) then 
        tprofs(ksmid+1) = T0
        do ks=ksmid+1, nstot-1

          zz(1) = zs(ks)
          call grav_accel('zsweep',ia, ja, zz, 1, gg)
          gravz = tune_zeq*gg(1)

          tmp =  tprofs(ks)
          call d_temp_dz(gravz,tmp,dtmpdz)
          tmp = tmp + dtmpdz*dzs/2.

          zz(1) =zs(ks) + dzs/2.
          call grav_accel('zsweep',ia, ja, zz, 1, gg)
          gravz = tune_zeq*gg(1)
          
          call d_temp_dz(gravz,tmp,dtmpdz)
          tprofs(ks+1) = tprofs(ks) + dtmpdz*dzs

!         write(*,"(i8,5(1x,e10.4))") ks, zs(ks), gravz,dtmpdz, tprofs(ks) 
        enddo
      endif
        
! Integration of the hydro-thermal equilibrium down to the bottom of the z-subgrid
      
      if(ksmid .gt. 1) then
        tprofs(ksmid) = T0
        do ks=ksmid, 2, -1
        
          zz(1) = zs(ks)
          call grav_accel('zsweep',ia, ja, zz, 1, gg)
          gravz = tune_zeq*gg(1)

          tmp =  tprofs(ks)
          call d_temp_dz(gravz,tmp,dtmpdz)
          tmp = tmp - dtmpdz*dzs/2.

          zz(1) =zs(ks) - dzs/2.
          call grav_accel('zsweep',ia, ja, zz, 1, gg)
          gravz = tune_zeq*gg(1)
          
          call d_temp_dz(gravz,tmp,dtmpdz)
          tprofs(ks-1) = tprofs(ks) - dtmpdz*dzs

!         write(*,"(i8,5(1x,e10.4))") ks, zs(ks), gravz,dtmpdz, tprofs(ks)
        enddo
      endif

! Computation of all the relavant quantities on the subgrid

      call cool(nstot, tprofs, cfuncs)
      dprofs(:) = G_sup1/(cfuncs(:) - G_uv1)
      cprofs(:) = sqrt(tprofs(:)*k_B / hydro_mass)
      eprofs(:) = cprofs(:)**2/(gamma-1.0) * dprofs(:)
      bprofs(:) = sqrt(2.*alpha*dprofs(:)*cprofs(:)**2) 

! Remapping of the subgrid quantities onto a basic z-grid 
            
      dprof(:) =0.0      
      eprof(:) =0.0      
      tprof(:) =0.0      
!     cprof(:) =0.0      
      bprof(:) =0.0      
      do k=1,nz
        do ks=1, nstot
          if(zs(ks) .gt. zl(k) .and. zs(ks) .lt. zr(k)) then
            dprof(k) = dprof(k) + dprofs(ks)/real(nsub)
            eprof(k) = eprof(k) + eprofs(ks)/real(nsub)
            tprof(k) = tprof(k) + tprofs(ks)/real(nsub)
!           cprof(k) = cprof(k) + cprofs(ks)/real(nsub)
            bprof(k) = bprof(k) + bprofs(ks)/real(nsub)

          endif
        enddo
      enddo

!      do k=1, nz
!       write(*,"(i8,10(1x,e10.4))") k, z(k), tprof(k),  dprof(k), eprof(k), bprof(k)
!      enddo
!      stop

    deallocate(zs,gprofs,tprofs,dprofs,cfuncs)

    return


999 format((1x,i4),10(1x,e10.4))

      end subroutine hydro_thermal_zeq

#endif GRAV
end module hydrostatic
