! $Id$
#include "piernik.def"
module thermal

! Written by R. Kosinki
! Adopted for this code by R. Kosinski & M. Hanasz, May 2006
!
! Calculates radiative energy loss. Based on Zeus3D pdvcool subroutine
!
!  PURPOSE:  This routine finds new energy in the source step
!       by solving the implicit eqn:
!
!       de / dt = - d**2 * COOL(e) + HEAT
!  where COOL is an empirical cooling function of e only, and HEAT is
!  an empirical heating function.

  use start, only     : nb, cool_model,heat_model,esrc_upper_lim,esrc_lower_lim, &
      gamma, alpha, G_uv1, G_sup1, L_C, h_coolheat_profile,C_heatcond, &
      K_heatcond, cfl_coolheat, dt, n_coolheat_profile, dt_coolheat

  use arrays, only    : dl, nzb, is,ie,js,je,ks,ke,xdim,ydim,zdim,nx,ny,nz, &
      u,b,coolheat_profile, z, idna, imxa, imya, imza, iena
!  use grid
  use constants, only : small, k_B, hydro_mass, cmps2, pc, m_H, chcf
  use mpi_setup

  real    eint_src_max, eint_src_min, dt_cool, dt_heat
  integer loc_dt_cool(3), loc_dt_heat(3)


contains

!--------------------------------------------------------------------------

  subroutine cool_heat (sweep, i1, i2, n, dens, eint,  esrc)

    implicit none
    character, intent(in) :: sweep*6
    integer,intent(in)    :: n, i1, i2
    real, dimension(n),intent(in) :: dens, eint
    real, dimension(n) :: cfunc, hfunc
    real, dimension(n),intent(out) :: esrc
    real, dimension(n) :: temp

      temp(:)  = 0.0
      cfunc(:) = 0.0
      hfunc(:) = 0.0
      esrc(:)  = 0.0

      temp = (gamma-1) * hydro_mass / k_B * eint / dens

      call cool(n, temp, cfunc)
      cfunc = cfunc / chcf

      call heat(n, dens, hfunc)
      hfunc = hfunc / chcf

      esrc =  hfunc - dens**2 * cfunc + small

      esrc = MIN(esrc, esrc_upper_lim * eint)
      esrc = MAX(esrc, esrc_lower_lim * eint)

      select case (sweep)
        case('xsweep')
          esrc = esrc * coolheat_profile(i2)
        case('ysweep')
          esrc = esrc * coolheat_profile(i1)
        case('zsweep')
          esrc = esrc * coolheat_profile
      end select



  end subroutine cool_heat

!--------------------------------------------------------------------------

  subroutine cool(n, temp, coolf)

    implicit none
    integer,           intent(in)  :: n
    real, dimension(n)             :: temp
    real, dimension(n),intent(out) :: coolf
    real con, skl

!    L_C = -1.23848d-55
    con = -0.012533
    skl =  2.8347d-10

      select case (cool_model)
        case ('c1')
          temp=max(temp,101.0)
          coolf = con * skl * (L_C * temp**10.73 ) &
                + 1.606d-24 * Exp(-40.d0/(temp-1.e2)**0.6)
        case ('null')
          return
        case default
          write(*,*) 'Cool model: ',cool_model,' not implemented'
      end select

  end subroutine cool

!--------------------------------------------------------------------------

  subroutine dcool (temp, dcoolf)

    implicit none

    real ::  temp
    real, parameter :: eps = 0.0001

    real, dimension(2) ::  tmp, cf
    real dtmp, dcf
    real  dcoolf

    tmp(1) = (1.0 - eps) * temp
    tmp(2) = (1.0 + eps) * temp

    call cool(2, tmp, cf)

    dcf   =  cf(2)  - cf(1)
    dtmp  = tmp(2) - tmp(1)

    if((abs(dcf).le.1.e-12).and.(abs(dtmp).le.1.e-12)) then
      dcoolf=0.0
      return
    endif

    dcoolf = dcf/dtmp

  end  subroutine dcool


!--------------------------------------------------------------------------

  subroutine heat(n, dens, heatf)

    implicit none
    integer,           intent(in)  :: n
    real, dimension(n),intent(in)  :: dens
    real, dimension(n),intent(out) :: heatf

      select case (heat_model)
        case ('sup1')
          heatf = dens**2 * G_uv1 + dens * G_sup1
        case ('null')
          return
        case default
          write(*,*) 'Heat model: ',heat_model,' not implemented'
      end select


  end subroutine heat

!--------------------------------------------------------------------------

    subroutine d_temp_dz(gz,temp,dtempdz)

      implicit NONE

        REAL gzpt
        REAL gz, temp, tmp(1), dtempdz
        REAL lambda(1), dlambdadt

        gzpt = m_H / ( k_B * (1.0+alpha) ) * gz/cmps2

        tmp(1) = temp
        call cool (1, tmp,  lambda)
        call dcool(   tmp(1), dlambdadt)

        select case (heat_model)
          case ('sup1')
            dtempdz =  -gzpt * (lambda(1) - G_uv1) / &
                              (temp * dlambdadt - lambda(1) + G_uv1)
            dtempdz = dtempdz*pc
          case ('null')
            return
          case default
            write(*,*) 'Heat model: ',heat_model,' not implemented'
        end select

    end subroutine d_temp_dz

!--------------------------------------------------------------------------

  subroutine cool_heat_profile

    implicit none

        if(n_coolheat_profile .ne. 0) then
          coolheat_profile = 1./cosh((z/(h_coolheat_profile))**n_coolheat_profile)
        else
          coolheat_profile(:) = 1.
        endif

!       write(*,*) coolheat_profile

  end subroutine  cool_heat_profile

!--------------------------------------------------------------------------

  subroutine timestep_coolheat

    implicit none


    real eint_src_min_sweep,eint_src_max_sweep
    real, dimension(nz) :: eint, dens, eint_src
    integer loc_dt_cool3(1), loc_dt_heat3(1)
    integer i,j

    real dt_coolheat_all, dt_coolheat_proc


    eint_src_max = 0.0
    eint_src_min = 0.0

    do j=js,je
      do i=is,ie
        dens = u(idna,i,j,:)
        eint = u(iena,i,j,:)-sum(u(imxa:imza,i,j,:)**2,1)/u(idna,i,j,:)/2-sum(b(:,i,j,:)**2,1)/2

        call cool_heat ('zsweep', i, j, nz, dens, eint,  eint_src)

        eint_src_min_sweep = MINVAL(eint_src(nb+1:nb+nzb)/eint(nb+1:nb+nzb))
        eint_src_max_sweep = MAXVAL(eint_src(nb+1:nb+nzb)/eint(nb+1:nb+nzb))


        if(eint_src_min_sweep .lt. eint_src_min) then
          eint_src_min = eint_src_min_sweep
          loc_dt_cool(1:2) = (/i,j/)
          loc_dt_cool3 = MINLOC(eint_src(nb+1:nb+nzb) / eint(nb+1:nb+nzb)) + nb
          loc_dt_cool(3) = loc_dt_cool3(1)
        endif

        if(eint_src_max_sweep .gt. eint_src_max) then
          eint_src_max = eint_src_max_sweep
          loc_dt_heat(1:2) = (/i,j/)
          loc_dt_heat3 = MAXLOC(eint_src(nb+1:nb+nzb) / eint(nb+1:nb+nzb)) + nb
          loc_dt_heat(3) = loc_dt_heat3(1)
        endif

      end do
    end do

    dt_heat = cfl_coolheat*abs(1./(eint_src_max+small))
    dt_cool = cfl_coolheat*abs(1./(eint_src_min+small))

    dt_coolheat_proc = MIN ( dt_cool, dt_heat)


    call MPI_REDUCE(dt_coolheat_proc, dt_coolheat_all, 1, MPI_DOUBLE_PRECISION, MPI_MIN, 0, comm, ierr)
    call MPI_BCAST(dt_coolheat_all, 1, MPI_DOUBLE_PRECISION, 0, comm, ierr)

    dt_coolheat = dt_coolheat_all

  end subroutine timestep_coolheat

!----------------------------------------------------------------------------

  subroutine heat_conduction

    implicit none

    real, allocatable   :: temp(:,:,:),dtemp(:,:,:)
    allocate (temp(nx,ny,nz),dtemp(nx,ny,nz))

!    temp = (gamma-1) * hydro_mass / k_B * eint / dens
!    dtemp = (eoshift(temp,1) - temp)/dx(1,:)
!    flux_heat = -1.5e-5*K_heatcond * dtemp
!    eint=eint-(flux_heat-eoshift(flux_heat,shift=-1,boundary=big))/dx(1,:)*dt

! Temperature
    temp = (gamma-1) * hydro_mass / k_B * (  &
           u(5,:,:,:)  - 0.5*( u(2,:,:,:)*u(2,:,:,:) &
                              +u(3,:,:,:)*u(3,:,:,:) &
                              +u(4,:,:,:)*u(4,:,:,:))/u(1,:,:,:) &
                       - 0.5*( b(1,:,:,:)*b(1,:,:,:) &
                              +b(2,:,:,:)*b(2,:,:,:) &
                              +b(3,:,:,:)*b(3,:,:,:))  )



    dtemp = (cshift(temp,shift=1,dim=1) - temp)/dl(xdim)

    u(5,:,:,:) = u(5,:,:,:)  + C_heatcond*K_heatcond * ( &
                (temp-cshift(temp,shift=-1,dim=1)) )/dl(xdim)*dt


    dtemp = (cshift(temp,shift=1,dim=2) - temp)/dl(ydim)

    u(5,:,:,:) = u(5,:,:,:)  + C_heatcond*K_heatcond * ( &
                (temp-cshift(temp,shift=-1,dim=2)) )/dl(ydim)*dt


    dtemp = (cshift(temp,shift=1,dim=3) - temp)/dl(zdim)

    u(5,:,:,:) = u(5,:,:,:)  + C_heatcond*K_heatcond * ( &
                (temp-cshift(temp,shift=-1,dim=3)) )/dl(zdim)*dt

   deallocate (temp,dtemp)

  end subroutine heat_conduction

!----------------------------------------------------------------------------

end module thermal

