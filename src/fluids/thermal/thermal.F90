#include "piernik.h"
module thermal
! pulled by THERM
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

use constants,      only: cbuff_len

implicit none
 private
 public ::  maxdeint, cool_heat, init_thermal, cool_model, heat_model, A_heat, A_cool, alpha_cool, thermal_active,  cfl_coolheat, beta_cool

   real    eint_src_max, eint_src_min, dt_cool, dt_heat
   integer loc_dt_cool(3), loc_dt_heat(3)
   logical                  :: thermal_active
   character(len=cbuff_len) :: cool_model, heat_model
   real                     :: A_cool, A_heat, alpha_cool, beta_cool, cfl_coolheat

contains

!--------------------------------------------------------------------------

#ifdef THERM
subroutine init_thermal

      use cg_list_global, only: all_cg
      use constants,      only: PIERNIK_INIT_MPI, cbuff_len, ndims
      use dataio_pub,     only: nh    ! QA_WARN required for diff_nml
      use dataio_pub,     only: printinfo, warn, die, code_progress
      use mpisetup,       only: ibuff, rbuff, cbuff, master, slave, lbuff, piernik_MPI_Bcast

      implicit none
      namelist /THERMAL/ thermal_active, cool_model, heat_model, A_cool, A_heat, alpha_cool, beta_cool, cfl_coolheat

      if (code_progress < PIERNIK_INIT_MPI) call die("[thermal:init_thermal] mpi not initialized.")

#ifdef VERBOSE
      if (master) call printinfo("[thermal:init_thermal] Commencing thermal module initialization")
#endif /* VERBOSE */

      thermal_active = .True.
      cool_model     = 'power_law'
      heat_model     = 'beta_coef'
      A_cool         = 1.0
      A_heat         = 1.0
      alpha_cool     = 1.0
      beta_cool      = 1.0
      cfl_coolheat   = 0.1

     if (master) then

         if (.not.nh%initialized) call nh%init()
         open(newunit=nh%lun, file=nh%tmp1, status="unknown")
         write(nh%lun,nml=THERMAL)
         close(nh%lun)
         open(newunit=nh%lun, file=nh%par_file)
         nh%errstr=''
         read(unit=nh%lun, nml=THERMAL, iostat=nh%ierrh, iomsg=nh%errstr)
         close(nh%lun)
         call nh%namelist_errh(nh%ierrh, "THERMAL")
         read(nh%cmdl_nml,nml=THERMAL, iostat=nh%ierrh)
         call nh%namelist_errh(nh%ierrh, "THERMAL", .true.)
         open(newunit=nh%lun, file=nh%tmp2, status="unknown")
         write(nh%lun,nml=THERMAL)
         close(nh%lun)
         call nh%compare_namelist()

         rbuff(1)   = A_cool
         rbuff(2)   = alpha_cool
         rbuff(3)   = A_heat
         rbuff(4)   = beta_cool
         rbuff(5)   = cfl_coolheat

         lbuff(1)   = thermal_active

         cbuff(1)   = cool_model
         cbuff(2)   = heat_model

      endif

      call piernik_MPI_Bcast(ibuff)
      call piernik_MPI_Bcast(lbuff)
      call piernik_MPI_Bcast(rbuff)
      call piernik_MPI_Bcast(cbuff, cbuff_len)

      if (slave) then

         cool_model     = cbuff(1)
         heat_model     = cbuff(2)


         thermal_active = lbuff(1)

         A_cool         = rbuff(1)
         alpha_cool     = rbuff(2)
         A_heat         = rbuff(3)
         beta_cool      = rbuff(4)
         cfl_coolheat   = rbuff(5)


      endif

   end subroutine init_thermal

  subroutine cool_heat(gamma, n, dens, eint, esrc)
  use units,      only: kboltz, mH

     implicit none
!    character, intent(in) :: sweep*6
     integer,intent(in)    :: n
     real, dimension(n),intent(in) :: dens, eint
     real, dimension(n) :: cfunc, hfunc
     real, dimension(n),intent(out) :: esrc
     real, dimension(n) :: temp
!    real, dimension(n) :: a,b,c
!    integer k
     real, intent(in)   :: gamma

      temp(:) = 0.0
      cfunc(:) = 0.0
      hfunc(:) = 0.0
      esrc(:) = 0.0

      temp = (gamma-1) * mH / kboltz * eint / dens

       call cool(n, temp, cfunc)

       call heat(n, dens, hfunc)
       esrc =  dens**2*cfunc + hfunc
!       write(*,*) "esrc", esrc
!      esrc = MIN(esrc, esrc_upper_lim * eint)
!      esrc = MAX(esrc, esrc_lower_lim * eint)

  end subroutine cool_heat

      subroutine maxdeint(cg, max_deint)
      use constants,        only: xdim, ydim, zdim, LO, HI
      use grid_cont,        only: grid_container
      use fluidtypes,       only: component_fluid
      use func,             only: emag, ekin
      use fluidindex,       only: flind

      implicit none
      integer                                                :: i, j, ifl
      type(grid_container), pointer, intent(in)              :: cg
      class(component_fluid), pointer                        :: pfl
      real, dimension(cg%ks:cg%ke)                           :: int_ener, kin_ener, mag_ener, eint_src
      real, dimension(cg%is:cg%ie, cg%js:cg%je, cg%ks:cg%ke) :: deint
      real, intent(out)                                      :: max_deint

do ifl = 1, flind%fluids
     pfl => flind%all_fluids(ifl)%fl
  do j=cg%js,cg%je
      do i=cg%is,cg%ie
            if (pfl%has_energy) then
               kin_ener = ekin(cg%u(pfl%imx,i,j,cg%ks:cg%ke), cg%u(pfl%imy,i,j,cg%ks:cg%ke), cg%u(pfl%imz,i,j,cg%ks:cg%ke), cg%u(pfl%idn,i,j,cg%ks:cg%ke))
               if (pfl%is_magnetized) then
                  mag_ener = emag(cg%b(xdim,i,j,cg%ks:cg%ke), cg%b(ydim,i,j,cg%ks:cg%ke), cg%b(zdim,i,j,cg%ks:cg%ke))
                  int_ener = cg%u(pfl%ien,i,j,cg%ks:cg%ke) - kin_ener - mag_ener
               else
               int_ener = cg%u(pfl%ien,i,j,cg%ks:cg%ke) - kin_ener
               endif
            endif
            call cool_heat(pfl%gam, cg%nzb, cg%u(pfl%idn,i,j,cg%ks:cg%ke), int_ener, eint_src)
            deint(i,j,:) = eint_src/int_ener
      enddo
   enddo
enddo
!print *,'cg%nxb, cg%nyb, cg%nzb', cg%nxb, cg%nyb, cg%nzb
!print *,'cg%is,cg%ie', cg%is,cg%ie
!print *,'cg%js,cg%je', cg%js,cg%je
!print *,'cg%ks,cg%ke', cg%ks,cg%ke
!stop
    max_deint = MAXVAL(ABS(deint))

!    print *, 'max_deint =', max_deint
      end subroutine maxdeint

subroutine cool(n, temp, coolf)

    implicit none
    integer,           intent(in)  :: n
    real, dimension(n)             :: temp
    real, dimension(n),intent(out) :: coolf


      select case (cool_model)
        case ('power_law')
          coolf = -A_cool * temp**(alpha_cool)
      case ('null')
          return
        case default
          write(*,*) 'Cool model: ',cool_model,' not implemented'
      end select

end subroutine cool

subroutine heat(n, dens, heatf)

    implicit none
    integer,           intent(in)  :: n
    real, dimension(n),intent(in)  :: dens
    real, dimension(n),intent(out) :: heatf

      select case (heat_model)
        case ('beta_coef')
          heatf =  A_heat * dens**beta_cool  !(trzeba dodac czytanie parametrow grzana z problem.par)
        case ('null')
          return
        case default
          write(*,*) 'Heat model: ',heat_model,' not implemented'
      end select

end subroutine heat

#endif /*THERM*/
!--------------------------------------------------------------------------
#ifdef IMPROPER
!--------------------------------------------------------------------------

subroutine timestep_coolheat

    implicit none


    real eint_src_min_sweep,eint_src_max_sweep
    real, dimension(nz) :: eint, dens, temp, eint_src
    integer loc_dt_cool3(1), loc_dt_heat3(1)
    integer i,j

    real dt_coolheat_all, dt_coolheat_proc


    eint_src_max = 0.0
    eint_src_min = 0.0

    do j=nb+1,nb+nyb
      do i=nb+1,nb+nxb
        dens = u(1,i,j,:)
        eint = u(5,i,j,:)-sum(u(2:4,i,j,:)**2,1)/u(1,i,j,:)/2-sum(b(:,i,j,:)**2,1)/2

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

  end subroutine timestep_cool_heat


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

subroutine cool(n, temp, coolf)

    implicit none
    integer,           intent(in)  :: n
    real, dimension(n)		   :: temp
    real, dimension(n),intent(out) :: coolf


      select case (cool_model)
        case ('power_law')
          coolf = -A_cool * temp**(alpha_cool)
	case ('null')
          return
        case default
          write(*,*) 'Cool model: ',cool_model,' not implemented'
      end select

end subroutine cool



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

!----------------------------------------------------------------------------

#endif /*IMPROPER*/
end module thermal

