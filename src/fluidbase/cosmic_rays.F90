! $Id$
#include "piernik.def"
#define SQR(var) var*var
#define SUM_SQR(x,y,z) ( SQR(x)+SQR(y)+SQR(z) )
#define RNG 2:n-1
#ifdef IONIZED
#define NUMBION IONIZED
#else /* IONIZED */
#define NUMBION 0
#endif /* IONIZED */

module cosmic_rays
   integer,dimension(COSM_RAYS)  :: iecr
   real  cr_active, gamma_cr, cr_eff, beta_cr, K_cr_paral, K_cr_perp, &
         cfl_cr, amp_cr, smallecr

   contains

   subroutine add_cr_index(nusf)
      implicit none
      integer :: iter, nusf
      do iter=1,COSM_RAYS
         iecr(iter)=nusf+iter
      enddo

   end subroutine add_cr_index

   subroutine specify_cosmrays
      use mpi_setup
      implicit none
      write(*,*) 'cosmic rays on as ', COSM_RAYS ,' bins'
      character par_file*(100), tmp_log_file*(100)
      namelist /COSMIC_RAYS/ cr_active, gamma_cr, cr_eff, beta_cr, &
                             K_cr_paral, K_cr_perp, amp_cr, cfl_cr, smallecr

      if(proc .eq. 0) then
         par_file = trim(cwd)//'/problem.par'
         tmp_log_file = trim(cwd)//'/tmp.log'
         open(101,file=par_file)
            read(unit=101,nml=COSMIC_RAYS)
         close(1)
         open(103, file=tmp_log_file, position='append')
         write(unit=103,nml=COSMIC_RAYS)
         close(103)
! export parameters to other procs
         rbuff(130) = cr_active
         rbuff(131) = gamma_cr
         rbuff(132) = cr_eff
         rbuff(133) = beta_cr
         rbuff(134) = K_cr_paral
         rbuff(135) = K_cr_perp
         rbuff(136) = amp_cr
         rbuff(137) = cfl_cr
         rbuff(138) = smallecr
      else
! import parameters from proc=0
         cr_active          = rbuff(130)
         gamma_cr           = rbuff(131)
         cr_eff             = rbuff(132)
         beta_cr            = rbuff(133)
         K_cr_paral         = rbuff(134)
         K_cr_perp          = rbuff(135)
         amp_cr             = rbuff(136)
         cfl_cr             = rbuff(137)
         smallecr           = rbuff(138)
      endif
   end subroutine specify_cosmrays

   subroutine flux_cr(flux,cfr,uu,bb,nu,n)
      use ionizeds, only : idni,imxi
      implicit none
      integer n, nfi, nu
!    real, dimension(COSM_RAYS,n) :: flux
      real, dimension(nu,n)        :: uu,cfr,flux
      real, dimension(3,n)         :: bb
      real, dimension(NUMBION,n)   :: vx

      nfi = NUMBION
!#ifdef COSM_RAYS
!    pcr = (gamma_cr-1)*uu(iecr,RNG)
!#endif COSM_RAYS
      vx(:,RNG)=uu(imxi(1:nfi),RNG)/uu(idni(1:nfi),RNG)
      flux(iecr,RNG)= uu(iecr,RNG)*spread(vx(1,RNG),1,COSM_RAYS)  ! niegotowe, vx to predkosc plynu zjonizowanego?
      cfr(iecr,:)=spread(cfr(idni(1),:),1,COSM_RAYS)

   end subroutine flux_cr

   subroutine src_cr(Duu,uu,sweep,i1,i2,n,dt,divvel,nfluid,nx,ny,nz,idna,imxa,iena,dx)
! pobieznie, do sprawdzenia i przetestowania
      use start,  only : gamma_cr, cr_active, smallecr
!    use arrays, only : divvel
      implicit none
      integer :: nx, ny, nz, nfluid, i1, i2, n
      real :: dt, dx
      character sweep*6
      integer, dimension(nfluid) :: idna, imxa, iena ! iena temporary - should have dimension of nadiab
      real, dimension(nfluid,nx,ny,nz) :: divvel
      real, dimension(nfluid,n) :: divv, Duu, uu
      real, dimension(n)    :: decr,gpcr,ecr,tmp

      select case (sweep)
         case('xsweep')
            divv = divvel(:,:,i1,i2)
         case('ysweep')
            divv = divvel(:,i2,:,i1)
         case('zsweep')
            divv = divvel(:,i1,i2,:)
      end select

      Duu(iecr,:) = Duu(iecr,:) - (gamma_cr-1.)*uu(iecr,:)*divv(1,:)*dt

!      vx  = uu(imxa,:)/uu(idna,:)     ! po co to tutaj?
      ecr = uu(iecr,:)

      gpcr(2:n-1) = cr_active*(gamma_cr -1.)*(ecr(3:n)-ecr(1:n-2))/(2.*dx) ; gpcr(1:2)=0.0 ; gpcr(n-1:n) = 0.0

#ifndef ISO
      Duu(iena,:) = Duu(iena,:) - uu(imxa,:)/uu(idna,:)*gpcr*dt
#endif /* ISO */
      Duu(imxa,:) = Duu(imxa,:) - gpcr*dt

   end subroutine src_cr

end module cosmic_rays
