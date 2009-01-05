! $Id$
#include "piernik.def"
#define NUPF 4

module dusts
   implicit none
   integer, dimension(DUST)       :: idnd, imxd, imyd, imzd
   integer, dimension(DUST*NUPF)  :: adust
   integer, parameter :: nudpf = NUPF

   contains

   subroutine add_dust_index(nusf)
      implicit none
      integer :: iter, nusf
      do iter=1,DUST
         idnd(iter)=nusf+(iter-1)*nudpf+1
         imxd(iter)=nusf+(iter-1)*nudpf+2
         imyd(iter)=nusf+(iter-1)*nudpf+3
         imzd(iter)=nusf+(iter-1)*nudpf+4
      enddo

   end subroutine add_dust_index

   subroutine specify_dust(nfl,ifl,idna,imxa,imya,imza,findex,gamma,fdust,idst)
      use mpi_setup
      implicit none
      integer :: nfl, ifl, iter, idst
      integer, dimension(nfl) :: idna, imxa, imya, imza, findex
      integer, dimension(DUST) :: fdust
      real, dimension(nfl) :: gamma
      character fluidtype*10, par_file*(100), tmp_log_file*(100)

!      namelist /DUST_GAS/

      if(proc .eq. 0) then
         par_file = trim(cwd)//'/problem.par'
         tmp_log_file = trim(cwd)//'/tmp.log'
         open(101,file=par_file)
!            read(unit=101,nml=DUST_GAS)
         close(1)
! export parameters to other procs
      else
! import parameters from proc=0
      endif
      do iter = 1, DUST
         idna(ifl)=idnd(iter)
         imxa(ifl)=imxd(iter)
         imya(ifl)=imyd(iter)
         imza(ifl)=imzd(iter)

         adust((iter-1)*nudpf+1)=idnd(iter)
         adust((iter-1)*nudpf+2)=imxd(iter)
         adust((iter-1)*nudpf+3)=imyd(iter)
         adust((iter-1)*nudpf+4)=imzd(iter)

         findex(ifl)=nudpf
         fdust(idst)=ifl
         fluidtype = 'isothermal'
         gamma(ifl) = 1.0

         write(*,'(a6,i2,a11,a10)') 'fluid ',ifl,' dust      ',fluidtype
         ifl=ifl+1
         idst=idst+1
      enddo

   end subroutine specify_dust

   subroutine flux_dust(dflux,dcfr,duu,nudo,n)
!      use arrays, only : nud
      implicit none
      integer n, nudo
      real,dimension(nudo,n) :: dflux,dcfr,duu

   end subroutine flux_dust

   subroutine src_dust
      implicit none
   end subroutine src_dust

end module dusts
