! $Id$
#include "piernik.def"
#ifdef ISO
#define NUPF 4
#else /* ISO */
#define NUPF 5
#endif /* ISO */

module ionizeds
   implicit none
   integer, dimension(IONIZED)       :: idni, imxi, imyi, imzi, ieni
   integer, dimension(IONIZED*NUPF)  :: aions
   integer, parameter :: nuipf = NUPF
   contains

   subroutine add_ion_index
      implicit none
      integer :: iter
      do iter=1,IONIZED
         idni(iter)=(iter-1)*nuipf+1
         imxi(iter)=(iter-1)*nuipf+2
         imyi(iter)=(iter-1)*nuipf+3
         imzi(iter)=(iter-1)*nuipf+4
#ifndef ISO
         ieni(iter)=(iter-1)*nuipf+5
#endif /* ISO */
      enddo

   end subroutine add_ion_index

   subroutine specify_ionized(nfl,nad,ifl,iadiab,ifion,imagn,idna,imxa,imya,imza,iena,findex,magn,fiond,fmagn,gamma,fadiab)
      use mpi_setup
      implicit none
      integer :: iadiab, ifion, imagn, nfl,nad
      integer :: ifl, iter
      integer, dimension(nfl) :: idna,imxa,imya,imza,findex,magn
      integer, dimension(nad) :: iena,fadiab
      integer, dimension(IONIZED) :: fiond,fmagn
      real, dimension(nfl) :: gamma
      character fluidtype*10, par_file*(100), tmp_log_file*(100)

!      namelist /IONIZED_GAS/

      if(proc .eq. 0) then
         par_file = trim(cwd)//'/problem.par'
         tmp_log_file = trim(cwd)//'/tmp.log'
         open(101,file=par_file)
!            read(unit=101,nml=IONIZED_GAS)
         close(1)
! export parameters to other procs
      else
! import parameters from proc=0
      endif
      do iter=1,IONIZED
         idna(ifl)=idni(iter)
         imxa(ifl)=imxi(iter)
         imya(ifl)=imyi(iter)
         imza(ifl)=imzi(iter)

         aions((iter-1)*nuipf+1)=idni(iter)
         aions((iter-1)*nuipf+2)=imxi(iter)
         aions((iter-1)*nuipf+3)=imyi(iter)
         aions((iter-1)*nuipf+4)=imzi(iter)

         findex(ifl)=nuipf
         magn(ifl)=1
         fiond(ifion)=ifl
         fmagn(imagn)=ifl
         ifion=ifion+1
         imagn=imagn+1

#ifdef ISO
         fluidtype = 'isothermal'
         gamma(ifl) = 1.0
#else /* ISO */
         iena(ifl)=ieni(iter)
         aions((iter-1)*nuipf+5)=ieni(iter)
         fadiab(iadiab)=ifl
         iadiab=iadiab+1
         fluidtype = 'adiabatic'
#endif /* ISO */
         write(*,'(a6,i2,a11,a10)') 'fluid ',ifl,' ionized   ',fluidtype
!#include "arrions.def"
         ifl=ifl+1
      enddo
   end subroutine specify_ionized

   subroutine flux_ionized(iflux,icfr,iuu,nuio,bb,n)
!use arra   ys, only : nui
      implicit none
      integer n, nuio
      real, dimension(nuio,n) :: iflux,icfr,iuu
      real, dimension(3,n)    :: bb


   end subroutine flux_ionized

   subroutine src_ionized
      implicit none
   end subroutine src_ionized

end module ionizeds
