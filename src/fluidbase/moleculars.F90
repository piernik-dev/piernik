! $Id$
#include "piernik.def"
#ifdef ISO
#define NUPF 4
#else /* ISO */
#define NUPF 5
#endif /* ISO */

module moleculars
   implicit none
   integer, dimension(MOLECULAR)      :: idnm, imxm, imym, imzm, ienm
   integer, dimension(MOLECULAR*NUPF) :: amols
   integer, parameter :: numpf = NUPF

   contains

   subroutine add_mol_index(nusf)
      implicit none
      integer :: iter, nusf
      do iter=1,MOLECULAR
         idnm(iter)=nusf+(iter-1)*numpf+1
         imxm(iter)=nusf+(iter-1)*numpf+2
         imym(iter)=nusf+(iter-1)*numpf+3
         imzm(iter)=nusf+(iter-1)*numpf+4
#ifndef ISO
         ienm(iter)=nusf+(iter-1)*numpf+5
#endif /* ISO */
      enddo

   end subroutine add_mol_index

   subroutine specify_molecular(nfl,nad,ifl,iadiab,idna,imxa,imya,imza,iena,findex,gamma,fadiab)
      implicit none
      integer :: iadiab, nfl, nad
      integer :: ifl, iter
      integer, dimension(nfl) :: idna, imxa, imya, imza, findex
      integer, dimension(nad) :: iena, fadiab
      real, dimension(nfl) :: gamma
      character fluidtype*10, par_file*(100), tmp_log_file*(100)

!      namelist /MOLECULAR_GAS/

      if(proc .eq. 0) then
         par_file = trim(cwd)//'/problem.par'
         tmp_log_file = trim(cwd)//'/tmp.log'
         open(101,file=par_file)
!            read(unit=101,nml=MOLECULAR_GAS)
         close(1)
! export parameters to other procs
      else
! import parameters from proc=0
      endif
      do iter = 1, MOLECULAR
         idna(ifl)=idnm(iter)
         imxa(ifl)=imxm(iter)
         imya(ifl)=imym(iter)
         imza(ifl)=imzm(iter)

         amols((iter-1)*numpf+1)=idnm(iter)
         amols((iter-1)*numpf+2)=imxm(iter)
         amols((iter-1)*numpf+3)=imym(iter)
         amols((iter-1)*numpf+4)=imzm(iter)

         findex(ifl)=numpf

#ifdef ISO
         fluidtype = 'isothermal'
         gamma(ifl) = 1.0
#else /* ISO */
         iena(ifl)=ienm(iter)
         amols((iter-1)*numpf+5)=ienm(iter)
         fadiab(iadiab)=ifl
         iadiab=iadiab+1
         fluidtype = 'adiabatic'
#endif /* ISO */
         write(*,'(a6,i2,a11,a10)') 'fluid ',ifl,' molecular ',fluidtype
         ifl=ifl+1
      enddo

   end subroutine specify_molecular

   subroutine flux_molecular(mflux,mcfr,muu,numo,n)
!      use arrays, only : num
      implicit none
      integer n, numo
      real,dimension(numo,n) :: mflux,mcfr,muu

   end subroutine flux_molecular

   subroutine src_molecular
      implicit none
   end subroutine src_molecular

end module moleculars
