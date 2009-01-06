! $Id$
#include "piernik.def"
#ifdef ISO
#define NUPF 4
#else /* ISO */
#define NUPF 5
#endif /* ISO */

module neutrals
   implicit none
   integer, dimension(NEUTRAL)       :: idnn, imxn, imyn, imzn, ienn
   integer, dimension(NEUTRAL*NUPF)  :: aneut
   integer, parameter :: nunpf = NUPF

   contains

   subroutine add_neut_index(nusf)
      implicit none
      integer :: iter, nusf
      do iter=1,NEUTRAL
         idnn(iter)=nusf+(iter-1)*nunpf+1
         imxn(iter)=nusf+(iter-1)*nunpf+2
         imyn(iter)=nusf+(iter-1)*nunpf+3
         imzn(iter)=nusf+(iter-1)*nunpf+4
#ifndef ISO
         ienn(iter)=nusf+(iter-1)*nunpf+5
#endif /* ISO */
      enddo

   end subroutine add_neut_index

   subroutine specify_neutral(nfl,nad,ifl,iadiab,idna,imxa,imya,imza,iena,findex,gamma,fadiab)
      use mpi_setup
      implicit none
      integer :: iadiab, nfl, nad
      integer :: ifl, iter
      integer, dimension(nfl) :: idna, imxa, imya, imza, findex
      integer, dimension(nad) :: iena, fadiab
      real, dimension(nfl) :: gamma
      character fluidtype*10, par_file*(100), tmp_log_file*(100)

!      namelist /NEUTRAL_GAS/

      if(proc .eq. 0) then
         par_file = trim(cwd)//'/problem.par'
         tmp_log_file = trim(cwd)//'/tmp.log'
         open(101,file=par_file)
!            read(unit=101,nml=NEUTRAL_GAS)
         close(1)
! export parameters to other procs
      else
! import parameters from proc=0
      endif
      do iter = 1, NEUTRAL
         idna(ifl)=idnn(iter)
         imxa(ifl)=imxn(iter)
         imya(ifl)=imyn(iter)
         imza(ifl)=imzn(iter)

         aneut((iter-1)*nunpf+1)=idnn(iter)
         aneut((iter-1)*nunpf+2)=imxn(iter)
         aneut((iter-1)*nunpf+3)=imyn(iter)
         aneut((iter-1)*nunpf+4)=imzn(iter)

         findex(ifl)=nunpf

#ifdef ISO
         fluidtype = 'isothermal'
         gamma(ifl) = 1.0
#else /* ISO */
         iena(ifl)=ienn(iter)
         aneut((iter-1)*nunpf+5)=ienn(iter)
         fadiab(iadiab)=ifl
         iadiab=iadiab+1
         fluidtype = 'adiabatic'
#endif /* ISO */
         write(*,'(a6,i2,a11,a10)') 'fluid ',ifl,' neutral   ',fluidtype
         ifl=ifl+1
      enddo
   end subroutine specify_neutral

   subroutine flux_neutral(nflux,ncfr,nuu,nuno,n)
!      use arrays, only : nun
      implicit none
      integer n, nuno
      real,dimension(nuno,n) :: nflux,ncfr,nuu

   end subroutine flux_neutral

   subroutine src_neutral
      implicit none
   end subroutine src_neutral

end module neutrals
