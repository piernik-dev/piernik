! $Id$
#include "piernik.def"

module fluidindex

   implicit none
   
   integer :: nvar,ncomponents,nfluid,nadiab
   
   integer,allocatable,  dimension(:) :: iarr_all_dn, iarr_all_mx, iarr_all_my, iarr_all_mz
#ifndef ISO 
   integer,allocatable,  dimension(:) :: iarr_all_en
#endif /* ISO */ 
#ifdef COSM_RAYS 
   integer,allocatable,  dimension(:) :: iarr_all_cr
#endif /* COSM_RAYS */ 
   integer, allocatable, dimension(:) :: iarr_all_swpx, iarr_all_swpy, iarr_all_swpz

#ifdef IONIZED
   integer :: nmag
   integer, allocatable, dimension(:) :: iarr_mag_swpx, iarr_mag_swpy, iarr_mag_swpz, iarr_all_mag
   integer, parameter  :: ibx=1, iby=2, ibz=3
#ifdef RESISTIVE
   integer, parameter  :: icx=1, icy=2, icz=3
#endif /RESISTIVE */
#endif /* IONIZED */

#ifdef IONIZED
   integer :: nvar_ion,beg_ion,end_ion,i_ion
#endif /* IONIZED */

#ifdef NEUTRAL
   integer :: nvar_neu,beg_neu,end_neu,i_neu
#endif /* NEUTRAL */

#ifdef DUST
   integer :: nvar_dust,beg_dust,end_dust,i_dust
#endif /* DUST */
   
#ifdef COSM_RAYS
   integer :: nvar_cr,beg_cr,end_cr,i_cr
#endif /* COSM_RAYS */


  contains
   
    subroutine fluid_index
   
#ifdef IONIZED
      use initionized, only 	: ionized_index
      use initionized, only 	: iarr_ion_swpx, iarr_ion_swpy, iarr_ion_swpz
      use initionized, only 	: idni,imxi,imyi,imzi
#ifndef ISO
      use initionized, only 	: ieni
#endif /* ISO */ 
#endif /* IONIZED */

#ifdef NEUTRAL
      use initneutral, only 	: neutral_index
      use initneutral, only 	: iarr_neu_swpx, iarr_neu_swpy, iarr_neu_swpz
      use initneutral, only 	: idnn,imxn,imyn,imzn
#ifndef ISO 
      use initneutral, only 	: ienn
#endif /* ISO */      
#endif /* NEUTRAL */

#ifdef DUST
      use initdust, only 	: dust_index
      use initdust, only 	: iarr_dst_swpx, iarr_dst_swpy, iarr_dst_swpz
      use initdust, only 	: idnd,imxd,imyd,imzd      
#endif /* DUST */

#ifdef COSM_RAYS
      use initcr, only 	: dust_index
      use initcr, only 	: iarr_crs_swpx, iarr_crs_swpy, iarr_crs_swpz
      use initcr, only 	: iecr           
#endif /* COSM_RAYS */

      implicit none


      nvar   = 0
      ncomponents  = 0
      nfluid = 0
      nadiab = 0
      
     
#ifdef IONIZED
      nmag = 3
      nvar_ion  = 0
      beg_ion   = nvar + 1
      call ionized_index(nvar,nvar_ion) 
      end_ion   = nvar  
      ncomponents  = ncomponents + 1
      nfluid = nfluid + 1  
      i_ion = ncomponents 
#ifndef ISO 
      nadiab = nadiab + 1
#endif /* ISO */      
#endif /* IONIZED */

#ifdef NEUTRAL
      nvar_neu = 0
      beg_neu = nvar + 1
      call neutral_index(nvar,nvar_neu) 
      end_neu = nvar    
      ncomponents  = ncomponents + 1
      nfluid = nfluid + 1   
      i_neu = ncomponents 
#ifndef ISO 
      nadiab = nadiab + 1
#endif /* ISO */      
#endif /* NEUTRAL */

#ifdef DUST
      nvar_dust = 0
      beg_dust = nvar + 1
      call dust_index(nvar,nvar_dust) 
      end_dust = nvar    
      ncomponents  = ncomponents + 1
      nfluid = nfluid + 1   
      i_dust = ncomponents 
#endif /* DUST */

#ifdef COSM_RAYS
      nvar_cr   = 0
      beg_cr = nvar + 1
      call cr_index(nvar,nvar_cr)  
      end_cr = nvar   
      ncomponents  = ncomponents + 1   
      i_cr = ncomponents                
#endif /* COSM_RAYS */     

#ifdef IONIZED
      allocate(iarr_mag_swpx(nmag),iarr_mag_swpy(nmag),iarr_mag_swpz(nmag),iarr_all_mag(nmag))
#endif /* IONIZED */
      allocate(iarr_all_swpx(nvar),iarr_all_swpy(nvar),iarr_all_swpz(nvar))
      allocate(iarr_all_dn(nfluid),iarr_all_mx(nfluid),iarr_all_my(nfluid),iarr_all_mz(nfluid))
#ifndef ISO      
      allocate(iarr_all_en(nadiab))  
#endif /* ISO */
#ifdef COSM_RAYS      
      allocate(iarr_all_cr(nvar_cr))  
#endif /* COSM_RAYS */

#ifdef IONIZED
      iarr_mag_swpx = [ibx,iby,ibz]
      iarr_mag_swpy = [iby,ibx,ibz]
      iarr_mag_swpz = [ibz,iby,ibx]
      iarr_all_mag  = [ibx,iby,ibz]
#endif /* IONIZED */


#ifdef IONIZED
      iarr_all_swpx(beg_ion:end_ion) = iarr_ion_swpx
      iarr_all_swpy(beg_ion:end_ion) = iarr_ion_swpy
      iarr_all_swpz(beg_ion:end_ion) = iarr_ion_swpz    

      iarr_all_dn(i_ion)     = idni
      iarr_all_mx(i_ion)     = imxi
      iarr_all_my(i_ion)     = imyi
      iarr_all_mz(i_ion)     = imzi
#ifndef ISO
      iarr_all_en(i_ion)     = ieni
#endif /* ISO */
#endif /* IONIZED */
      
#ifdef NEUTRAL
      iarr_all_swpx(beg_neu:end_neu) = iarr_neu_swpx
      iarr_all_swpy(beg_neu:end_neu) = iarr_neu_swpy
      iarr_all_swpz(beg_neu:end_neu) = iarr_neu_swpz    

      iarr_all_dn(i_neu)    = idnn
      iarr_all_mx(i_neu)    = imxn
      iarr_all_my(i_neu)    = imyn
      iarr_all_mz(i_neu)    = imzn
#ifndef ISO
      iarr_all_en(i_neu)    = ienn
#endif /* ISO */
#endif /* NEUTRAL */
      
#ifdef DUST
      iarr_all_swpx(beg_dust:end_dust) = iarr_dust_swpx
      iarr_all_swpy(beg_dust:end_dust) = iarr_dust_swpy
      iarr_all_swpz(beg_dust:end_dust) = iarr_dust_swpz    

      iarr_all_dn(i_dust)    = idnd
      iarr_all_mx(i_dust)    = imxd
      iarr_all_my(i_dust)    = imyd
      iarr_all_mz(i_dust)    = imzd
#endif /* DUST */
      
#ifdef COSM_RAYS      
      iarr_all_swpx(beg_cr:end_cr) = iarr_cr_swpx
      iarr_all_swpy(beg_cr:end_cr) = iarr_cr_swpy
      iarr_all_swpz(beg_cr:end_cr) = iarr_cr_swpz    

      iarr_all_cr(1:nvar_cr) = iecr 
#endif /* COSM_RAYS */
   
!      write(*,*) 'fluid_index', iarr_ion_swpx
!      write(*,*) 'fluid_index', iarr_all_swpx

   end subroutine fluid_index


end module fluidindex

