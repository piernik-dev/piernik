! $Id$
!
! PIERNIK Code Copyright (C) 2006 Michal Hanasz
!
!    This file is part of PIERNIK code.
!
!    PIERNIK is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    PIERNIK is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with PIERNIK.  If not, see <http://www.gnu.org/licenses/>.
!
!    Initial implemetation of PIERNIK code was based on TVD split MHD code by
!    Ue-Li Pen 
!        see: Pen, Arras & Wong (2003) for algorithm and
!             http://www.cita.utoronto.ca/~pen/MHD 
!             for original source code "mhd.f90" 
!   
!    For full list of developers see $PIERNIK_HOME/license/pdt.txt
!
#include "piernik.def"

module fluidindex

   implicit none
   
   integer, parameter  :: nmag=3
   integer, parameter  :: ibx=1,iby=2,ibz=3
   integer, parameter  :: idn=1,imx=2,imy=3,imz=4,ien=5,icr=1
   integer             :: nvar,ncomponents,nfluid,nadiab

   
   integer,allocatable,  dimension(:) :: iarr_all_dn, iarr_all_mx, iarr_all_my, iarr_all_mz
#ifndef ISO 
   integer,allocatable,  dimension(:) :: iarr_all_en
#endif /* ISO */ 
#ifdef COSM_RAYS 
   integer,allocatable,  dimension(:) :: iarr_all_cr
#endif /* COSM_RAYS */ 
   integer, allocatable, dimension(:) :: iarr_all_swpx, iarr_all_swpy, iarr_all_swpz

#ifdef IONIZED
   integer, allocatable, dimension(:) :: iarr_mag_swpx, iarr_mag_swpy, iarr_mag_swpz, iarr_all_mag
#ifdef RESISTIVE
   integer, parameter  :: icx=1, icy=2, icz=3
#endif /* RESISTIVE */
#endif /* IONIZED */

#ifdef IONIZED
   integer :: nvar_ion,beg_ion,end_ion,i_ion
#endif /* IONIZED */

#ifdef NEUTRAL
   integer :: nvar_neu,beg_neu,end_neu,i_neu
#endif /* NEUTRAL */

#ifdef DUST
   integer :: nvar_dst,beg_dst,end_dst,i_dst
#endif /* DUST */
   
#ifdef COSM_RAYS
   integer :: nvar_crs,beg_crs,end_crs,i_crs
#endif /* COSM_RAYS */


  contains
   
    subroutine fluid_index
   
#ifdef IONIZED
      use initionized,    only : ionized_index
      use initionized,    only : iarr_ion_swpx, iarr_ion_swpy, iarr_ion_swpz
      use initionized,    only : idni,imxi,imyi,imzi
#ifndef ISO
      use initionized,    only : ieni
#endif /* ISO */ 
#endif /* IONIZED */

#ifdef NEUTRAL
      use initneutral,    only : neutral_index
      use initneutral,    only : iarr_neu_swpx, iarr_neu_swpy, iarr_neu_swpz
      use initneutral,    only : idnn,imxn,imyn,imzn
#ifndef ISO 
      use initneutral,    only : ienn
#endif /* ISO */      
#endif /* NEUTRAL */

#ifdef DUST
      use initdust,       only : dust_index
      use initdust,       only : iarr_dst_swpx, iarr_dst_swpy, iarr_dst_swpz
      use initdust,       only : idnd,imxd,imyd,imzd      
#endif /* DUST */

#ifdef COSM_RAYS
      use initcosmicrays, only : cosmicray_index
      use initcosmicrays, only : iarr_crs_swpx, iarr_crs_swpy, iarr_crs_swpz
      use initcosmicrays, only : iecr           
#endif /* COSM_RAYS */

      implicit none


      nvar   = 0
      ncomponents  = 0
      nfluid = 0
      nadiab = 0
      
     
#ifdef IONIZED
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
      nvar_dst = 0
      beg_dst = nvar + 1
      call dust_index(nvar,nvar_dst) 
      end_dst = nvar    
      ncomponents  = ncomponents + 1
      nfluid = nfluid + 1   
      i_dst = ncomponents 
#endif /* DUST */

#ifdef COSM_RAYS
      nvar_crs   = 0
      beg_crs = nvar + 1
      call cosmicray_index(nvar,nvar_crs)  
      end_crs = nvar   
      ncomponents  = ncomponents + 1   
      i_crs = ncomponents                
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
      allocate(iarr_all_cr(nvar_crs))  
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

      iarr_all_dn(i_ion)      = idni
      iarr_all_mx(i_ion)      = imxi
      iarr_all_my(i_ion)      = imyi
      iarr_all_mz(i_ion)      = imzi
#ifndef ISO
      iarr_all_en(i_ion)      = ieni
#endif /* ISO */
#endif /* IONIZED */
      
#ifdef NEUTRAL
      iarr_all_swpx(beg_neu:end_neu) = iarr_neu_swpx
      iarr_all_swpy(beg_neu:end_neu) = iarr_neu_swpy
      iarr_all_swpz(beg_neu:end_neu) = iarr_neu_swpz    

      iarr_all_dn(i_neu)      = idnn
      iarr_all_mx(i_neu)      = imxn
      iarr_all_my(i_neu)      = imyn
      iarr_all_mz(i_neu)      = imzn
#ifndef ISO
      iarr_all_en(i_neu)      = ienn
#endif /* ISO */
#endif /* NEUTRAL */
      
#ifdef DUST
      iarr_all_swpx(beg_dst:end_dst) = iarr_dst_swpx
      iarr_all_swpy(beg_dst:end_dst) = iarr_dst_swpy
      iarr_all_swpz(beg_dst:end_dst) = iarr_dst_swpz    

      iarr_all_dn(i_dst)      = idnd
      iarr_all_mx(i_dst)      = imxd
      iarr_all_my(i_dst)      = imyd
      iarr_all_mz(i_dst)      = imzd
#endif /* DUST */
      
#ifdef COSM_RAYS      
      iarr_all_swpx(beg_crs:end_crs) = iarr_crs_swpx
      iarr_all_swpy(beg_crs:end_crs) = iarr_crs_swpy
      iarr_all_swpz(beg_crs:end_crs) = iarr_crs_swpz    

      iarr_all_cr(1:nvar_crs) = iecr 
#endif /* COSM_RAYS */
   
!      write(*,*) 'fluid_index', iarr_ion_swpx
!      write(*,*) 'fluid_index', iarr_all_swpx

   end subroutine fluid_index


end module fluidindex

