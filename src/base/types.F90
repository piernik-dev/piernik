!$Id$
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
module types
   type :: indx
      integer :: dnd = -1, dnn = -1, dni = -1
      integer :: mxd = -1, mxn = -1, mxi = -1
      integer :: myd = -1, myn = -1, myi = -1
      integer :: mzd = -1, mzn = -1, mzi = -1
      integer :: enn = -1, eni = -1
      integer :: ecr = -1
      integer :: bx = -1, by = -1, bz = -1
   end type indx

   type :: hdf
      integer :: nhdf, ntsl, nres, nlog, step_hdf, step_res, log_lun, &
         nstep, nrestart
      real    :: last_hdf_time
      character(len=128) :: log_file
      character(len=16)  :: domain
      character(len=3)   :: new_id
   end type hdf

   type :: value
      real    :: val
      integer, dimension(3) :: loc
      integer :: proc
   end type value

   type :: tsl_container
#ifdef NEUTRAL 
      real :: denn_min, denn_max, vxn_max, vyn_max, vzn_max, &
              pren_min, pren_max, temn_min, temn_max, csn_max
#endif /* NEUTRAL */

#ifdef DUST
      real :: dend_min, dend_max, vxd_max, vyd_max, vzd_max
#endif /* DUST */

#ifdef COSM_RAYS
      real :: encr_min, encr_max
#endif /* COSM_RAYS */

#ifdef RESISTIVE
      real :: etamax
#endif /* RESISTIVE */

#ifdef MAGNETIC
      real :: b_min, b_max, divb_max
#endif /* MAGNETIC */ 

#ifdef IONIZED 
      real :: deni_min, deni_max, vxi_max, vyi_max, vzi_max, &
              prei_min, prei_max, temi_min, temi_max, vai_max, csi_max      
#endif /* IONIZED */
   end type tsl_container

end module types
