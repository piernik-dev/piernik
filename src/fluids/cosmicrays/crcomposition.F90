! $Id: types.F90 1403 2009-12-06 14:08:48Z xarth $
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
!>
!! \brief (MH)
!<
module crcomposition

   type cr_component
      character(len=8)  :: isotope     !< isotope name, eg. Be10
      integer  :: index=-1             !< relative index (with respect to crn_beg)
      real     :: abund=0.0            !< initial abundance relative to H
   end type cr_component

   integer     :: icr_H1   = 1
   integer     :: icr_C12  = 2
   integer     :: icr_N14  = 3
   integer     :: icr_O16  = 4
   integer     :: icr_Li7  = 5         !< from decay of Be7 with tau of 0.3 years
   integer     :: icr_Be9  = 6
   integer     :: icr_Be10 = 7    !!! BEWARE: ncrn should be set up gerater than maximum isotope numeber,
                                  !!! which should be smaller than ncr_max (<10 currently)


!<====Cross sections for spallation from Garcia-Munoz 1987 (see also Longair)====>

   real,parameter   :: mbarn=1e-27 !cm2   !!! BEWARE: this line breaks unit consistency, move it to constants.F90 and use scaling

   real   :: sigma_C12_Li7  = 10*mbarn
   real   :: sigma_C12_Be9  = 6*mbarn
   real   :: sigma_C12_Be10 = 3.5*mbarn

   real   :: sigma_N14_Li7  = 9.5*mbarn

   real   :: sigma_O16_Li7  = 9.5*mbarn
   real   :: sigma_O16_Be9  = 4.5*mbarn
   real   :: sigma_O16_Be10 = 2*mbarn

!<====Decay half live times from Garcia-Munoz 1987====>

   real,parameter :: Myear=1d6*365*24*60*60 !s !!! BEWARE: this line breaks unit consistency, move it to constants.F90 and use scaling

   real   :: tau_Be10 = 1.6 !Myr !!! BEWARE: this line breaks unit consistency, move it to constants.F90 and use scaling

!<Initial source abundances (in numer density) relative to hydrogen (compare e.g. Longair)>

   real   :: primary_C12  =  4.5e-3
   real   :: primary_N14  =  1e-3
   real   :: primary_O16  =  4e-3

end module crcomposition
