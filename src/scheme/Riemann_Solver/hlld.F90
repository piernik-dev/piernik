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
!    Initial implementation of PIERNIK code was based on TVD split MHD code by
!    Ue-Li Pen
!        see: Pen, Arras & Wong (2003) for algorithm and
!             http://www.cita.utoronto.ca/~pen/MHD
!             for original source code "mhd.f90"
!
!    For full list of developers see $PIERNIK_HOME/license/pdt.txt
!
!  References:
!
!  A multi-state HLLD approximate Riemann solver for ideal magnetohydrodynamics.
!  Takahiro Miyoshi, Kanya Kusano
!  Journal of Computational Physics 208 (2005) 315-344
! 
!  ->Solve one dimensional Riemann problem using adiabatic HLLD scheme 
!
!  Varadarajan Parthasarathy, CAMK, Warszawa. 2015.
! 
!---------------------------------------------------------------------------------------------------------------------------


#include "piernik.def"

module hlld
! pulled by RIEMANN

  implicit none 
  
  private
  public :: riemann_hlld

contains

  subroutine riemann_hlld(q, f, n, b)


    use constants,  only: zero, one, half, idn, imx, imy, imz, ien

    implicit none

    real, dimension(:,:), intent(in)              :: b

    integer,              intent(in)              :: n

    real, dimension(:,:), pointer, intent(in)     :: q

    real, dimension(:,:), pointer, intent(inout)  :: f

    ! local variables 

    real    :: sl, sr, sm, sml, smr
    real    :: srl, srml, slmv, srmv, slmm, srmm, smvl, smvr
    real    :: dn, dnl, dnr, dlsq, drsq
    real    :: mx, dv, fc, b2, bs, ds, vbl, vbr, vb1l, vb1r, vb2
    real    :: pml, pmr, ptl, ptr, pt, pm

    

    ! Follow the references and Kowal code as an example for computing fluxes.


  end subroutine riemann_hlld

end module hlld
