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
#include "macros.h"

!>
!! \brief (MH) Initialization of the neutral fluid
!!
!!
!!
!! In this module following namelist of parameters is specified:
!! \copydetails initneutral::init_neutral
!<

module initneutral

  implicit none

    real                  :: gamma_neu             !< adiabatic index for the neutral gas
    real                  :: cs_iso_neu            !< isothermal sound speed (p = cs_iso_neu<sup>2</sup>\f$\rho\f$), active only if neutral gas is \ref isothermal
    real                  :: cs_iso_neu2
    real                  :: global_gradP_neu
    real                  :: eta_gas_neu
    real                  :: csvk
    logical               :: selfgrav_neu

    integer               :: idnn, imxn, imyn, imzn
#ifndef ISO
    integer               :: ienn
#endif /* !ISO */

    integer, allocatable, dimension(:)  :: iarr_neu
    integer, allocatable, dimension(:)  :: iarr_neu_swpx, iarr_neu_swpy, iarr_neu_swpz

  contains

!>
!! \brief Routine to set parameters values from namelist FLUID_NEUTRAL
!!
!! \n \n
!! @b FLUID_NEUTRAL
!! \n \n
!! <table border="+1">
!! <tr><td width="150pt"><b>parameter</b></td><td width="135pt"><b>default value</b></td><td width="200pt"><b>possible values</b></td><td width="315pt"> <b>description</b></td></tr>
!! <tr><td>gamma_neu  </td><td>1.66666666</td><td>real value</td><td>\copydoc initneutral::gamma_neu  </td></tr>
!! <tr><td>cs_iso_neu </td><td>1.0       </td><td>real value</td><td>\copydoc initneutral::cs_iso_neu </td></tr>
!! <tr><td>eta_gas_neu</td><td>          </td><td>real value</td><td>\copydoc initneutral::eta_gas_neu</td></tr>
!! <tr><td>csvk       </td><td>          </td><td>real value</td><td>\copydoc initneutral::csvk       </td></tr>
!! </table>
!! \n \n
!<
  subroutine init_neutral

    use mpisetup, only: proc, ierr, comm, rbuff, lbuff, buffer_dim, MPI_LOGICAL, MPI_DOUBLE_PRECISION
    use errh,     only: namelist_errh
    use dataio_public, only: par_file, ierrh
    use func,        only : compare_namelist

#ifdef SHEAR
    use shear,    only: omega
#endif /* SHEAR */

    implicit none


    namelist /FLUID_NEUTRAL/ gamma_neu, cs_iso_neu, eta_gas_neu, csvk, selfgrav_neu

    gamma_neu    = 1.66666666
    cs_iso_neu   = 1.0
    selfgrav_neu = .false.

    if (proc == 0) then

       diff_nml(FLUID_NEUTRAL)

       lbuff(1)   = selfgrav_neu

       rbuff(1)   = gamma_neu
       rbuff(2)   = cs_iso_neu
       rbuff(3)   = eta_gas_neu
       rbuff(4)   = csvk

    end if

    call MPI_Bcast(rbuff,    buffer_dim, MPI_DOUBLE_PRECISION, 0, comm, ierr)
    call MPI_Bcast(lbuff,    buffer_dim, MPI_LOGICAL,          0, comm, ierr)

    if (proc /= 0) then

      selfgrav_neu = lbuff(1)

      gamma_neu   = rbuff(1)
      cs_iso_neu  = rbuff(2)
      eta_gas_neu = rbuff(3)
      csvk        = rbuff(4)

    endif

    cs_iso_neu2      = cs_iso_neu**2
#ifdef SHEAR
    global_gradP_neu = 2.0*omega*eta_gas_neu * cs_iso_neu / csvk
#else /* SHEAR */
    global_gradP_neu = 0.0
#endif /* SHEAR */

  end subroutine init_neutral



  subroutine neutral_index(nvar,nvar_neu)

    implicit none
    integer :: nvar, nvar_neu


      idnn = nvar + 1
      imxn = nvar + 2
      imyn = nvar + 3
      imzn = nvar + 4

#ifdef ISO
      nvar_neu      = 4
      nvar          = imzn

      allocate(iarr_neu(nvar_neu),iarr_neu_swpx(nvar_neu), iarr_neu_swpy(nvar_neu), iarr_neu_swpz(nvar_neu))

      iarr_neu      = [idnn,imxn,imyn,imzn]
      iarr_neu_swpx = [idnn,imxn,imyn,imzn]
      iarr_neu_swpy = [idnn,imyn,imxn,imzn]
      iarr_neu_swpz = [idnn,imzn,imyn,imxn]
#else /* ISO */
      ienn          = nvar + 5
      nvar_neu      = 5
      nvar          = ienn

      allocate(iarr_neu(nvar_neu),iarr_neu_swpx(nvar_neu), iarr_neu_swpy(nvar_neu), iarr_neu_swpz(nvar_neu))

      iarr_neu      = [idnn,imxn,imyn,imzn,ienn]
      iarr_neu_swpx = [idnn,imxn,imyn,imzn,ienn]
      iarr_neu_swpy = [idnn,imyn,imxn,imzn,ienn]
      iarr_neu_swpz = [idnn,imzn,imyn,imxn,ienn]
#endif /* ISO */

   end subroutine neutral_index

end module initneutral
