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
!    Initial implementation of PIERNIK code was based on TVD split MHD code by
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

    use dataio_public,  only: par_file, ierrh, namelist_errh, compare_namelist
    use mpisetup,       only: proc, ierr, comm, rbuff, lbuff, buffer_dim, MPI_LOGICAL, MPI_DOUBLE_PRECISION
#ifdef SHEAR
    use shear,          only: omega
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

    endif

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

   subroutine neutral_index(nvar)
      use fluidtypes,   only: var_numbers
      use diagnostics,  only: my_allocate
      implicit none
      type(var_numbers), intent(inout) :: nvar

      nvar%neu%beg    = nvar%all + 1

      idnn = nvar%all + 1
      imxn = nvar%all + 2
      imyn = nvar%all + 3
      imzn = nvar%all + 4

      nvar%neu%idn = idnn
      nvar%neu%imx = imxn
      nvar%neu%imy = imyn
      nvar%neu%imz = imzn

      nvar%neu%all  = 4
      nvar%all      = imzn
#ifndef ISO
      ienn          = imzn + 1
      nvar%neu%ien  = ienn
      nvar%all      = nvar%all + 1
      nvar%neu%all  = nvar%neu%all +1
#endif /* !ISO */

      call my_allocate(nvar%iarr_neu,       [nvar%neu%all], "iarr_neu")
      call my_allocate(nvar%iarr_neu_swpx,  [nvar%neu%all], "iarr_neu_swpx")
      call my_allocate(nvar%iarr_neu_swpy,  [nvar%neu%all], "iarr_neu_swpy")
      call my_allocate(nvar%iarr_neu_swpz,  [nvar%neu%all], "iarr_neu_swpz")

      nvar%iarr_neu(1:4)      = [idnn,imxn,imyn,imzn]
      nvar%iarr_neu_swpx(1:4) = [idnn,imxn,imyn,imzn]
      nvar%iarr_neu_swpy(1:4) = [idnn,imyn,imxn,imzn]
      nvar%iarr_neu_swpz(1:4) = [idnn,imzn,imyn,imxn]

#ifndef ISO
      nvar%iarr_neu(5)      = ienn
      nvar%iarr_neu_swpx(5) = ienn
      nvar%iarr_neu_swpy(5) = ienn
      nvar%iarr_neu_swpz(5) = ienn

      nvar%adiab = nvar%adiab + 1
#endif /* ISO */

      nvar%neu%end    = nvar%all
      nvar%components = nvar%components + 1
      nvar%fluids     = nvar%fluids + 1
      nvar%neu%pos    = nvar%components
      if (selfgrav_neu)  nvar%fluids_sg = nvar%fluids_sg + 1

   end subroutine neutral_index

end module initneutral
