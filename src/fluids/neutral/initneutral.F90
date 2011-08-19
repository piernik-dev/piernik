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
#include "piernik.h"
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
! pulled by NEUTRAL
   implicit none

   public ! QA_WARN no secrets are kept here

   real                  :: gamma_neu             !< adiabatic index for the neutral gas component
   real                  :: cs_iso_neu            !< isothermal sound speed (p = cs_iso_neu<sup>2</sup>\f$\rho\f$), active only if neutral gas is \ref isothermal
   real                  :: cs_iso_neu2
   logical               :: selfgrav_neu          !< true if neutral gas is selfgravitating
   integer(kind=4)       :: idnn, imxn, imyn, imzn
#ifndef ISO
   integer(kind=4)       :: ienn
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
!! <tr><td>gamma_neu     </td><td>5./3.   </td><td>real value </td><td>\copydoc initneutral::gamma_neu    </td></tr>
!! <tr><td>cs_iso_neu    </td><td>1.0     </td><td>real value </td><td>\copydoc initneutral::cs_iso_neu   </td></tr>
!! <tr><td>selfgrav_neu  </td><td>.false. </td><td>logical    </td><td>\copydoc initneutral::selfgrav_neu </td></tr>
!! </table>
!! \n \n
!<
   subroutine init_neutral

      use dataio_pub,    only: par_file, ierrh, namelist_errh, compare_namelist, cmdl_nml ! QA_WARN required for diff_nml
      use mpisetup,      only: rbuff, lbuff, comm, ierr, buffer_dim, master, slave, FIRST
      use mpi,           only: MPI_DOUBLE_PRECISION, MPI_LOGICAL

      implicit none

      namelist /FLUID_NEUTRAL/ gamma_neu, cs_iso_neu, selfgrav_neu

      gamma_neu    = 5./3.
      cs_iso_neu   = 1.0
      selfgrav_neu = .false.

      if (master) then

         diff_nml(FLUID_NEUTRAL)

         lbuff(1)   = selfgrav_neu

         rbuff(1)   = gamma_neu
         rbuff(2)   = cs_iso_neu

      endif

      call MPI_Bcast(rbuff,    buffer_dim, MPI_DOUBLE_PRECISION, FIRST, comm, ierr)
      call MPI_Bcast(lbuff,    buffer_dim, MPI_LOGICAL,          FIRST, comm, ierr)

      if (slave) then

         selfgrav_neu = lbuff(1)

         gamma_neu   = rbuff(1)
         cs_iso_neu  = rbuff(2)

      endif

      cs_iso_neu2      = cs_iso_neu**2

   end subroutine init_neutral

   subroutine neutral_index(flind)
      use constants,    only: NEU, INT4, xdim, ydim, zdim, ndims
      use diagnostics,  only: ma1d, ma2d, my_allocate
      use fluidtypes,   only: var_numbers

      implicit none

      type(var_numbers), intent(inout) :: flind

      flind%neu%beg    = flind%all + 1_INT4

      idnn = flind%all + 1_INT4
      imxn = flind%all + 2_INT4
      imyn = flind%all + 3_INT4
      imzn = flind%all + 4_INT4

      flind%neu%idn = idnn
      flind%neu%imx = imxn
      flind%neu%imy = imyn
      flind%neu%imz = imzn

      flind%neu%all  = 4
      flind%all      = imzn
#ifndef ISO
      ienn          = imzn + 1_INT4
      flind%all      = flind%all + 1_INT4
      flind%neu%all  = flind%neu%all +1_INT4
      flind%neu%ien  = ienn
#endif /* !ISO */

      ma1d = [flind%neu%all]
      call my_allocate(flind%neu%iarr,       ma1d)
      ma2d = [ndims, flind%neu%all]
      call my_allocate(flind%neu%iarr_swp,   ma2d)

      !\deprecated repeated magic integers
      flind%neu%iarr(1:4)           = [idnn,imxn,imyn,imzn]
      flind%neu%iarr_swp(xdim, 1:4) = [idnn,imxn,imyn,imzn]
      flind%neu%iarr_swp(ydim, 1:4) = [idnn,imyn,imxn,imzn]
      flind%neu%iarr_swp(zdim, 1:4) = [idnn,imzn,imyn,imxn]

#ifndef ISO
      flind%neu%iarr(5)       = ienn
      flind%neu%iarr_swp(:,5) = ienn
      flind%neu%has_energy    = .true.

      flind%energ = flind%energ + 1_INT4
#endif /* !ISO */

      flind%neu%end    = flind%all
      flind%components = flind%components + 1_INT4
      flind%fluids     = flind%fluids + 1_INT4
      flind%neu%pos    = flind%components
      if (selfgrav_neu)  flind%fluids_sg = flind%fluids_sg + 1_INT4

      flind%neu%gam   = gamma_neu
      flind%neu%gam_1 = gamma_neu-1.0
      flind%neu%cs    = cs_iso_neu
      flind%neu%cs2   = cs_iso_neu**2
      flind%neu%tag   = NEU

      flind%neu%is_selfgrav   = selfgrav_neu
      flind%neu%is_magnetized = .false.
#ifndef ISO
      flind%neu%has_energy    = .true.
#endif /* !ISO */

   end subroutine neutral_index

   subroutine cleanup_neutral

      implicit none

   end subroutine cleanup_neutral

end module initneutral
