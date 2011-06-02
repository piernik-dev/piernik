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
!! \brief (MH) Initialization of the ionized fluid
!!
!!
!!
!! In this module following namelist of parameters is specified:
!! \copydetails initionized::init_ionized
!<

module initionized
! pulled by IONIZED
   implicit none

   public ! QA_WARN no secrets are kept here

   real                  :: gamma_ion       !< adiabatic index for the ionized gas component
   real                  :: cs_iso_ion      !< isothermal sound speed (p = cs_iso_ion<sup>2</sup>\f$\rho\f$), active only if ionized gas is \ref isothermal
   real                  :: cs_iso_ion2
   real                  :: cs_ion
   logical               :: selfgrav_ion    !< true if ionized gas is selfgravitating
   integer               :: idni, imxi, imyi, imzi
#ifndef ISO
   integer               :: ieni
#endif /* !ISO */

contains

!>
!! \brief Routine to set parameters values from namelist FLUID_IONIZED
!!
!! \n \n
!! @b FLUID_IONIZED
!! \n \n
!! <table border="+1">
!! <tr><td width="150pt"><b>parameter</b></td><td width="135pt"><b>default value</b></td><td width="200pt"><b>possible values</b></td><td width="315pt"> <b>description</b></td></tr>
!! <tr><td>gamma_ion     </td><td>1.66666666</td><td>real value </td><td>\copydoc initionized::gamma_ion </td></tr>
!! <tr><td>cs_iso_ion    </td><td>1.0        </td><td>real value</td><td>\copydoc initionized::cs_iso_ion</td></tr>
!! <tr><td>cs_ion        </td><td>           </td><td>real value</td><td>\copydoc initionized::cs_ion    </td></tr>
!! <tr><td>selfgrav_ion  </td><td>.false.    </td><td>logical   </td><td>\copydoc initionized::selfgrav_ion </td></tr>
!! </table>
!! \n \n
!<
   subroutine init_ionized

      use dataio_pub,    only: par_file, ierrh, namelist_errh, compare_namelist, cmdl_nml ! QA_WARN required for diff_nml
      use mpisetup,      only: rbuff, lbuff, comm, ierr, buffer_dim, master, slave
      use mpi,           only: MPI_DOUBLE_PRECISION, MPI_LOGICAL

      implicit none

      namelist /FLUID_IONIZED/ gamma_ion, cs_iso_ion, cs_ion, selfgrav_ion

      gamma_ion     = 1.66666666
      cs_iso_ion    = 1.0
      selfgrav_ion  = .false.

      if (master) then

         diff_nml(FLUID_IONIZED)

         lbuff(1)   = selfgrav_ion

         rbuff(1)   = gamma_ion
         rbuff(2)   = cs_iso_ion
         rbuff(3)   = cs_ion

      endif

      call MPI_Bcast(rbuff,    buffer_dim, MPI_DOUBLE_PRECISION, 0, comm, ierr)
      call MPI_Bcast(lbuff,    buffer_dim, MPI_LOGICAL,          0, comm, ierr)

      if (slave) then

         selfgrav_ion = lbuff(1)

         gamma_ion    = rbuff(1)
         cs_iso_ion   = rbuff(2)
         cs_ion       = rbuff(3)

      endif

      cs_iso_ion2  = cs_iso_ion**2

   end subroutine init_ionized

   subroutine ionized_index(flind)

      use constants,    only: ION
      use diagnostics,  only: ma1d, my_allocate
      use fluidtypes,   only: var_numbers

      implicit none

      type(var_numbers), intent(inout) :: flind

      flind%ion%beg  = flind%all + 1

      idni = flind%all + 1
      imxi = flind%all + 2
      imyi = flind%all + 3
      imzi = flind%all + 4

      flind%ion%idn  = idni
      flind%ion%imx  = imxi
      flind%ion%imy  = imyi
      flind%ion%imz  = imzi

      flind%ion%all  = 4
      flind%all      = imzi
#ifndef ISO
      ieni          = imzi + 1
      flind%all      = flind%all + 1
      flind%ion%all  = flind%ion%all + 1
      flind%ion%ien  = ieni
#endif /* !ISO */

      ma1d = [flind%ion%all]
      call my_allocate(flind%ion%iarr,       ma1d)
      call my_allocate(flind%ion%iarr_swpx,  ma1d)
      call my_allocate(flind%ion%iarr_swpy,  ma1d)
      call my_allocate(flind%ion%iarr_swpz,  ma1d)

      !\deprecated repeated magic integers (multifile: initneutral)
      flind%ion%iarr(1:4)      = [idni,imxi,imyi,imzi]
      flind%ion%iarr_swpx(1:4) = [idni,imxi,imyi,imzi]
      flind%ion%iarr_swpy(1:4) = [idni,imyi,imxi,imzi]
      flind%ion%iarr_swpz(1:4) = [idni,imzi,imyi,imxi]

#ifndef ISO
      flind%ion%iarr(5)      = ieni
      flind%ion%iarr_swpx(5) = ieni
      flind%ion%iarr_swpy(5) = ieni
      flind%ion%iarr_swpz(5) = ieni
      flind%ion%has_energy   = .true.

      flind%energ = flind%energ + 1
#endif /* !ISO */

      flind%ion%end    = flind%all
      flind%components = flind%components + 1
      flind%fluids     = flind%fluids + 1
      flind%ion%pos    = flind%components
      if (selfgrav_ion)  flind%fluids_sg = flind%fluids_sg + 1

      flind%ion%gam   = gamma_ion
      flind%ion%gam_1 = gamma_ion-1.0
      flind%ion%cs    = cs_iso_ion
      flind%ion%cs2   = cs_iso_ion**2
      flind%ion%tag   = ION
      flind%ion%is_selfgrav = selfgrav_ion

   end subroutine ionized_index

   subroutine cleanup_ionized

      implicit none

   end subroutine cleanup_ionized

end module initionized
