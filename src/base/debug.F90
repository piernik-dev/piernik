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
!! \brief Module providing few parameter for debugging and developing the code.
!!
!! \details Keep this module as free of dependences as possible to allow its use early after start and eliminate risk of cyclic dependences.
!<
module piernikdebug
! pulled by DEBUG

   use constants, only: cbuff_len

   implicit none

   private
   public :: init_piernikdebug, has_const_dt, constant_dt, aux_R, aux_I, aux_L, aux_S

   ! Auxiliary input parameters for debugging, quick tweaks and tests of new features.
   ! Their purpose is to avoid messing up existing namelists until it becomes clear that certain parameter is really useful.
   ! There is no reason to give them protected attribute.
   integer, parameter                        :: naux = 5 !< number of auxiliary variables of each kind
   real,                     dimension(naux) :: aux_R    !< real auxiliary parameter
   integer(kind=4),          dimension(naux) :: aux_I    !< integer auxiliary parameter
   logical,                  dimension(naux) :: aux_L    !< boolean auxiliary parameter
   character(len=cbuff_len), dimension(naux) :: aux_S    !< string auxiliary parameter

   real,    protected :: constant_dt               !< value of timestep regardless of fluid state
   logical, protected :: has_const_dt              !< true if piernikdebug::constant_dt > 0

   namelist /PIERNIK_DEBUG/ constant_dt, aux_R, aux_I, aux_L, aux_S

contains

!>
!! \brief Routine to set debug auxiliary parameters
!!
!! \n \n
!! @b PIERNIK_DEBUG
!! \n \n
!! <table border="+1">
!!   <tr><td width="150pt"><b>parameter</b></td><td width="135pt"><b>default value</b></td><td width="200pt"><b>possible values</b></td><td width="315pt"> <b>description</b></td></tr>
!!   <tr><td>constant_dt      </td><td>0.0    </td><td>real value             </td><td>\copydoc piernikdebug::constant_dt      </td></tr>
!!   <tr><td>aux_R            </td><td>0.0    </td><td>5-element real array   </td><td>\copydoc piernikdebug::aux_r            </td></tr>
!!   <tr><td>aux_I            </td><td>0      </td><td>5-element integer array</td><td>\copydoc piernikdebug::aux_i            </td></tr>
!!   <tr><td>aux_L            </td><td>.false.</td><td>5-element logical array</td><td>\copydoc piernikdebug::aux_l            </td></tr>
!!   <tr><td>aux_S            </td><td>""     </td><td>5-element string array </td><td>\copydoc piernikdebug::aux_s            </td></tr>
!! </table>
!! The list is active while DEBUG is defined.
!! \n \n
!<

   subroutine init_piernikdebug

      use constants,             only: PIERNIK_INIT_MPI
      use dataio_pub,            only: par_file, ierrh, namelist_errh, compare_namelist, cmdl_nml, lun  ! QA_WARN required for diff_nml
      use dataio_pub,            only: code_progress, die
      use mpisetup,              only: master, slave, rbuff, lbuff, cbuff, ibuff, buffer_dim, piernik_MPI_Bcast

      implicit none

      if (code_progress < PIERNIK_INIT_MPI) call die("[debug:init_piernikdebug] MPI not initialized.")

      constant_dt = 0.0
      aux_R(:) = 0.
      aux_I(:) = 0
      aux_L(:) = .false.
      aux_S(:) = ''

      if (master) then
         diff_nml(PIERNIK_DEBUG)

         rbuff(1) = constant_dt

         rbuff(buffer_dim-naux+1:buffer_dim) = aux_R(:)
         ibuff(buffer_dim-naux+1:buffer_dim) = aux_I(:)
         lbuff(buffer_dim-naux+1:buffer_dim) = aux_L(:)
         cbuff(buffer_dim-naux+1:buffer_dim) = aux_S(:)

      endif

      call piernik_MPI_Bcast(cbuff, cbuff_len)
      call piernik_MPI_Bcast(ibuff)
      call piernik_MPI_Bcast(rbuff)
      call piernik_MPI_Bcast(lbuff)

      if (slave) then
         constant_dt       = rbuff(1)

         aux_R(:) = rbuff(buffer_dim-naux+1:buffer_dim)
         aux_I(:) = ibuff(buffer_dim-naux+1:buffer_dim)
         aux_L(:) = lbuff(buffer_dim-naux+1:buffer_dim)
         aux_S(:) = cbuff(buffer_dim-naux+1:buffer_dim)

      endif

      has_const_dt = (constant_dt > 0.0)

   end subroutine init_piernikdebug

end module piernikdebug
