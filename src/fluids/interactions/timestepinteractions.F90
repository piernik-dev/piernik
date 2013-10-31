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
!>
!! \brief Module containing a routine to compute upper limit of %timestep due to fluids %interactions.
!<
module timestepinteractions
! pulled by ANY
   implicit none
   private
   public :: timestep_interactions

contains
!>
!! \brief Routine that computes upper limit of %timestep due to fluids %interactions.
!! \warning works only with neutrals and dust case !!!!
!! \deprecated BEWARE: works only with neu+dust!!!!
!! \todo check if subtraction of momenta is really the case and rewrite for all fluids
!<
   real function timestep_interactions(cg) result(dt)

      use constants,    only: small, I_ONE
      use fluidindex,   only: flind
      use func,         only: L2norm
      use grid_cont,    only: grid_container
      use interactions, only: collfaq, cfl_interact, has_interactions
      use mpi,          only: MPI_MIN, MPI_DOUBLE_PRECISION
      use mpisetup,     only: comm, mpi_err, FIRST, piernik_MPI_Bcast

      implicit none

      real :: dt_interact_proc        !< timestep due to %interactions for the current process (MPI block) only
      real :: dt_interact_all         !< timestep due to %interactions for all MPI blocks
      real :: val                     !< variable used to store the maximum value of relative momentum

      type(grid_container), pointer, intent(in) :: cg

      !    dt_interact_proc = 1.0 / (maxval(collfaq)+small) / maxval(cg%u(iarr_all_dn,:,:,:))

      if (has_interactions) then
         !       val = maxval (  sqrt( (cg%u(flind%dst%imx,:,:,:)-cg%u(flind%neu%imx,:,:,:))**2 + (cg%u(flind%dst%imy,:,:,:)-cg%u(flind%neu%imy,:,:,:))**2 + &
         !                             (cg%u(flind%dst%imz,:,:,:)-cg%u(flind%neu%imz,:,:,:))**2   ) * cg%u(flind%dst%idn,:,:,:) )
         val = maxval ( L2norm(cg%u(flind%dst%imx,:,:,:),cg%u(flind%dst%imy,:,:,:),cg%u(flind%dst%imz,:,:,:), &
                            &  cg%u(flind%neu%imx,:,:,:),cg%u(flind%neu%imy,:,:,:),cg%u(flind%neu%imz,:,:,:) ) * cg%u(flind%dst%idn,:,:,:) )
         dt_interact_proc = flind%neu%cs / (maxval(collfaq) * val + small)

         call MPI_Reduce(dt_interact_proc, dt_interact_all, I_ONE, MPI_DOUBLE_PRECISION, MPI_MIN, FIRST, comm, mpi_err)
         call piernik_MPI_Bcast(dt_interact_all)
         dt = cfl_interact*dt_interact_all
      else
         dt = huge(1.)
      endif

   end function timestep_interactions
!-------------------------------------------------------------------------------
end module timestepinteractions
