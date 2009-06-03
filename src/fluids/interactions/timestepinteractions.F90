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

module timestepinteractions

   real :: dt_interact

 contains

  subroutine timestep_interactions
    use mpisetup
    use arrays,         only : u,b
    use constants,      only : small
    use fluidindex,     only : iarr_all_dn
    use initdust,       only : imxd,imyd,imzd,idnd
    use initneutral,    only : imxn,imyn,imzn,cs_iso_neu
    use interactions,   only : collfaq, cfl_interact

    implicit none

    real :: dt_interact_proc, dt_interact_all
    real :: val


!    dt_interact_proc = 1.0 / (maxval(collfaq)+small) / maxval(u(iarr_all_dn,:,:,:))

    !!!!BEWARE: works only with neu+dust!!!!

    val = maxval (  sqrt( (u(imxd,:,:,:)-u(imxn,:,:,:))**2 + &
                    (u(imyd,:,:,:)-u(imyn,:,:,:))**2 + &  
                    (u(imzd,:,:,:)-u(imzn,:,:,:))**2   ) * u(idnd,:,:,:) )
    dt_interact_proc = cs_iso_neu / (maxval(collfaq) * val)


    call MPI_REDUCE(dt_interact_proc, dt_interact_all, 1, MPI_DOUBLE_PRECISION, MPI_MIN, 0, comm, ierr)
    call MPI_BCAST(dt_interact_all, 1, MPI_DOUBLE_PRECISION, 0, comm, ierr)
    dt_interact = cfl_interact*dt_interact_all

  end subroutine timestep_interactions

!-------------------------------------------------------------------------------
end module timestepinteractions

