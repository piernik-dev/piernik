! $Id: initdsttral.F90 594 2009-01-21 12:33:44Z xarth $
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
!! \brief (MH) Initialization of the dust fluid
!!
!!
!!
!! In this module following namelist of parameters is specified:
!! \copydetails initdust::init_dust
!<

module initdust

  implicit none

    integer               :: idnd, imxd, imyd, imzd
    real                  :: dragc_gas_dust, taus, dalpha
    logical               :: selfgrav_dst

    integer, allocatable, dimension(:)  :: iarr_dst
    integer, allocatable, dimension(:)  :: iarr_dst_swpx, iarr_dst_swpy, iarr_dst_swpz

  contains

!>
!! \brief Routine to set parameter values from namelist FLUID_DUST
!!
!! \n \n
!! @b FLUID_DUST
!! \n \n
!! <table border="+1">
!! <tr><td width="150pt"><b>parameter</b></td><td width="135pt"><b>default value</b></td><td width="200pt"><b>possible values</b></td><td width="315pt"> <b>description</b></td></tr>
!! <tr><td>dragc_gas_dust</td><td>1.0  </td><td>real value</td><td>\copydoc initdust::dragc_gas_dust</td></tr>
!! <tr><td>dalpha        </td><td>1.0  </td><td>real value</td><td>\copydoc initdust::dalpha        </td></tr>
!! </table>
!! \n \n
!<
  subroutine init_dust

    use mpisetup,       only: proc, rbuff, lbuff, MPI_LOGICAL, MPI_DOUBLE_PRECISION, buffer_dim, comm, ierr
    use dataio_public,  only: par_file, ierrh, namelist_errh
    use func,           only: compare_namelist

    implicit none


    namelist /FLUID_DUST/ dragc_gas_dust, dalpha, selfgrav_dst

    dragc_gas_dust  = 1.0
    dalpha = 1.0
    selfgrav_dst = .false.

    if (proc == 0) then
       diff_nml(FLUID_DUST)

       rbuff(1)   = dragc_gas_dust
       rbuff(2)   = dalpha

       lbuff(1)   = selfgrav_dst
    endif

    call MPI_Bcast(rbuff,    buffer_dim, MPI_DOUBLE_PRECISION, 0, comm, ierr)
    call MPI_Bcast(lbuff,    buffer_dim, MPI_LOGICAL,          0, comm, ierr)

    if (proc /= 0) then

      selfgrav_dst    = lbuff(1)

      dragc_gas_dust  = rbuff(1)
      dalpha          = rbuff(2)

    endif
    taus = 1. / dragc_gas_dust

  end subroutine init_dust


  subroutine dust_index(nvar,nvar_dst)

    implicit none
    integer :: nvar, nvar_dst


      idnd = nvar + 1
      imxd = nvar + 2
      imyd = nvar + 3
      imzd = nvar + 4

      nvar_dst      = 4
      nvar          = imzd

      allocate(iarr_dst(nvar_dst),iarr_dst_swpx(nvar_dst), iarr_dst_swpy(nvar_dst), iarr_dst_swpz(nvar_dst))

      iarr_dst      = [idnd,imxd,imyd,imzd]
      iarr_dst_swpx = [idnd,imxd,imyd,imzd]
      iarr_dst_swpy = [idnd,imyd,imxd,imzd]
      iarr_dst_swpz = [idnd,imzd,imyd,imxd]


   end subroutine dust_index

end module initdust
