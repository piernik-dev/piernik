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
module timestepthermal
! pulled by THERM
   implicit none
   private
   public :: timestep_thermal

contains
!>
!! \brief Routine that computes upper limit of %timestep due to fluids %interactions.
!! \warning works only with neutrals and dust case !!!!
!! \deprecated BEWARE: works only with neu+dust!!!!
!! \todo check if subtraction of momenta is really the case and rewrite for all fluids
!<
   real function timestep_thermal(cg) result(dt)
      use grid_cont,    only: grid_container
      use constants,    only: small, I_ONE
      use fluidindex,   only: flind
      use func,         only: L2norm
      use grid_cont,    only: grid_container
      use interactions, only: collfaq, cfl_interact, has_interactions
      use mpi,          only: MPI_MIN, MPI_DOUBLE_PRECISION
      use mpisetup,     only: comm, mpi_err, FIRST, piernik_MPI_Bcast
      use thermal,      only: maxdeint, cfl_coolheat, thermal_active

      implicit none

      real :: dt_coolheat_proc        !< timestep due to %interactions for the current process (MPI block) only
      real :: dt_coolheat_all         !< timestep due to %interactions for all MPI blocks
!      real :: val                     !< variable used to store the maximum value of relative momentum
      real :: mxdeint

      type(grid_container), pointer, intent(in) :: cg

      !    dt_interact_proc = 1.0 / (maxval(collfaq)+small) / maxval(cg%u(iarr_all_dn,:,:,:))

      if (thermal_active) then
         !       val = maxval (  sqrt( (cg%u(flind%dst%imx,:,:,:)-cg%u(flind%neu%imx,:,:,:))**2 + (cg%u(flind%dst%imy,:,:,:)-cg%u(flind%neu%imy,:,:,:))**2 + &
         !                             (cg%u(flind%dst%imz,:,:,:)-cg%u(flind%neu%imz,:,:,:))**2   ) * cg%u(flind%dst%idn,:,:,:) )
         !val = maxval ( L2norm(cg%u(flind%dst%imx,:,:,:),cg%u(flind%dst%imy,:,:,:),cg%u(flind%dst%imz,:,:,:), &
         !                   &  cg%u(flind%neu%imx,:,:,:),cg%u(flind%neu%imy,:,:,:),cg%u(flind%neu%imz,:,:,:) ) * cg%u(flind%dst%idn,:,:,:) )
         call maxdeint(cg, mxdeint)
         dt_coolheat_proc = ABS(1./(mxdeint+small))
         call MPI_Reduce(dt_coolheat_proc, dt_coolheat_all, I_ONE, MPI_DOUBLE_PRECISION, MPI_MIN, FIRST, comm, mpi_err)
         call piernik_MPI_Bcast(dt_coolheat_all)
         dt = cfl_coolheat*dt_coolheat_all
      else
         dt = huge(1.)
      endif
    end function timestep_thermal 

#ifdef IMPROPER
!-------------------------------------------------------------------------------
   subroutine timestep_coolheat 
  
    implicit none
    
    
    real eint_src_min_sweep,eint_src_max_sweep
    real, dimension(nz) :: eint, dens, temp, eint_src
    integer loc_dt_cool3(1), loc_dt_heat3(1)
    integer i,j

    real dt_coolheat_all, dt_coolheat_proc

    
    eint_src_max = 0.0
    eint_src_min = 0.0

    do j=nb+1,nb+nyb
      do i=nb+1,nb+nxb
        dens = u(1,i,j,:)
        eint = u(5,i,j,:)-sum(u(2:4,i,j,:)**2,1)/u(1,i,j,:)/2-sum(b(:,i,j,:)**2,1)/2
	   
        call cool_heat ('zsweep', i, j, nz, dens, eint,  eint_src)

        eint_src_min_sweep = MINVAL(eint_src(nb+1:nb+nzb)/eint(nb+1:nb+nzb))
        eint_src_max_sweep = MAXVAL(eint_src(nb+1:nb+nzb)/eint(nb+1:nb+nzb))


        if(eint_src_min_sweep .lt. eint_src_min) then	
	  eint_src_min = eint_src_min_sweep
	  loc_dt_cool(1:2) = (/i,j/)
	  loc_dt_cool3 = MINLOC(eint_src(nb+1:nb+nzb) / eint(nb+1:nb+nzb)) + nb
	  loc_dt_cool(3) = loc_dt_cool3(1)	     
	endif					       

        if(eint_src_max_sweep .gt. eint_src_max) then	
	  eint_src_max = eint_src_max_sweep
	  loc_dt_heat(1:2) = (/i,j/)
	  loc_dt_heat3 = MAXLOC(eint_src(nb+1:nb+nzb) / eint(nb+1:nb+nzb)) + nb
	  loc_dt_heat(3) = loc_dt_heat3(1)	     
	endif					       

      end do
    end do
    
    dt_heat = cfl_coolheat*abs(1./(eint_src_max+small))
    dt_cool = cfl_coolheat*abs(1./(eint_src_min+small))
    
    dt_coolheat_proc = MIN ( dt_cool, dt_heat)
    
 
    call MPI_REDUCE(dt_coolheat_proc, dt_coolheat_all, 1, MPI_DOUBLE_PRECISION, MPI_MIN, 0, comm, ierr)
    call MPI_BCAST(dt_coolheat_all, 1, MPI_DOUBLE_PRECISION, 0, comm, ierr)
    
    dt_coolheat = dt_coolheat_all
 
  end subroutine timestep_coolheat
#endif /*IMPROPER*/
end module timestepthermal
