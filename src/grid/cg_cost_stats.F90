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

!> \brief This module provides a structure for calculating basic statistics of computational costs of a bunch of cg

module cg_cost_stats

   use cg_cost_data, only: cg_cost_data_t

   implicit none

   private
   public :: cg_stats_t, stat_labels, I_MIN, I_MAX, I_AVG, I_SIGMA, I_SUM, I_SUM2

   type :: cg_stats_t
      type(cg_cost_data_t), private :: min     !< element-wise minimum
      type(cg_cost_data_t), private :: max     !< element-wise maximum
      type(cg_cost_data_t), private :: w_sum   !< sum of elements
      type(cg_cost_data_t), private :: w_sum2  !< sum of squares
      ! type(cg_cost_data_t), private :: bias    !< bias for improved accuracy of the mean
      integer,              private :: n       !< number of elements
   contains
      procedure :: reset        !< initialize
      procedure :: add          !< add sample
      procedure :: get_minimum  !< return maxima
      procedure :: get_maximum  !< return minima
      procedure :: get_average  !< return average
      procedure :: get_sigma    !< return standard deviation
      procedure :: get          !< return minimum, maximum mean, sigma and sum for selected accumulator
      procedure :: get_sum      !< return total accumulated time
   end type cg_stats_t

   enum, bind(C)
      enumerator :: I_MIN = 1, I_MAX, I_AVG, I_SIGMA, I_SUM, I_SUM2
   end enum

   character(len=*), dimension(I_MIN:I_SUM2), parameter :: stat_labels = &
        [ "minimum  ", &
        & "maximum  ", &
        & "average  ", &
        & "deviation", &
        & "sum      ", &
        & "sum^2    " ]

contains

!> \brief Initialize

   subroutine reset(this)

      implicit none

      class(cg_stats_t), intent(inout) :: this

      this%min%wtime =  huge(1.)
      this%max%wtime = -huge(1.)
      this%w_sum%wtime  = 0.
      this%w_sum2%wtime = 0.
      this%n = 0

   end subroutine reset

!> \brief Add sample

   subroutine add(this, data)

      implicit none

      class(cg_stats_t),     intent(inout) :: this
      class(cg_cost_data_t), intent(in)    :: data

      if (any(data%wtime > 0.)) then
         this%min%wtime = min(this%min%wtime, data%wtime)
         this%max%wtime = max(this%max%wtime, data%wtime)
         this%w_sum%wtime = this%w_sum%wtime + data%wtime
         this%w_sum2%wtime = this%w_sum2%wtime + data%wtime**2
         this%n = this%n + 1
      endif

   end subroutine add

!> \brief Return minima
   type(cg_cost_data_t) function get_minimum(this)

      implicit none

      class(cg_stats_t), intent(in) :: this

      get_minimum = this%min

   end function get_minimum

!> \brief Return maxima

   type(cg_cost_data_t) function get_maximum(this)

      implicit none

      class(cg_stats_t), intent(in) :: this

      get_maximum = this%max

   end function get_maximum

!> \brief Return average

   type(cg_cost_data_t) function get_average(this)

      implicit none

      class(cg_stats_t), intent(in) :: this

      if (this%n /= 0) then
         get_average%wtime = this%w_sum%wtime / this%n
      else
         get_average%wtime = 0.  ! this should be safe value
      endif

   end function get_average

!> \brief Return standard deviation

   type(cg_cost_data_t) function get_sigma(this)

      implicit none

      class(cg_stats_t), intent(in) :: this

      if (this%n /= 0) then
         get_sigma%wtime = sqrt(this%w_sum2%wtime / this%n - (this%w_sum%wtime / this%n)**2)
      else
         get_sigma%wtime = 0.  ! this should be safe value
      endif

   end function get_sigma

!> \brief Return minimum, maximum mean and sigma for selected accumulator

   function get(this, ind)

      use dataio_pub, only: die

      implicit none

      class(cg_stats_t), intent(in) :: this
      integer(kind=4),   intent(in) :: ind

      real, dimension(lbound(stat_labels,1):ubound(stat_labels,1)) :: get

      if (ind >= lbound(this%min%wtime, 1) .and. ind <= ubound(this%min%wtime, 1)) then
         if (this%n /= 0) then
            ! Beware: formulas repeated from this%get_average and this%get_sigma
            get = [ this%min%wtime(ind), &
                 &  this%max%wtime(ind), &
                 &  this%w_sum%wtime(ind) / this%n, &
                 &  sqrt(this%w_sum2%wtime(ind) / this%n - (this%w_sum%wtime(ind) / this%n)**2), &
                 &  this%w_sum%wtime(ind) , &
                 &  this%w_sum2%wtime(ind) ]
         else
            get = 0.
         endif
      else
         call die("[cg_cost_stats:get] invalid index")
         get = 0.
      endif

   end function get

!> \brief Return total accumulated time

   real function get_sum(this)

      implicit none

      class(cg_stats_t), intent(in) :: this

      get_sum = sum(this%w_sum%wtime)

   end function get_sum

end module cg_cost_stats
