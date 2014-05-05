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
module initproblem

   use constants, only: dsetnamelen

   implicit none

   private
   public  :: read_problem_par, problem_initial_conditions, problem_pointers

   integer(kind=4) :: n_sn
   real            :: d0, p0, bx0, by0, bz0, Eexpl, x0, y0, z0, r0, dt_sn, r, t_sn

   character(len=dsetnamelen), parameter :: ngp_n = "ngp"
   character(len=dsetnamelen), parameter :: cic_n = "cic"
   character(len=dsetnamelen), parameter :: tsc_n = "tsc"

   namelist /PROBLEM_CONTROL/ d0, p0, bx0, by0, bz0, Eexpl, x0, y0, z0, r0, n_sn, dt_sn

contains
!-----------------------------------------------------------------------------
   subroutine problem_pointers

      use dataio_user, only: user_vars_hdf5

      implicit none

      user_vars_hdf5 => map_vars_hdf5

   end subroutine problem_pointers
!-----------------------------------------------------------------------------
   subroutine read_problem_par

      implicit none

   end subroutine read_problem_par
!-----------------------------------------------------------------------------
   subroutine map_vars_hdf5(var, tab, ierrh, cg)

      use grid_cont,        only: grid_container
      use particle_pub,     only: pset
      use named_array_list, only: qna

      implicit none

      character(len=*),               intent(in)    :: var
      real(kind=4), dimension(:,:,:), intent(inout) :: tab
      integer,                        intent(inout) :: ierrh
      type(grid_container), pointer,  intent(in)    :: cg

      ierrh = 0
      select case (trim(var))
         case ("ngp")
            call pset%map_ngp(qna%ind(ngp_n), 1.0)
            tab(:,:,:) = real(cg%q(qna%ind(ngp_n))%span(cg%ijkse), kind=4)
         case ("cic")
            call pset%map_cic(qna%ind(cic_n), 1.0)
            tab(:,:,:) = real(cg%q(qna%ind(cic_n))%span(cg%ijkse), kind=4)
         case ("tsc")
            call pset%map_tsc(qna%ind(tsc_n), 1.0)
            tab(:,:,:) = real(cg%q(qna%ind(tsc_n))%span(cg%ijkse), kind=4)
         case default
            ierrh = -1
      end select

      if (.true. .or. cg%grid_id > 0) return ! suppress compiler warnings

   end subroutine map_vars_hdf5
!-----------------------------------------------------------------------------
   subroutine problem_initial_conditions

      use cg_leaves,        only: leaves
      use cg_list,          only: cg_list_element
      use cg_list_global,   only: all_cg
      use constants,        only: xdim, ydim, zdim, LO, HI
      use fluidindex,       only: flind
      use named_array_list, only: qna
      use particle_pub,     only: pset

      implicit none

      integer                         :: i, j, k, p
      type(cg_list_element),  pointer :: cgl
      integer :: ncg_i, cic_i, tsc_i
      logical, save :: frun = .true.

      real, parameter, dimension(185) :: px = [ &
         -4.264752, -4.264106, -4.263515, -4.263066, -4.252018, -4.251456, -4.281534, -4.260255, -4.249320, -4.229333, -4.250500, &
         -4.269165, -4.248252, -4.257781, -4.172945, -4.122122, -3.979266, -3.856818, -3.683379, -3.509434, -3.478597, -3.447704, &
         -3.416867, -3.396262, -3.365257, -3.395475, -3.507354, -3.670280, -3.823171, -4.016962, -4.118973, -3.148528, -3.148528, &
         -3.148528, -3.158479, -3.198705, -3.198058, -3.197299, -3.166406, -3.084521, -2.992489, -2.941525, -2.870182, -3.159744, &
         -3.200251, -3.047304, -3.057676, -2.249595, -2.198547, -2.198547, -2.045488, -1.892597, -1.750079, -1.679044, -1.689698, &
         -1.832779, -2.006162, -2.138588, -2.240403, -2.311437, -2.331086, -2.259209, -2.126108, -1.993175, -1.819454, -1.748027, &
         -1.625607, -1.503272, -1.098205, -1.026834, -0.873774, -0.812073, -0.780702, -0.769992, -0.800042, -0.809402, -0.799030, &
         -0.788882, -0.788882, -0.728361, -0.667784, -0.607094, -0.536088, -0.454906, -0.332711, -0.312359, -0.312359, -0.240622, &
         -0.148674, -0.016163,  0.270222,  0.260131,  0.260131,  0.300356,  0.442847,  0.544942,  0.596412,  0.617748,  0.638746, &
          0.639308,  0.639673,  0.619603,  0.619603,  0.619603,  0.588316,  0.617888,  0.678522,  0.739183,  0.820590,  0.891877, &
          0.932693,  0.932693,  1.024473,  1.146976,  1.208313,  1.270352,  1.290985,  1.322102,  1.373966,  1.384535,  1.405112, &
          1.456216,  1.456216,  1.558227,  2.045094,  2.045094,  2.045094,  2.025277,  2.046247,  2.067189,  2.088131,  2.098728, &
          2.180557,  2.303004,  2.415107,  2.033682,  2.013358,  2.013358,  2.064462,  2.064462,  2.115341,  2.094765,  2.867793, &
          2.867793,  2.867793,  2.878531,  2.889100,  2.869339,  2.869929,  2.870435,  2.871110,  2.881567,  2.881960,  2.882578, &
          2.893204,  2.903774,  2.904392,  2.904814,  2.905348,  2.915692,  2.915692,  2.915833,  2.915833,  2.915861,  2.915861, &
          2.884434,  2.894300,  2.924378,  2.985152,  3.045983,  3.127052,  3.248966,  3.401688,  3.513847,  3.697828,  3.107459, &
          3.117579,  3.117579,  3.373493,  3.394435,  3.507241,  3.711686,  3.874809,  3.945956,  3.230272 ]

      real, parameter, dimension(185) :: py = [ &
         3.523095, 3.333015, 3.159463, 3.027233, 2.779336, 2.614048, 2.456924, 2.200796, 1.985956, 2.109989, 2.333060, 1.820601, &
         1.671910, 1.473532, 3.531663, 3.589682, 3.590155, 3.590560, 3.599397, 3.459476, 3.393463, 3.310920, 3.244906, 3.187123, &
         3.071523, 2.955721, 2.847913, 2.748201, 2.698109, 2.672675, 2.664073, 2.353232, 2.353232, 2.353232, 2.278818, 2.105132, &
         1.915051, 1.691913, 1.609371, 1.535261, 1.477714, 1.494412, 1.519441, 2.650715, 2.559672, 2.593236, 2.642788, 2.066947, &
         2.058851, 2.058851, 2.059357, 2.109449, 2.209094, 2.325031, 2.457228, 2.522870, 2.497504, 2.430950, 2.364498, 2.248561, &
         2.025355, 1.893362, 1.761570, 1.679365, 1.605559, 1.605795, 1.614464, 1.647927, 2.558357, 2.575121, 2.575627, 2.435336, &
         2.212299, 2.063574, 1.898185, 1.650220, 1.600668, 1.617230, 1.617230, 1.824042, 2.014325, 2.171550, 2.295752, 2.428252, &
         2.503036, 2.519632, 2.519632, 2.428960, 2.396206, 2.437967, 2.240566, 2.207475, 2.207475, 2.381162, 2.489071, 2.472879, &
         2.340818, 2.068161, 1.894677, 1.729390, 1.621953, 1.522713, 1.522713, 1.522713, 1.720957, 2.026840, 2.200594, 2.366083, &
         2.432468, 2.474026, 2.474161, 2.474161, 2.490994, 2.474870, 2.442014, 2.202550, 2.136503, 1.987845, 1.740083, 1.632680, &
         1.583161, 1.558536, 1.558536, 1.567138, 2.428252, 2.428252, 2.428252, 2.254633, 2.089413, 1.932457, 1.775501, 1.659834, &
         1.602253, 1.602658, 1.644351, 2.783586, 2.758725, 2.758725, 2.734101, 2.734101, 2.775591, 2.825110, 3.554938, 3.554938, &
         3.554938, 3.397949, 3.290546, 3.100398, 2.926846, 2.778088, 2.579743, 2.505397, 2.389696, 2.207880, 2.083948, 1.976545, &
         1.794729, 1.670763, 1.513740, 1.472452, 1.472452, 1.431130, 1.431130, 1.422866, 1.422866, 1.662431, 1.761638, 1.918762, &
         2.051194, 2.167098, 2.332655, 2.490083, 2.589761, 2.614925, 2.524624, 2.092921, 2.117748, 2.117748, 1.878924, 1.721969, &
         1.557052, 1.450290, 1.492152, 1.575031, 1.985889 ]

      if (frun) then
         call all_cg%reg_var(ngp_n)
         call all_cg%reg_var(cic_n)
         call all_cg%reg_var(tsc_n)

         do i = lbound(px, 1), ubound(px, 1)
            call pset%add(1.0, [px(i), py(i), 0.0], [0.0, 0.0, 0.0])
         enddo

         frun = .false.
      endif

      ncg_i = qna%ind(ngp_n)
      cic_i = qna%ind(cic_n)
      tsc_i = qna%ind(tsc_n)

      do p = lbound(flind%all_fluids, dim=1), ubound(flind%all_fluids, dim=1)
         cgl => leaves%first
         do while (associated(cgl))
            associate(cg => cgl%cg)
               do k = cg%lhn(zdim,LO), cg%lhn(zdim,HI)
                  do j = cg%lhn(ydim,LO), cg%lhn(ydim,HI)
                     do i = cg%lhn(xdim,LO), cg%lhn(xdim,HI)
                        associate(fl => flind%all_fluids(p)%fl)
                           cg%u(fl%idn,i,j,k) = 1.0
                           cg%u(fl%imx,i,j,k) = 0.0
                           cg%u(fl%imy,i,j,k) = 0.0
                           cg%u(fl%imz,i,j,k) = 0.0
                        end associate
                        cg%q(ncg_i)%arr(i, j, k) = 0.0
                        cg%q(cic_i)%arr(i, j, k) = 0.0
                        cg%q(tsc_i)%arr(i, j, k) = 0.0
                     enddo
                  enddo
               enddo
            end associate
            cgl => cgl%nxt
         enddo
      enddo


   end subroutine problem_initial_conditions
!-----------------------------------------------------------------------------
end module initproblem
