! $Id$
!
! PIERNIK Code Copyright (C) 2006-2012 Michal Hanasz
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
module gdf

   implicit none
   private
   public :: gdf_create_root_datasets, gdf_create_simulation_parameters, gdf_create_format_stamp, gdf_create_field_types, gdf_field_type, fmax

   integer, parameter :: fmax = 60

   type :: gdf_field_type
      real(kind=8) :: f2cgs
      integer(kind=8) :: stag
      character(len=fmax) :: fu
      character(len=fmax) :: fn
   end type gdf_field_type


contains

   subroutine gdf_create_root_datasets(file, cg_all_n_b, cg_all_off, cg_all_rl, cg_all_parents, cg_all_particles)

      use helpers_hdf5, only: create_dataset
      use hdf5,         only: HID_T

      implicit none

      integer(HID_T), intent(in) :: file !> File identifier
      integer(kind=4),  dimension(:,:), pointer, intent(in) :: cg_all_n_b       !> sizes of all cg
      integer(kind=4),  dimension(:,:), pointer, intent(in) :: cg_all_rl        !> refinement levels of all cgs
      integer(kind=8),  dimension(:,:), pointer, intent(in) :: cg_all_off       !> offsets of all cgs
      integer(kind=8),  dimension(:), pointer, intent(in) :: cg_all_parents     !> parents IDs of all cgs
      integer(kind=4),  dimension(:,:), pointer, intent(in) :: cg_all_particles !> particles counts in all cgs

      call create_dataset(file, 'grid_dimensions', cg_all_n_b)
      call create_dataset(file, 'grid_left_index', cg_all_off)
      call create_dataset(file, 'grid_level', cg_all_rl)
      call create_dataset(file, 'grid_parent_id', cg_all_parents)
      call create_dataset(file, 'grid_particle_count', cg_all_particles)
   end subroutine gdf_create_root_datasets

   subroutine gdf_create_simulation_parameters(file_id, t, n_d, nb, edge, domain_dump)

      use hdf5,         only: HID_T, h5gcreate_f, h5gclose_f
      use helpers_hdf5, only: create_attribute

      implicit none

      integer(HID_T), intent(in) :: file_id
      real, intent(in)              :: t    !< time in physical units
      integer(kind=4), dimension(:) :: n_d  !< number of grid cells in physical domain in x-, y- and z-direction
      integer(kind=4)               :: nb   !< number of boundary cells surrounding the physical domain, same for all directions
      real, dimension(:, :)         :: edge !< physical domain boundary positions
      character(len=*), intent(in)  :: domain_dump

      integer(HID_T) :: g_id
      integer(kind=4) :: error
      integer(kind=8), dimension(:), pointer :: ibuf
      integer, parameter :: uniqid_len = 12
      character(len=uniqid_len), target :: uniq_id = "ala123"
      character(len=uniqid_len), pointer :: p_str
      integer(kind=4) :: ndims

      ndims = size(n_d, kind=4)

      call h5gcreate_f(file_id, 'simulation_parameters', g_id, error)

      allocate(ibuf(1))
      ibuf = int(2,8)    !!!! constant for now
      call create_attribute(g_id, 'refine_by', ibuf)

      ibuf = int(ndims,8)
      call create_attribute(g_id, 'dimensionality', ibuf)

      ibuf = int(0,8)
      call create_attribute(g_id, 'cosmological_simulation', ibuf)

      select case (trim(domain_dump))
         case ('phys_domain')
            ibuf = int(0,8)
         case ('full_domain')
            ibuf = int(nb,8)
      end select
      call create_attribute(g_id, 'num_ghost_zones', ibuf)
      deallocate(ibuf)

      allocate(ibuf(ndims))
      select case (trim(domain_dump))
         case ('phys_domain')
            ibuf = n_d            !???
         case ('full_domain')
            ibuf = n_d + 2*nb !???
      end select
      call create_attribute(g_id, 'domain_dimensions', ibuf)
      deallocate(ibuf)

      call create_attribute(g_id, 'domain_left_edge', edge(:, lbound(edge, dim=2)))
      call create_attribute(g_id, 'domain_right_edge', edge(:, ubound(edge, dim=2)))
      call create_attribute(g_id, 'current_time', [t])
      call create_attribute(g_id, 'field_ordering', [int(1, kind=4)])
      p_str => uniq_id
      call create_attribute(g_id, 'unique_identifier', p_str)
      call create_attribute(g_id, 'boundary_conditions', int([0,0,0,0,0,0], kind=4))
      call h5gclose_f(g_id, error)

   end subroutine gdf_create_simulation_parameters

   subroutine gdf_create_format_stamp(file)

      use hdf5,   only: HID_T, h5gcreate_f, h5gclose_f
      use h5lt,   only: h5ltset_attribute_string_f

      implicit none

      integer(HID_T), intent(in) :: file
      integer(HID_T) :: g_id
      integer(kind=4)                :: error

      character(len=*), parameter :: gname = 'gridded_data_format'
      character(len=*), parameter :: gname2 = '/gridded_data_format'

      call h5gcreate_f(file, gname, g_id, error)
      call h5ltset_attribute_string_f(g_id, gname2, 'data_software', 'piernik', error )
      call h5ltset_attribute_string_f(g_id, gname2, 'data_software_version', '1.0', error )
      call h5gclose_f(g_id, error)

   end subroutine gdf_create_format_stamp

   subroutine gdf_create_field_types(filename, o_func)

      use hdf5,    only: HID_T, h5gcreate_f, h5gclose_f, h5fopen_f, h5fclose_f, H5F_ACC_RDWR_F, h5open_f, h5close_f

      implicit none

      character(len=*), intent(in) :: filename
      interface
         subroutine o_func(group_id)
            use hdf5, only: HID_T
            implicit none
            integer(HID_T), intent(in) :: group_id
         end subroutine o_func
      end interface

      integer(HID_T) :: g_id, file_id
      integer(kind=4) :: error
      character(len=*), parameter :: gname = 'field_types'

      call h5open_f(error)
      call h5fopen_f(filename, H5F_ACC_RDWR_F, file_id, error)
      call h5gcreate_f(file_id, gname, g_id, error)
      call o_func(g_id)
      call h5gclose_f(g_id, error)
      call h5fclose_f(file_id, error)
      call h5close_f(error)

   end subroutine gdf_create_field_types

end module gdf
