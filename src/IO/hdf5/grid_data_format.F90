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
!! \brief Implementation of Grid Data Format
!<
module gdf

   implicit none
   private
   public :: gdf_create_root_datasets, gdf_create_simulation_parameters, gdf_create_format_stamp, gdf_create_field_types, gdf_field_type, fmax
   public :: gdf_parameters_T, gdf_root_datasets_T

   integer, parameter :: fmax = 60

   type :: gdf_field_type
      real(kind=8)        :: f2cgs
      integer(kind=8)     :: stag
      character(len=fmax) :: fu
      character(len=fmax) :: fn
   end type gdf_field_type

   integer, parameter :: uniqid_len = 12

   type :: gdf_parameters_T
      integer(kind=8), dimension(:), pointer :: refine_by                !< relative global refinement
      integer(kind=8), dimension(:), pointer :: dimensionality           !< 1-, 2- or 3-D data
      integer(kind=8), dimension(:), pointer :: cosmological_simulation  !< 0 or 1 == True or False
      integer(kind=8), dimension(:), pointer :: num_ghost_zones          !< number of ghost zones
      integer(kind=4), dimension(:), pointer :: field_ordering           !< integer: 0 for C, 1 for Fortran
      integer(kind=8), dimension(:), pointer :: domain_dimensions        !< dimensions in the top grid
      real(kind=8),    dimension(:), pointer :: domain_left_edge         !< the left edge of the domain, in cm
      real(kind=8),    dimension(:), pointer :: domain_right_edge        !< the right edge of the domain, in cm
      real(kind=8),    dimension(:), pointer :: current_time             !< current time in simulation, in seconds, from "start" of simulation
      character(len=uniqid_len),     pointer :: unique_identifier        !< regarded as a string, but can be anything
      integer(kind=4), dimension(:), pointer :: boundary_conditions      !< 0 for periodic, 1 for mirrored, 2 for outflow.  Needs one for each face
                                                                         !< of the cube.  Any past the dimensionality should be set to -1.  The order of specification
                                                                         !< goes left in 0th dimension, right in 0th dimension, left in 1st dimension,
                                                                         !< right in 1st dimensions, left in 2nd dimension, right in 2nd dimension.
                                                                         !< Note also that yt does not currently support non-periodic boundary conditions, and that
                                                                         !< the assumption of periodicity shows up primarily in plots and covering grids.
   contains
      procedure :: init => gdf_parameters_init
      procedure :: cleanup => gdf_parameters_finalize
   end type gdf_parameters_T

   type :: gdf_root_datasets_T
      integer(kind=8),  dimension(:),   pointer :: grid_parent_id       !< optional, may only reference a single parent
      integer(kind=8),  dimension(:,:), pointer :: grid_left_index      !< global, relative to current level, and only the active region
      integer(kind=4),  dimension(:,:), pointer :: grid_dimensions      !< only the active regions
      integer(kind=4),  dimension(:,:), pointer :: grid_level           !< level, indexed by zero
      integer(kind=4),  dimension(:,:), pointer :: grid_particle_count  !< total number of particles.  (May change in subsequent versions.)
   contains
      procedure :: gdf_root_datasets_init_existing
      procedure :: gdf_root_datasets_init_new
      generic, public :: init => gdf_root_datasets_init_new, gdf_root_datasets_init_existing
      procedure :: cleanup => gdf_root_datasets_cleanup
      ! finalize :: gdf_root_datasets_cleanup
   end type gdf_root_datasets_T

contains

   subroutine gdf_create_root_datasets(file, rd)

      use hdf5,         only: HID_T
      use helpers_hdf5, only: create_dataset

      implicit none

      integer(HID_T),            intent(in) :: file !< File identifier
      type(gdf_root_datasets_T), intent(in) :: rd   !< cointainer for root datasets

      call create_dataset(file, 'grid_dimensions', rd%grid_dimensions)
      call create_dataset(file, 'grid_left_index', rd%grid_left_index)
      call create_dataset(file, 'grid_level', rd%grid_level)
      call create_dataset(file, 'grid_parent_id', rd%grid_parent_id)
      call create_dataset(file, 'grid_particle_count', rd%grid_particle_count)
   end subroutine gdf_create_root_datasets

   subroutine gdf_create_simulation_parameters(file_id, sp)

      use hdf5,         only: HID_T, h5gcreate_f, h5gclose_f
      use helpers_hdf5, only: create_attribute

      implicit none

      integer(HID_T),         intent(in) :: file_id    !< hdf5 file identification number
      type(gdf_parameters_T), intent(in) :: sp         !< container for simulation parameters

      integer(HID_T)                     :: g_id
      integer(kind=4)                    :: error

      call h5gcreate_f(file_id, 'simulation_parameters', g_id, error)
      call create_attribute(g_id, 'refine_by', sp%refine_by)
      call create_attribute(g_id, 'dimensionality', sp%dimensionality)
      call create_attribute(g_id, 'cosmological_simulation', sp%cosmological_simulation)
      call create_attribute(g_id, 'num_ghost_zones', sp%num_ghost_zones)
      call create_attribute(g_id, 'domain_dimensions', sp%domain_dimensions)
      call create_attribute(g_id, 'domain_left_edge', sp%domain_left_edge)
      call create_attribute(g_id, 'domain_right_edge', sp%domain_right_edge)
      call create_attribute(g_id, 'current_time', sp%current_time)
      call create_attribute(g_id, 'field_ordering', sp%field_ordering)
      call create_attribute(g_id, 'unique_identifier', sp%unique_identifier)
      call create_attribute(g_id, 'boundary_conditions', sp%boundary_conditions)
      call h5gclose_f(g_id, error)

   end subroutine gdf_create_simulation_parameters

   subroutine gdf_create_format_stamp(file)

      use hdf5, only: HID_T, h5gcreate_f, h5gclose_f
      use h5lt, only: h5ltset_attribute_string_f

      implicit none

      integer(HID_T), intent(in) :: file
      integer(HID_T)             :: g_id
      integer(kind=4)            :: error

      character(len=*), parameter :: gname = 'gridded_data_format'
      character(len=*), parameter :: gname2 = '/gridded_data_format'

      call h5gcreate_f(file, gname, g_id, error)
      call h5ltset_attribute_string_f(g_id, gname2, 'data_software', 'piernik', error )
      call h5ltset_attribute_string_f(g_id, gname2, 'data_software_version', '1.0', error )
      call h5gclose_f(g_id, error)

   end subroutine gdf_create_format_stamp

   subroutine gdf_create_field_types(filename, o_func)

      use hdf5, only: HID_T, h5gcreate_f, h5gclose_f, h5fopen_f, h5fclose_f, H5F_ACC_RDWR_F, h5open_f, h5close_f

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

   subroutine gdf_root_datasets_init_existing(this, cg_all_n_b, cg_all_rl, cg_all_off, cg_all_parents, cg_all_particles)
      implicit none
      class(gdf_root_datasets_T),                intent(inout) :: this
      integer(kind=4),  dimension(:,:), pointer, intent(in)    :: cg_all_n_b       !> sizes of all cg
      integer(kind=4),  dimension(:,:), pointer, intent(in)    :: cg_all_rl        !> refinement levels of all cgs
      integer(kind=8),  dimension(:,:), pointer, intent(in)    :: cg_all_off       !> offsets of all cgs
      integer(kind=8),  dimension(:),   pointer, intent(in)    :: cg_all_parents   !> parents IDs of all cgs
      integer(kind=4),  dimension(:,:), pointer, intent(in)    :: cg_all_particles !> particles counts in all cgs

      allocate(this%grid_dimensions(size(cg_all_n_b,1), size(cg_all_n_b, 2)), source=cg_all_n_b)
      allocate(this%grid_left_index(size(cg_all_off,1), size(cg_all_off, 2)), source=cg_all_off)
      allocate(this%grid_level(size(cg_all_rl,1), size(cg_all_rl, 2)), source=cg_all_rl)
      allocate(this%grid_particle_count(size(cg_all_particles,1), size(cg_all_particles, 2)), source=cg_all_particles)
      allocate(this%grid_parent_id(size(cg_all_parents,1)), source=cg_all_parents)
   end subroutine gdf_root_datasets_init_existing

   subroutine gdf_root_datasets_init_new(this, n)
      implicit none
      class(gdf_root_datasets_T), intent(inout) :: this
      integer(kind=4),            intent(in)    :: n

      allocate(this%grid_dimensions(3, n))
      allocate(this%grid_left_index(3, n))
      allocate(this%grid_level(1, n))
      allocate(this%grid_particle_count(1, n))
      allocate(this%grid_parent_id(n))
   end subroutine gdf_root_datasets_init_new

   subroutine gdf_root_datasets_cleanup(this)
      implicit none
      class(gdf_root_datasets_T), intent(inout) :: this

      if (associated(this%grid_dimensions)) deallocate(this%grid_dimensions)
      if (associated(this%grid_left_index)) deallocate(this%grid_left_index)
      if (associated(this%grid_level)) deallocate(this%grid_level)
      if (associated(this%grid_parent_id)) deallocate(this%grid_parent_id)
      if (associated(this%grid_particle_count)) deallocate(this%grid_particle_count)
   end subroutine gdf_root_datasets_cleanup

   subroutine gdf_parameters_init(this)
      implicit none
      class(gdf_parameters_T), intent(inout) :: this

      allocate(this%refine_by(1))
      allocate(this%dimensionality(1))
      allocate(this%cosmological_simulation(1))
      allocate(this%num_ghost_zones(1))
      allocate(this%field_ordering(1))
      allocate(this%domain_dimensions(3))
      allocate(this%domain_left_edge(3))
      allocate(this%domain_right_edge(3))
      allocate(this%current_time(1))
      allocate(this%boundary_conditions(6))
      allocate(this%unique_identifier)
   end subroutine gdf_parameters_init

   subroutine gdf_parameters_finalize(this)
      implicit none
      class(gdf_parameters_T), intent(inout) :: this

      deallocate(this%refine_by)
      deallocate(this%dimensionality)
      deallocate(this%cosmological_simulation)
      deallocate(this%num_ghost_zones)
      deallocate(this%field_ordering)
      deallocate(this%domain_dimensions)
      deallocate(this%domain_left_edge)
      deallocate(this%domain_right_edge)
      deallocate(this%current_time)
      deallocate(this%boundary_conditions)
      deallocate(this%unique_identifier)
   end subroutine gdf_parameters_finalize
end module gdf
