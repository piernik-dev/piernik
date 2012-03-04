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
!! \brief Module responsible for integration of OpenCL in Piernik
!!
!!
!! In this module a namelist of parameters is specified:
!! \copydetails piernikcl::init_opencl
!<
module piernikcl
! pulled PIERNIK_OPENCL
   use cl,           only: cl_context, cl_command_queue

   implicit none
   private
   public  :: init_opencl, cleanup_opencl, context, command_queue

   integer :: ierr
   type(cl_context)       :: context
   type(cl_command_queue) :: command_queue

contains

   subroutine init_opencl

      use cl,        only: clCreateContext, clCreateCommandQueue, CL_QUEUE_PROFILING_ENABLE, cl_platform_id, cl_device_id

      implicit none

      type(cl_platform_id)   :: platform
      type(cl_device_id)     :: device

      call piernikCL_get_suitable_gpu(platform,device)

      ! create the context and the command queue
      context = clCreateContext(platform, device, ierr)
      command_queue = clCreateCommandQueue(context, device, CL_QUEUE_PROFILING_ENABLE, ierr)

   end subroutine init_opencl

   subroutine cleanup_opencl

      use cl,         only: clReleaseCommandQueue, clReleaseContext

      implicit none

      integer(kind=4) :: ierr

      call clReleaseCommandQueue(command_queue, ierr)
      call clReleaseContext(context, ierr)

   end subroutine cleanup_opencl

   subroutine clerr(ierr,msg)
      use cl,         only: CL_SUCCESS
      use dataio_pub, only: die

      implicit none

      integer(kind=4), intent(in) :: ierr
      character(len=*), intent(in) :: msg

      if (ierr /= CL_SUCCESS) call die(msg)

   end subroutine clerr

   subroutine piernikCL_get_suitable_gpu(platform,device)

      use cl,         only: cl_platform_id, cl_device_id, clGetPlatformIDs, clGetPlatformIDs, clGetDeviceIDs, &
                       &   CL_DEVICE_TYPE_GPU, clGetPlatformInfo, CL_PLATFORM_NAME, CL_PLATFORM_VERSION, &
                       &   CL_DEVICE_NAME, clGetDeviceInfo
      use dataio_pub, only: printinfo, msg, msglen

      implicit none

      type(cl_platform_id), intent(out) :: platform
      type(cl_device_id), intent(out) :: device

      type(cl_platform_id), dimension(:), allocatable :: platforms
      type(cl_device_id), dimension(:), allocatable :: devices
      integer(kind=4) :: num_platforms, num_devices, ierr, iplat, idev
      integer(kind=8) :: val
      character(len=msglen) :: clinfo

      ! get the number of platforms
      call clGetPlatformIDs(num_platforms, ierr)
      call clerr(ierr,"Cannot get CL platforms number.")
      allocate(platforms(num_platforms))

      ! get an array of platforms
      call clGetPlatformIDs(platforms, num_platforms, ierr)
      call clerr(ierr,"Cannot get CL platforms.")
      do iplat = 1, num_platforms
         ! get the device ID
         call clGetDeviceIDs(platforms(iplat), CL_DEVICE_TYPE_GPU, num_devices, ierr)
         if (num_devices == 0) then
            cycle ! skip platform without gpus
         else
            exit ! on first gpu that is found
         endif
         call clerr(ierr,"I've hit wall while getting devices on CL platform")
      enddo

      allocate(devices(num_devices))
      call clGetDeviceIDs(platforms(iplat), CL_DEVICE_TYPE_GPU, devices, num_devices, ierr)

      platform = platforms(iplat)
      device   = devices(1)   ! settle for first available gpu for now
      deallocate(devices, platforms)
      ! get the device name and print it
      call clGetPlatformInfo(platform, CL_PLATFORM_NAME, clinfo, ierr)
      write(msg, '(2a)') 'Name      : ', trim(clinfo)
      call printinfo(msg)
      call clGetPlatformInfo(platform, CL_PLATFORM_VERSION, clinfo, ierr)
      write(msg, '(2a)') 'Version   : ', trim(clinfo)
      call printinfo(msg)
      call clGetDeviceInfo(device, CL_DEVICE_NAME, clinfo, ierr)
      write(msg, '(2a)') "CL device : ", trim(clinfo)
      call printinfo(msg)

      return
  end subroutine piernikCL_get_suitable_gpu

end module piernikcl
