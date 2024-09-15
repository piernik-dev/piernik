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

! Initial condition for blob test
! Blob test by Agertz et al., 2007, MNRAS, 380, 963.
!
! For description and Initial Condition files look here:
! http://www.astrosim.net/code/doku.php?id=home:codetest:hydrotest:wengen:wengen3

   use constants, only: cwdlen

   implicit none

   private
   public :: read_problem_par, problem_initial_conditions, problem_pointers

   integer, parameter :: expected_type_size = 4
   real(kind=expected_type_size), dimension(:,:,:,:), allocatable :: data

   ! namelist parameters
   real   :: chi, rblob, blobxc, blobyc, blobzc, Mext, denv, venv
   character(len=cwdlen) :: ICfile

   namelist /PROBLEM_CONTROL/  chi, rblob, blobxc, blobyc, blobzc, Mext, denv, venv, ICfile

contains

!-----------------------------------------------------------------------------

   subroutine problem_pointers

      use dataio_user, only: user_tsl
      use user_hooks,  only: problem_post_IC

      implicit none

      problem_post_IC => deallocate_h5IC
      user_tsl        => clump_mass

   end subroutine problem_pointers

!-----------------------------------------------------------------------------

   subroutine read_problem_par

      use bcast,      only: piernik_MPI_Bcast
      use constants,  only: cwdlen
      use dataio_pub, only: nh
      use mpisetup,   only: rbuff, master, slave

      implicit none

      character(len=cwdlen) :: lcbuff

      chi     =   10.0
      rblob   =  197.0
      blobxc  =    0.0
      blobyc  =    0.0
      blobzc  =    0.0
      Mext    =    2.7
      denv    =    3.13e-8
      venv    = 1000.0

      ICfile = ""

      if (master) then

         if (.not.nh%initialized) call nh%init()
         open(newunit=nh%lun, file=nh%tmp1, status="unknown")
         write(nh%lun,nml=PROBLEM_CONTROL)
         close(nh%lun)
         open(newunit=nh%lun, file=nh%par_file)
         nh%errstr=""
         read(unit=nh%lun, nml=PROBLEM_CONTROL, iostat=nh%ierrh, iomsg=nh%errstr)
         close(nh%lun)
         call nh%namelist_errh(nh%ierrh, "PROBLEM_CONTROL")
         read(nh%cmdl_nml,nml=PROBLEM_CONTROL, iostat=nh%ierrh)
         call nh%namelist_errh(nh%ierrh, "PROBLEM_CONTROL", .true.)
         open(newunit=nh%lun, file=nh%tmp2, status="unknown")
         write(nh%lun,nml=PROBLEM_CONTROL)
         close(nh%lun)
         call nh%compare_namelist()

         rbuff(1) = chi
         rbuff(2) = rblob
         rbuff(3) = blobxc
         rbuff(4) = blobyc
         rbuff(5) = blobzc
         rbuff(6) = Mext
         rbuff(7) = denv
         rbuff(8) = venv

         lcbuff   = ICfile

      endif

      call piernik_MPI_Bcast(rbuff)
      call piernik_MPI_Bcast(lcbuff, cwdlen)

      if (slave) then

         chi      = rbuff(1)
         rblob    = rbuff(2)
         blobxc   = rbuff(3)
         blobyc   = rbuff(4)
         blobzc   = rbuff(5)
         Mext     = rbuff(6)
         denv     = rbuff(7)
         venv     = rbuff(8)

         ICfile   = lcbuff

      endif

   end subroutine read_problem_par

!-----------------------------------------------------------------------------

   subroutine problem_initial_conditions

      use dataio_pub, only: msg, printinfo
      use mpisetup,   only: master

      implicit none

      if (len(trim(ICfile)) > 0) then
         write(msg, *) "[initproblem:problem_initial_conditions] Reading Initial Conditions from '", trim(ICfile), "'"
         if (master) call printinfo(msg)
         call problem_initial_conditions_original
      else
         if (master) call printinfo("[initproblem:problem_initial_conditions] Using analytical Initial Conditions")
         call problem_initial_conditions_analytical
      endif

   end subroutine problem_initial_conditions

!-----------------------------------------------------------------------------

   subroutine problem_initial_conditions_original

      use cg_leaves,  only: leaves
      use cg_list,    only: cg_list_element
      use dataio_pub, only: msg, die
      use domain,     only: dom
      use fluidindex, only: flind
      use fluidtypes, only: component_fluid
      use func,       only: ekin
      use grid_cont,  only: grid_container
      use units,      only: km, sek
      implicit none

      logical :: firstcall = .true.
      integer :: f
      class(component_fluid), pointer :: fl
      type(cg_list_element),  pointer :: cgl
      type(grid_container),   pointer :: cg

      if (firstcall) then
#ifdef HDF5
         call problem_initial_conditions_readh5
#else /* !HDF5 */
         call die("[initproblem:problem_initial_conditions_original] Without HDF5 available try to use analytical initial conditions (leave empty ICfile variable)")
#endif /* !HDF5 */
         firstcall = .false.
      endif

      ! Most naive: put the data 1:1, expect correct sizes, don't try to be friendly (ToDo: add more flexibility for AMR)

      if (any(dom%n_d /= shape(data(1, :, :, :)))) then
         write(msg, *)"[initproblem:problem_initial_conditions_original] domain doesn't match IC data: ",dom%n_d," /= ", shape(data(1, :, :, :)), " (read it in Z-Y-X order)"
         call die(msg)
      endif

      fl => flind%neu
      cgl => leaves%first
      do while (associated(cgl))
         cg => cgl%cg
#define RNG1 1+cg%is:1+cg%ie, 1+cg%js:1+cg%je, 1+cg%ks:1+cg%ke
         cg%u(fl%idn, RNG) = data(1, RNG1)
         do f = fl%imx, fl%imz
            cg%u(f, RNG) = sek/km * data(2+f-fl%imx, RNG1) * cg%u(fl%idn, RNG)
         enddo
         cg%u(fl%ien, RNG) = data(5, RNG1) * cg%u(fl%idn, RNG)+ekin(cg%u(fl%imx, RNG), cg%u(fl%imy, RNG), cg%u(fl%imz, RNG), cg%u(fl%idn, RNG))
#undef RNG1
         cgl => cgl%nxt
      enddo

   end subroutine problem_initial_conditions_original

#ifdef HDF5
   subroutine problem_initial_conditions_readh5

      use constants,  only: cbuff_len
      use dataio_pub, only: msg, die
      use hdf5,       only: H5F_ACC_RDONLY_F, HID_T, HSIZE_T, SIZE_T, H5T_NATIVE_REAL, &
           &                h5open_f, h5close_f, h5fopen_f, h5fclose_f
      use h5lt,       only: h5ltfind_dataset_f, h5ltget_dataset_ndims_f, h5ltget_dataset_info_f, &
           &                h5ltread_dataset_f
      use mpisetup,   only: master, slave

      implicit none

      logical :: file_exist
      character(len=cbuff_len), dimension(5), parameter :: dsets = ["Density    ", "x-velocity ", "y-velocity ", "z-velocity ", "TotalEnergy" ]
      integer(kind=4) :: error, rank, type_class
      integer(HID_T) :: file_id
      integer :: d
      integer(SIZE_T) :: type_size
      integer(HSIZE_T), allocatable, dimension(:) :: dims

      inquire(file = trim(ICfile), exist = file_exist)
      if (.not. file_exist) then
         write(msg,'(3a)') '[initproblem:problem_initial_conditions_readh5] IC file: ', trim(ICfile),' does not exist'
         call die(msg)
      endif

      call h5open_f(error)
      if (master .or. slave) then  ! ToDo do somewhat more intelligent reading as this may easily clutter memory on lagre IC file
         call h5fopen_f(trim(ICfile), H5F_ACC_RDONLY_F, file_id, error)
         if (error /= 0) call die("[initproblem:problem_initial_conditions_readh5] error opening IC file")
         do d = lbound(dsets, 1), ubound(dsets, 1)
            if (h5ltfind_dataset_f(file_id, trim(dsets(d))) == 0) then
               write(msg, *)"[initproblem:problem_initial_conditions_readh5] Cannot find dataset '",trim(dsets(d)),"'"
               call die(msg)
            endif
            call h5ltget_dataset_ndims_f(file_id, dsets(d), rank, error)
            if (rank /= 3) call die("[initproblem:problem_initial_conditions_readh5] Wrong dataset rank")
            allocate(dims(rank+1))
            dims(1) = size(dsets)
            call h5ltget_dataset_info_f(file_id, dsets(d), dims(2:), type_class, type_size, error)
            if (type_size /= expected_type_size) call die("[initproblem:problem_initial_conditions_readh5] Wrong type size")
            if (d==1) then
               allocate(data(dims(1), dims(2), dims(3), dims(4)))
            else
               if (any(dims /= shape(data))) call die("[initproblem:problem_initial_conditions_readh5] datasets differ in shape")
            endif
            call h5ltread_dataset_f(file_id, dsets(d), H5T_NATIVE_REAL, data(d, :, :, :), dims, error)
            deallocate(dims)
         enddo
         call h5fclose_f(file_id, error)
      endif
      call h5close_f(error)

   end subroutine problem_initial_conditions_readh5
#endif /* HDF5 */

   subroutine deallocate_h5IC

      use dataio_pub, only: printinfo, msg
      use mpisetup,   only: master
      use units,      only: newtong

      implicit none

      if (master) then
         write(msg, *) "Newton's constant, G = ", newtong
         call printinfo(msg)
      endif

      if (allocated(data)) deallocate(data)

   end subroutine deallocate_h5IC

!-----------------------------------------------------------------------------

   subroutine problem_initial_conditions_analytical

      use cg_leaves,  only: leaves
      use cg_list,    only: cg_list_element
      use constants,  only: xdim, ydim, zdim, LO, HI
      use domain,     only: dom
      use fluidindex, only: flind
      use fluidtypes, only: component_fluid
      use func,       only: ekin
      use grid_cont,  only: grid_container

      implicit none

      class(component_fluid), pointer :: fl
      real                            :: uenv, rcx, rcy, rrel, rblob2
      integer                         :: i, j, k
      type(cg_list_element),  pointer :: cgl
      type(grid_container),   pointer :: cg

      fl => flind%neu
      rcx = 0.
      rcy = 0.
      rblob2 = rblob**2

      uenv = (venv / Mext)**2 * denv / (fl%gam * fl%gam_1)

      cgl => leaves%first
      do while (associated(cgl))
         cg => cgl%cg

         cg%u(fl%imx, :, :, :) = 0.0
         cg%u(fl%imy, :, :, :) = 0.0

         do i = cg%lhn(xdim,LO), cg%lhn(xdim,HI)
            if (dom%has_dir(xdim)) rcx = (cg%x(i)-blobxc)**2
            do j = cg%lhn(ydim,LO), cg%lhn(ydim,HI)
               if (dom%has_dir(ydim)) rcy = (cg%y(j)-blobyc)**2
               do k = cg%lhn(zdim,LO), cg%lhn(zdim,HI)
                  rrel = rcx + rcy
                  if (dom%has_dir(zdim)) rrel = rrel + (cg%z(k)-blobzc)**2

                  if (rblob2 >= rrel) then         ! inside clump
                     cg%u(fl%idn,i,j,k) = chi*denv
                     cg%u(fl%imz,i,j,k) = 0.0
                  else                             ! outsize clump
                     cg%u(fl%idn,i,j,k) = denv
                     cg%u(fl%imz,i,j,k) = cg%u(fl%idn,i,j,k)*venv
                  endif
               enddo
            enddo
         enddo

         cg%u(fl%ien, :, :, :) = uenv + ekin(cg%u(fl%imx, :, :, :), cg%u(fl%imy, :, :, :), cg%u(fl%imz, :, :, :), cg%u(fl%idn, :, :, :))

         cgl => cgl%nxt
      enddo

   end subroutine problem_initial_conditions_analytical

!------------------------------------------------------------------------------------------

   subroutine clump_mass(user_vars, tsl_names)

      use allreduce,   only: piernik_MPI_Allreduce
      use cg_leaves,   only: leaves
      use cg_list,     only: cg_list_element
      use constants,   only: pSUM
      use diagnostics, only: pop_vector
      use fluidindex,  only: flind
      use grid_cont,   only: grid_container
      use mpisetup,    only: master

      implicit none

      real, dimension(:), intent(inout), allocatable                       :: user_vars
      character(len=*), dimension(:), intent(inout), allocatable, optional :: tsl_names

      real :: m_clump
      type(cg_list_element), pointer :: cgl
      type(grid_container),  pointer :: cg
      real, parameter :: rho_thr = 0.64 * 3.13e-7, T_thr = 0.9 *123900. ! ToDo: get more exact eatimate of T_ext

      if (present(tsl_names)) then
         call pop_vector(tsl_names, len(tsl_names(1)), ["m_clump"])      !   add to header
      else

         m_clump = 0.
         cgl => leaves%first
         do while (associated(cgl))
            cg => cgl%cg
            m_clump = m_clump + cg%dvol * sum(cg%u(flind%neu%idn, RNG), &
                 &                    mask = (cg%leafmap(         RNG) .and. &
                 &                            cg%u(flind%neu%idn, RNG) > rho_thr .and. &
                 &                            cg%u(flind%neu%ien, RNG) < T_thr))
            cgl => cgl%nxt
         enddo

         call piernik_MPI_Allreduce(m_clump, pSUM)
         if (master) call pop_vector(user_vars,[m_clump])                 !   pop value
      endif


   end subroutine clump_mass

end module initproblem
