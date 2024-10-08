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

   use constants, only: dsetnamelen, ndims

   implicit none

   private
   public :: read_problem_par, problem_initial_conditions, problem_pointers

   ! namelist parameters
   integer(kind=4) :: ord_prolong   !< Prolongation order used for working field
   logical :: point                 !< If .true. then IC contains just one non-zero value, some more complex pattern is used otherwise
   real, dimension(ndims) :: coords !< Coordinates of the point mentioned above
   integer(kind=4) :: checkerboard  !< If > 0 and .not. point then paint a global checkerboard pattern of this length
   namelist /PROBLEM_CONTROL/ ord_prolong, point, coords, checkerboard

   ! other private data
   character(len=dsetnamelen), parameter :: fld_n = "fld"

contains

!-----------------------------------------------------------------------------

   subroutine problem_pointers

      use dataio_user, only: user_reg_var_restart

      implicit none

      user_reg_var_restart => fld_reg

   end subroutine problem_pointers

!-----------------------------------------------------------------------------

   subroutine read_problem_par

      use bcast,      only: piernik_MPI_Bcast
      use constants,  only: PIERNIK_INIT_MPI, xdim, zdim
      use dataio_pub, only: code_progress, die, nh
      use mpisetup,   only: master, slave, ibuff, lbuff, rbuff

      implicit none

      if (code_progress < PIERNIK_INIT_MPI) call die("[initproblem:read_problem_par] MPI not initialized.")

      ord_prolong  = 0
      checkerboard = 0
      point        = .false.
      coords(:)    = 0.

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

         ibuff(1) = ord_prolong
         ibuff(2) = checkerboard

         lbuff(1) = point

         rbuff(xdim:zdim) = coords(:)
      endif

      call piernik_MPI_Bcast(ibuff)
      call piernik_MPI_Bcast(lbuff)
      call piernik_MPI_Bcast(rbuff)

      if (slave) then

         ord_prolong  = ibuff(1)
         checkerboard = ibuff(2)

         point        = lbuff(1)

         coords(:)    = rbuff(xdim:zdim)

      endif

   end subroutine read_problem_par

!-----------------------------------------------------------------------------

   subroutine problem_initial_conditions

      use constants,  only: base_level_id
      use dataio_pub, only: die
      use mpisetup,   only: master

      implicit none

      integer(kind=4) :: i
      integer, parameter :: lmax = 3

      call fld_reg

      do i = lmax, 1, -1
         call prepare_fld(base_level_id-i, 0)
      enddo
      do i = lmax, 0, -1
         call prepare_fld(base_level_id, int(i))
      enddo

      if (master) call die("[initproblem:problem_initial_conditions] End of test")

   end subroutine problem_initial_conditions

!-----------------------------------------------------------------------------

   subroutine prepare_fld(lev, n_c)

      use cg_level_base,      only: base
      use cg_level_connected, only: cg_level_connected_t
      use dataio_pub,    only: printinfo, warn, msg
      use mpisetup,      only: master, proc
      use named_array_list, only: qna

      implicit none

      integer(kind=4), intent(in) :: lev !< which level to initialize
      integer,         intent(in) :: n_c !< how many times to coarsen the data

      type(cg_level_connected_t), pointer :: curl
      integer :: c

      if (master) then
         write(msg,'(2(a,i4))')"[initproblem:prepare_fld] ^",lev," coarsen ", n_c
         call printinfo(msg)
      endif

      call clear_fld(huge(1_4))
      call set_up_fld(lev)
      if (point) then
         write(msg,*)"ip:pf set ^",lev
         call printinfo(msg)
         call find_non_0_or_write_hdf5
      endif

      curl => base%level
      do while (associated(curl) .and. curl%l%id /= lev)
         call clear_lev(curl)
         curl => curl%coarser
      enddo
      if (.not. associated(curl)) then
         write(msg,'(a,i4,a)')"[initproblem:prepare_fld] ^",lev," not found"
         if (master) call warn(msg)
         return
      endif

      do c = 1, n_c
         if (.not. associated(curl%coarser)) then
            write(msg,'(2(a,i4))')"[initproblem:prepare_fld] no coarser level ^", lev, ", n_c", n_c
            if (master) call warn(msg)
            return
         endif
         call curl%restrict_1var(qna%ind(fld_n))
         call clear_lev(curl)
         if (point) then
            write(msg,*)"ip:pf restricted ^",curl%l%id
            call printinfo(msg)
            call find_non_0_or_write_hdf5
         endif
         curl => curl%coarser
      enddo

      do while (associated(curl))
         if (associated(curl%finer)) then
            call curl%prolong_1var(qna%ind(fld_n))
            call clear_lev(curl)
            if (point) then
               write(msg,*)"ip:pf prolonged ^",curl%l%id, " @", proc
               call printinfo(msg)
               call find_non_0_or_write_hdf5
            endif
         endif
         curl => curl%finer
      enddo

      if (master) call printinfo("ip:pf finished")
      call find_non_0_or_write_hdf5

   end subroutine prepare_fld

!-----------------------------------------------------------------------------

   subroutine set_up_fld(lev)

      use cg_list,            only: cg_list_element
      use cg_level_base,      only: base
      use cg_level_connected, only: cg_level_connected_t
      use constants,          only: xdim, ydim, zdim, ndims, LO, HI
      use dataio_pub,         only: msg, warn
      use domain,             only: dom
      use fluidindex,         only: iarr_all_dn
      use grid_cont,          only: grid_container
      use mpisetup,           only: master
      use named_array_list,   only: qna, wna

      implicit none

      integer(kind=4), intent(in) :: lev

      type(cg_level_connected_t), pointer :: glev
      type(cg_list_element), pointer :: cgl
      type(grid_container), pointer :: cg
      integer, dimension(ndims) :: ijk

      call clear_fld(lev)

      if (point .or. checkerboard > 0) return

      glev => base%level
      do while (associated(glev))
         if (glev%l%id == lev) exit
         glev => glev%coarser
      enddo

      if (.not. associated(glev)) then
         write(msg, '(a,i3)')"[initproblem:set_up_fld] Can't find level# ", lev
         if (master) call warn(msg)
         return
      endif

      cgl => glev%first
      do while (associated(cgl))
         cg => cgl%cg

         cg%q(qna%ind(fld_n))%arr(cg%is, cg%js, cg%ks) = cg%grid_id * 1.1

         if (all(cg%n_b > 3 .or. .not. dom%has_dir(:))) then
            ijk = min( cg%ijkse(:, LO) + 3_4*dom%D_(:), cg%ijkse(:, HI) )
            cg%q(qna%ind(fld_n))%arr(ijk(xdim), ijk(ydim), ijk(zdim)) = cg%grid_id * 1.2
         endif

         if (all(cg%n_b > 7 .or. .not. dom%has_dir(:))) then
            ijk = min( cg%ijkse(:, LO) + 7_4*dom%D_(:), cg%ijkse(:, HI) )
            cg%q(qna%ind(fld_n))%arr(ijk(xdim), ijk(ydim), cg%ks:cg%ke) = cg%grid_id * 1.4
         endif

         if (all(cg%n_b > 11 .or. .not. dom%has_dir(:))) then
            ijk = min( cg%ijkse(:, LO) + 11_4*dom%D_(:), cg%ijkse(:, HI) )
            cg%q(qna%ind(fld_n))%arr(ijk(xdim), cg%js:cg%je, cg%ks:cg%ke) = cg%grid_id * 1.8
         endif

         if (all(cg%n_b > 3 .or. .not. dom%has_dir(:))) then
            cg%q(qna%ind(fld_n))%arr(cg%ie-3*dom%D_x:cg%ie, cg%je-3*dom%D_y:cg%je, cg%ke-3*dom%D_z:cg%ke) = cg%grid_id * 2.6
         endif

         cgl => cgl%nxt
      enddo

      if (associated(glev, base%level)) call glev%qw_copy(qna%ind(fld_n), wna%fi, iarr_all_dn(1))

   end subroutine set_up_fld

!-----------------------------------------------------------------------------

   subroutine clear_fld(lev)

      use cg_list,          only: cg_list_element
      use cg_list_global,   only: all_cg
      use constants,        only: xdim, ydim, zdim, LEFT, RIGHT
      use grid_cont,        only: grid_container
      use named_array_list, only: qna

      implicit none

      integer(kind=4), intent(in) :: lev

      type(cg_list_element), pointer :: cgl
      type(grid_container), pointer :: cg
      integer :: i, j, k, i_fld

      i_fld = qna%ind(fld_n)

      cgl => all_cg%first
      do while (associated(cgl))
         cg => cgl%cg
         cg%q(i_fld)%arr(:,:,:) = 0.
         if (cg%l%id == lev) then
            do i = cg%is, cg%ie
               do j = cg%js, cg%je
                  do k = cg%ks, cg%ke
                     if (point) then
                        if ( cg%coord(LEFT, xdim)%r(i) < coords(xdim) .and. cg%coord(RIGHT, xdim)%r(i) >= coords(xdim) .and. &
                             cg%coord(LEFT, ydim)%r(j) < coords(ydim) .and. cg%coord(RIGHT, ydim)%r(j) >= coords(ydim) .and. &
                             cg%coord(LEFT, zdim)%r(k) < coords(zdim) .and. cg%coord(RIGHT, zdim)%r(k) >= coords(zdim)) then
                           cg%q(i_fld)%arr(i, j, k) = 1.
                        else
                           cg%q(i_fld)%arr(i, j, k) = 0.
                        endif
                     else if (checkerboard > 0) then
                        cg%q(i_fld)%arr(i, j, k) = sum(mod([ i, j, k ], int(checkerboard)))
                     else
                        cg%q(i_fld)%arr(i, j, k) = cg%grid_id * (1. + xab(i, cg%is, cg%ie) + xab(j, cg%js, cg%je) + xab(k, cg%ks, cg%ke))
                     endif
                  enddo
               enddo
            enddo
         endif
         cgl => cgl%nxt
      enddo

   end subroutine clear_fld

!-----------------------------------------------------------------------------

   subroutine find_non_0_or_write_hdf5

      use cg_list,          only: cg_list_element
      use cg_list_global,   only: all_cg
      use constants,        only: stdout
      use dataio_pub,       only: msg, printinfo
      use func,             only: operator(.notequals.)
      use grid_cont,        only: grid_container
      use mpisetup,         only: proc, err_mpi
      use MPIF,             only: MPI_COMM_WORLD
      use named_array_list, only: qna
#ifdef HDF5
      use data_hdf5,        only: write_hdf5
#else /* !HDF5 */
      use dataio_pub,       only: warn
#endif /* !HDF5 */

      implicit none

      type(cg_list_element), pointer :: cgl
      type(grid_container), pointer :: cg
      integer :: i, j, k, i_fld
      integer, save :: n = 0

      i_fld = qna%ind(fld_n)

      if (point) then
         cgl => all_cg%first
         do while (associated(cgl))
            cg => cgl%cg
            do i = cg%is, cg%ie
               do j = cg%js, cg%je
                  do k = cg%ks, cg%ke
                     if (cg%q(i_fld)%arr(i, j, k) .notequals. 0.) then
                        write(msg,'(a,i5.5,2(a,i5),a,3f10.5,a,i3,a,g15.7)')"fn",n," @",proc," #",cg%grid_id," [",cg%x(i),cg%y(j),cg%z(k),"]^",cg%l%id," =",cg%q(i_fld)%arr(i, j, k)
                        call printinfo(msg)
                     endif
                  enddo
               enddo
            enddo
            cgl => cgl%nxt
         enddo
         n = n + 1
      else
#ifdef HDF5
         call write_hdf5(.false.)
#else /* !HDF5 */
         call warn("[initproblem:find_non_0_or_write_hdf5] no HDF5 available")
#endif /* !HDF5 */
      endif

      call flush(stdout)
      call MPI_Barrier(MPI_COMM_WORLD, err_mpi)

   end subroutine find_non_0_or_write_hdf5

!-----------------------------------------------------------------------------

   real function xab(x, a, b)

      implicit none

      integer,         intent(in) :: x
      integer(kind=4), intent(in) :: a, b

      if (a /= b) then
         xab = real(2*x-a-b)/real(a-b)
         xab = xab**2
      else
         xab = 1.
      endif

   end function xab

!-----------------------------------------------------------------------------

   subroutine clear_lev(lev)

      use cg_level_connected, only: cg_level_connected_t
      use dataio_pub  ,  only: warn
      use named_array_list, only: qna

      implicit none

      type(cg_level_connected_t), pointer, intent(in) :: lev

      if (.not. associated(lev)) then
         call warn("cl: null level")
         return
      endif

      call lev%set_q_value(qna%ind(fld_n), 0.)

   end subroutine clear_lev

!-----------------------------------------------------------------------------

   subroutine fld_reg

      use cg_list_global, only: all_cg
      use constants,      only: AT_OUT_B

      implicit none

      call all_cg%reg_var(fld_n, restart_mode = AT_OUT_B, multigrid=.true., ord_prolong = ord_prolong)

   end subroutine fld_reg

!-----------------------------------------------------------------------------

end module initproblem
