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

!>
!! \brief Initial condition for a Propagation of Circularly polarized Alfvén Waves (CPAW)
!!
!! \details See section 4.1 in
!!   Andrea Mignone, Petros Tzeferacos
!!   A Second-Order Unsplit Godunov Scheme for Cell-Centered MHD: the CTU-GLM scheme
!!   https://arxiv.org/pdf/0911.3410
!<

   use constants, only: dsetnamelen

   implicit none

   private
   public :: read_problem_par, problem_initial_conditions, problem_pointers

   ! namelist parameters
   real :: d0   !! uniform density
   real :: p0   !! uniform pressure
   real :: vx0  !! x-velocity in the domain
   real :: vy0  !! y-velocity in the domain
   real :: vz0  !! z-velocity in the domain
   real :: A    !! perturbation wave amplitude
   real :: vA   !! fraction of the Alfvén speed for perturbed wave
   real :: kx   !! x-wavenumber
   real :: ky   !! y-wavenumber
   real :: kz   !! z-wavenumber
   namelist /PROBLEM_CONTROL/  d0, p0, vx0, vy0, vz0, A, vA, kx, ky, kz

   ! other private data
   character(len=dsetnamelen), parameter :: ini_B_n = "ini_B"

contains

!> \brief Set up custom pointers to tweak the code execution according to our needs

   subroutine problem_pointers

      use user_hooks, only: finalize_problem

      implicit none

      finalize_problem => verify_test

   end subroutine problem_pointers

!> \brief Read the runtime parameters specified in the namelist

   subroutine read_problem_par

      use bcast,          only: piernik_MPI_Bcast
      use cg_list_global, only: all_cg
      use constants,      only: AT_NO_B, ndims
      use dataio_pub,     only: nh
      use mpisetup,       only: rbuff, master, slave

      implicit none

      ! the default values
      d0  = 1.0
      p0  = 0.1
      vx0 = 1.0
      vy0 = 0.0
      vz0 = 0.0
      A   = 0.1
      vA  = 1.0
      kx  = 1.0
      ky  = 0.0
      kz  = 0.0

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

         rbuff(1)  = d0
         rbuff(2)  = p0
         rbuff(3)  = vx0
         rbuff(4)  = vy0
         rbuff(5)  = vz0
         rbuff(6)  = A
         rbuff(7)  = vA
         rbuff(8)  = kx
         rbuff(9)  = ky
         rbuff(10) = kz
      endif

      call piernik_MPI_Bcast(rbuff)

      if (slave) then

         d0  = rbuff(1)
         p0  = rbuff(2)
         vx0 = rbuff(3)
         vy0 = rbuff(4)
         vz0 = rbuff(5)
         A   = rbuff(6)
         vA  = rbuff(7)
         kx  = rbuff(8)
         ky  = rbuff(9)
         kz  = rbuff(10)

      endif

      call all_cg%reg_var(ini_B_n, restart_mode = AT_NO_B, dim4=ndims)

   end subroutine read_problem_par

!> \brief Set up the initial conditions. Note that this routine can be called multiple times during initial iterations of refinement structure

   subroutine problem_initial_conditions

      use cg_leaves,        only: leaves
      use cg_list,          only: cg_list_element
      use constants,        only: xdim, ydim, zdim, LO, HI, ndims
      use dataio_pub,       only: die
      use fluidindex,       only: flind
      use fluidtypes,       only: component_fluid
      use func,             only: ekin, emag
      use grid_cont,        only: grid_container
      use named_array_list, only: wna

      implicit none

      class(component_fluid), pointer :: fl
      integer                         :: i, j, k
      real                            :: xi, yj, zk, phi
      type(cg_list_element),  pointer :: cgl
      type(grid_container),   pointer :: cg
      integer                         :: p

      real, dimension(ndims, ndims)     :: rot_matrix, rot_matrix_inv
      real, dimension(ndims)            :: v_prime, b_prime, v_final, b_final
      real, dimension(ndims)            :: pos_vec, pos_vec_prime
      real                              :: aa, g, k_prime
      real, dimension(:,:,:,:), pointer :: ini_B

      aa = atan(ky/kx)
      g = atan(cos(aa) * (kz/kx))
      k_prime = sqrt(kx**2 + ky**2 + kz**2)

      rot_matrix = reshape([cos(aa)*cos(g),  sin(aa)*cos(g), sin(g), &
                           -sin(aa),         cos(aa),        0.0,    &
                           -cos(aa)*sin(g), -sin(aa)*sin(g), cos(g)], [ndims, ndims])

      rot_matrix_inv = transpose(rot_matrix)

      do p = 1, flind%fluids

         fl => flind%all_fluids(p)%fl

         cgl => leaves%first
         do while (associated(cgl))
            cg => cgl%cg

            ini_B => cg%w(wna%ind(ini_B_n))%arr
            if (.not. associated(ini_B)) call die("[initproblem:problem_initial_conditions] ini_B not found")

            do j = cg%lhn(ydim,LO), cg%lhn(ydim,HI)
               yj = cg%y(j)
               do i = cg%lhn(xdim,LO), cg%lhn(xdim,HI)
                  xi = cg%x(i)
                  do k = cg%lhn(zdim,LO), cg%lhn(zdim,HI)
                     zk = cg%z(k)

                     pos_vec = [xi, yj, zk]
                     pos_vec_prime = matmul(rot_matrix_inv, pos_vec)

                     phi = k_prime * pos_vec_prime(1)

                     v_prime = [vx0, vy0 + A * sin(phi), vz0 + A * cos(phi)]
                     b_prime = [vA * sqrt(d0), -sqrt(d0) * A * sin(phi), -sqrt(d0) * A * cos(phi)]

                     v_final = matmul(rot_matrix, v_prime)
                     b_final = matmul(rot_matrix, b_prime)

                     cg%u(fl%idn, i, j, k) = d0
                     cg%u(fl%imx, i, j, k) = v_final(xdim)
                     cg%u(fl%imy, i, j, k) = v_final(ydim)
                     cg%u(fl%imz, i, j, k) = v_final(zdim)

                     if (fl%has_energy) then
                        cg%u(fl%ien, i, j, k) = p0/(fl%gam_1) + ekin(v_final(xdim), v_final(ydim), v_final(zdim), d0)

                        if (fl%is_magnetized) then
                           cg%b(xdim:zdim, i, j, k) = b_final(:)
                           cg%u(fl%ien, i, j, k) = cg%u(fl%ien, i, j, k) + emag(b_final(xdim), b_final(ydim), b_final(zdim))
                        endif
                     endif

                     ini_B(:, i, j, k) = b_final(:)
                  enddo
               enddo
            enddo

            cgl => cgl%nxt
         enddo
      enddo

   end subroutine problem_initial_conditions

!>
!! \brief Final routine called after the endo of the simulation
!!
!! \ToDo For long runs maybe call it every period?
!<

   subroutine verify_test

      use allreduce,        only: piernik_MPI_Allreduce
      use cg_list,          only: cg_list_element
      use cg_leaves,        only: leaves
      use constants,        only: V_ESSENTIAL, V_INFO, xdim, zdim, pSUM
      use dataio_pub,       only: msg, printinfo, die
      use domain,           only: dom
      use grid_cont,        only: grid_container
      use mpisetup,         only: master
      use named_array_list, only: wna

      implicit none

      enum, bind(C)
         enumerator :: N_D, N_2
      end enum
      real, dimension(N_D:N_2) :: norm
      type(cg_list_element),  pointer   :: cgl
      type(grid_container),   pointer   :: cg
      real, dimension(:,:,:,:), pointer :: ini_B
      integer :: d

      if (master) then
         call printinfo('Copy compare_slices.py from problem folder to the run folder.', V_INFO)
         call printinfo('Make sure that the cpaw_tst_0000.h5 file is the intial file and cpaw_tst_0001.h5 is the last file.', V_INFO)
         call printinfo('Then run python compare_slices.py to see the comparison.', V_INFO)
      endif

      norm = 0.
      cgl => leaves%first
      do while (associated(cgl))
         cg => cgl%cg

         ini_B => cg%w(wna%ind(ini_B_n))%arr
         if (.not. associated(ini_B)) call die("[initproblem:verify_test] ini_B not found")

         do d = xdim, zdim
            if (dom%has_dir(d)) then
               cg%wa(RNG) = ini_B(d, RNG) - cg%b(d, RNG)
               norm(N_D) = norm(N_D) + sum(cg%wa(RNG)**2, mask=cg%leafmap)
               norm(N_2) = norm(N_2) + sum(ini_B(d, RNG)**2, mask=cg%leafmap)
            endif
         enddo

         cgl => cgl%nxt
      enddo

      do d = N_D, N_2
         call piernik_MPI_Allreduce(norm(d), pSUM)
      enddo

      if (master) then
         write(msg,'(a,f12.8)')"[initproblem:verify_test] L2 error norm od the B field = ", sqrt(norm(N_D)/norm(N_2))
         call printinfo(msg, V_ESSENTIAL)
      endif

   end subroutine verify_test

end module initproblem
