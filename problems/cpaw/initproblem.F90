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
! ----------------------------------------- !
! Initial condition for a CPAW              !
! See: Andrea Mignone, Petros Tzeferacosa   !
! A Second-Order Unsplit Godunov Scheme for !
! Cell-Centered MHD: the CTU-GLM scheme     !
! arXiv:0911.3410v1                         !
! ----------------------------------------- !

   implicit none

   private
   public :: read_problem_par, problem_initial_conditions, problem_pointers

   real                   :: d0, p0, vx0, vy0, vz0, A, vA, kx, ky, kz
   namelist /PROBLEM_CONTROL/  d0, p0, vx0, vy0, vz0, A, vA, kx, ky, kz

contains

   subroutine problem_pointers

      use user_hooks, only: finalize_problem

      implicit none

      finalize_problem => verify_test

   end subroutine problem_pointers

   subroutine read_problem_par

      use bcast,      only: piernik_MPI_Bcast
      use constants,  only: cbuff_len
      use dataio_pub, only: nh
      use mpisetup,   only: rbuff, cbuff, master, slave

      implicit none

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
      call piernik_MPI_Bcast(cbuff, cbuff_len)

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

   end subroutine read_problem_par
!-----------------------------------------------------------------------------
subroutine problem_initial_conditions

      use cg_leaves,   only: leaves
      use cg_list,     only: cg_list_element
      use constants,   only: xdim, ydim, zdim, LO, HI, pi
      use dataio_pub,  only: die
      use fluidindex,  only: flind
      use fluidtypes,  only: component_fluid
      use func,        only: ekin, emag
      use grid_cont,   only: grid_container
#ifndef ISO
!      use global,      only: smallei
#endif /* !ISO */
      implicit none

      class(component_fluid), pointer :: fl
      integer                         :: i, j, k
      real                            :: xi, yj, zk, phi
      type(cg_list_element),  pointer :: cgl
      type(grid_container),   pointer :: cg
      integer                         :: p
      
      real, dimension(3, 3)           :: rot_matrix, rot_matrix_inv
      real, dimension(3)              :: v_prime, b_prime, v_final, b_final
      real, dimension(3)              :: pos_vec, pos_vec_prime
      real                            :: aa, g, k_prime

      aa = atan(ky/kx)
      g = atan(cos(aa) * (kz/kx))
      k_prime = sqrt(kx**2 + ky**2 + kz**2)

      rot_matrix = reshape([cos(aa)*cos(g), sin(aa)*cos(g), sin(g), &
                           -sin(aa),         cos(aa),         0.0,   &
                           -cos(aa)*sin(g), -sin(aa)*sin(g), cos(g)], [3, 3])

      rot_matrix_inv = transpose(rot_matrix)

      do p = 1, flind%fluids

         fl => flind%all_fluids(p)%fl

         cgl => leaves%first
         do while (associated(cgl))
            cg => cgl%cg
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

                     cg%u(fl%idn,i,j,k) = d0
                     cg%u(fl%imx,i,j,k) = v_final(1)
                     cg%u(fl%imy,i,j,k) = v_final(2)
                     cg%u(fl%imz,i,j,k) = v_final(3)

                     if (fl%has_energy) then

                        cg%u(fl%ien,i,j,k) = p0/(fl%gam_1) + ekin(v_final(1), v_final(2), v_final(3), d0)

                        if (fl%is_magnetized) then

                           cg%b(xdim,i,j,k)   =  b_final(1)
                           cg%b(ydim,i,j,k)   =  b_final(2)
                           cg%b(zdim,i,j,k)   =  b_final(3)
                           cg%u(fl%ien,i,j,k) = cg%u(fl%ien,i,j,k) + emag(b_final(1), b_final(2), b_final(3))
                        endif

                     endif
                  enddo
               enddo
            enddo
            cgl => cgl%nxt
         enddo
      enddo

   end subroutine problem_initial_conditions

   subroutine verify_test

      use mpisetup,   only: master
      use dataio_pub, only: msg, printinfo
      use constants,  only: V_INFO

      implicit none

      if (master) then
         call printinfo("", V_INFO, .true.)
         write(msg, *)  'Copy compare_slices.py from problem folder to the run folder and make sure '
         write(msg, *)  'the cpaw_tst_0000.h5 file is the intial file '
         write(msg, *)  'cpaw_tst_0001.h5 is the last file .  Then run python compare_slices.py'
      endif
   end subroutine verify_test
   
end module initproblem
