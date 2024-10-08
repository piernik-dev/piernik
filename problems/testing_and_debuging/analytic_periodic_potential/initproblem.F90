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

!> \brief Implementation of some analytic, periodic potentials

module initproblem

   use constants, only: cbuff_len, dsetnamelen, ndims

   implicit none

   private
   public :: read_problem_par, problem_initial_conditions, problem_pointers

   ! namelist parameters
   character(len=cbuff_len) :: type !< type of potential
   real                     :: a    !< proportionality constant
   real, dimension(ndims)   :: kx   !< wave number
   integer(kind=4)          :: n    !< exponent
   real :: ref_thr                  !< refinement threshold

   namelist /PROBLEM_CONTROL/ type, a, n, kx, ref_thr

   ! private data
   character(len=dsetnamelen), parameter :: apot_n = "apot" !< name of the analytical potential field
   character(len=dsetnamelen), parameter :: asrc_n = "asrc" !< name of the source field used for "ares" calculation (auxiliary space)
   character(len=dsetnamelen), parameter :: ares_n = "ares" !< name of the numerical residuum with respect to analytical potential field

   !< recognized types of potential
   enum, bind(c)
      enumerator :: SIN_MUL
      enumerator :: SIN_ADD
   end enum
   integer(kind=4) :: itype  !< type of potential encoded in integers

contains

!> \brief Set up custom pointers to tweak the code execution according to our needs

   subroutine problem_pointers

#ifdef HDF5
      use dataio_user, only: user_vars_hdf5
#endif /* HDF5 */
      use user_hooks,  only: finalize_problem

      implicit none

      finalize_problem => finalize_problem_app
#ifdef HDF5
      user_vars_hdf5   => app_error_vars
#endif /* HDF5 */

   end subroutine problem_pointers

!> \brief Read the runtime parameters specified in the namelist and set up some auxiliary values

   subroutine read_problem_par

      use bcast,            only: piernik_MPI_Bcast
      use cg_list_global,   only: all_cg
      use constants,        only: INVALID, cbuff_len, pi
      use dataio_pub,       only: die, msg, nh
      use mpisetup,         only: rbuff, ibuff, cbuff, master, slave
      use multigridvars,    only: ord_prolong
      use named_array_list, only: qna
      use unified_ref_crit_list, only: urc_list

      implicit none

      ! namelist default parameter values
      type       = "none"
      a          = 0.
      kx(:)      = pi
      n          = 0
      ref_thr    = 3e-2     !< Refine if density difference is greater than this value

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

         rbuff(1) = a
         rbuff(2:4) = kx(:)
         rbuff(5) = ref_thr

         ibuff(1) = n

         cbuff(1) = type

      endif

      call piernik_MPI_Bcast(ibuff)
      call piernik_MPI_Bcast(rbuff)
      call piernik_MPI_Bcast(cbuff, cbuff_len)

      if (slave) then

         a         = rbuff(1)
         kx(:)     = rbuff(2:4)
         ref_thr   = rbuff(5)

         n    = ibuff(1)

         type = cbuff(1)

      endif

      if (n<=0) call die("[initproblem:read_problem_par] n must be positive")
      if (n==1) call die("[initproblem:read_problem_par] n == 1 not properly implemented") ! requires Dirac deltas along the domain boundary

      itype = INVALID
      select case (type)
         case ("none")
            call die("[initproblem:read_problem_par] Nothing is set by default")
         case ("sin*")
            itype = SIN_MUL
         case ("sin+")
            itype = SIN_ADD
         case default
            write(msg, '(3a)')"[initproblem:read_problem_par] unrecognized potential type '", type, "'."
            call die(msg)
      end select

      call all_cg%reg_var(apot_n, ord_prolong = ord_prolong)
      call all_cg%reg_var(ares_n)
      call all_cg%reg_var(asrc_n)

      call urc_list%add_user_urcv(qna%ind(apot_n), INVALID, ref_thr, 0., "grad", .true.)

   end subroutine read_problem_par

!> \brief Set up the initial conditions. Note that this routine can be called multiple times during initial iterations of refinement structure

   subroutine problem_initial_conditions

      use cg_list,           only: cg_list_element
      use cg_leaves,         only: leaves
      use constants,         only: ndims, xdim, ydim, zdim
      use dataio_pub,        only: die, msg, printinfo
      use domain,            only: dom
      use fluidindex,        only: iarr_all_dn, iarr_all_mx, iarr_all_my, iarr_all_mz
      use global,            only: dirty_debug, no_dirty_checks
      use grid_cont,         only: grid_container
      use named_array_list,  only: qna
      use mpisetup,          only: master
      use multigrid_Laplace, only: residual
      use units,             only: fpiG

      implicit none

      integer                        :: i, j, k
      real                           :: dens
      real, dimension(ndims)         :: s
      type(cg_list_element), pointer :: cgl
      type(grid_container),  pointer :: cg

      call compute_app_potential

      cgl => leaves%first
      do while (associated(cgl))
         cg => cgl%cg

         do k = cg%ks, cg%ke
            do j = cg%js, cg%je
               do i = cg%is, cg%ie

                  dens = 0.
                  select case (itype)
                     case (SIN_MUL)
                        s(:) = 1.
                        where (dom%has_dir(:)) s(:) = sin([kx(xdim)*cg%x(i), kx(ydim)*cg%y(j), kx(zdim)*cg%z(k)])**n
                        if (dom%has_dir(xdim)) dens = dens + kx(xdim)**2*s(ydim) * s(zdim) * fsinadd(kx(xdim)*cg%x(i))
                        if (dom%has_dir(ydim)) dens = dens + kx(ydim)**2*s(xdim) * s(zdim) * fsinadd(kx(ydim)*cg%y(j))
                        if (dom%has_dir(zdim)) dens = dens + kx(zdim)**2*s(xdim) * s(ydim) * fsinadd(kx(zdim)*cg%z(k))
                     case (SIN_ADD)
                        if (dom%has_dir(xdim)) dens = dens + kx(xdim)**2*fsinadd(kx(xdim)*cg%x(i))
                        if (dom%has_dir(ydim)) dens = dens + kx(ydim)**2*fsinadd(kx(ydim)*cg%y(j))
                        if (dom%has_dir(zdim)) dens = dens + kx(zdim)**2*fsinadd(kx(zdim)*cg%z(k))
                     case default
                        write(msg, '(3a)')"[initproblem:problem_initial_conditions] unrecognized potential type '", type, "'."
                        call die(msg)
                  end select

                  cg%u(iarr_all_dn, i, j, k) = dens/fpiG

               enddo
            enddo
         enddo

         cg%u(iarr_all_mx, RNG) = 0.0
         cg%u(iarr_all_my, RNG) = 0.0
         cg%u(iarr_all_mz, RNG) = 0.0

#ifdef MAGNETIC
         call cg%set_constant_b_field([0., 0., 0.])
#endif /* MAGNETIC */

         cgl => cgl%nxt
      enddo

      ! Compute residuum for analytical solution and print its norm
      cgl => leaves%first
      do while (associated(cgl))
         cg => cgl%cg
         if (dirty_debug .and. .not. no_dirty_checks) then
            cg%wa(RNG) = cg%q(qna%ind(apot_n))%arr(RNG)
            cg%q(qna%ind(apot_n))%arr = huge(1.)
            cg%q(qna%ind(apot_n))%arr(RNG) = cg%wa(RNG)
            cg%q(qna%ind(asrc_n))%arr = huge(1.)
         endif
         cg%q(qna%ind(asrc_n))%arr(RNG) = fpiG * sum(cg%u(iarr_all_dn, RNG), dim=1)
         cgl => cgl%nxt
      enddo
      call residual(leaves, qna%ind(asrc_n), qna%ind(apot_n), qna%ind(ares_n))
      call leaves%check_dirty(qna%ind(ares_n), "a-residual")

      write(msg,'(a,f13.10)')"[initproblem:problem_initial_conditions] Analytical norm residual/source= ",leaves%norm_sq(qna%ind(ares_n))/leaves%norm_sq(qna%ind(asrc_n))
      if (master) call printinfo(msg)

   end subroutine problem_initial_conditions

   real function fsinadd(x)

      implicit none

      real, intent(in) :: x

      real :: sx

      sx = sin(x)
      fsinadd = - sx**n
      if (n > 1) fsinadd = fsinadd + (n-1)*sx**(n-2)* (1. - sx**2)
      fsinadd = a * n * fsinadd

   end function fsinadd

!>
!! \brief Calculate analytical potential for given spheroid.
!!
!! \details We assume that the spheroid is fully contained in the computational domain here.
!! Oblate potential formula from Ricker_2008ApJS..176..293R
!! Prolate potential formula from http://scienceworld.wolfram.com/physics/ProlateSpheroidGravitationalPotential.html
!<
   subroutine compute_app_potential

      use cg_list,          only: cg_list_element
      use cg_leaves,        only: leaves
      use constants,        only: xdim, ydim, zdim
      use dataio_pub,       only: die, msg
      use domain,           only: dom
      use grid_cont,        only: grid_container
      use named_array_list, only: qna

      implicit none

      integer                        :: i, j, k
      integer(kind=4)                :: apot_i
      real                           :: pot
      type(cg_list_element), pointer :: cgl
      type(grid_container),  pointer :: cg

      apot_i = qna%ind(apot_n)

      cgl => leaves%first
      do while (associated(cgl))
         cg => cgl%cg

         do k = cg%ks, cg%ke
            do j = cg%js, cg%je
               do i = cg%is, cg%ie
                  pot = 0.
                  select case (itype)
                     case (SIN_MUL)
                        pot = 1.
                        if (dom%has_dir(xdim)) pot = pot * sin(kx(xdim)*cg%x(i))**n
                        if (dom%has_dir(ydim)) pot = pot * sin(kx(ydim)*cg%y(j))**n
                        if (dom%has_dir(zdim)) pot = pot * sin(kx(zdim)*cg%z(k))**n
                     case (SIN_ADD)
                        if (dom%has_dir(xdim)) pot = pot + sin(kx(xdim)*cg%x(i))**n
                        if (dom%has_dir(ydim)) pot = pot + sin(kx(ydim)*cg%y(j))**n
                        if (dom%has_dir(zdim)) pot = pot + sin(kx(zdim)*cg%z(k))**n
                     case default
                        write(msg, '(3a)')"[initproblem:problem_initial_conditions] unrecognized potential type '", type, "'."
                        call die(msg)
                  end select
                  cg%q(apot_i)%arr(i, j, k) = pot
               enddo
            enddo
         enddo

         cgl => cgl%nxt
      enddo

      ! for sure this can be computed analytically :-)
      call leaves%subtract_average(apot_i)

   end subroutine compute_app_potential

!> \brief Compute the L2 error norm of the error of computed potential with respect to the analytical solution

   subroutine finalize_problem_app

      use allreduce,        only: piernik_MPI_Allreduce
      use cg_list,          only: cg_list_element
      use cg_leaves,        only: leaves
      use constants,        only: GEO_RPZ, pSUM, pMIN, pMAX
      use dataio_pub,       only: msg, printinfo, warn
      use domain,           only: dom
      use grid_cont,        only: grid_container
      use mpisetup,         only: master
      use named_array_list, only: qna

      implicit none

      integer                        :: i, j, k, apot_i
      real, dimension(2)             :: norm, dev
      real                           :: potential, fac
      type(cg_list_element), pointer :: cgl
      type(grid_container),  pointer :: cg

      fac = 1.
      norm(:) = 0.
      dev(1) = huge(1.0)
      dev(2) = -dev(1)
      apot_i = qna%ind(apot_n)

      cgl => leaves%first
      do while (associated(cgl))
         cg => cgl%cg

         if (.not. qna%exists(apot_n)) then
            if (master) call warn("[initproblem:finalize_problem_app] Cannot compare results with the analytical potential.")
            return
         endif

         do k = cg%ks, cg%ke
            do j = cg%js, cg%je
               do i = cg%is, cg%ie
                  potential = cg%q(apot_i)%arr(i, j, k)
                  if (dom%geometry_type == GEO_RPZ) fac = cg%x(i)
                  norm(1) = norm(1) + (potential - cg%sgp(i, j, k))**2 * fac
                  norm(2) = norm(2) + potential**2 * fac
                  dev(1) = min(dev(1), (potential - cg%sgp(i, j, k))/potential)
                  dev(2) = max(dev(2), (potential - cg%sgp(i, j, k))/potential)
               enddo
            enddo
         enddo
         cgl => cgl%nxt
      enddo

      call piernik_MPI_Allreduce(norm,   pSUM)
      call piernik_MPI_Allreduce(dev(1), pMIN)
      call piernik_MPI_Allreduce(dev(2), pMAX)

      if (master) then
         write(msg,'(a,f12.6,a,2f12.6)')"[initproblem:finalize_problem_app] L2 error norm = ", sqrt(norm(1)/norm(2)), ", min and max error = ", dev(1:2)
         call printinfo(msg)
      endif

   end subroutine finalize_problem_app

!>
!! \brief This routine provides the "apot" and "errp" variablesvalues to be dumped to the .h5 file
!!
!! \details
!! * "apot" is the analytical potential solution for cell centers
!! * "errp" is the difference between analytical potential and computed potential
!<

   subroutine app_error_vars(var, tab, ierrh, cg)

      use dataio_pub,       only: die
      use grid_cont,        only: grid_container
      use func,             only: operator(.notequals.)
      use named_array_list, only: qna

      implicit none

      character(len=*),               intent(in)    :: var
      real, dimension(:,:,:),         intent(inout) :: tab
      integer,                        intent(inout) :: ierrh
      type(grid_container), pointer,  intent(in)    :: cg

      if (.not. qna%exists(apot_n)) call die("[initproblem:app_error_vars] Cannot find apot_n")

      ierrh = 0
      select case (trim(var))
         case ("errp")
            tab(:,:,:) = cg%q(qna%ind(apot_n))%span(cg%ijkse) - cg%sgp(RNG)
         case ("relerr")
            where (cg%q(qna%ind(apot_n))%span(cg%ijkse) .notequals. 0.)
               tab(:,:,:) = cg%sgp(RNG)/cg%q(qna%ind(apot_n))%span(cg%ijkse) -1.
            elsewhere
               tab(:,:,:) = 0.
            endwhere
         case default
            ierrh = -1
      end select

   end subroutine app_error_vars

end module initproblem
