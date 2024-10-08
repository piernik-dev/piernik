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

! Initial condition for Sedov-Taylor explosion
! Written by: M. Hanasz, March 2006

   implicit none

   private
   public  :: read_problem_par, problem_initial_conditions, problem_pointers

   integer(kind=4) :: n_sn
   real            :: d0, p0, bx0, by0, bz0, Eexpl, x0, y0, z0, r0, smooth, dt_sn, r, dtrig
   real :: ref_thr   !< refinement threshold
   real :: ref_eps   !< smoother filter

   namelist /PROBLEM_CONTROL/ d0, p0, bx0, by0, bz0, Eexpl, x0, y0, z0, r0, smooth, n_sn, dt_sn, ref_thr, ref_eps, dtrig

contains

!-----------------------------------------------------------------------------

   subroutine problem_pointers

      use dataio_user, only: user_tsl
      use user_hooks,  only: problem_domain_update, late_initial_conditions

      implicit none

      user_tsl                => sedov_tsl
      problem_domain_update   => sedov_dist_to_edge
      late_initial_conditions => sedov_late_init

   end subroutine problem_pointers

!-----------------------------------------------------------------------------

   subroutine read_problem_par

      use bcast,                 only: piernik_MPI_Bcast
      use constants,             only: DST
      use dataio_pub,            only: msg, printinfo, die, nh
      use domain,                only: dom
      use fluidindex,            only: flind
      use mpisetup,              only: ibuff, rbuff, master, slave
      use named_array_list,      only: wna
      use unified_ref_crit_list, only: urc_list
      use user_hooks,            only: problem_domain_update

      implicit none

      integer :: p, id

      d0      = 1.0
      dtrig   = -1.5
      p0      = 1.e-3
      bx0     =   0.
      by0     =   0.
      bz0     =   0.
      Eexpl   = 1.e5
      x0      = 0.0
      y0      = 0.0
      z0      = 0.0
      r0      = minval(dom%L_(:)/dom%n_d(:), mask=dom%has_dir(:))/2.
      smooth  = 0.0 ! smoothing relative to r0, smooth == 1 => profile like cos(r)**2, without uniform core
      n_sn    = 1
      dt_sn   = 0.0
      ref_thr   = 0.5 ! Lower this value if you want to track the shock wave as it gets weaker
      ref_eps   = 0.01

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

         rbuff(1) = d0
         rbuff(2) = p0
         rbuff(3) = bx0
         rbuff(4) = by0
         rbuff(5) = bz0
         rbuff(6) = Eexpl
         rbuff(7) = x0
         rbuff(8) = y0
         rbuff(9) = z0
         rbuff(10)= r0
         rbuff(11)= dt_sn
         rbuff(12)= ref_thr
         rbuff(14)= dtrig
         rbuff(15)= smooth
         rbuff(16)= ref_eps

         ibuff(1) = n_sn

      endif

      call piernik_MPI_Bcast(ibuff)
      call piernik_MPI_Bcast(rbuff)

      if (slave) then

         d0           = rbuff(1)
         p0           = rbuff(2)
         bx0          = rbuff(3)
         by0          = rbuff(4)
         bz0          = rbuff(5)
         Eexpl        = rbuff(6)
         x0           = rbuff(7)
         y0           = rbuff(8)
         z0           = rbuff(9)
         r0           = rbuff(10)
         dt_sn        = rbuff(11)
         ref_thr      = rbuff(12)
         dtrig        = rbuff(14)
         smooth       = rbuff(15)
         ref_eps      = rbuff(16)

         n_sn         = ibuff(1)

      endif

      if (flind%energ < flind%fluids) call die("[initproblem:problem_initial_conditions] Not all fluids are adiabatic!")
      do p = 1, flind%energ
         if (flind%all_fluids(p)%fl%tag == DST) call die("[initproblem:problem_initial_conditions] This setup is not suitable for dust!")
         write(msg, '(a,i2)')"Working with fluid#", flind%all_fluids(p)%fl%tag
         if (master) call printinfo(msg)
      enddo

      if (ref_thr > 0.) then
         do id = 1, flind%energ
            call urc_list%add_user_urcv(wna%fi, flind%all_fluids(id)%fl%ien, ref_thr, ref_eps, "Loechner", .true.)
         enddo
      else
         if (master) call printinfo("[initproblem:problem_initial_conditions] automatic refinement disabled by negative ref_thr")
      endif

      if (dtrig < 0.) nullify(problem_domain_update)

   end subroutine read_problem_par
!-----------------------------------------------------------------------------
   subroutine problem_initial_conditions

      use cg_cost_data, only: I_IC
      use cg_leaves,    only: leaves
      use cg_list,      only: cg_list_element
      use constants,    only: ION, xdim, ydim, zdim, LO, HI, pi, ndims
      use domain,       only: dom
      use fluidindex,   only: flind
      use func,         only: operator(.notequals.)
      use grid_cont,    only: grid_container

      implicit none

      integer, parameter              :: isub = 4
      integer                         :: i, j, k, p, ii, jj, kk
      type(cg_list_element),  pointer :: cgl
      type(grid_container),   pointer :: cg
      real :: x, y, z, s

      ! BEWARE:
      !  3 triple loop are completely unnecessary here, but this problem serves
      !  as an educational example

      do p = 1, flind%energ
         associate(fl => flind%all_fluids(p)%fl)

! Uniform equilibrium state

         cgl => leaves%first
         do while (associated(cgl))
            cg => cgl%cg
            call cg%costs%start

            do k = cg%lhn(zdim,LO), cg%lhn(zdim,HI)
               do j = cg%lhn(ydim,LO), cg%lhn(ydim,HI)
                  do i = cg%lhn(xdim,LO), cg%lhn(xdim,HI)
                     cg%u(fl%idn,i,j,k) = d0
                     cg%u(fl%imx,i,j,k) = 0.0
                     cg%u(fl%imy,i,j,k) = 0.0
                     cg%u(fl%imz,i,j,k) = 0.0
                     cg%u(fl%ien,i,j,k) = p0/(fl%gam_1)
                     cg%u(fl%ien,i,j,k) = cg%u(fl%ien,i,j,k) + 0.5*(cg%u(fl%imx,i,j,k)**2 +cg%u(fl%imy,i,j,k)**2 + cg%u(fl%imz,i,j,k)**2)/cg%u(fl%idn,i,j,k)
                  enddo
               enddo
            enddo

! Explosion

            do k = cg%lhn(zdim,LO), cg%lhn(zdim,HI)
               do j = cg%lhn(ydim,LO), cg%lhn(ydim,HI)
                  do i = cg%lhn(xdim,LO), cg%lhn(xdim,HI)
                     r = (cg%x(i)-x0)**2 + (cg%y(j)-y0)**2 + (cg%z(k)-z0)**2
                     if ( r < ((1 + smooth) * (r0 + maxval(cg%dl, mask=dom%has_dir) ))**2) then
                        do ii = 1, isub
                           x = cg%x(i)-x0
                           if (dom%has_dir(xdim)) x = x + cg%dx*(2*ii -isub - 1)/real(2*isub)
                           do jj = 1, isub
                              y = cg%y(j)-y0
                              if (dom%has_dir(ydim)) y = y + cg%dy*(2*jj -isub - 1)/real(2*isub)
                              do kk = 1, isub
                                 z = cg%z(k)-z0
                                 if (dom%has_dir(zdim)) z = z + cg%dx*(2*kk -isub - 1)/real(2*isub)
                                 r = sqrt(x*x+y*y+z*z)/r0 - 1.
                                 if (r < -smooth) then
                                    s = 1.
                                 else if (r > smooth) then
                                    s = 0.
                                 else
                                    s = 0.5
                                    if (smooth .notequals. 0.) s = 0.5 * (1. - sin(pi / 2. * r/smooth))
                                 endif
                                 if (s > 0.) cg%u(fl%ien,i,j,k) = cg%u(fl%ien,i,j,k) + Eexpl/isub**ndims * s
                              enddo
                           enddo
                        enddo
                     endif
                  enddo
               enddo
            enddo

            if (fl%tag == ION) then
               call cg%set_constant_b_field([bx0, by0, bz0])
               do k = cg%lhn(zdim,LO), cg%lhn(zdim,HI)
                  do j = cg%lhn(ydim,LO), cg%lhn(ydim,HI)
                     do i = cg%lhn(xdim,LO), cg%lhn(xdim,HI)
                        cg%u(fl%ien,i,j,k) = cg%u(fl%ien,i,j,k) + 0.5*(cg%b(xdim,i,j,k)**2 + cg%b(ydim,i,j,k)**2 + cg%b(zdim,i,j,k)**2)
                     enddo
                  enddo
               enddo
            endif

            call cg%costs%stop(I_IC)
            cgl => cgl%nxt
         enddo

#ifdef TRACER
#error Check this code and move into loop over cg above
         do k = cg%lhn(zdim,LO), cg%lhn(zdim,HI)
            do j = cg%lhn(ydim,LO), cg%lhn(ydim,HI)
               do i = cg%lhn(xdim,LO), cg%lhn(xdim,HI)
                  r = sqrt( (cg%x(i)-x0)**2 + (cg%y(j)-y0)**2 + (cg%z(k)-z0)**2 )
                  if ( r**2 < r0**2) then
                     cg%u(flind%trc%pos, i, j, k)   = 1.0
                  else
                     cg%u(flind%trc%pos, i, j, k)   = 0.0
                  endif
               enddo
            enddo
         enddo
#endif /* TRACER */
         end associate
      enddo

   end subroutine problem_initial_conditions
!-----------------------------------------------------------------------------
   subroutine sedov_tsl(user_vars, tsl_names)

      use allreduce,   only: piernik_MPI_Allreduce
      use constants,   only: pSUM
      use diagnostics, only: pop_vector
      use mpisetup,    only: proc, master

      implicit none

      real, dimension(:), intent(inout), allocatable                       :: user_vars
      character(len=*), dimension(:), intent(inout), allocatable, optional :: tsl_names
      real :: output

      if (present(tsl_names)) then
         call pop_vector(tsl_names, len(tsl_names(1)), ["foobar_sedov"])    !   add to header
      else
         ! do mpi stuff here...
         output = real(proc,8)
         call piernik_MPI_Allreduce(output, pSUM)
         if (master) call pop_vector(user_vars,[output])                 !   pop value
      endif

   end subroutine sedov_tsl

!> \brief Find how close is the shockwave to the external edges and call expansion routine if necessary

   subroutine sedov_dist_to_edge

      use cg_leaves,      only: leaves
      use cg_expand_base, only: expand_base
      use cg_list,        only: cg_list_element
      use constants,      only: xdim, ydim, zdim, LO, HI
      use domain,         only: dom
      use fluidindex,     only: iarr_all_dn

      implicit none

      type(cg_list_element),  pointer :: cgl
      real, dimension(xdim:zdim, LO:HI) :: ddist
      integer :: i
      integer, parameter :: iprox = 2

      if (dtrig < 0.) return

      ddist = huge(1.)
      cgl => leaves%first
      do while (associated(cgl))
         if (any(cgl%cg%ext_bnd)) then
            !> \todo roll it to a nested loop
            if (dom%has_dir(xdim)) then
               if (cgl%cg%ext_bnd(xdim, LO)) then
                  do i = cgl%cg%is, cgl%cg%ie
                     if (any(cgl%cg%u(iarr_all_dn, i, :, :) > d0*dtrig)) then
                        ddist(xdim, LO) = min(ddist(xdim, LO), (cgl%cg%x(i) - cgl%cg%fbnd(xdim, LO))/cgl%cg%dx)
                        exit
                     endif
                  enddo
               endif
               if (cgl%cg%ext_bnd(xdim, HI)) then
                  do i = cgl%cg%ie, cgl%cg%is, -1
                     if (any(cgl%cg%u(iarr_all_dn, i, :, :) > d0*dtrig)) then
                        ddist(xdim, HI) = min(ddist(xdim, HI), (cgl%cg%fbnd(xdim, HI) - cgl%cg%x(i))/cgl%cg%dx)
                        exit
                     endif
                  enddo
               endif
            endif

            if (dom%has_dir(ydim)) then
               if (cgl%cg%ext_bnd(ydim, LO)) then
                  do i = cgl%cg%js, cgl%cg%je
                     if (any(cgl%cg%u(iarr_all_dn, :, i, :) > d0*dtrig)) then
                        ddist(ydim, LO) = min(ddist(ydim, LO), (cgl%cg%y(i) - cgl%cg%fbnd(ydim, LO))/cgl%cg%dy)
                        exit
                     endif
                  enddo
               endif
               if (cgl%cg%ext_bnd(ydim, HI)) then
                  do i = cgl%cg%je, cgl%cg%js, -1
                     if (any(cgl%cg%u(iarr_all_dn, :, i, :) > d0*dtrig)) then
                        ddist(ydim, HI) = min(ddist(ydim, HI), (cgl%cg%fbnd(ydim, HI) - cgl%cg%y(i))/cgl%cg%dy)
                        exit
                     endif
                  enddo
               endif
            endif

            if (dom%has_dir(zdim)) then
               if (cgl%cg%ext_bnd(zdim, LO)) then
                  do i = cgl%cg%ks, cgl%cg%ke
                     if (any(cgl%cg%u(iarr_all_dn, :, :, i) > d0*dtrig)) then
                        ddist(zdim, LO) = min(ddist(zdim, LO), (cgl%cg%z(i) - cgl%cg%fbnd(zdim, LO))/cgl%cg%dz)
                        exit
                     endif
                  enddo
               endif
               if (cgl%cg%ext_bnd(zdim, HI)) then
                  do i = cgl%cg%ke, cgl%cg%ks, -1
                     if (any(cgl%cg%u(iarr_all_dn, :, :, i) > d0*dtrig)) then
                        ddist(zdim, HI) = min(ddist(zdim, HI), (cgl%cg%fbnd(zdim, HI) - cgl%cg%z(i))/cgl%cg%dz)
                        exit
                     endif
                  enddo
               endif
            endif
         endif
         cgl => cgl%nxt
      enddo

      call expand_base(ddist(:,:) < iprox)

   end subroutine sedov_dist_to_edge

!> \brief Performa late initialization of the cg added after domain expansion

   subroutine sedov_late_init

      use cg_list,          only: cg_list_element
      use cg_list_dataop,   only: expanded_domain
      use constants,        only: xdim, ydim, zdim, ION, psi_n
      use fluidindex,       only: flind
      use func,             only: ekin, emag
      use named_array_list, only: qna

      implicit none

      type(cg_list_element),  pointer :: cgl

      integer :: p

      cgl => expanded_domain%first
      do while (associated(cgl))
         if (.not. cgl%cg%is_old) then
            do p = 1, flind%energ
               associate (fl => flind%all_fluids(p)%fl)
                  cgl%cg%u(fl%idn, :, :, :) = d0
                  cgl%cg%u(fl%imx, :, :, :) = 0.
                  cgl%cg%u(fl%imy, :, :, :) = 0.
                  cgl%cg%u(fl%imz, :, :, :) = 0.
                  cgl%cg%u(fl%ien, :, :, :) = p0/(fl%gam_1) + ekin(cgl%cg%u(fl%imx, :, :, :), cgl%cg%u(fl%imy, :, :, :), cgl%cg%u(fl%imz, :, :, :), cgl%cg%u(fl%idn, :, :, :))
                  if (fl%tag == ION) then
                     call cgl%cg%set_constant_b_field([bx0, by0, bz0])
                     cgl%cg%u(fl%ien, :, :, :) = cgl%cg%u(fl%ien, :, :, :) + emag(cgl%cg%b(xdim, :, :, :), cgl%cg%b(ydim, :, :, :), cgl%cg%b(zdim, :, :, :)**2)
                  endif
               end associate
            enddo
         endif
         cgl => cgl%nxt
      enddo

      ! ToDo: move this call to some other place where it would be automagically called
      if (qna%exists(psi_n)) call expanded_domain%set_q_value(qna%ind(psi_n), 0.)

   end subroutine sedov_late_init

end module initproblem
