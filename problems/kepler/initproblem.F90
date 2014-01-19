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

! Initial condition for Keplerian disk
! Written by: M. Hanasz, March 2006

   use constants,    only: cbuff_len

   implicit none

   private
   public :: read_problem_par, problem_initial_conditions, problem_pointers

   real                     :: d0, r_max, dout, alpha, r_in, r_out, f_in, f_out
   character(len=cbuff_len) :: mag_field_orient

   namelist /PROBLEM_CONTROL/  alpha, d0, dout, r_max, mag_field_orient, r_in, r_out, f_in, f_out

contains
!-----------------------------------------------------------------------------
   subroutine problem_pointers

      implicit none

   end subroutine problem_pointers
!-----------------------------------------------------------------------------
   subroutine read_problem_par

      use dataio_pub,            only: nh      ! QA_WARN required for diff_nml
      use mpisetup,              only: cbuff, rbuff, master, slave, piernik_MPI_Bcast

      implicit none

      d0               = 1.0
      dout             = 1.0e-4
      r_max            = 1.0
      mag_field_orient = 'none'
      alpha            = 1.0

      r_in             = 0.5
      f_in             = 10.0
      r_out            = 2.1
      f_out            = 0.0

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

         cbuff(1) = mag_field_orient

         rbuff(1) = d0
         rbuff(2) = dout
         rbuff(3) = r_max
         rbuff(4) = alpha
         rbuff(5) = r_in
         rbuff(6) = r_out
         rbuff(7) = f_in
         rbuff(8) = f_out

      endif

      call piernik_MPI_Bcast(cbuff, cbuff_len)
      call piernik_MPI_Bcast(rbuff)

      if (slave) then

         mag_field_orient = cbuff(1)

         d0               = rbuff(1)
         dout             = rbuff(2)
         r_max            = rbuff(3)
         alpha            = rbuff(4)
         r_in             = rbuff(5)
         r_out            = rbuff(6)
         f_in             = rbuff(7)
         f_out            = rbuff(8)

      endif

   end subroutine read_problem_par
!-----------------------------------------------------------------------------
   subroutine problem_initial_conditions

      use cg_leaves,          only: leaves
      use cg_list,            only: cg_list_element
      use constants,          only: xdim, ydim, zdim, GEO_XYZ, LO, HI
      use dataio_pub,         only: die
      use domain,             only: dom, is_multicg
      use fluidindex,         only: flind
      use fluidtypes,         only: component_fluid
      use gravity,            only: r_smooth, r_grav, n_gravr, ptmass
      use grid_cont,          only: grid_container
      use hydrostatic,        only: hydrostatic_zeq_densmid, set_default_hsparams, dprof
      use units,              only: newtong

      implicit none

      integer                         :: i, j, k, kmid
      real                            :: xi, yj, rc, vx, vy, vz, b0, sqr_gm
      real                            :: csim2

      class(component_fluid), pointer :: fl
      type(cg_list_element),  pointer :: cgl
      type(grid_container),   pointer :: cg

!   Secondary parameters
      cgl => leaves%first
      do while (associated(cgl))
         cg => cgl%cg

         if (is_multicg) call die("[initproblem:problem_initial_conditions] multiple grid pieces per procesor not implemented yet") !nontrivial kmid, allocate

         sqr_gm = sqrt(newtong*ptmass)
         do k = cg%lhn(zdim,LO), cg%lhn(zdim,HI)
            if (cg%z(k) < 0.0) kmid = k       ! the midplane is in between ksmid and ksmid+1
         enddo

         if (associated(flind%ion) .and. dom%geometry_type == GEO_XYZ) then
            fl => flind%ion
            csim2 = fl%cs2*(1.0+alpha)
            b0    = sqrt(2.*alpha*d0*fl%cs2)

            if (dom%has_dir(zdim)) call set_default_hsparams(cg)

            do j = cg%lhn(ydim,LO), cg%lhn(ydim,HI)
               yj = cg%y(j)
               do i = cg%lhn(xdim,LO), cg%lhn(xdim,HI)
                  xi = cg%x(i)
                  rc = sqrt(xi**2+yj**2)

                  if (dom%has_dir(zdim)) call hydrostatic_zeq_densmid(i, j, d0, csim2)

                  do k = cg%lhn(zdim,LO), cg%lhn(zdim,HI)

                     vx = sqr_gm * (-yj)/(rc**2+r_smooth**2)**0.75
                     vy = sqr_gm * ( xi)/(rc**2+r_smooth**2)**0.75
                     vz = 0.0

                     cg%u(fl%idn,i,j,k) = min((rc/r_grav)**n_gravr,100.0)

                     if (dom%has_dir(zdim)) then
                        cg%u(fl%idn,i,j,k) = dprof(k)/cosh(cg%u(fl%idn,i,j,k))
                        cg%u(fl%idn,i,j,k) = max(cg%u(fl%idn,i,j,k), dout)
                     else
                        cg%u(fl%idn,i,j,k) = dout + (d0 - dout)/cosh(cg%u(fl%idn,i,j,k))
                     endif
                     cg%u(fl%idn,i,j,k) = cg%u(fl%idn,i,j,k)
                     cg%u(fl%imx,i,j,k) = vx*cg%u(fl%idn,i,j,k)
                     cg%u(fl%imy,i,j,k) = vy*cg%u(fl%idn,i,j,k)
                     cg%u(fl%imz,i,j,k) = vz*cg%u(fl%idn,i,j,k)
                     if (fl%ien > 0) then
                        cg%u(fl%ien,i,j,k) = fl%cs2/(fl%gam_1)*cg%u(fl%idn,i,j,k)
!                     cg%u(fl%ien,i,j,k) = max(cg%u(fl%ien,i,j,k), smallei)
                        cg%u(fl%ien,i,j,k) = cg%u(fl%ien,i,j,k) + 0.5 * (vx**2 + vy**2 + vz**2) * cg%u(fl%idn,i,j,k)
                     endif
                     if (trim(mag_field_orient) == 'toroidal') then
                        cg%b(xdim,i,j,k)   = -b0*sqrt(cg%u(fl%idn,i,j,k)/d0)*yj/rc
                        cg%b(ydim,i,j,k)   =  b0*sqrt(cg%u(fl%idn,i,j,k)/d0)*xi/rc
                        cg%b(zdim,i,j,k)   =  0.0
                     else if (trim(mag_field_orient) == 'vertical') then
                        cg%b(xdim,i,j,k)   =  0.0
                        cg%b(ydim,i,j,k)   =  0.0
                        cg%b(zdim,i,j,k)   =  b0
                     else if (trim(mag_field_orient) == 'none') then
                        cg%b(:,i,j,k)     =  0.0
                     endif

                     if (fl%ien > 0) cg%u(fl%ien,i,j,k)   = cg%u(fl%ien,i,j,k) + 0.5*sum(cg%b(:,i,j,k)**2,1)
                  enddo
               enddo
            enddo
         else
            call die("[initproblem:problem_initial_conditions] I don't know what to do... :/")
         endif

         cgl => cgl%nxt
      enddo

   end subroutine problem_initial_conditions

end module initproblem
! vim: set tw=120:
