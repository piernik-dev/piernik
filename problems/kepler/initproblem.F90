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
!    Initial implemetation of PIERNIK code was based on TVD split MHD code by
!    Ue-Li Pen
!        see: Pen, Arras & Wong (2003) for algorithm and
!             http://www.cita.utoronto.ca/~pen/MHD
!             for original source code "mhd.f90"
!
!    For full list of developers see $PIERNIK_HOME/license/pdt.txt
!
#include "piernik.def"

module initproblem

! Initial condition for Keplerian disk
! Written by: M. Hanasz, March 2006

   use mpisetup
   use problem_pub, only: problem_name, run_id

   real :: d0, r_max, dout, alpha
   character(len=32) :: mag_field_orient

   namelist /PROBLEM_CONTROL/  problem_name, run_id, alpha, &
                               d0,dout,r_max,mag_field_orient

   contains

!-----------------------------------------------------------------------------

   subroutine read_problem_par

      implicit none

      problem_name = 'aaa'
      run_id  = 'aa'
      d0      = 1.0
      dout    = 1.0e-4
      r_max   = 1.0
      mag_field_orient = 'none'
      alpha   = 1.0

      if(proc .eq. 0) then
         open(1,file='problem.par')
         read(unit=1,nml=PROBLEM_CONTROL)
         write(*,nml=PROBLEM_CONTROL)
         close(1)
         open(3, file='tmp.log', position='append')
         write(3,nml=PROBLEM_CONTROL)
         write(3,*)
         close(3)
      endif

      if(proc .eq. 0) then

         cbuff(1) =  problem_name
         cbuff(2) =  run_id

         rbuff(1) = d0
         rbuff(2) = dout
         rbuff(3) = r_max
         rbuff(4) = alpha

         call MPI_BCAST(cbuff, 32*buffer_dim, MPI_CHARACTER,        0, comm, ierr)
         call MPI_BCAST(ibuff,    buffer_dim, MPI_INTEGER,          0, comm, ierr)
         call MPI_BCAST(rbuff,    buffer_dim, MPI_DOUBLE_PRECISION, 0, comm, ierr)

      else

         call MPI_BCAST(cbuff, 32*buffer_dim, MPI_CHARACTER,        0, comm, ierr)
         call MPI_BCAST(ibuff,    buffer_dim, MPI_INTEGER,          0, comm, ierr)
         call MPI_BCAST(rbuff,    buffer_dim, MPI_DOUBLE_PRECISION, 0, comm, ierr)

         problem_name = cbuff(1)
         run_id       = cbuff(2)

         d0           = rbuff(1)
         dout         = rbuff(2)
         r_max        = rbuff(3)
         alpha        = rbuff(4)

      endif

   end subroutine read_problem_par

!-----------------------------------------------------------------------------

   subroutine init_prob
      use arrays,      only : u,b, dprof
      use grid,        only : x,y,z,nx,ny,nz,nzd
      use constants,   only : newtong
      use hydrostatic, only : hydrostatic_zeq
      use gravity,     only : r_smooth,r_grav,n_gravr,ptmass
      use initionized, only : idni, imxi, imyi, imzi
      use fluidindex,  only : ibx, iby, ibz
      use initfluids,  only : gamma,cs_iso
#ifndef ISO
      use initionized, only : ieni, gamma_ion, cs_ion
#endif /* !ISO */
      implicit none

      integer :: i,j,k,kmid
      real :: xi,yj,zk, rc, rs, vx, vy, vz, h2, dgdz, b0, sqr_gm, v_phi
      real :: csim2

!   Secondary parameters

      csim2 = cs_iso**2*(1.0+alpha)
#ifndef ISO
      gamma(1) = gamma_ion
#endif /* !ISO */

      sqr_gm = sqrt(newtong*ptmass)

      b0 = sqrt(2.*alpha*d0*cs_iso**2)

      do k=1, nz
         if(z(k) .lt. 0.0) kmid = k       ! the midplane is in between
      enddo                                  ! ksmid and ksmid+1

      do j = 1,ny
         yj = y(j)
         do i = 1,nx
            xi = x(i)
            rc = sqrt(xi**2+yj**2)

            if(nzd /= 1) then
               call hydrostatic_zeq(i, j, d0, csim2, dprof)
            endif

            do k = 1,nz

               vx = sqr_gm * (-yj)/(rc**2+r_smooth**2)**0.75
               vy = sqr_gm * ( xi)/(rc**2+r_smooth**2)**0.75
               vz = 0.0

               u(idni,i,j,k) = min((rc/r_grav)**n_gravr,100.0)

               if(nzd /= 1) then
                  u(idni,i,j,k) = dout + (dprof(k)-dout)/cosh(u(idni,i,j,k))
               else
                  u(idni,i,j,k) = dout + (d0 - dout)/cosh(u(idni,i,j,k))
               endif
               u(idni,i,j,k) = max(u(idni,i,j,k), smalld)
               u(imxi,i,j,k) = vx*u(idni,i,j,k)
               u(imyi,i,j,k) = vy*u(idni,i,j,k)
               u(imzi,i,j,k) = vz*u(idni,i,j,k)
#ifndef ISO
               u(ieni,i,j,k) = cs_ion**2/(gamma_ion-1.0)*u(idni,i,j,k)
               u(ieni,i,j,k) = max(u(ieni,i,j,k), smallei)
               u(ieni,i,j,k) = u(ieni,i,j,k) +0.5*(vx**2+vy**2+vz**2)*u(idni,i,j,k)
#endif /* ISO */
               if(trim(mag_field_orient) .eq. 'toroidal') then
                  b(ibx,i,j,k)   = -b0*sqrt(u(idni,i,j,k)/d0)*yj/rc
                  b(iby,i,j,k)   =  b0*sqrt(u(idni,i,j,k)/d0)*xi/rc
                  b(ibz,i,j,k)   =  0.0
               else if(trim(mag_field_orient) .eq. 'vertical') then
                  b(ibx,i,j,k)   =  0.0
                  b(iby,i,j,k)   =  0.0
                  b(ibz,i,j,k)   =  b0
               else if(trim(mag_field_orient) .eq. 'none') then
                  b(:,i,j,k)   =  0.0
               endif

#ifndef ISO
               u(ieni,i,j,k)   = u(ieni,i,j,k) +0.5*sum(b(:,i,j,k)**2,1)
#endif /* ISO */
            enddo
         enddo
      enddo

      return
   end subroutine init_prob

end module initproblem

