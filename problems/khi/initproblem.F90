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

! Initial condition for fluid flows for Kelvin-Helmholtz Instability
! based on Agertz et al. 2008
! Written by: D. Woltanski, February 2008

   use mpisetup

   character :: problem_name*32,run_id*3
   real :: chi,dbot,lpert,Mtop,Mbot,dpert,tkh,vtransf

   namelist /PROBLEM_CONTROL/  problem_name, run_id, &
                               chi,dbot,lpert,Mtop,Mbot,dpert,tkh,vtransf

   contains

!-----------------------------------------------------------------------------

   subroutine read_problem_par
      use errh
      implicit none

      character(len=100) :: par_file, tmp_log_file
      integer :: cwd_status, ierrh

      par_file = trim(cwd)//'/problem.par'
      tmp_log_file = trim(cwd)//'/tmp.log'

      problem_name = 'aaa'
      run_id  = 'aa'
      chi     = 8.0
      dbot    = 1.0
      lpert   = 0.05
      Mtop    = 0.11
      Mbot    = 0.34
      dpert   = 80.0
      tkh     = 1.70
      vtransf = 0.0

      if(proc .eq. 0) then
         open(1,file=par_file)
         read(unit=1,nml=PROBLEM_CONTROL,iostat=ierrh)
         call namelist_errh(ierrh,'PROBLEM_CONTROL')
         close(1)
         open(3, file=tmp_log_file, position='append')
         write(3,nml=PROBLEM_CONTROL)
         write(3,*)
         close(3)
      endif

      if(proc .eq. 0) then

         cbuff(1) =  problem_name
         cbuff(2) =  run_id

         rbuff(1) = chi
         rbuff(2) = dbot
         rbuff(3) = lpert
         rbuff(4) = Mtop
         rbuff(5) = Mbot
         rbuff(6) = dpert
         rbuff(7) = tkh
         rbuff(8) = vtransf

         call MPI_BCAST(cbuff, 32*buffer_dim, MPI_CHARACTER,        0, comm, ierr)
         call MPI_BCAST(ibuff,    buffer_dim, MPI_INTEGER,          0, comm, ierr)
         call MPI_BCAST(rbuff,    buffer_dim, MPI_DOUBLE_PRECISION, 0, comm, ierr)

      else

         call MPI_BCAST(cbuff, 32*buffer_dim, MPI_CHARACTER,        0, comm, ierr)
         call MPI_BCAST(ibuff,    buffer_dim, MPI_INTEGER,          0, comm, ierr)
         call MPI_BCAST(rbuff,    buffer_dim, MPI_DOUBLE_PRECISION, 0, comm, ierr)

         problem_name = cbuff(1)
         run_id       = cbuff(2)

         chi          = rbuff(1)
         dbot         = rbuff(2)
         lpert        = rbuff(3)
         Mtop         = rbuff(4)
         Mbot         = rbuff(5)
         dpert        = rbuff(6)
         tkh          = rbuff(7)
         vtransf      = rbuff(8)

      endif

   end subroutine read_problem_par

!-----------------------------------------------------------------------------

   subroutine init_prob
      use arrays,    only  : u
      use grid,      only  : ymin,ymax,x,y,nx,ny,nzd
      use constants, only  : dpi
      use initneutral
      implicit none

      real :: dtop,lambda,boxlen,p0,vtop,vbot,k0,vp,rcx,rcy,rc
      integer :: i,j

      dtop = dbot/chi
      lambda = 1./6.
      boxlen = ymax-ymin

      p0 = lambda**2 * (1.+chi)**2 /(chi*tkh**2) *dbot &
            / ( (Mtop*dsqrt(chi)+Mbot)**2 * gamma_neu)
      vtop  =  1.*Mtop*dsqrt(gamma_neu*p0/dtop)
      vbot  = -1.*Mbot*dsqrt(gamma_neu*p0/dbot)
      k0    = dpi/lambda
      vp    = (Mtop*dsqrt(chi)+Mbot)*dsqrt(gamma_neu*p0/dbot)/dpert

      do i = 1,nx
         rcx = x(i)
         do j = 1,ny
            rcy = y(j)
            rc=rcy-0.5*boxlen
            if(rc .gt. 0.0) then
               u(idnn,i,j,:) = dtop
               u(imxn,i,j,:) = vtop*dtop
            endif
            if(rc .le. 0.0) then
               u(idnn,i,j,:) = dbot
               u(imxn,i,j,:) = vbot*dbot
            endif
            if(abs(rc) .lt. lpert) then
               u(imyn,i,j,:) = vp*sin(k0*rcx)*u(idnn,i,j,:)
            endif
            if(nzd /= 1) then
               u(imzn,i,j,:) = vtransf*u(1,i,j,:)
            else
               u(imzn,i,j,:) = 0.0
            endif
            u(ienn,i,j,:) = p0/(gamma_neu-1.0)
         enddo
      enddo

      return
   end subroutine init_prob

end module initproblem
