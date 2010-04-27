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

  real :: d0,p0,bx0,by0,bz0,Eexpl, x0,y0,z0,r0, lamx,lamy,lamz,amp
  character(len=32) :: problem_name
  character(len=3)  :: run_id

  namelist /PROBLEM_CONTROL/  problem_name, run_id, &
                              d0,p0,lamx,lamy,lamz,amp

contains

!-----------------------------------------------------------------------------

   subroutine read_problem_par
      use errh, only : namelist_errh
      use mpisetup, only : cwd, ierr, rbuff, cbuff, ibuff, proc, &
         MPI_CHARACTER, MPI_DOUBLE_PRECISION, MPI_INTEGER, &
         buffer_dim, comm

      implicit none
      integer :: ierrh
      character(LEN=100) :: par_file, tmp_log_file

      problem_name = 'aaa'
      run_id  = 'aaa'
      d0      = 1.0
      p0      = 1.e-3
      lamx    = 0.0
      lamy    = 0.0
      lamz    = 0.0
      amp     = 0.0

      if(proc == 0) then
         par_file = trim(cwd)//'/problem.par'
         tmp_log_file = trim(cwd)//'/tmp.log'

         open(1,file=par_file)
            read(unit=1,nml=PROBLEM_CONTROL,iostat=ierrh)
            call namelist_errh(ierrh,'PROBLEM_CONTROL')
            write(*,nml=PROBLEM_CONTROL)
         close(1)
         open(3, file=tmp_log_file, position='append')
            write(3,nml=PROBLEM_CONTROL)
            write(3,*)
         close(3)
      endif


      if(proc == 0) then

         cbuff(1) =  problem_name
         cbuff(2) =  run_id

         rbuff(1) = d0
         rbuff(2) = p0
         rbuff(3) = lamx
         rbuff(4) = lamy
         rbuff(5) = lamz
         rbuff(6) = amp


         call MPI_BCAST(cbuff, 32*buffer_dim, MPI_CHARACTER,        0, comm, ierr)
         call MPI_BCAST(ibuff,    buffer_dim, MPI_INTEGER,          0, comm, ierr)
         call MPI_BCAST(rbuff,    buffer_dim, MPI_DOUBLE_PRECISION, 0, comm, ierr)

      else

         call MPI_BCAST(cbuff, 32*buffer_dim, MPI_CHARACTER,        0, comm, ierr)
         call MPI_BCAST(ibuff,    buffer_dim, MPI_INTEGER,          0, comm, ierr)
         call MPI_BCAST(rbuff,    buffer_dim, MPI_DOUBLE_PRECISION, 0, comm, ierr)

         problem_name = cbuff(1)
         run_id       = cbuff(2)(1:3)

         d0           = rbuff(1)
         p0           = rbuff(2)
         lamx         = rbuff(3)
         lamy         = rbuff(4)
         lamz         = rbuff(5)
         amp          = rbuff(6)

      endif

   end subroutine read_problem_par

!-----------------------------------------------------------------------------

   subroutine init_prob
      use mpisetup, only : proc
      use arrays, only : u,b
      use constants, only: fpiG,dpi,pi,newtong
      use grid, only : x,y,z,nx,ny,nz
      use initionized, only : gamma_ion, idni, imxi, imyi, imzi, ieni
      implicit none

      integer :: i,j,k
      real :: xi,yj,zk,kn,kx,ky,kz,Lbox,Tamp,pres,Tamp_rounded, Tamp_aux
      real :: cs0,omg,kJ

      kx = dpi/lamx
      ky = dpi/lamy
      kz = dpi/lamz

      cs0 = sqrt(gamma_ion*p0/d0)
      kn  = sqrt(kx**2+ky**2+kz**2)
      omg = sqrt(cs0**2 * kn**2 - fpiG*d0)
      kJ  = sqrt(fpiG*d0)/cs0
      Lbox = 0.5*sqrt(pi*gamma_ion*p0/newtong / d0**2)
      Tamp = (d0*amp**2*omg**2*Lbox**2)/(8.0*kn**2)

      write(*,*) 'Unperturbed adiabatic sound speed = ', cs0
      write(*,*) 'Gravitational constant * 4pi      = ', fpiG
      write(*,*) 'Lbox                              = ', Lbox
      write(*,*) 'gamma                             = ', gamma_ion
      write(*,*) 'Perturbation wavenumber           = ', kn
      write(*,*) 'Jeans wavenumber                  = ', kJ
      write(*,*) 'characteristic frequency          = ', omg
      write(*,*) '------------------------------------'
      write(*,*) 'T(t) = ',Tamp,'*[1 - cos(',2.0*omg,'t)]'
      write(*,*) '------------------------------------'
      write(*,*) 'Divide T(t) for .tsl by L to get proper amplitude !'
! Uniform equilibrium state

      do k = 1,nz
         zk = z(k)
         do j = 1,ny
            yj = y(j)
            do i = 1,nx
               xi = x(i)
               u(idni,i,j,k)   = d0 * (1.d0 + amp* dcos(kx*xi+ky*yj+kz*zk))
               pres            = p0 * (1.d0 + gamma_ion*amp*dcos(kx*xi+ky*yj+kz*zk))

               u(imxi:imzi,i,j,k)   = 0.0
#ifndef ISO
               u(ieni,i,j,k)   = pres/(gamma_ion-1.0) + &
                            0.5*sum(u(imxi:imzi,i,j,k)**2,1)/u(idni,i,j,k)

               b(:,i,j,k)      = 0.0
               u(ieni,i,j,k)   = u(ieni,i,j,k) + 0.5*sum(b(:,i,j,k)**2,1)
#endif /* ISO */
        enddo
      enddo
    enddo

    Tamp_aux = 10**int(log(Tamp)/log(10.))
    Tamp_rounded = (int(Tamp/Tamp_aux)+1)*Tamp_aux
    if(proc == 0) then
      write(*,*) '======================================'
      write(*,*) 'Run:'
      write(*,*) ' % gnuplot verify.gpl; display jeans-*.png'
      write(*,*) 'to verify results'
      write(*,*) '======================================'
      open(137,file="verify.gpl",status="unknown")
         write(137,*) "set sample 1000"
         write(137,*) "set term png #font luximr"
#ifdef SELF_GRAV
         write(137,*) "set output 'jeans-fft.png'"
         write(137,*) 'se tit "Jeans oscillations (FFT)"'
#else
         write(137,*) "set output 'jeans-mg.png'"
         write(137,*) 'se tit "Jeans oscillations (multigrid)"'
#endif
         write(137,*) "a = ", Tamp
         write(137,*) "b = ", 2.0*omg
         write(137,*) "L = ", Lbox
         write(137,*) "T = 2*pi/b"
         write(137,*) "y(x) = a *(1-cos(b*x))"
         write(137,'(2(a,/))') 'se xla "time [periods]"', 'se yla "E_int"'
         write(137,*) 'se ytics ',Tamp_rounded
         write(137,'(3(a,/))') 'se mytics 2', 'se xtics 1', 'se mxtics 2'
         write(137,*) 'plot [0:int(1./T)][0:',2*Tamp_rounded,'] "jeans_ts1_000.tsl" u ($2/T):($11/L) w p t "calculated", "" u ($2/T):($11/L) smoo cspl t "" w l 1, y(x*T) t "analytical"'
      close(137)
    endif
    return
  end subroutine init_prob

end module initproblem

