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

  real              :: d0, p0, amp, kx, ky, kz
  integer           :: ix, iy, iz, mode
  character(len=32) :: problem_name
  character(len=3)  :: run_id

  namelist /PROBLEM_CONTROL/  problem_name, run_id, d0, p0, ix, iy, iz, amp, mode

contains

!-----------------------------------------------------------------------------

   subroutine read_problem_par
      use grid,     only : xmin, xmax, ymin, ymax, zmin, zmax
      use errh,     only : namelist_errh, die
      use mpisetup, only : cwd, ierr, rbuff, cbuff, ibuff, proc, buffer_dim, comm, &
           &               MPI_CHARACTER, MPI_DOUBLE_PRECISION, MPI_INTEGER
      use constants, only: pi

      implicit none

      integer :: ierrh
      character(LEN=100) :: par_file, tmp_log_file

      ! namelist default parameter values
      problem_name = 'Jeans oscillations'  !< The default problem name
      run_id       = '_'                   !< Auxiliary run identifier
      d0           = 1.0                   !< Average density of the medium (density bias required for correct EOS evaluation)
      p0           = 1.e-3                 !< Average pressure of the medium (for calculating sound speed ot temperature)
      ix           = 2                     !< Number of perturbation waves in the x direction
      iy           = 0                     !< Number of perturbation waves in the y direction
      iz           = 0                     !< Number of perturbation waves in the z direction
      amp          = 0.0                   !< Perturbation relative amplitude
      mode         = 0                     !< Variant of the test. 0: cos(kx *x + ky*y + kz*z), 1: cos(kx *x) * cos(ky*y) * cos(kz*z)

      if(proc == 0) then
         par_file = trim(cwd)//'/problem.par'
         tmp_log_file = trim(cwd)//'/tmp.log'

         open(1,file=par_file)
            read(unit=1,nml=PROBLEM_CONTROL,iostat=ierrh)
            call namelist_errh(ierrh,'PROBLEM_CONTROL')
         close(1)
         open(3, file=tmp_log_file, position='append')
            write(3,nml=PROBLEM_CONTROL)
            write(3,*)
         close(3)
      endif

      if (proc == 0) then

         cbuff(1) =  problem_name
         cbuff(2) =  run_id

         rbuff(1) = d0
         rbuff(2) = p0
         rbuff(3) = amp

         ibuff(1) = ix
         ibuff(2) = iy
         ibuff(3) = iz
         ibuff(4) = mode

      end if

      call MPI_BCAST(cbuff, 32*buffer_dim, MPI_CHARACTER,        0, comm, ierr)
      call MPI_BCAST(ibuff,    buffer_dim, MPI_INTEGER,          0, comm, ierr)
      call MPI_BCAST(rbuff,    buffer_dim, MPI_DOUBLE_PRECISION, 0, comm, ierr)

      if (proc /= 0) then

         problem_name = cbuff(1)
         run_id       = cbuff(2)(1:3)

         d0           = rbuff(1)
         p0           = rbuff(2)
         amp          = rbuff(3)

         ix           = ibuff(1)
         iy           = ibuff(2)
         iz           = ibuff(3)
         mode         = ibuff(4)

      endif

#if !(defined(MULTIGRID) || defined(SELF_GRAV))
      call die("You must define either MULTIGRID or SELF_GRAV for this problem")
#endif

      if (mode < 0 .or. mode > 1)     call die("[initproblem:read_problem_par] Invalid mode.")
      if (d0 < 0. .or. abs(amp) > 1.) call die("[initproblem:read_problem_par] Negative average density or amplitude too high.")
      if (p0 < 0.)                    call die("[initproblem:read_problem_par] Negative average pressure.")
      if (ix<0 .or. iy<0 .or. iz<0 .or. amp<0) &
           write(*,'(a,g12.4,a,3i6,a)')        "[initproblem:read_problem_par] Warning: suspicious values for some parameters were detected: amp=", &
           amp," (ix iy iz) = (", ix, iy, iz, ")."

      kx = 2. * pi * ix / (xmax - xmin)
      ky = 2. * pi * iy / (ymax - ymin)
      kz = 2. * pi * iz / (zmax - zmin)
      if (mode == 1) then
         kx = kx / 2.
         ky = ky / 2.
         kz = kz / 2.
      end if

   end subroutine read_problem_par

!-----------------------------------------------------------------------------

   subroutine init_prob

      use mpisetup,     only: proc
      use arrays,       only: u, b
      use constants,    only: fpiG, pi, newtong
      use grid,         only: x, y, z, nx, ny, nz, xmin, ymin, zmin, dx, dy, dz
      use initionized,  only: gamma_ion, idni, imxi, imzi, ieni

      implicit none

      integer :: i, j, k
      real    :: xi, yj, zk, kn, Lbox, Tamp, pres, Tamp_rounded, Tamp_aux
      real    :: cs0, omg, omg2, kJ

      cs0  = sqrt(gamma_ion * p0 / d0)
      kn   = sqrt(kx**2 + ky**2 + kz**2)
      omg2 = cs0**2 * kn**2 - fpiG * d0
      omg  = sqrt(abs(omg2))
      kJ   = sqrt(fpiG * d0) / cs0
      Lbox = 0.5 * sqrt(pi * gamma_ion * p0 / newtong / d0**2)
      Tamp = (d0 * amp**2 * omg2 * Lbox**2)/(8.0 * kn**2)
      if (mode == 1) Tamp = Tamp / 4.

      if (proc == 0) then
         write(*,*) 'Unperturbed adiabatic sound speed = ', cs0
         write(*,*) 'Gravitational constant * 4pi      = ', fpiG
         write(*,*) 'Lbox                              = ', Lbox
         write(*,*) 'gamma                             = ', gamma_ion
         write(*,*) 'Perturbation wavenumber           = ', kn
         write(*,*) 'Jeans wavenumber                  = ', kJ
         if (omg2>0) then
            write(*,*) 'characteristic frequency          = ', omg
            write(*,*)
            write(*,*) 'T(t) = ',Tamp,'*[1 - cos(',2.0*omg,'t)]'
         else
            write(*,*) 'characteristic timescale          = ', omg
            !write(*,*) 'Something like T(t) = ',Tamp,'* exp(',omg,'t)]' !BEWARE the formula is not completed
         end if
         write(*,'(/,a)') 'Divide T(t) for .tsl by L to get proper amplitude !'
      end if
! Uniform equilibrium state

      do k = 1,nz
         zk = z(k)-zmin
         do j = 1,ny
            yj = y(j)-ymin
            do i = 1,nx
               xi = x(i)-xmin
               select case (mode)
               case (0)
                  u(idni,i,j,k)   = d0 * (1.d0 +             amp * dsin(kx*xi + ky*yj + kz*zk))
                  pres            = p0 * (1.d0 + gamma_ion * amp * dsin(kx*xi + ky*yj + kz*zk))
               case (1)
                  u(idni,i,j,k)   = d0 * (1.d0 +             amp * dsin(kx*xi) * dsin(ky*yj) * dsin(kz*zk))
                  pres            = p0 * (1.d0 + gamma_ion * amp * dsin(kx*xi) * dsin(ky*yj) * dsin(kz*zk))
               case default ! should not happen
                  u(idni,i,j,k)   = d0
                  pres            = p0
               end select

               u(imxi:imzi,i,j,k) = 0.0
#ifndef ISO
               u(ieni,i,j,k)      = pres/(gamma_ion-1.0) + 0.5*sum(u(imxi:imzi,i,j,k)**2,1) / u(idni,i,j,k)

               b(:,i,j,k)         = 0.0
               u(ieni,i,j,k)      = u(ieni,i,j,k) + 0.5*sum(b(:,i,j,k)**2,1)
#endif /* ISO */
        enddo
      enddo
    enddo

    Tamp_aux = 10**int(log(Tamp)/log(10.))
    Tamp_rounded = (int(Tamp/Tamp_aux)+1)*Tamp_aux
    if (proc == 0) then
      write(*,'(/,a)') 'Run:'
#ifdef MULTIGRID
      write(*,'(a)')   ' % gnuplot verify.gpl; display jeans-mg.png'
#else
      write(*,'(a)')   ' % gnuplot verify.gpl; display jeans-fft.png'
#endif
      write(*,'(a,/)') 'to verify results'

      open(137,file="verify.gpl",status="unknown")
         write(137,'(a)') "set sample 1000"
         write(137,'(a)') "set term png #font luximr"
#ifdef MULTIGRID
         write(137,'(a)') "set output 'jeans-mg.png'"
         write(137,'(a)') 'set title "Jeans oscillations (multigrid)"'
#else
         write(137,'(a)') "set output 'jeans-fft.png'"
         write(137,'(a)') 'set title "Jeans oscillations (FFT)"'
#endif
         write(137,'(a,g13.5)') "a = ", Tamp
         write(137,'(a,g13.5)') "b = ", 2.0*omg
         write(137,'(a,g13.5)') "L = ", Lbox
         write(137,'(a,g13.5)') "T = 2*pi/b"
         write(137,'(a,g13.5)') "y(x) = a *(1-cos(b*x))"
         write(137,'(a,/,a)') 'set xlabel "time [periods]"', 'set ylabel "E_int"'
         write(137,'(a)') "set key left Left reverse bottom"
         write(137,'(3(a,/),a,g10.2)') 'set xtics 1', 'set mxtics 2', 'set mytics 2', 'set ytics ',Tamp_rounded
         write(137,'(2(a,g10.2),a)') 'plot [0:int(1./T)][',Tamp_rounded/(-2.),':',2*Tamp_rounded,'] "jeans_ts1_000.tsl" u ($2/T):($11/L) w p t "calculated", "" u ($2/T):($11/L) smoo cspl t "" w l 1, y(x*T) t "analytical", "" u ($2/T):(10*(y($2)-$11/L)) t "10 * difference" w lp, 0 t "" w l 0'
      close(137)
    endif
    return
  end subroutine init_prob

end module initproblem
