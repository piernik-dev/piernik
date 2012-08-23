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
!    Initial implementation of PIERNIK code was based on TVD split MHD code by
!    Ue-Li Pen
!        see: Pen, Arras & Wong (2003) for algorithm and
!             http://www.cita.utoronto.ca/~pen/MHD
!             for original source code "mhd.f90"
!
!    For full list of developers see $PIERNIK_HOME/license/pdt.txt
!
#include "piernik.h"
#include "macros.h"

module initproblem

   implicit none

   private
   public :: read_problem_par, init_prob, problem_pointers
   public :: d0, mode

   ! Private variables
   real :: kx, ky, kz

   ! Namelist variables
   real            :: d0, p0, amp
   integer(kind=4) :: ix, iy, iz, mode

   namelist /PROBLEM_CONTROL/  d0, p0, ix, iy, iz, amp, mode

contains

!-----------------------------------------------------------------------------

   subroutine problem_pointers

      implicit none

   end subroutine problem_pointers

!-----------------------------------------------------------------------------

   subroutine read_problem_par

      use constants,   only: xdim, ydim, zdim, pi
      use dataio_pub,  only: ierrh, par_file, namelist_errh, compare_namelist, cmdl_nml, lun    ! QA_WARN required for diff_nml
      use dataio_pub,  only: tend, msg, die, warn, printinfo
      use domain,      only: dom
      use fluidindex,  only: flind
      use fluidtypes,  only: component_fluid
      use mpisetup,    only: rbuff, ibuff, master, slave, piernik_MPI_Bcast
      use problem_pub, only: jeans_d0, jeans_mode
      use units,       only: fpiG, newtong

      implicit none

      class(component_fluid), pointer :: fl
      real                            :: kn, cs0, omg2, kJ, Tamp_rounded, Tamp_aux, Tamp, omg
      integer, parameter              :: g_lun = 137

      ! namelist default parameter values
      d0          = 1.0                   !< Average density of the medium (density bias required for correct EOS evaluation)
      p0          = 1.e-3                 !< Average pressure of the medium (for calculating sound speed ot temperature)
      ix          = 2                     !< Number of perturbation waves in the x direction
      iy          = 0                     !< Number of perturbation waves in the y direction
      iz          = 0                     !< Number of perturbation waves in the z direction
      amp         = 0.0                   !< Perturbation relative amplitude
      mode        = 0                     !< Variant of the test. 0: cos(kx *x + ky*y + kz*z), 1: cos(kx *x) * cos(ky*y) * cos(kz*z)

      if (master) then

         diff_nml(PROBLEM_CONTROL)

         rbuff(1) = d0
         rbuff(2) = p0
         rbuff(3) = amp

         ibuff(1) = ix
         ibuff(2) = iy
         ibuff(3) = iz
         ibuff(4) = mode

      endif

      call piernik_MPI_Bcast(ibuff)
      call piernik_MPI_Bcast(rbuff)

      if (slave) then

         d0       = rbuff(1)
         p0       = rbuff(2)
         amp      = rbuff(3)

         ix       = ibuff(1)
         iy       = ibuff(2)
         iz       = ibuff(3)
         mode     = ibuff(4)

      endif

#if !(defined(MULTIGRID) || defined(POISSON_FFT))
      call die("You must define either MULTIGRID or POISSON_FFT for this problem")
#endif /* !(MULTIGRID || POISSON_FFT) */

      if (mode < 0 .or. mode > 1)     call die("[initproblem:read_problem_par] Invalid mode.")
      if (d0 < 0. .or. abs(amp) > 1.) call die("[initproblem:read_problem_par] Negative average density or amplitude too high.")
      if (p0 < 0.)                    call die("[initproblem:read_problem_par] Negative average pressure.")
      if (ix<0 .or. iy<0 .or. iz<0 .or. amp<0) then
         write(msg, '(a,g12.4,a,3i6,a)')     "[initproblem:read_problem_par] Suspicious values for some parameters were detected: amp=", &
              &                               amp," (ix iy iz) = (", ix, iy, iz, ")."
         call warn(msg)
      endif

      ! suppress waves in nonexistent directions
      if (.not. dom%has_dir(xdim)) ix = 0
      if (.not. dom%has_dir(ydim)) iy = 0
      if (.not. dom%has_dir(zdim)) iz = 0

      kx = 2. * pi * ix / dom%L_(xdim)
      ky = 2. * pi * iy / dom%L_(ydim)
      kz = 2. * pi * iz / dom%L_(zdim)
      if (mode == 1) then
         kx = kx / 2.
         ky = ky / 2.
         kz = kz / 2.
      endif

      ! export variables required for a test in multigrid_gravity
      jeans_d0 = d0
      jeans_mode = mode

      fl => flind%ion

      cs0  = sqrt(fl%gam * p0 / d0)
      kn   = sqrt(kx**2 + ky**2 + kz**2)
      omg2 = cs0**2 * kn**2 - fpiG * d0
      omg  = sqrt(abs(omg2))
      kJ   = sqrt(fpiG * d0) / cs0
      if (kn > 0) then
         Tamp = (d0 * amp**2 * omg2 * dom%Vol)/(4.0 * kn**2)
      else
         Tamp = 0.
         if (master) call warn("[initproblem:init_prob] No waves (kn == 0)")
      endif
      if (mode == 1) Tamp = Tamp / 4.
      if (Tamp > 0) then
         Tamp_aux = 10**int(log(Tamp)/log(10.))
         Tamp_rounded = (int(1.05*Tamp/Tamp_aux)+1)*Tamp_aux
      else
         Tamp_rounded = 0.
      endif

      if (master) then
         write(msg, *) 'Unperturbed adiabatic sound speed = ', cs0
         call printinfo(msg, .true.)
         write(msg, *) 'Gravitational constant * 4pi      = ', fpiG
         call printinfo(msg, .true.)
         write(msg, *) 'L_critical                        = ', sqrt(pi * fl%gam * p0 / newtong / d0**2)
         call printinfo(msg, .true.)
         write(msg, *) 'gamma                             = ', fl%gam
         call printinfo(msg, .true.)
         write(msg, *) 'Perturbation wavenumber           = ', kn
         call printinfo(msg, .true.)
         write(msg, *) 'Jeans wavenumber                  = ', kJ
         call printinfo(msg, .true.)
         if (omg2>0) then
            write(msg, *) 'characteristic frequency          = ', omg
            call printinfo(msg, .true.)
            call printinfo("", .true.)
            write(msg, *) 'T(t) = ',Tamp,'*[1 - cos(',omg,'t)]'
            call printinfo(msg, .true.)
         else
            write(msg, *) 'characteristic timescale          = ', omg
            call printinfo(msg, .true.)
            !write(msg, *) 'Something like T(t) = ',Tamp,'* exp(',omg,'t)]' !BEWARE the formula is not completed
            !call printinfo(msg, .true.)
         endif
         call printinfo('Divide T(t) for .tsl by L to get proper amplitude !', .true.)
         call printinfo('', .true.)
         call printinfo('To verify results, run:', .true.)
#ifdef MULTIGRID
         call printinfo(' % gnuplot verify.gpl; display jeans-mg.png', .true.)
#else /* !MULTIGRID */
         call printinfo(' % gnuplot verify.gpl; display jeans-fft.png', .true.)
#endif /* !MULTIGRID */
         call printinfo('', .true.)

         open(g_lun,file="verify.gpl",status="unknown")
         write(g_lun,'(a)') "set sample 1000"
         write(g_lun,'(a)') "set term png #font luximr"
#ifdef MULTIGRID
         write(g_lun,'(a)') "set output 'jeans-mg.png'"
         write(g_lun,'(a)') 'set title "Jeans oscillations (multigrid)"'
#else /* !MULTIGRID */
         write(g_lun,'(a)') "set output 'jeans-fft.png'"
         write(g_lun,'(a)') 'set title "Jeans oscillations (FFT)"'
#endif /* !MULTIGRID */
         write(g_lun,'(3(a,/),a)') 'set ylabel "E_int"', 'set xtics 1', 'set mxtics 2', 'set mytics 2'
         if (Tamp_rounded /= 0 .and. Tamp >0) then
            write(g_lun,'(a,g11.3)')'set ytics ',Tamp_rounded/2.
            write(g_lun,'(2(a,g11.3),a)')'set yrange [ ',Tamp_rounded/(-4.),':',Tamp_rounded,']'
         else
            write(g_lun,'(a)')'set yrange [ * : * ]'
         endif
         if (Tamp >0) then
            write(g_lun,'(a)') "set key left Left reverse bottom"
            write(g_lun,'(a,g13.5)') "a = ", Tamp
            write(g_lun,'(a,g13.5)') "b = ", omg
            write(g_lun,'(a,g13.5)') "T = 2*pi/b"
            write(g_lun,'(a,g13.5)') "y(x) = a * sin(b*x)**2"
            write(g_lun,'(a)') 'set xlabel "time [periods]"'
            if (tend > pi/omg) then
               write(g_lun,'(a,g11.3,a)')'set xrange [ 0 : int(',tend,'/T)]'
            else
               write(g_lun,'(a)')'set xrange [ * : * ]'
            endif
            write(g_lun,'(a)') 'plot "jeans_ts1_000.tsl" u ($2/T):($11) w p t "calculated", "" u ($2/T):($11) smoo cspl t "" w l lt 1, y(x*T) t "analytical", "" u ($2/T):(10*(y($2)-$11)) t "10 * difference" w lp, 0 t "" w l lt 0'
         else

            write(g_lun,'(a,g13.5)') "a = ", amp**2 * omg**2 * 800000. !BEWARE: stronger dependence on omg, magic number 800000
            write(g_lun,'(a,g13.5)') "b = ", 2.0*omg
            write(g_lun,'(a,g13.5)') "T = 2*pi/b"
            write(g_lun,'(a,g13.5)') "y(x) = a * exp(b*x)"
            write(g_lun,'(3(a,/),a)') 'set key left Left reverse top', 'set log y', 'set xlabel "time"', 'set xrange [ * : * ]'
            write(g_lun,'(a)') 'plot "jeans_ts1_000.tsl" u ($2):($11) w p t "calculated", y(x) t "exp(2 om T)"'
         endif
         close(g_lun)
      endif

   end subroutine read_problem_par

!-----------------------------------------------------------------------------

   subroutine init_prob

      use cg_list,     only: cg_list_element
      use cg_list_bnd, only: leaves
      use constants,   only: xdim, ydim, zdim, LO
      use domain,      only: dom
      use fluidindex,  only: flind
      use fluidtypes,  only: component_fluid
      use func,        only: ekin, emag
      use grid_cont,   only: grid_container

      implicit none

      class(component_fluid), pointer :: fl
      integer                         :: i, j, k
      real                            :: xi, yj, zk, pres
      type(cg_list_element),  pointer :: cgl
      type(grid_container),   pointer :: cg

! Uniform equilibrium state
      fl => flind%ion
      cgl => leaves%first
      do while (associated(cgl))
         cg => cgl%cg
         do k = cg%ks, cg%ke
            zk = cg%z(k)-dom%edge(zdim, LO)
            do j = cg%js, cg%je
               yj = cg%y(j)-dom%edge(ydim, LO)
               do i = cg%is, cg%ie
                  xi = cg%x(i)-dom%edge(xdim, LO)
                  select case (mode)
                     case (0)
                        cg%u(fl%idn,i,j,k)  = d0 * (1. +          amp * sin(kx*xi + ky*yj + kz*zk))
                        pres                = p0 * (1. + fl%gam * amp * sin(kx*xi + ky*yj + kz*zk))
                     case (1)
                        cg%u(fl%idn,i,j,k)  = d0 * (1. +          amp * sin(kx*xi) * sin(ky*yj) * sin(kz*zk))
                        pres                = p0 * (1. + fl%gam * amp * sin(kx*xi) * sin(ky*yj) * sin(kz*zk))
                     case default ! should not happen
                        cg%u(fl%idn,i,j,k)  = d0
                        pres                = p0
                  end select

                  cg%u(fl%imx:fl%imz,i,j,k) = 0.0
#ifndef ISO
                  cg%u(fl%ien,i,j,k)        = pres/fl%gam_1 + ekin(cg%u(fl%imx,i,j,k), cg%u(fl%imy,i,j,k), cg%u(fl%imz,i,j,k), cg%u(fl%idn,i,j,k))

#ifdef MAGNETIC
                  cg%b(:,i,j,k)             = 0.0
                  cg%u(fl%ien,i,j,k)        = cg%u(fl%ien,i,j,k) + emag(cg%b(xdim,i,j,k)), cg%b(ydim,i,j,k), cg%b(zdim,i,j,k)
#endif /* MAGNETIC */
#endif /* !ISO */
               enddo
            enddo
         enddo
         cgl => cgl%nxt
      enddo

   end subroutine init_prob

end module initproblem
