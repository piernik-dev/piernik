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
#define RNG is:ie, js:je, ks:ke
module initproblem

! Initial condition for Sedov-Taylor explosion
! Written by: M. Hanasz, March 2006

   implicit none

   private
   public  :: read_problem_par, problem_initial_conditions, problem_pointers

   integer(kind=4) :: n_sn
   real            :: d0, p0, bx0, by0, bz0, Eexpl, x0, y0, z0, r0, dt_sn, r, t_sn
   

   namelist /PROBLEM_CONTROL/ d0, p0, bx0, by0, bz0, Eexpl, x0, y0, z0, r0, n_sn, dt_sn

contains
!-----------------------------------------------------------------------------
   subroutine problem_pointers

      implicit none

   end subroutine problem_pointers
!-----------------------------------------------------------------------------
   subroutine read_problem_par

      implicit none
      !call print_essential_units
      !call print_essential_units

   end subroutine read_problem_par
!-----------------------------------------------------------------------------
   subroutine problem_initial_conditions

      use cg_leaves,    only: leaves
      use cg_list,      only: cg_list_element
      use constants,    only: xdim, ydim, zdim, LO, HI
      use dataio_pub,   only: printinfo
      use fluidindex,   only: flind
      use particle_pub, only: pset
      !use gravity,      only: ptmass
      !use units,        only: newtong
      !use particles_io_hdf5

      implicit none

      integer                          :: i, j, k, p, n_particles
      !real                              :: dtheta
      real(kind=4)                     :: e
      !real,dimension(3)                :: pos_init, vel_init
      !real,parameter                   :: pi2=6.283185307
      logical,save                     :: first_run = .true.
      type(cg_list_element),  pointer  :: cgl

      
      do p = lbound(flind%all_fluids, dim=1), ubound(flind%all_fluids, dim=1)
         cgl => leaves%first
         do while (associated(cgl))
            associate(cg => cgl%cg)
               do k = cg%lhn(zdim,LO), cg%lhn(zdim,HI)
                  do j = cg%lhn(ydim,LO), cg%lhn(ydim,HI)
                     do i = cg%lhn(xdim,LO), cg%lhn(xdim,HI)
                        associate( fl => flind%all_fluids(p)%fl )
                           cg%u(fl%idn,i,j,k) = 1.0!e-6
                           cg%u(fl%imx,i,j,k) = 0.0
                           cg%u(fl%imy,i,j,k) = 0.0
                           cg%u(fl%imz,i,j,k) = 0.0
                        end associate
                     enddo
                  enddo
               enddo
            end associate
            cgl => cgl%nxt
         enddo
         
      enddo


      e = 0.0

      n_particles = 1
      !write(*,*) "int(3/3)=",int(3.0/3.0)
      
      call orbits(n_particles, e, first_run)
      !call relax_time(n_particles, first_run)
      !call read_buildgal
      
      


      contains

      function positions(dtheta, pos_init)
         implicit none
            real, dimension(3) :: positions, pos_init
            real :: dtheta
               positions = rotate(dtheta, pos_init)
      end function positions


      function velocities(pos_init, e)
         implicit none
            real, dimension(3) :: pos_init, velocities
            real :: a, r
            real(kind=4) :: e
            real, parameter :: mu = 1.0,G=1.0, M=1.0, zero = 0.0
            if((e<0.0) .or. (e>=1.0)) then
               write(*,*) "Bledna wartosc mimosrodu! Zatrzymano!"
               stop
            else
               r = sqrt(pos_init(1)**2+pos_init(2)**2+pos_init(3)**2)
               if (e==zero) then
                  velocities(2) = sqrt(G*M/r)
                  write(*,*) "Orbita kolowa"
               else

                  a = r/(1.0 + e)
                  velocities(2) = sqrt(mu*(2.0/r - 1.0/a))
                  write(*,'(A11,F4.2,A3,F5.3,A3,F5.3)') "#Elipsa: e=", e, " a=",a, " b=", a*sqrt(1.0 - e**2)
                  !lenght = pi*( (x0/(1.0+e))*(1.5*(1.0+sqrt(1.0-e**2)) - (1.0-e**2)**0.25))
               endif
            endif!
            velocities(1) = 0.0
            !velocities(2) = vy
            velocities(3) = 0.0
      end function velocities
      
      function rotate (theta, vector)
         implicit none
            real, dimension(3) :: vector, rotate
            real :: theta
               !x-axis
               !rotate(1) = vector(1)
               !rotate(2) = vector(2)*cos(theta) - vector(3)*sin(theta)
               !rotate(3) = vector(2)*sin(theta) + vector(3)*cos(theta)
               !y-axis
               !rotate(1) = vector(1)*cos(theta) - vector(3)*sin(theta)
               !rotate(2) = vector(2)
               !rotate(3) = vector(1)*sin(theta) + vector(3)*cos(theta)
               !z-axis
               rotate(1) = vector(1)*cos(theta) - vector(2)*sin(theta)
               rotate(2) = vector(1)*sin(theta) + vector(2)*cos(theta)
               rotate(3) = vector(3)
      end function rotate
   
   
      subroutine orbits(n_particles, e, first_run)
         use particle_pub, only: pset
         implicit none
         integer,intent(in)            :: n_particles
         real(kind=4),intent(in)       :: e
         real,dimension(3)             :: pos_init, vel_init
         real                           :: dtheta,zero = 0.0
         real,parameter                :: pi2=6.283185307
         logical,intent(inout)         :: first_run
         
         write(*,*) "Number of particles: ", n_particles
         
         dtheta = pi2/n_particles
         !write(*,*) "1/0: ", 1.0/zero

         pos_init(1) = 2.0
         pos_init(2) = 0.0
         pos_init(3) = 0.0
         
         vel_init = velocities(pos_init, e)
         
         if(first_run) then
            !call pset%add(1.1, [ 0.9700436, -0.24308753, 0.0], [ 0.466203685, 0.43236573, 0.0], 0.0)
            !call pset%add(1.1, [-0.9700436, 0.24308753, 0.0], [ 0.466203685, 0.43236573, 0.0],0.0)
            !call pset%add(1.1, [ 0.0, 0.0, 0.0], [-0.932407370, -0.86473146, 0.0], 0.0 )
            !do i = 1, n_particles, 1
               !call pset%add(1.0, pos_init, vel_init,0.0 ) !orbita eliptyczna
               !call pset%add(1.0, [4.54545454545455, 0.0, 0.0],[-0.909090909090909, 0.0, 0.0], 0.0)
               
               !call pset%add(1.0, [5.0, 0.0, 0.0],[-0.909090909090909, 0.0, 0.0], 0.0) !10 na komorke
               !call pset%add(1.0, [4.54545454545455, 0.0, 0.0],[-9.09090909090909, 0.0, 0.0], 0.0) !srodki
               !call pset%add(1.0, [4.54545454545455, 0.454545454545455, 0.0],[-0.909090909090909, 0.0, 0.0], 0.0) !j=1/2 sciana
               !call pset%add(1.0, [4.54545454545455, 3.0*0.909090909090909, 0.0],[-0.909090909090909, 0.0, 0.0], 0.0) !j=1
               
               !pos_init = positions(dtheta, pos_init)
               !vel_init = rotate(dtheta, vel_init)
            !enddo

            call pset%add(10.0, [0.0,0.0,0.0],[0.0,0.0,0.0],0.0)
            call pset%add(1.0, [2.0,2.0,0.0],[0.0,0.0,0.0],0.0)
           
            first_run = .false.
            
            write(*,*) "Obliczono pozycje czastek "
         endif
      end subroutine orbits
      
      subroutine relax_time(n_particles, first_run)
         use particle_pub, only: pset
         use particles_io_hdf5!, only: write_hdf5, read_hdf5
         implicit none
         integer :: i, j
         integer,parameter :: seed = 86437
         integer,intent(in) :: n_particles
         logical,intent(inout) :: first_run
         real, dimension(n_particles, 3) :: pos_init, pos_init2!,vel_init
         real,dimension(3,2) :: domain
         real :: tmp, factor, x, y, z, r, r_dom
         real, parameter :: onesixth = 1.0/6.0
         real, dimension(n_particles,3) :: vel2hdf5,pos2hdf5
         
         !allocate(pos2hdf5(n_particles,3), vel2hdf5(n_particles,3))
         
         domain(1,1) = -5.0
         domain(2,1) = -5.0
         domain(3,1) = -5.0
         domain(1,2) = 5.0
         domain(2,2) = 5.0
         domain(3,2) = 5.0
         
         write(*,*) "Number of particles: ", n_particles
         call srand(seed)
         r_dom = onesixth*sqrt(domain(1,2)**2 + domain(2,2)**2 + domain(3,2)**2)
         

         if(first_run) then
            !open(unit=37, file="particles.dat")
            do i = 1, n_particles
               r = r_dom
               do while ((r>=r_dom))
                  do j = 1, 3
                     factor = rand(0)
                     pos_init(i, j) = sign(rand(0)*domain(j, 2),factor-0.5)
                  enddo
                  r = sqrt(pos_init(i,1)**2 + pos_init(i,2)**2 + pos_init(i,3)**2)
                  !write(*,*) i, r
               enddo
               !write(37, *) i, pos_init(i,:)
               pos2hdf5(i,:) = pos_init(i,:)
               call pset%add(1.0, pos_init(i,:), [0.0,0.0,0.0],0.0 )
            enddo
            !close(37)
            first_run = .false.
            write(*,*) "Obliczono pozycje czastek"
            !write(*,*) pos_init
            call write_hdf5(pos_init, n_particles)
            call read_hdf5(pos_init2, n_particles)
            !write(*,*) pos_init
            !write(*,*)
            !write(*,*)
            !write(*,*) pos_init2
            !write(*,*)
            !write(*,*) pos_init - pos_init2
         endif
         
         !stop

      end subroutine relax_time

   end subroutine problem_initial_conditions
   
      subroutine print_essential_units

      use dataio_pub, only: msg, printinfo
      use mpisetup,   only: master
      use units,      only: units_set, s_len_u, s_time_u, s_mass_u, miu0, km, au, lyr, mH, gmu, erg, eV, clight, Gs, mGs, cm, ppcm2, ppcm3

      implicit none

      if (master) then

         call printinfo('VINE to PIERNIK convertion of units done:')
         write(msg,'(a,a28       )') units_set, ' is used in this simulation.'                         ; call printinfo(msg)
         write(msg,'(a16,3a      )') 'It is based on: ', trim(s_len_u), trim(s_time_u), trim(s_mass_u) ; call printinfo(msg)
         write(msg,'(a30,e22.15  )') '                       miu0 = ', miu0                            ; call printinfo(msg)
         call printinfo('Some essential units for this set:')
         write(msg,'(a30,e22.15,a)') '                         km = ', km, s_len_u                     ; call printinfo(msg)
         write(msg,'(a30,e22.15,a)') '                         AU = ', au, s_len_u                     ; call printinfo(msg)
         write(msg,'(a30,e22.15,a)') '                        lyr = ', lyr, s_len_u                    ; call printinfo(msg)
         write(msg,'(a30,e22.15,a)') 'hydrogen atom mass:      mH = ', mH, s_mass_u                    ; call printinfo(msg)
         write(msg,'(a30,e22.15,a)') 'galactic mass unit:     gmu = ', gmu, s_mass_u                   ; call printinfo(msg)
         write(msg,'(a30,e22.15  )') '    hydrogen atoms / cm2    = ', mH/cm**2                        ; call printinfo(msg)
         write(msg,'(a30,e22.15  )') '    hydrogen atoms / cm3    = ', mH/cm**3                        ; call printinfo(msg)
         write(msg,'(a30,e22.15  )') 'averaged particles / cm2    = ', ppcm2                           ; call printinfo(msg)
         write(msg,'(a30,e22.15  )') 'averaged particles / cm3    = ', ppcm3                           ; call printinfo(msg)
         write(msg,'(a30,e22.15  )') '                        erg = ', erg                             ; call printinfo(msg)
         write(msg,'(a30,e22.15  )') '                         eV = ', eV                              ; call printinfo(msg)
         write(msg,'(a30,e22.15  )') '                  erg / cm3 = ', erg/cm**3                       ; call printinfo(msg)
         write(msg,'(a30,e22.15  )') '                   eV / cm3 = ', eV/cm**3                        ; call printinfo(msg)
         write(msg,'(a30,e22.15  )') 'speed of light in vacuum: c = ', clight                          ; call printinfo(msg)
         write(msg,'(a30,e22.15  )') 'magnetic induction unit: Gs = ', Gs                              ; call printinfo(msg)
         write(msg,'(a30,e22.15  )') '                    microGs = ', mGs                            ; call printinfo(msg)
         
         write(msg,*) '                    rho = ',  cm**3/mH                           ; call printinfo(msg)
      endif

   end subroutine print_essential_units


   subroutine read_buildgal
      use particle_pub, only: pset
      implicit none
      integer :: i, j, nbodies,dims=3
      integer :: galfile=1
      character(len=6) :: galname="SPIRAL"
      real,dimension(:,:), allocatable :: pos, vel
      real,dimension(:), allocatable :: mass

      open(unit=galfile,file=galname,action="read",status="old")
         read(galfile,*) nbodies
         write(*,*) "nbodies=",nbodies

         allocate(mass(nbodies),pos(nbodies,3),vel(nbodies,3))

         read(galfile,*) (mass(i),i=1,nbodies), &
         ((pos(i,j),j=1,dims),i=1,nbodies),&
         ((vel(i,j),j=1,dims),i=1,nbodies)

      close(galfile)
      
      open(unit=2,file='galtest.dat')
         do i=1,nbodies
         if (modulo(i,1000) .eq. 0.0) then
            write(*,*) i!, mass(i), pos(i,:), vel(i,:)
         endif
            call pset%add(mass(i), pos(i,:), vel(i,:),0.0)
         enddo
      close(2)
      !!pset%add(mass(i), [pos(i,:)], [vel(i,:)],0.0)

   end subroutine read_buildgal

!-----------------------------------------------------------------------------
end module initproblem
