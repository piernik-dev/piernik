! $Id$

#include "piernik.h"
#include "macros.h"

module initproblem

   implicit none

   private
   public :: read_problem_par, init_prob, problem_pointers

   integer           :: n_sn
   real              :: d0, Mrms, t_sn, c_si

   namelist /PROBLEM_CONTROL/  d0, c_si, Mrms

contains

!-----------------------------------------------------------------------------

   subroutine problem_pointers

      implicit none

   end subroutine problem_pointers

!-----------------------------------------------------------------------------

   subroutine read_problem_par

      use dataio_pub,    only: ierrh, par_file, namelist_errh, compare_namelist, cmdl_nml, lun      ! QA_WARN required for diff_nml
      use mpi,           only: MPI_DOUBLE_PRECISION
      use mpisetup,      only: rbuff, buffer_dim, comm, mpi_err, master, slave, FIRST

      implicit none

      t_sn = 0.0

      d0   = 1.0
      c_si = 0.1
      Mrms = 5.0

      if (master) then

         diff_nml(PROBLEM_CONTROL)

         rbuff(1) = d0
         rbuff(2) = Mrms
         rbuff(3) = c_si

      endif

      call MPI_Bcast(rbuff, buffer_dim, MPI_DOUBLE_PRECISION, FIRST, comm, mpi_err)

      if (slave) then

         d0       = rbuff(1)
         Mrms     = rbuff(2)
         c_si     = rbuff(3)

      endif

   end subroutine read_problem_par

!-----------------------------------------------------------------------------

   subroutine init_prob

      use constants,  only: xdim, ydim, zdim, ndims
      use dataio_pub, only: msg, printinfo
      use domain,     only: dom
      use fluidindex, only: flind
      use fluidtypes, only: component_fluid
      use func,       only: resample_gauss, ekin, emag
      use gc_list,    only: cg_list_element
      use grid,       only: leaves
      use grid_cont,  only: grid_container
      use mpisetup,   only: proc

      implicit none

      class(component_fluid), pointer       :: fl
      integer                               :: i, j, k, m, n, l
      integer, parameter                    :: kp = 8
      real, dimension(6)                    :: mn
      real, dimension(ndims)                :: deltav
      real, dimension(:,:,:,:), allocatable :: dv
      real                                  :: rms, cma
      type(cg_list_element), pointer        :: cgl
      type(grid_container),  pointer        :: cg
!      real :: somx,somy,somz

! Uniform equilibrium state

      call random_seed()

      fl  => flind%neu
      cgl => leaves%first
      do while (associated(cgl))
         cg => cgl%cg

         allocate(dv(ndims, cg%n_(xdim), cg%n_(ydim), cg%n_(zdim)))

         do k = 1, cg%n_(zdim)
            do j = 1, cg%n_(ydim)
               do i = 1, cg%n_(xdim)
                  deltav(:) = 0.0
                  do m=-kp,kp
                     do n = -kp,kp
                        call random_number(mn)
!                     somx = dpi*(float(n)*y(j) + float(m)*z(k)) / dom%L_(xdim)
!                     somy = dpi*(float(n)*x(i) + float(m)*z(k)) / dom%L_(ydim)
!                     somz = dpi*(float(n)*x(i) + float(m)*y(j)) / dom%L_(zdim)
!                     deltav(xdim) = deltav(xdim) + mn(1)*sin(somx) + mn(2)*cos(somx)
!                     deltav(ydim) = deltav(ydim) + mn(3)*sin(somy) + mn(4)*cos(somy)
!                     deltav(zdim) = deltav(zdim) + mn(5)*sin(somz) + mn(6)*cos(somz)
                        deltav(xdim) = deltav(xdim) + mn(1)*sin(float(m)*cg%z(k)) + mn(2)*cos(float(n)*cg%y(j))
                        deltav(ydim) = deltav(ydim) + mn(3)*sin(float(m)*cg%x(i)) + mn(4)*cos(float(n)*cg%z(k))
                        deltav(zdim) = deltav(zdim) + mn(5)*sin(float(m)*cg%y(j)) + mn(6)*cos(float(n)*cg%x(i))
                     enddo
                  enddo
                  dv(:,i,j,k) = deltav(:)
               enddo
            enddo
         enddo
         rms = sqrt( sum(dv**2) / real(product(cg%n_)) )

         cma = 1.0
         if ( rms/c_si < 0.1) then
            cma = rms/c_si
         else
            l=1
            do while (cma >= 0.1 .and. l <=10)
               cma = rms / c_si * (0.1)**l
               l=l+1
            enddo
         endif

         write(msg,'(2(a,g12.5),a,i4)')   "[initproblem:init_prob] cma = ", cma, " rms = ", rms, " on ", proc
         call printinfo(msg, .true.)
         write(msg,'(2(a,g12.5),a,i4)')   "[initproblem:init_prob] c_si = ", c_si, " l = ", l, " on ", proc
         call printinfo(msg, .true.)
         do k = 1, cg%n_(zdim)
            do j = 1, cg%n_(ydim)
               do i = 1, cg%n_(xdim)
                  cg%u(fl%idn,i,j,k) = d0
                  cg%u(fl%imx,i,j,k) = cg%u(fl%idn,i,j,k) * dv(xdim,i,j,k) * cma
                  cg%u(fl%imy,i,j,k) = cg%u(fl%idn,i,j,k) * dv(ydim,i,j,k) * cma
                  cg%u(fl%imz,i,j,k) = cg%u(fl%idn,i,j,k) * dv(zdim,i,j,k) * cma
                  cg%u(fl%ien,i,j,k) = c_si**2*d0/(fl%gam*fl%gam_1)
                  cg%u(fl%ien,i,j,k) = cg%u(fl%ien,i,j,k) + ekin(cg%u(fl%imx,i,j,k), cg%u(fl%imy,i,j,k), cg%u(fl%imz,i,j,k),cg%u(fl%idn,i,j,k))
                  cg%b(:,i,j,k)      = 0.0
                  cg%u(fl%ien,i,j,k) = cg%u(fl%ien,i,j,k) + emag(cg%b(xdim,i,j,k), cg%b(ydim,i,j,k), cg%b(zdim,i,j,k))
                  cg%u(flind%trc%beg:flind%trc%end, i, j, k)   = &
                        resample_gauss( cg%x(i), cg%y(j), cg%z(k), cg%dl(xdim), cg%dl(ydim), cg%dl(zdim), 0.05*dom%L_(xdim), 0.05*dom%L_(ydim), 0.05*dom%L_(zdim), 10)
               enddo
            enddo
         enddo

! Explosions

         deallocate(dv)

         cgl => cgl%nxt
      enddo

   end subroutine init_prob

end module initproblem
