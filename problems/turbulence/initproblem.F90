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
      use mpisetup,      only: rbuff, buffer_dim, comm, ierr, master, slave, FIRST
      use mpi,           only: MPI_DOUBLE_PRECISION

      implicit none

      t_sn = 0.0

      d0      = 1.0
      c_si    = 0.1
      Mrms    = 5.0

      if (master) then

         diff_nml(PROBLEM_CONTROL)

         rbuff(1) = d0
         rbuff(2) = Mrms
         rbuff(3) = c_si

      endif

      call MPI_Bcast(rbuff, buffer_dim, MPI_DOUBLE_PRECISION, FIRST, comm, ierr)

      if (slave) then

         d0           = rbuff(1)
         Mrms         = rbuff(2)
         c_si         = rbuff(3)

      endif

   end subroutine read_problem_par

!-----------------------------------------------------------------------------

   subroutine init_prob

      use constants,   only: xdim, ydim, zdim
      use domain,      only: dom
      use dataio_pub,  only: msg, printinfo
      use grid,        only: leaves
      use gc_list,     only: cg_list_element
      use grid_cont,   only: grid_container
      use initneutral, only: idnn,imxn,imyn,imzn,ienn, gamma_neu
      use fluidindex,  only: flind
      use mpisetup,    only: proc
      use func,        only: resample_gauss

      implicit none

      integer            :: i,j,k,m,n,l
      integer, parameter :: kp = 8
      real, dimension(6) :: mn
      real, dimension(3) :: deltav
      real, dimension(:,:,:,:), allocatable :: dv
      real :: rms, cma, vol
      type(cg_list_element), pointer :: cgl
      type(grid_container), pointer :: cg
!      real :: somx,somy,somz

! Uniform equilibrium state

      call random_seed()

      cgl => leaves%first
      do while (associated(cgl))
         cg => cgl%cg

         allocate(dv(3, cg%n_(xdim), cg%n_(ydim), cg%n_(zdim)))

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
!                     deltav(1) = deltav(1) + mn(1)*sin(somx) + mn(2)*cos(somx)
!                     deltav(2) = deltav(2) + mn(3)*sin(somy) + mn(4)*cos(somy)
!                     deltav(3) = deltav(3) + mn(5)*sin(somz) + mn(6)*cos(somz)
                        deltav(1) = deltav(1) + mn(1)*sin(float(m)*cg%z(k)) + mn(2)*cos(float(n)*cg%y(j))
                        deltav(2) = deltav(2) + mn(3)*sin(float(m)*cg%x(i)) + mn(4)*cos(float(n)*cg%z(k))
                        deltav(3) = deltav(3) + mn(5)*sin(float(m)*cg%y(j)) + mn(6)*cos(float(n)*cg%x(i))
                     enddo
                  enddo
                  dv(:,i,j,k) = deltav(:)
               enddo
            enddo
         enddo
         vol = cg%n_(xdim)*cg%n_(ydim)*cg%n_(zdim)
         rms = sqrt( sum(dv**2) / vol )

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
                  cg%u(idnn,i,j,k) = d0
                  cg%u(imxn,i,j,k) = cg%u(idnn,i,j,k) * dv(1,i,j,k) * cma
                  cg%u(imyn,i,j,k) = cg%u(idnn,i,j,k) * dv(2,i,j,k) * cma
                  cg%u(imzn,i,j,k) = cg%u(idnn,i,j,k) * dv(3,i,j,k) * cma
                  cg%u(ienn,i,j,k) = c_si**2*d0/(gamma_neu*(gamma_neu-1.0))
                  cg%u(ienn,i,j,k) = cg%u(ienn,i,j,k) + 0.5*(cg%u(imxn,i,j,k)**2 + &
                       &                 cg%u(imyn,i,j,k)**2 + cg%u(imzn,i,j,k) )/cg%u(idnn,i,j,k)
                  cg%b(:,i,j,k)    = 0.0
                  cg%u(ienn,i,j,k) = cg%u(ienn,i,j,k) + 0.5*sum(cg%b(:,i,j,k)**2,1)
                  cg%u(flind%trc%beg:flind%trc%end, i, j, k)   = &
                        resample_gauss( cg%x(i), cg%y(j), cg%z(k), cg%dl(xdim), cg%dl(ydim), cg%dl(zdim), &
                                        0.05*dom%L_(xdim), 0.05*dom%L_(ydim), 0.05*dom%L_(zdim), 10)
               enddo
            enddo
         enddo

! Explosions

         deallocate(dv)

         cgl => cgl%nxt
      enddo

   end subroutine init_prob

end module initproblem
