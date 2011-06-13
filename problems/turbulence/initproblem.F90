! $Id$

#include "piernik.h"
#include "macros.h"

#define RNG nb+1:nx-nb, nb+1:ny-nb, nb+1:nz-nb
module initproblem

! Initial condition for Sedov-Taylor explosion
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

      use dataio_pub,    only: ierrh, par_file, namelist_errh, compare_namelist, cmdl_nml      ! QA_WARN required for diff_nml
      use mpisetup,      only: rbuff, buffer_dim, comm, ierr, master, slave
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

      call MPI_Bcast(rbuff, buffer_dim, MPI_DOUBLE_PRECISION, 0, comm, ierr)

      if (slave) then

         d0           = rbuff(1)
         Mrms         = rbuff(2)
         c_si         = rbuff(3)

      endif

   end subroutine read_problem_par

!-----------------------------------------------------------------------------

   subroutine init_prob

      use dataio_pub,  only: msg, printinfo
      use grid,        only: cga
      use grid_cont,   only: cg_list_element, grid_container
      use initneutral, only: idnn,imxn,imyn,imzn,ienn, gamma_neu
      use mpisetup,    only: proc

      implicit none

      integer            :: i,j,k,m,n,l
      integer, parameter :: kp = 8
      real, dimension(6) :: mn
      real, dimension(3) :: deltav
      real, dimension(:,:,:,:), allocatable :: dv
      real :: rms,cma, vol
      type(cg_list_element), pointer :: cgl
      type(grid_container), pointer :: cg
!      real :: somx,somy,somz

! Uniform equilibrium state

      call random_seed()

      cgl => cga%cg_leafs%cg_l(1)
      do while (associated(cgl))
         cg => cgl%cg

         allocate(dv(3, cg%nx, cg%ny, cg%nz))

         do k = 1, cg%nz
            do j = 1, cg%ny
               do i = 1, cg%nx
                  deltav(:) = 0.0
                  do m=-kp,kp
                     do n = -kp,kp
                        call random_number(mn)
!                     somx = dpi*(float(n)*y(j) + float(m)*z(k)) / dom%Lx
!                     somy = dpi*(float(n)*x(i) + float(m)*z(k)) / dom%Ly
!                     somz = dpi*(float(n)*x(i) + float(m)*y(j)) / dom%Lz
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
         vol = cg%nx*cg%ny*cg%nz
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
         do k = 1, cg%nz
            do j = 1, cg%ny
               do i = 1, cg%nx
                  cg%u%arr(idnn,i,j,k) = d0
                  cg%u%arr(imxn,i,j,k) = cg%u%arr(idnn,i,j,k) * dv(1,i,j,k) * cma
                  cg%u%arr(imyn,i,j,k) = cg%u%arr(idnn,i,j,k) * dv(2,i,j,k) * cma
                  cg%u%arr(imzn,i,j,k) = cg%u%arr(idnn,i,j,k) * dv(3,i,j,k) * cma
                  cg%u%arr(ienn,i,j,k) = c_si**2*d0/(gamma_neu*(gamma_neu-1.0))
                  cg%u%arr(ienn,i,j,k) = cg%u%arr(ienn,i,j,k) + 0.5*(cg%u%arr(imxn,i,j,k)**2 + &
                       &                 cg%u%arr(imyn,i,j,k)**2 + cg%u%arr(imzn,i,j,k) )/cg%u%arr(idnn,i,j,k)
                  cg%b%arr(:,i,j,k)    = 0.0
                  cg%u%arr(ienn,i,j,k) = cg%u%arr(ienn,i,j,k) + 0.5*sum(cg%b%arr(:,i,j,k)**2,1)
               enddo
            enddo
         enddo

! Explosions

         deallocate(dv)

         cgl => cgl%nxt
      enddo

   end subroutine init_prob

end module initproblem
