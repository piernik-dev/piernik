! $Id$

#include "piernik.h"
#include "macros.h"

#define RNG nb+1:nx-nb, nb+1:ny-nb, nb+1:nz-nb
module initproblem

! Initial condition for Sedov-Taylor explosion
   implicit none

   private
   public :: read_problem_par, init_prob

   integer           :: n_sn
   real              :: d0, Mrms, t_sn, c_si

   namelist /PROBLEM_CONTROL/  d0, c_si, Mrms

contains

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

      use arrays,        only: u,b
      use dataio_pub,    only: msg, printinfo
      use grid,          only: cg
      use initneutral,   only: idnn,imxn,imyn,imzn,ienn, gamma_neu
      use mpisetup,      only: proc !, dom

      implicit none

      integer            :: i,j,k,m,n,l
      integer, parameter :: kp = 8
      real, dimension(6) :: mn
      real, dimension(3) :: deltav
      real, dimension(3, cg%nx, cg%ny, cg%nz) :: dv
      real :: rms,cma, vol
!      real :: somx,somy,somz

! Uniform equilibrium state

      call random_seed()
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
!                     deltav(1) = deltav(1) + mn(1)*dsin(somx) + mn(2)*dcos(somx)
!                     deltav(2) = deltav(2) + mn(3)*dsin(somy) + mn(4)*dcos(somy)
!                     deltav(3) = deltav(3) + mn(5)*dsin(somz) + mn(6)*dcos(somz)
                     deltav(1) = deltav(1) + mn(1)*dsin(float(m)*cg%z(k)) &
                                           + mn(2)*dcos(float(n)*cg%y(j))
                     deltav(2) = deltav(2) + mn(3)*dsin(float(m)*cg%x(i)) &
                                           + mn(4)*dcos(float(n)*cg%z(k))
                     deltav(3) = deltav(3) + mn(5)*dsin(float(m)*cg%y(j)) &
                                           + mn(6)*dcos(float(n)*cg%x(i))
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
               u(idnn,i,j,k) = d0
               u(imxn,i,j,k) = u(idnn,i,j,k) * dv(1,i,j,k) * cma
               u(imyn,i,j,k) = u(idnn,i,j,k) * dv(2,i,j,k) * cma
               u(imzn,i,j,k) = u(idnn,i,j,k) * dv(3,i,j,k) * cma
               u(ienn,i,j,k) = c_si**2*d0/(gamma_neu*(gamma_neu-1.0))
               u(ienn,i,j,k) = u(ienn,i,j,k) + 0.5*(u(imxn,i,j,k)**2 &
                             + u(imyn,i,j,k)**2 + u(imzn,i,j,k) )/u(idnn,i,j,k)
               b(:,i,j,k)   = 0.0
               u(ienn,i,j,k)   = u(ienn,i,j,k) + 0.5*sum(b(:,i,j,k)**2,1)
            enddo
         enddo
      enddo

! Explosions

   end subroutine init_prob

end module initproblem
