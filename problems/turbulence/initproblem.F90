
#include "piernik.h"
#include "macros.h"

module initproblem

   implicit none

   private
   public :: read_problem_par, problem_initial_conditions, problem_pointers

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

      use dataio_pub, only: nh      ! QA_WARN required for diff_nml
      use mpisetup,   only: rbuff, master, slave, piernik_MPI_Bcast

      implicit none

      t_sn = 0.0

      d0   = 1.0
      c_si = 0.1
      Mrms = 5.0

      if (master) then

         if (.not.nh%initialized) call nh%init()
         open(newunit=nh%lun, file=nh%tmp1, status="unknown")
         write(nh%lun,nml=PROBLEM_CONTROL)
         close(nh%lun)
         open(newunit=nh%lun, file=nh%par_file)
         nh%errstr=""
         read(unit=nh%lun, nml=PROBLEM_CONTROL, iostat=nh%ierrh, iomsg=nh%errstr)
         close(nh%lun)
         call nh%namelist_errh(nh%ierrh, "PROBLEM_CONTROL")
         read(nh%cmdl_nml,nml=PROBLEM_CONTROL, iostat=nh%ierrh)
         call nh%namelist_errh(nh%ierrh, "PROBLEM_CONTROL", .true.)
         open(newunit=nh%lun, file=nh%tmp2, status="unknown")
         write(nh%lun,nml=PROBLEM_CONTROL)
         close(nh%lun)
         call nh%compare_namelist()

         rbuff(1) = d0
         rbuff(2) = Mrms
         rbuff(3) = c_si

      endif

      call piernik_MPI_Bcast(rbuff)

      if (slave) then

         d0       = rbuff(1)
         Mrms     = rbuff(2)
         c_si     = rbuff(3)

      endif

   end subroutine read_problem_par

!-----------------------------------------------------------------------------

   subroutine problem_initial_conditions

      use cg_leaves,  only: leaves
      use cg_list,    only: cg_list_element
      use constants,  only: xdim, ydim, zdim, ndims, LO, HI
      use dataio_pub, only: msg, printinfo
      use domain,     only: dom
      use fluidindex, only: flind
      use fluidtypes, only: component_fluid
      use func,       only: resample_gauss, ekin, emag
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

         allocate(dv(ndims, cg%lhn(xdim,LO):cg%lhn(xdim,HI), cg%lhn(ydim,LO):cg%lhn(ydim,HI), cg%lhn(zdim,LO):cg%lhn(zdim,HI)))

         do k = cg%lhn(zdim,LO), cg%lhn(zdim,HI)
            do j = cg%lhn(ydim,LO), cg%lhn(ydim,HI)
               do i = cg%lhn(xdim,LO), cg%lhn(xdim,HI)
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

         write(msg,'(2(a,g12.5),a,i4)')   "[initproblem:problem_initial_conditions] cma = ", cma, " rms = ", rms, " on ", proc
         call printinfo(msg, .true.)
         write(msg,'(2(a,g12.5),a,i4)')   "[initproblem:problem_initial_conditions] c_si = ", c_si, " l = ", l, " on ", proc
         call printinfo(msg, .true.)
         do k = cg%lhn(zdim,LO), cg%lhn(zdim,HI)
            do j = cg%lhn(ydim,LO), cg%lhn(ydim,HI)
               do i = cg%lhn(xdim,LO), cg%lhn(xdim,HI)
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

   end subroutine problem_initial_conditions

end module initproblem
