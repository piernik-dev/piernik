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

module initproblem

   use constants,    only: dsetnamelen

   implicit none

   private
   public :: read_problem_par, problem_initial_conditions, problem_pointers

   enum, bind(C)
      enumerator :: ADDED = 1, NOT_ADDED
   end enum
   enum, bind(C)
      enumerator :: RANDOM = 1, SINE
   end enum

   real                     :: d0, r_in, r_out, f_in, f_out
   real                     :: dens_exp      !< exponent in profile density \f$\rho(R) = \rho_0 R^{-k}\f$
   real                     :: eps           !< dust to gas ratio
   real                     :: amplify       !< relative amplitude of radial bar (build test only)
   integer(kind=4)          :: cutoff_ncells !< width of cut-off profile
   real, save               :: T_inner = 0.0 !< Orbital period at the inner boundary
   real, save               :: max_vy = -HUGE(1.0) !< Maximum tangential dust velocity
   integer(kind=4), save    :: noise_added = NOT_ADDED !< whether noise has been already added
   !>
   !! \f$\tau\f$ in \f$\frac{Du}{Dt} = - \frac{u-u_0}{\tau}f(R)
   !! when initproblem::problem_customize_solution is used
   !<
   real                     :: dumping_coeff, amp_noise
   logical                  :: use_inner_orbital_period  !< use 1./T_inner as dumping_coeff
   integer(kind=4)          :: amp_func  !< 1 - random, 2 - sine

   integer(kind=4), parameter :: ngauss = 4
   real, dimension(ngauss)  :: gauss
   character(len=dsetnamelen), parameter :: inid_n = "u_0"

   namelist /PROBLEM_CONTROL/  d0, r_in, r_out, f_in, f_out, dens_exp, eps, dumping_coeff, use_inner_orbital_period, &
      & amp_noise, amp_func, gauss, amplify

contains
!-----------------------------------------------------------------------------------------------------------------------
   subroutine problem_pointers

      use dataio_user,           only: user_attrs_wr, user_attrs_rd
      use user_hooks,            only: problem_customize_solution, problem_grace_passed, problem_post_restart
      use gravity,               only: grav_pot_3d
#ifdef HDF5
      use dataio_user,           only: user_vars_hdf5
#endif /* HDF5 */

      implicit none

      user_attrs_wr => my_attrs_wr
      user_attrs_rd => my_attrs_rd
      problem_customize_solution => problem_customize_solution_kepler
      problem_grace_passed => si_grace_passed
      problem_post_restart => kepler_problem_post_restart
      grav_pot_3d => my_grav_pot_3d
#ifdef HDF5
      user_vars_hdf5 => prob_vars_hdf5
#endif /* HDF5 */

   end subroutine problem_pointers
!-----------------------------------------------------------------------------------------------------------------------
   subroutine register_user_var

      use cg_list_global,   only: all_cg
      use constants,        only: AT_NO_B
      use named_array_list, only: wna

      implicit none

      integer(kind=4) :: dim4 !< BEWARE: workaround for gcc-4.6 bug

      dim4 = wna%lst(wna%fi)%dim4
      call all_cg%reg_var(inid_n, restart_mode = AT_NO_B, dim4 = dim4)

   end subroutine register_user_var
!-----------------------------------------------------------------------------------------------------------------------
   subroutine read_problem_par

      use dataio_pub,            only: nh      ! QA_WARN required for diff_nml
      use mpisetup,              only: rbuff, ibuff, lbuff, master, slave, piernik_MPI_Bcast

      implicit none

      d0               = 1.0
      amp_noise        = 1.e-6
      amp_func         = RANDOM
      gauss            = 0.0

      r_in             = 0.5
      f_in             = 10.0
      r_out            = 2.1
      f_out            = 0.0

      dens_exp         = 0.0
      eps              = 1.0

      dumping_coeff    = 1.0

      use_inner_orbital_period = .false.

      amplify = 0.0

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

         ibuff(1) = amp_func

         lbuff(1) = use_inner_orbital_period

         rbuff(1) = d0
         rbuff(2) = r_in
         rbuff(3) = r_out
         rbuff(4) = f_in
         rbuff(5) = f_out
         rbuff(6) = dens_exp
         rbuff(7) = eps
         rbuff(8) = dumping_coeff
         rbuff(9) = amp_noise
         rbuff(10:13) = gauss
         rbuff(14) = amplify

      endif

      call piernik_MPI_Bcast(rbuff)
      call piernik_MPI_Bcast(ibuff)
      call piernik_MPI_Bcast(lbuff)

      if (slave) then

         amp_func         = ibuff(1)

         use_inner_orbital_period = lbuff(1)

         d0               = rbuff(1)
         r_in             = rbuff(2)
         r_out            = rbuff(3)
         f_in             = rbuff(4)
         f_out            = rbuff(5)
         dens_exp         = rbuff(6)
         eps              = rbuff(7)
         dumping_coeff    = rbuff(8)
         amp_noise        = rbuff(9)
         gauss            = rbuff(10:13)
         amplify          = rbuff(14)

      endif

      call register_user_var

   end subroutine read_problem_par
!-----------------------------------------------------------------------------------------------------------------------
   subroutine si_grace_passed

      implicit none

      select case (amp_func)
         case (SINE)
            call add_sine
         case default
            call add_random_noise
      end select

   end subroutine si_grace_passed
!-----------------------------------------------------------------------------------------------------------------------
   subroutine add_sine

      use cg_leaves,  only: leaves
      use cg_list,    only: cg_list_element
      use constants,  only: dpi, xdim, zdim
      use dataio_pub, only: printinfo, warn
      use domain,     only: dom
      use fluidindex, only: flind
      use grid_cont,  only: grid_container
      use mpisetup,   only: master

      implicit none

      real                           :: kx, kz
      integer                        :: i, k
      type(cg_list_element), pointer :: cgl
      type(grid_container),  pointer :: cg

      if (amp_noise <= 0.0 .or. noise_added == ADDED) then
         if (master) call printinfo("[initproblem:add_sine]: Noise has been added in previous run. Bailing out!")
         return
      endif
      if (.not. associated(flind%dst)) then
         if (master) call warn("[initproblem:add_sine]: Cannot add sine wave perturbation to dust, because there is no dust.")
         return
      endif
      if (master) call printinfo("[initproblem:add_sine]: adding sine wave perturbation to dust")

      kx = 15.*dpi/dom%L_(xdim)
      kz =  5.*dpi/dom%L_(zdim)

      cgl => leaves%first
      do while (associated(cgl))
         cg => cgl%cg

         do i = cg%is, cg%ie
            do k = cg%ks, cg%ke
               cg%u(flind%dst%imx,i,:,k) = cg%u(flind%dst%imx,i,:,k) + amp_noise*sin(kx*cg%x(i) + kz*cg%z(k)) * cg%u(flind%dst%idn,i,:,k)
               cg%u(flind%dst%imy,i,:,k) = cg%u(flind%dst%imy,i,:,k) + amp_noise*sin(kx*cg%x(i) + kz*cg%z(k)) * cg%u(flind%dst%idn,i,:,k)
               cg%u(flind%dst%imz,i,:,k) = cg%u(flind%dst%imz,i,:,k) + amp_noise*sin(kx*cg%x(i) + kz*cg%z(k)) * cg%u(flind%dst%idn,i,:,k)
            enddo
         enddo

         cgl => cgl%nxt
      enddo

      noise_added = ADDED

   end subroutine add_sine
!-----------------------------------------------------------------------------------------------------------------------
   subroutine add_random_noise

      use cg_leaves,  only: leaves
      use cg_list,    only: cg_list_element
      use constants,  only: xdim, ydim, zdim, LO, HI
      use dataio_pub, only: printinfo
      use grid_cont,  only: grid_container
      use fluidindex, only: flind
      use mpisetup,   only: proc, master

      implicit none

      integer, dimension(:), allocatable    :: seed
      integer                               :: n, clock, i
      real, dimension(:,:,:,:), allocatable :: noise
      type(cg_list_element), pointer        :: cgl
      type(grid_container),  pointer        :: cg

      if (amp_noise <= 0.0 .or. noise_added == ADDED) then
         if (master) call printinfo("[initproblem:add_random_noise]: Noise has been added in previous run. Bailing out!")
         return
      endif
      if (master) call printinfo("[initproblem:add_random_noise]: adding random noise to dust")
      call random_seed(size=n)
      allocate(seed(n))
      call system_clock(count=clock)
      seed = clock*proc + 37 * [ (i-1, i = 1, n) ]
      call random_seed(put=seed)
      deallocate(seed)

      cgl => leaves%first
      do while (associated(cgl))
         cg => cgl%cg

         allocate(noise(xdim:zdim, cg%lhn(xdim, LO):cg%lhn(xdim, HI), &
            & cg%lhn(ydim, LO):cg%lhn(ydim, HI), cg%lhn(zdim, LO):cg%lhn(zdim, HI)))
         call random_number(noise)
         cg%u(flind%dst%imx,:,:,:) = cg%u(flind%dst%imx,:,:,:) +amp_noise*(1.0-2.0*noise(xdim,:,:,:)) * cg%u(flind%dst%idn,:,:,:)
         cg%u(flind%dst%imy,:,:,:) = cg%u(flind%dst%imy,:,:,:) +amp_noise*(1.0-2.0*noise(ydim,:,:,:)) * cg%u(flind%dst%idn,:,:,:)
         cg%u(flind%dst%imz,:,:,:) = cg%u(flind%dst%imz,:,:,:) +amp_noise*(1.0-2.0*noise(zdim,:,:,:)) * cg%u(flind%dst%idn,:,:,:)
         deallocate(noise)

         cgl => cgl%nxt
      enddo

      noise_added = ADDED

   end subroutine add_random_noise
!-----------------------------------------------------------------------------------------------------------------------
   subroutine problem_initial_conditions

      use cg_leaves,          only: leaves
      use cg_list,            only: cg_list_element
      use constants,          only: dpi, xdim, ydim, zdim, GEO_RPZ, DST, LO, HI, pMAX
      use dataio_pub,         only: msg, printinfo, die
      use domain,             only: dom, is_multicg
      use fluidindex,         only: flind
      use fluidtypes,         only: component_fluid
      use gravity,            only: r_smooth, ptmass
      use grid_cont,          only: grid_container
      use mpisetup,           only: master, piernik_MPI_Allreduce
      use named_array_list,   only: wna
      use units,              only: newtong, gram, cm, kboltz, mH

      implicit none

      integer                         :: i, j, k, kmid, p
      real                            :: xi, yj, zk, rc, vz, sqr_gm, vr, vphi
      real                            :: gprim, H2

      real, dimension(:), allocatable :: grav, dens_prof, ln_dens_der
      class(component_fluid), pointer :: fl
      type(cg_list_element),  pointer :: cgl
      type(grid_container),   pointer :: cg

      integer :: xl, xr

!   Secondary parameters
      allocate(ln_dens_der(0)) ! suppress compiler warnings
      cgl => leaves%first
      do while (associated(cgl))
         cg => cgl%cg

         if (is_multicg) call die("[initproblem:problem_initial_conditions] multiple grid pieces per procesor not implemented yet") !nontrivial kmid, allocate

         sqr_gm = sqrt(newtong*ptmass)
         do k = cg%lhn(zdim, LO), cg%lhn(zdim, HI)
            if (cg%z(k) < 0.0) kmid = k       ! the midplane is in between ksmid and ksmid+1
         enddo

         if (dom%geometry_type == GEO_RPZ) then
            if (master) then
               call printinfo("------------------------------------------------------------------")
               call printinfo(" Assuming temperature profile for MMSN ")
               call printinfo(" T(R) = 150 ( R / 1 AU )^(-0.429) K")
               write(msg,'(A,F5.1,A)') " T(xmin) = ", mmsn_T(dom%edge(xdim, LO))," K"
               call printinfo(msg)
               write(msg,'(A,F5.1,A)') " T(xmax) = ", mmsn_T(dom%edge(xdim, HI))," K"
               call printinfo(msg)
               write(msg,'(A,F5.1,A)') " T_mean  = ", 0.5*(mmsn_T(dom%edge(xdim, LO))+mmsn_T(dom%edge(xdim, HI)))," K"
               call printinfo(msg)
               write(msg,'(A,F9.5)') " cs2(T_mean) = ", kboltz * 0.5*(mmsn_T(dom%edge(xdim, LO))+mmsn_T(dom%edge(xdim, HI))) / mH
               call printinfo(msg)
               if (associated(flind%neu)) then
                  write(msg,'(A,F12.3,A)') " T_real(cs2) = ", flind%neu%cs2*mH/kboltz, " K"
                  call printinfo(msg)
               endif
               call printinfo("------------------------------------------------------------------")
            endif

            xl = cg%lhn(xdim, LO)
            xr = cg%lhn(xdim, HI)
            if (.not.allocated(grav)) allocate(grav(xl:xr))
            if (size(ln_dens_der) /= xr-xl+1) deallocate(ln_dens_der)
            allocate(ln_dens_der(xl:xr))
            if (.not.allocated(dens_prof)) allocate(dens_prof(xl:xr))

            grav = compute_gravaccelR(cg)
            dens_prof(:) = d0 * cg%x(:)**(-dens_exp)  * gram / cm**2

            !! \f$ v_\phi = \sqrt{R\left(c_s^2 \partial_R \ln\rho + \partial_R \Phi \right)} \f$
            ln_dens_der  = log(dens_prof)
            ln_dens_der(xl+1:xr)  = ( ln_dens_der(xl+1:xr) - ln_dens_der(xl:xr-1) ) / cg%dx
            ln_dens_der(xl)       = ln_dens_der(xl+1)
            T_inner               = dpi*cg%x(cg%is) / sqrt( abs(grav(cg%is)) * cg%x(cg%is) )
            write(msg,*) "T_inner = ", T_inner
            if (master) call printinfo(msg)
            write(msg,*) "III Kepler Law gives T = ", sqr_gm/dpi , " yr at 1 AU"
            if (master) call printinfo(msg)

            do p = 1, flind%fluids
               fl => flind%all_fluids(p)%fl
               if (fl%tag /= DST .and. master) then
                  write(msg,'(A,F9.5)') "[initproblem:initprob] cs2 used = ", fl%cs2
                  call printinfo(msg)
               endif

               do j = cg%lhn(ydim, LO), cg%lhn(ydim, HI)
                  yj = cg%y(j)
                  do i = cg%lhn(xdim, LO), cg%lhn(xdim, HI)
                     xi = cg%x(i)
                     rc = xi + r_smooth

                     gprim = newtong*ptmass / xi**3
                     if (fl%cs > 0) then
                        H2 = 2.0*fl%cs2/gprim
                     else
                        H2 = 1.0
                     endif

                     vphi = 0.
                     do k = cg%lhn(zdim, LO), cg%lhn(zdim, HI)
                        zk = cg%z(k)
                        ! Next few lines are build test related
                        if (abs(yj - dom%C_(ydim)) < 0.5*xi*cg%dl(ydim)) then
                           cg%u(fl%idn,i,j,k) = dens_prof(i) * amplify
                        else
                           cg%u(fl%idn,i,j,k) = dens_prof(i)
                        endif
                        if (fl%tag == DST) cg%u(fl%idn,i,j,k) = eps * cg%u(fl%idn,i,j,k)

                        vr   = 0.0
                        ! that condition is not necessary since cs2 == 0.0 for dust
                        if (fl%tag /= DST) then
                           vphi = sqrt( max(cg%x(i)*(fl%cs2*ln_dens_der(i) + abs(grav(i))),0.0) )
                        else
                           vphi = sqrt( max(abs(grav(i)) * rc, 0.0))
                        endif
                        vz   = 0.0

                        cg%u(fl%imx,i,j,k) = vr   * cg%u(fl%idn,i,j,k)
                        cg%u(fl%imy,i,j,k) = vphi * cg%u(fl%idn,i,j,k)
                        cg%u(fl%imz,i,j,k) = vz   * cg%u(fl%idn,i,j,k)
                        if (fl%ien > 0) then
                           cg%u(fl%ien,i,j,k) = fl%cs2/(fl%gam_1)*cg%u(fl%idn,i,j,k)
                           cg%u(fl%ien,i,j,k) = cg%u(fl%ien,i,j,k) + 0.5*(vr**2+vphi**2+vz**2)*cg%u(fl%idn,i,j,k)
                        endif
                     enddo
                  enddo
               enddo

            enddo
            cg%w(wna%ind(inid_n))%arr(:,:,:,:) = cg%u(:,:,:,:)
            cg%b(:,:,:,:) = 0.0
            if (allocated(grav)) deallocate(grav)
            if (allocated(dens_prof)) deallocate(dens_prof)
            if (allocated(ln_dens_der)) deallocate(ln_dens_der)
         else
            call die("[initproblem:problem_initial_conditions] I don't know what to do... :/")
         endif
         max_vy = max(max_vy, maxval(abs(cg%u(flind%dst%imy,:,:,:))/cg%u(flind%dst%idn,:,:,:)) )
         cgl => cgl%nxt
      enddo

      call piernik_MPI_Allreduce(max_vy, pMAX)

   end subroutine problem_initial_conditions
!-----------------------------------------------------------------------------------------------------------------------
   real function mmsn_T(r)
      implicit none
      real, intent(in) :: r         ! [AU]
      real, parameter  :: T_0 = 150 ! [K]
      real, parameter  :: k   = 0.429

      mmsn_T = T_0 * r**(-k)
   end function mmsn_T
!-----------------------------------------------------------------------------------------------------------------------
   subroutine kepler_problem_post_restart

      use cg_leaves,        only: leaves
      use cg_list,          only: cg_list_element
      use constants,        only: b0_n, fluid_n
      use all_boundaries,   only: all_fluid_boundaries
      use named_array_list, only: wna
      use grid_cont,        only: grid_container
#ifdef TRACER
      use constants,        only: xdim, ydim, zdim
      use func,             only: resample_gauss
      use fluidindex,       only: flind
#endif /* TRACER */

      implicit none

      type(cg_list_element), pointer :: cgl
      type(grid_container),  pointer :: cg
#ifdef TRACER
      integer                        :: i, j, k
#endif /* TRACER */

      call my_grav_pot_3d   ! reset gp, to get right values on the boundaries

      cgl => leaves%first
      do while (associated(cgl))
         cg => cgl%cg
         cg%u  => cg%w(wna%ind(inid_n))%arr  ! BEWARE: Don't do things like that without parental supervision
         cg%b = 0.0
         cg%w(wna%ind(b0_n))%arr = 0.0
         cgl => cgl%nxt
      enddo

      wna%fi = wna%ind(inid_n)    ! BEWARE: or things like that...
      call all_fluid_boundaries   ! all_fluid_boundaries properly set boundaries for %u pointer and %wna%fi
      wna%fi = wna%ind(fluid_n)   ! revert first ugly hack

      cgl => leaves%first
      do while (associated(cgl))
         cg => cgl%cg
         cg%u  => cg%w(wna%fi)%arr ! Quick! Revert to sane state before anyone notices (2nd ugly hack)
         cgl => cgl%nxt
      enddo

#ifdef TRACER
      if (gauss(4) > 0.0) then
         cgl => leaves%first
         do while (associated(cgl))
            cg => cgl%cg

            do k = cg%lhn(zdim, LO), cg%lhn(zdim, HI)
               do j = cg%lhn(ydim, LO), cg%lhn(ydim, HI)
                  do i = cg%lhn(xdim, LO), cg%lhn(xdim, HI)

                     cg%u(flind%trc%beg:flind%trc%end, i, j, k)   = &
                          resample_gauss( cg%x(i) - gauss(1), cg%y(j) - gauss(2), cg%z(k) - gauss(3), &
                                          cg%dl(xdim), cg%dl(ydim), cg%dl(zdim), &
                                          gauss(4), gauss(4), gauss(4), 10)

                  enddo
               enddo
            enddo
            cgl => cgl%nxt
         enddo
      endif

      call all_fluid_boundaries
#endif /* TRACER */

   end subroutine kepler_problem_post_restart
!-----------------------------------------------------------------------------------------------------------------------
   subroutine problem_customize_solution_kepler(forward)

      use cg_leaves,        only: leaves
      use cg_list,          only: cg_list_element
      use constants,        only: xdim, ydim, zdim, LO, HI, pMAX
      use dataio_pub,       only: die!, warn, msg
      use domain,           only: is_multicg
      use global,           only: dt, smalld
      use grid_cont,        only: grid_container
      use all_boundaries,   only: all_fluid_boundaries
      use fluidindex,       only: flind!, iarr_all_mz, iarr_all_dn
      use mpisetup,         only: piernik_MPI_Allreduce
      use named_array_list, only: wna

      implicit none

      logical, intent(in)                     :: forward
      integer                                 :: i, j, k
      logical, save                           :: frun = .true.
      real, dimension(:,:), allocatable, save :: funcR
      logical, dimension(:,:,:), allocatable  :: adjust
      real, dimension(:,:,:), allocatable     :: vx_sign, vz_sign
      real                                    :: max_vx, mean_vy
      type(cg_list_element), pointer          :: cgl
      type(grid_container),  pointer          :: cg

      if (is_multicg) call die("[initproblem:problem_customize_solution_kepler] multiple grid pieces per procesor not implemented yet") !nontrivial

      max_vx = -HUGE(1.0)
      cgl => leaves%first
      do while (associated(cgl))
         cg => cgl%cg

         if (frun) then
            allocate(funcR(size(cg%u,dim=1), cg%lhn(xdim, LO):cg%lhn(xdim, HI)) )

            funcR(1,:) = -tanh((cg%x(:)-r_in+1.0)**f_in) + 1.0 + max( tanh((cg%x(:)-r_out+1.0)**f_out), 0.0)

            if (use_inner_orbital_period) then
               funcR(1,:) = funcR(1,:) / T_inner
            else
               funcR(1,:) = funcR(1,:) * dumping_coeff
            endif
            frun = .false.
            funcR(:,:) = spread(funcR(1,:),1,size(cg%u,dim=1))

         endif

         do j = cg%lhn(ydim, LO), cg%lhn(ydim, HI)
            do k = cg%lhn(zdim, LO), cg%lhn(zdim, HI)
               cg%u(:,:,j,k) = cg%u(:,:,j,k) - dt*(cg%u(:,:,j,k) - cg%w(wna%ind(inid_n))%arr(:,:,j,k))*funcR(:,:)
            enddo
         enddo

!         where ( cg%u(iarr_all_dn,:,:,:) < 2.0*smalld )
!            cg%u(iarr_all_mz,:,:,:) = cg%u(iarr_all_mz,:,:,:)*0.1
!         endwhere

         max_vx = max(max_vx, maxval( abs(cg%u(flind%neu%imx,:,:,:))/cg%u(flind%neu%idn,:,:,:) ))
         cgl => cgl%nxt
      enddo

      call piernik_MPI_Allreduce(max_vx, pMAX)

      cgl => leaves%first
      do while (associated(cgl))
         cg => cgl%cg
         if (.not.allocated(adjust)) &
            & allocate(adjust(cg%lhn(xdim, LO):cg%lhn(xdim, HI), &
               & cg%lhn(ydim, LO):cg%lhn(ydim, HI), cg%lhn(zdim, LO):cg%lhn(zdim, HI)))
         if (.not.allocated(vx_sign)) &
            & allocate(vx_sign(cg%lhn(xdim, LO):cg%lhn(xdim, HI), &
               & cg%lhn(ydim, LO):cg%lhn(ydim, HI), cg%lhn(zdim, LO):cg%lhn(zdim, HI)))
         if (.not.allocated(vz_sign)) &
            & allocate(vz_sign(cg%lhn(xdim, LO):cg%lhn(xdim, HI), &
               & cg%lhn(ydim, LO):cg%lhn(ydim, HI), cg%lhn(zdim, LO):cg%lhn(zdim, HI)))

         where (abs(cg%u(flind%dst%imx,:,:,:))/cg%u(flind%dst%idn,:,:,:) >= max_vx)
            adjust = .true.
         elsewhere
            adjust = .false.
         endwhere
         if ( any(adjust) ) then
            where (adjust)
               vx_sign = signum(cg%u(flind%dst%imx,:,:,:))
               vz_sign = signum(cg%u(flind%dst%imz,:,:,:))
            endwhere
            where (adjust)
               cg%u(flind%dst%idn,:,:,:) = max(cg%u(flind%dst%idn,:,:,:), 1.1*smalld)
               cg%u(flind%dst%imx,:,:,:) = vx_sign * cg%u(flind%neu%imx,:,:,:)/cg%u(flind%neu%idn,:,:,:) * cg%u(flind%dst%idn,:,:,:)
               cg%u(flind%dst%imz,:,:,:) = vz_sign * cg%u(flind%neu%imz,:,:,:)/cg%u(flind%neu%idn,:,:,:) * cg%u(flind%dst%idn,:,:,:)
            endwhere
         endif

         do i = cg%lhn(xdim, LO), cg%lhn(xdim, HI)
            do j = cg%lhn(ydim, LO), cg%lhn(ydim, HI)
               mean_vy = sum( cg%u(flind%dst%imy, i, j, :) /  cg%u(flind%dst%idn, i, j, :) ) / cg%n_(zdim)
               do k = cg%lhn(zdim, LO), cg%lhn(zdim, HI)
                  if ( (abs(cg%u(flind%dst%imy, i, j, k) - mean_vy * cg%u(flind%dst%idn, i, j, k)) &
                     / (mean_vy * cg%u(flind%dst%idn, i, j, k)) >= 0.1) .and. &
                         cg%u(flind%dst%idn, i, j, k) < 10.0*smalld ) then
                     cg%u(flind%dst%idn, i ,j ,k) = max(cg%u(flind%dst%idn, i, j, k), 1.1*smalld)
                     cg%u(flind%dst%imy, i, j, k) = mean_vy * cg%u(flind%dst%idn, i, j, k)
                  endif
               enddo
            enddo
         enddo

!        where ( cg%u(flind%dst%imy,:,:,:) > max_vy*cg%u(flind%dst%idn,:,:,:) )
!           cg%u(flind%dst%idn,:,:,:) = max(cg%u(flind%dst%idn,:,:,:), 1.1*smalld)
!           cg%u(flind%dst%imy,:,:,:) = cg%u(flind%neu%imy,:,:,:)/cg%u(flind%neu%idn,:,:,:) * cg%u(flind%dst%idn,:,:,:)
!        endwhere

         if (allocated(adjust))  deallocate(adjust)
         if (allocated(vx_sign)) deallocate(vx_sign)
         if (allocated(vz_sign)) deallocate(vz_sign)

         cgl => cgl%nxt
      enddo

      call all_fluid_boundaries

      return
      if (forward) i = j ! suppress compiler warnings on unused arguments

   end subroutine problem_customize_solution_kepler
!-----------------------------------------------------------------------------------------------------------------------
   subroutine my_grav_pot_3d

      use cg_leaves, only: leaves
      use cg_list,   only: cg_list_element
      use constants, only: xdim, zdim, LO, HI
      use gravity,   only: ptmass
      use grid_cont, only: grid_container
      use units,     only: newtong

      implicit none

      real                           :: r2
      integer                        :: i, k
      type(cg_list_element), pointer :: cgl
      type(grid_container),  pointer :: cg

      cgl => leaves%first
      do while (associated(cgl))
         cg => cgl%cg

         if (.not. cg%is_old) then
            do i = cg%lhn(xdim, LO), cg%lhn(xdim, HI)
               do k = cg%lhn(zdim, LO), cg%lhn(zdim, HI)
                  r2 = cg%x(i)**2! + cg%z(k)**2
                  cg%gp(i,:,k) = -newtong*ptmass / sqrt(r2)
               enddo
            enddo
         endif

         cgl => cgl%nxt
      enddo

   end subroutine my_grav_pot_3d
!-----------------------------------------------------------------------------------------------------------------------
   !>
   !! This function is a redundant code that does exactly what grav_pot2accel
   !! is supposed to. However, it allows to get rid of chicken-egg problem of
   !! density and gpot and fugly magic that used to do it in init_prob.
   !<
   function compute_gravaccelR(cg) result (grav)

      use constants,        only: xdim, LO, HI, half
      use gravity,          only: ptmass
      use grid_cont,        only: grid_container
      use units,            only: newtong

      implicit none

      type(grid_container),  pointer, intent(in)          :: cg
      real, dimension(cg%lhn(xdim, LO):cg%lhn(xdim, HI))  :: grav, gpot

      gpot = - newtong * ptmass / sqrt(cg%x(:)**2)
      grav(cg%lhn(xdim, LO)+1:cg%lhn(xdim, HI)-1) = &
         half*(gpot(cg%lhn(xdim, LO):cg%lhn(xdim, HI)-2) - gpot(cg%lhn(xdim, LO)+2:cg%lhn(xdim, HI))) / cg%dl(xdim)
      grav(cg%lhn(xdim, LO)) = grav(cg%lhn(xdim, LO)+1)
      grav(cg%lhn(xdim, HI)) = grav(cg%lhn(xdim, HI)-1)
   end function compute_gravaccelR
!-----------------------------------------------------------------------------------------------------------------------
   subroutine prob_vars_hdf5(var,tab, ierrh, cg)

      use fluidindex,   only: flind
      use grid_cont,    only: grid_container
      use interactions, only: epstein_factor

      implicit none

      character(len=*),               intent(in)    :: var
      real(kind=4), dimension(:,:,:), intent(inout) :: tab
      integer,                        intent(inout) :: ierrh
      type(grid_container), pointer,  intent(in)    :: cg

      ierrh = 0
      select case (trim(var))
         case ("tauf")
            tab(:,:,:) = real(epstein_factor(flind%neu%pos) / cg%u(flind%neu%idn,cg%is:cg%ie,cg%js:cg%je,cg%ks:cg%ke), 4)
         case default
            ierrh = -1
      end select

   end subroutine prob_vars_hdf5
!-----------------------------------------------------------------------------------------------------------------------
   subroutine my_attrs_wr(file_id)
      use hdf5, only: HID_T, SIZE_T
      use h5lt, only: h5ltset_attribute_double_f, h5ltset_attribute_int_f
      implicit none
      integer(HID_T), intent(in) :: file_id
      integer(SIZE_T), parameter :: bufsize = 1
      integer(kind=4)            :: error

      call h5ltset_attribute_double_f(file_id, "/", "T_inner", [T_inner], bufsize, error)
      call h5ltset_attribute_double_f(file_id, "/", "max_vy", [max_vy], bufsize, error)
      call h5ltset_attribute_int_f(file_id, "/", "noise_added", [noise_added], bufsize, error)

   end subroutine my_attrs_wr
!-----------------------------------------------------------------------------------------------------------------------
   subroutine my_attrs_rd(file_id)
      use constants, only: I_ONE
      use hdf5,      only: HID_T
      use h5lt,      only: h5ltget_attribute_double_f, h5ltget_attribute_int_f
      implicit none
      integer(HID_T), intent(in) :: file_id
      integer(kind=4)            :: error
      real, dimension(I_ONE)     :: rbuff
      integer(kind=4), dimension(I_ONE) :: ibuff

      call h5ltget_attribute_double_f(file_id, "/", "T_inner", rbuff, error)
      T_inner = rbuff(1)
      call h5ltget_attribute_double_f(file_id, "/", "max_vy", rbuff, error)
      max_vy = rbuff(1)

      call h5ltget_attribute_int_f(file_id, "/", "noise_added", ibuff, error)
      noise_added = ibuff(1)

   end subroutine my_attrs_rd
!-----------------------------------------------------------------------------------------------------------------------
   elemental function signum(a) result (b)
      implicit none
      real, intent(in) :: a
      real             :: b
      if (a > 0.0) then
         b = 1.0
      else
         b = -1.0
      endif
   end function signum
!-----------------------------------------------------------------------------------------------------------------------
end module initproblem
! vim: set tw=120:
