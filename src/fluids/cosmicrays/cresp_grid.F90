module cresp_grid
! pulled by COSM_RAY_ELECTRONS

! This module contains routines necessary to initialize, compute timestep for cre bins and to update spectrum in the whole domain
! as the crspectrum module operates on a single grid cell.
      use initcosmicrays, only: iarr_cre_e, iarr_cre_n !, iarr_cre_pl
      use initcrspectrum, only: ncre
      use global,         only: dt, t

      
      public        dt_cre      
      real(kind=8)                    :: dt_cre
contains

 subroutine cresp_update_grid
  use cg_leaves,      only: leaves
  use cg_list,        only: cg_list_element
  use constants,      only: xdim, ydim, zdim, zero !, pMAX, 
  use grid_cont,      only: grid_container
  use cresp_crspectrum, only:cresp_update_cell, printer
  use units,           only: s_len_u
  use initcrspectrum, only: spec_mod_trms, virtual_e, virtual_n !, p_fix, crel
  use named_array_list, only: qna
  use crhelpers,        only: divv_n
  use func,             only: emag, ekin, operator(.equals.), operator(.notequals.)
  use fluidindex, only: iarr_all_dn
  implicit none
    integer                         :: i, j, k 
    type(cg_list_element),  pointer :: cgl
    type(grid_container),   pointer :: cg
    real(kind=8), dimension(1:ncre)  :: n_cell, e_cell
    type(spec_mod_trms)  :: sptab
        i = 0; j = 0;  k = 0
        cgl => leaves%first

        n_cell = zero
        e_cell = zero
        do while (associated(cgl))
            cg => cgl%cg
            do k = cg%ks, cg%ke
                do j = cg%js, cg%je
                    do i = cg%is, cg%ie
                        sptab%ud = 0.0 ; sptab%ub = 0.0 ; sptab%ucmb = 0.0
                        n_cell    = cg%u(iarr_cre_n, i, j, k)
                        e_cell    = cg%u(iarr_cre_e, i, j, k)
                        sptab%ub = emag(cg%b(xdim,i,j,k), cg%b(ydim,i,j,k), cg%b(zdim,i,j,k))
                        sptab%ud = cg%q(qna%ind(divv_n))%point([i,j,k])
#ifdef VERBOSE
                        print *, 'Output of cosmic ray electrons module for grid cell with coordinates i,j,k:', i, j, k
#endif /* VERBOSE */    
                        call cresp_update_cell(2*dt, n_cell, e_cell, sptab, virtual_n(:,i,j,k), virtual_e(:,i,j,k))
                        cg%u(iarr_cre_n, i, j, k) = n_cell
                        cg%u(iarr_cre_e, i, j, k) = e_cell
!              diagnostic:
                        if (i.eq.12.and.j.eq.12.and.k.eq.0) then
                            call printer(t)      
                        endif
                    enddo
                enddo
            enddo
            cgl=>cgl%nxt
        enddo

  end subroutine cresp_update_grid
!----------------------------------------------------------------------------------------------------      
 subroutine cresp_init_grid
  use cg_leaves,      only: leaves
  use cg_list,        only: cg_list_element
  use constants,      only: I_ONE, I_FOUR, fpi, LO, HI, xdim, ydim, zdim, zero
  use grid_cont,      only: grid_container
  use initcrspectrum, only: ncre, f_init, p_up_init, p_lo_init, q_init, cre_eff, spec_mod_trms, initial_condition, bump_amp, &
                            virtual_e, virtual_n
  use initcosmicrays, only: iarr_crn !, iarr_cre
  use cresp_crspectrum, only: cresp_init_state, allocate_all_allocatable, printer
  use units,          only: clight
  implicit none
    integer                         :: i, j, k
    type(cg_list_element),  pointer :: cgl
    type(grid_container),   pointer :: cg
    real(kind=8), allocatable, dimension(:)  :: cresp_arguments
    real(kind=8)                             :: max_amp_cr, f_amplitude
    real(kind=8), dimension(1:ncre) :: n_cell, e_cell
    type(spec_mod_trms)  :: sptab
    logical, save :: first_run = .true.
      if (first_run .eqv. .true.) then
        call allocate_all_allocatable

        i = 0; j = 0; k = 0
        cgl => leaves%first
        do while (associated(cgl))
            cg => cgl%cg
            cg%u(iarr_cre_e,:,:,:) = zero
            cg%u(iarr_cre_n,:,:,:) = zero
            max_amp_cr = maxval(cg%u(iarr_crn(1),:,:,:))
        
        if (.not. allocated(virtual_e)) allocate(virtual_e(1:2,cg%lhn(xdim,LO):cg%lhn(xdim,HI),cg%lhn(ydim,LO):cg%lhn(ydim,HI),&
                        cg%lhn(zdim,LO):cg%lhn(zdim,HI)))
        if (.not. allocated(virtual_n)) allocate(virtual_n(1:2,cg%lhn(xdim,LO):cg%lhn(xdim,HI),cg%lhn(ydim,LO):cg%lhn(ydim,HI),&
                        cg%lhn(zdim,LO):cg%lhn(zdim,HI)))
        virtual_e = zero
        virtual_n = zero
        
        n_cell = zero        
        e_cell = zero
            do k = cg%lhn(zdim,LO), cg%lhn(zdim,HI)
                do j = cg%lhn(ydim,LO), cg%lhn(ydim,HI)
                    do i = cg%lhn(xdim,LO), cg%lhn(xdim,HI)
                        cresp_arguments = 0.0
                   ! Every initial condition should be normalized before initializing Cosmic Ray Electron SPectrum module
                        if ( initial_condition == "powl") then
                            f_amplitude = 1/(fpi*clight*(p_lo_init**(I_FOUR))*(((p_up_init/p_lo_init)**(I_FOUR-q_init))-I_ONE)/ &
                                                                                                            (I_FOUR-q_init))
                     ! amplitude and distribution of electron energy density is inherited after those of nucleons, see crspectrum.pdf, eq. 29
                            f_amplitude = f_init*cg%u(iarr_crn(1),i,j,k)*cre_eff
                        endif
                        if (initial_condition == "bump") then
                            f_amplitude = cre_eff * cg%u(iarr_crn(1),i,j,k) ! for gaussian distribution & inheritance of spatial energy/number density after nucleons
!                             f_amplitude = bump_amp ! * clight 
                        endif
                        call cresp_init_state(n_cell, e_cell, f_amplitude, sptab)
#ifdef VERBOSE
                        print *, 'Output of cosmic ray electrons module for grid cell with coordinates i,j,k:', i, j, k
#endif /* VERBOSE */
                        cg%u(iarr_cre_n, i, j, k) = n_cell
                        cg%u(iarr_cre_e, i, j, k) = e_cell
                        if (i.eq.12.and.j.eq.12.and.k.eq.0) then ! diagnostics
                            call printer(t)
                        endif
                    enddo
                enddo
            enddo
            cgl=>cgl%nxt
        enddo
        first_run = .false. ! FIXME uncommenting results inf SIGFPE for some reason; whole subroutine is called twice.
      endif
   end subroutine cresp_init_grid

! ------------------------------------------------

 subroutine grid_cresp_timestep
  use cg_leaves,        only: leaves
  use cg_list,          only: cg_list_element
  use crhelpers,        only: divv_n
  use func,             only: emag !, operator(.equals.), operator(.notequals.)
  use grid_cont,        only: grid_container
  use constants,        only: xdim, ydim, zdim, LO, HI !, pMAX, 
  use named_array_list, only: qna
  use constants,        only: one
  use initcrspectrum,   only: spec_mod_trms
  use timestep_cresp,  only: cresp_timestep
  implicit none
    integer                         :: i, j, k 
    type(grid_container), pointer   :: cg
    type(cg_list_element), pointer  :: cgl
    real(kind=8)                    :: dt_cre_tmp
    type(spec_mod_trms)             :: sptab
        dt_cre = huge(one)
        dt_cre_tmp = huge(one)
        cgl => leaves%first
        do while (associated(cgl))
            cg => cgl%cg
            do k = cg%ks, cg%ke
                do j = cg%js, cg%je
                    do i = cg%is, cg%ie
                        sptab%ud = 0.0 ; sptab%ub = 0.0 ; sptab%ucmb = 0.0
                        sptab%ub = emag(cg%b(xdim,i,j,k), cg%b(ydim,i,j,k), cg%b(zdim,i,j,k))
                        sptab%ud = cg%q(qna%ind(divv_n))%point([i,j,k])
                        call cresp_timestep(dt_cre_tmp, sptab, cg%u(iarr_cre_n, i, j, k), cg%u(iarr_cre_e, i, j, k))
                        dt_cre = min(dt_cre, dt_cre_tmp)
                    enddo
                enddo
            enddo
            cgl=>cgl%nxt
        enddo
        dt_cre = 0.5 * dt_cre ! dt comes in to cresp_crspectrum with factor 2
  end subroutine grid_cresp_timestep
!----------------------------------------------------------------------------------------------------
 subroutine append_dissipative_terms(i,j,k) ! To be fixed
  use initcrspectrum, only: spec_mod_trms
  use named_array_list, only: qna
  use crhelpers,        only: divv_n
  use func,             only: emag
  use grid_cont,        only: grid_container
  use cg_list,          only: cg_list_element
  use constants,        only: xdim, ydim, zdim
  use cg_leaves,        only: leaves
  implicit none
    type(spec_mod_trms)  :: sptab
    type(grid_container), pointer :: cg
    type(cg_list_element), pointer:: cgl
    integer :: i,j,k
!Below - magnetic energy density and velocity divergence values are passed to sptab   
        cgl => leaves%first
        do while (associated(cgl))
            cg => cgl%cg
            sptab%ub = emag(cg%b(xdim,i,j,k), cg%b(ydim,i,j,k), cg%b(zdim,i,j,k))
            sptab%ud = cg%q(qna%ind(divv_n))%point([i,j,k])
            sptab%ucmb = 0.0
            cgl =>cgl%nxt
        enddo
 end subroutine append_dissipative_terms
   
end module cresp_grid
