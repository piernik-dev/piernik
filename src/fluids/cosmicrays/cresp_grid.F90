module cresp_grid
! pulled by COSM_RAY_ELECTRONS

contains


 subroutine grid_cresp_update
      use cg_leaves,      only: leaves
      use cg_list,        only: cg_list_element
      use constants,      only: xdim, ydim, zdim, LO, HI, pMAX
      use domain,         only: dom, is_multicg
      use func,           only: ekin, emag, operator(.equals.), operator(.notequals.)
      use grid_cont,      only: grid_container
      use cresp_variables, only: ind_p_lo, ind_p_up, cresp_taylor_order, taylor_coeff_2nd, taylor_coeff_3rd, &
                                ind_e_beg, ind_e_end, ind_n_beg, ind_n_end

      implicit none

!       class(component_fluid), pointer :: fl
      integer                         :: i, j, k !, icr, ipm, jpm, kpm
!       real                            :: cs_iso, xsn, ysn, zsn, r2, maxv
      type(cg_list_element),  pointer :: cgl
      type(grid_container),   pointer :: cg
      
      
!       logical, save :: frun = .true.
!       integer       :: cr_id         ! maybe we should make this variable global in the module and do not pass it as an argument?

!       i,j,k    = 0
   cgl => leaves%first
   cg => cgl%cg

!       ts =  set_timer(tmr_mgd, .true.)
!       call all_dirty

!       if (diff_explicit .or. (allow_explicit .and. dt/diff_dt_crs_orig<1)) then
! 
!          if (frun) then
!             if (master .and. diff_explicit) call warn("[cresp_multigrid_update:multigrid_solve_diff] Multigrid was initialized but is not used")
!             frun = .false.
!          endif
        print *,'PASS - cresp_grid'
!         print *, 'in domain cell(-2,-2,0) p_lo_init = cg%u(:, -24, -24, 0) = ',cg%u(ind_n_beg:ind_p_up, -2, -2, 0)  ! just some check, to be removed
!       do k = cg%lhn(zdim,LO), cg%lhn(zdim,HI)
!          do j = cg%lhn(ydim,LO), cg%lhn(ydim,HI)
!             do i = cg%lhn(xdim,LO), cg%lhn(xdim,HI)
              
!               cresp_vector%cresp_ind(1:ncre) = cg%u()
!             u_b = emag = pmag(bx(i,j,k), by(i,j,k), bz(i,j,k))
!             u_d = div_v
#ifdef VERBOSE
              print *, 'Output of cosmic ray electrons module for grid cell with coordinates i,j,k:', i, j, k
#endif /* VERBOSE */
!             call cresp_crs_update(2*cresp_dt, cg%u(cr_table(cren)), cg%u(cr_table(cree)), cg%u(cr_table(crepl), &
!             cg%u(cr_table(crepu), div_v ,emag(cg%b(i), cg%b(j), cg%b(k)) ! most likely 2 sweeps in one dt are executed in one step, we will test whether it's true or not
          
!           enddo
!         enddo
!       enddo

      
   end subroutine grid_cresp_update
   
   subroutine grid_cresp_initialization
   
!        do k = cg%lhn(zdim,LO), cg%lhn(zdim,HI)
!            do j = cg%lhn(ydim,LO), cg%lhn(ydim,HI)
!               do i = cg%lhn(xdim,LO), cg%lhn(xdim,HI)
              
!             u_b = emag = pmag(bx(i,j,k), by(i,j,k), bz(i,j,k))
!             u_d = div_v
#ifdef VERBOSE
              print *, 'Output of cosmic ray electrons module for grid cell with coordinates i,j,k:', i, j, k
#endif /* VERBOSE */
!               call cresp_init_state(two*dt,)
          
!              enddo
!           enddo
!        enddo
   
   end subroutine grid_cresp_initialization
   
end module cresp_grid