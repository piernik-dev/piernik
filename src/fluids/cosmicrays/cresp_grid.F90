module cresp_grid
! pulled by COSM_RAY_ELECTRONS
      use initcosmicrays, only: ncre, iarr_cre
      use global,         only: dt, t

contains


 subroutine grid_cresp_update
      use cg_leaves,      only: leaves
      use cg_list,        only: cg_list_element
      use constants,      only: I_ONE, I_TWO, I_FOUR, xdim, ydim, zdim!, LO, HI, pMAX, 
      use domain,         only: dom!, is_multicg
      use func,           only: ekin, emag, operator(.equals.), operator(.notequals.)
      use grid_cont,      only: grid_container
      use cresp_variables, only: ind_p_lo, ind_p_up, cresp_taylor_order, taylor_coeff_2nd, taylor_coeff_3rd, &
                                ind_e_beg, ind_e_end, ind_n_beg, ind_n_end
      use cresp_crspectrum, only:cresp_crs_update, printer
      

      implicit none

!       class(component_fluid), pointer :: fl
      integer                         :: i, j, k !, icr, ipm, jpm, kpm
!       real                            :: cs_iso, xsn, ysn, zsn, r2, maxv
      type(cg_list_element),  pointer :: cgl
      type(grid_container),   pointer :: cg
      real(kind=8), allocatable, dimension(:)  :: cresp_arguments
      real(kind=8)                           :: dt_cre_tmp
      
      allocate(cresp_arguments(I_ONE:I_TWO*ncre+I_FOUR))
      !       logical, save :: frun = .true.
!       integer       :: cr_id         ! maybe we should make this variable global in the module and do not pass it as an argument?

   i = 0
   j = 0
   k = 0
   
   dt_cre_tmp = 10.0
   dt_cre = dt_cre_tmp
   
   cgl => leaves%first
   do while (associated(cgl))
     cg => cgl%cg
     cresp_arguments = 0.0
   !       ts =  set_timer(tmr_mgd, .true.)
!       call all_dirty

!       if (diff_explicit .or. (allow_explicit .and. dt/diff_dt_crs_orig<1)) then
! 
!          if (frun) then
!             if (master .and. diff_explicit) call warn("[cresp_multigrid_update:multigrid_solve_diff] Multigrid was initialized but is not used")
!             frun = .false.
!          endif
        print *,'PASS - cresp_grid'
!           print *, 'in domain cell(-2,-2,0) p_lo_init = cg%u(:, -2, -2, 0) = ',cg%u(ind_n_beg:ind_p_up, -2, -2, 0)  ! just some check, to be removed
        do k = cg%ks, cg%ke
    !            print *,'entering k', k
!              print *, 'emag = ',emag(cg%b(xdim,64,64,k), cg%b(ydim,64,64,k), cg%b(zdim,64,64,k))*1.0e-6
           do j = cg%js, cg%je
!              print *,'entering j', j
              do i = cg%is, cg%ie
!                 print *,'entering i', i
                  cresp_arguments(I_ONE:I_TWO*ncre+I_TWO)    = cg%u(iarr_cre, i, j, k)
!                   cresp_arguments(ncre+1:2*ncre) = cg%u(ind_e_beg:ind_e_end, i, j, k)
!                   cresp_arguments(2*ncre+1) = cg%u(ind_p_lo, i, j, k)
!                   cresp_arguments(2*ncre+2) = cg%u(ind_p_up, i, j, k)
                  
                  cresp_arguments(2*ncre+3) = 2.5e-7
                  cresp_arguments(2*ncre+4) = emag(cg%b(xdim,i,j,k), cg%b(ydim,i,j,k), cg%b(zdim,i,j,k))*1.0e-6
!                   print *,'emag = ', cresp_arguments(2*ncre+4)
!                   cresp_arguments(2*ncre+4) = 5.0e-5


! !                   print *, 'c_args = ', cresp_arguments
#ifdef VERBOSE
              print *, 'Output of cosmic ray electrons module for grid cell with coordinates i,j,k:', i, j, k
#endif /* VERBOSE */
<<<<<<< HEAD
                 call cresp_crs_update(2*dt, cresp_arguments, dt_cre_tmp) ! <- in this module n, e, p_lo and p_up are updated
                 cg%u(iarr_cre, i, j, k) = cresp_arguments(I_ONE:I_TWO*ncre+I_TWO)
                 
                 if (dt_cre .ge. dt_cre_tmp) then
                      dt_cre = dt_cre_tmp
                  endif
                 
=======
              call cresp_crs_update(2*dt, cresp_arguments) !cg%u(cr_table(cren)), cg%u(cr_table(cree)), cg%u(cr_table(crepl), &
              cg%u(iarr_cre, i, j, k) = cresp_arguments(I_ONE:I_TWO*ncre+I_TWO)
>>>>>>> parent of 6ad4962... !
!              diagnostic:
              if(i.eq.1.and.j.eq.1.and.k.eq.0) call printer(t)          
           enddo
         enddo
       enddo
      cgl=>cgl%nxt
      enddo
<<<<<<< HEAD
      
      print *,'dt_cre grid update = ', dt_cre, ' ==0.5!'
      dt_cre = 0.5
     
=======

      
      
>>>>>>> parent of 6ad4962... !
   end subroutine grid_cresp_update
   
   subroutine grid_cresp_initialization
   
      use cg_leaves,      only: leaves
      use cg_list,        only: cg_list_element
      use constants,      only: xdim, ydim, zdim, LO, HI, pMAX, I_ONE, I_TWO, I_FOUR
      use domain,         only: dom, is_multicg
      use func,           only: ekin, emag, operator(.equals.), operator(.notequals.)
      use grid_cont,      only: grid_container
      use cresp_variables, only: ind_p_lo, ind_p_up, cresp_taylor_order, taylor_coeff_2nd, taylor_coeff_3rd, &
                                ind_e_beg, ind_e_end, ind_n_beg, ind_n_end
      use initcosmicrays, only: ncre, iarr_cre
      use cresp_crspectrum, only: cresp_init_state
      implicit none

!       class(component_fluid), pointer :: fl
      integer                         :: i, j, k !, icr, ipm, jpm, kpm
!       real                            :: cs_iso, xsn, ysn, zsn, r2, maxv
      type(cg_list_element),  pointer :: cgl
      type(grid_container),   pointer :: cg
      real(kind=8), allocatable, dimension(:)  :: cresp_arguments
      
      allocate(cresp_arguments(I_ONE:I_TWO*ncre+I_FOUR))
      !       logical, save :: frun = .true.
!       integer       :: cr_id         ! maybe we should make this variable global in the module and do not pass it as an argument?

   i = 0
   j = 0
   k = 0
   cgl => leaves%first
   do while (associated(cgl))
     cg => cgl%cg
     cresp_arguments = 0.0
     
         do k = cg%ks, cg%ke
!             print *,'entering k', k
           do j = cg%js, cg%je
!              print *,'entering j', j
              do i = cg%is, cg%ie
!                 print *,'entering i', i
                  cresp_arguments(I_ONE:I_TWO*ncre+I_TWO)    = cg%u(iarr_cre, i, j, k)
!                   cresp_arguments(ncre+1:2*ncre) = cg%u(ind_e_beg:ind_e_end, i, j, k)
!                   cresp_arguments(2*ncre+1) = cg%u(ind_p_lo, i, j, k)
!                   cresp_arguments(2*ncre+2) = cg%u(ind_p_up, i, j, k)
!                   cresp_arguments(2*ncre+3) = 2.5e-7
!                   cresp_arguments(2*ncre+4) = 5.0e-5
   
                  call cresp_init_state(dt, cresp_arguments)
#ifdef VERBOSE
              print *, 'Output of cosmic ray electrons module for grid cell with coordinates i,j,k:', i, j, k
#endif /* VERBOSE */

                  cg%u(iarr_cre, i, j, k) = cresp_arguments(I_ONE:I_TWO*ncre+I_TWO)
<<<<<<< HEAD
                  dt_cre = 0.5
              enddo
=======

>>>>>>> parent of 6ad4962... !
           enddo
         enddo
       enddo
      cgl=>cgl%nxt
      enddo

   
   end subroutine grid_cresp_initialization
   
end module cresp_grid