module cresp_grid
! pulled by COSM_RAY_ELECTRONS
      use initcosmicrays, only: ncre, iarr_cre
      use global,         only: dt, t
!       implicit none
  
      public        dt_cre      
      
      real(kind=8)                    :: dt_cre
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
!           print *, 'in domain cell(-2,-2,0) p_lo_init = cg%u(:, -2, -2, 0) = ',cg%u(ind_n_beg:ind_p_up, -2, -2, 0)  ! just some check, to be removed
        do k = cg%ks, cg%ke
!            print *,'entering k', k
!              print *, 'emag = ',emag(cg%b(xdim,64,64,k), cg%b(ydim,64,64,k), cg%b(zdim,64,64,k))*1.0e-6
           do j = cg%js, cg%je
!              print *,'entering j', j
              do i = cg%is, cg%ie
!                 print *,'entering i', i
                  cresp_arguments(I_ONE:I_TWO*ncre+I_TWO)    = cg%u(iarr_cre, i, j, k)

                  cresp_arguments(2*ncre+3) = 2.5e-7
                  cresp_arguments(2*ncre+4) = emag(cg%b(xdim,i,j,k), cg%b(ydim,i,j,k), cg%b(zdim,i,j,k))/cg%dvol
#ifdef VERBOSE
              print *, 'Output of cosmic ray electrons module for grid cell with coordinates i,j,k:', i, j, k
#endif /* VERBOSE */
                 call cresp_crs_update(2*dt, cresp_arguments, dt_cre_tmp) ! <- in this module n, e, p_lo and p_up are updated
                 cg%u(iarr_cre, i, j, k) = cresp_arguments(I_ONE:I_TWO*ncre+I_TWO)
                 
                 if (dt_cre .ge. dt_cre_tmp) then
                      dt_cre = dt_cre_tmp
                  endif
                 
!              diagnostic:
              if(i.eq.1.and.j.eq.1.and.k.eq.0) call printer(t)          
           enddo
         enddo
       enddo
      cgl=>cgl%nxt
      enddo
      
      print *,'dt_cre grid update = ', dt_cre, ' ==0.5!'
      dt_cre = 0.5
     
   end subroutine grid_cresp_update
   
!-------------------------------------------

   subroutine grid_cresp_initialization
   
      use cg_leaves,      only: leaves
      use cg_list,        only: cg_list_element
      use constants,      only: one, xdim, ydim, zdim, LO, HI, pMAX, I_ZERO, I_ONE, I_TWO, I_FOUR
      use domain,         only: dom, is_multicg
      use func,           only: ekin, emag, operator(.equals.), operator(.notequals.)
      use grid_cont,      only: grid_container
      use cresp_variables, only: ind_p_lo, ind_p_up, cresp_taylor_order, taylor_coeff_2nd, taylor_coeff_3rd, &
                                 ind_e_beg, ind_e_end, ind_n_beg, ind_n_end
      use initcosmicrays, only: ncre, iarr_cre
      use cresp_crspectrum, only: cresp_init_state
      use crhelpers,      only: divv_n
      use named_array_list, only: qna, wna
      implicit none

      integer                                :: i, j, k !, icr, ipm, jpm, kpm
      type(cg_list_element),  pointer        :: cgl
      type(grid_container),   pointer        :: cg
      real(kind=8), allocatable, dimension(:):: cresp_arguments
      real(kind=8)                           :: dt_cre_tmp
      real(kind=8)    :: divvel
  
  allocate(cresp_arguments(I_ONE:I_TWO*ncre+I_FOUR))

   i = 0
   j = 0
   k = 0

   cgl => leaves%first
   
   dt_cre_tmp = 10.0
   dt_cre = dt_cre_tmp
   do while (associated(cgl))
     cg => cgl%cg
     cresp_arguments = 0.0
     
        do k = cg%ks, cg%ke
           do j = cg%js, cg%je
              do i = cg%is, cg%ie
                  cresp_arguments(I_ONE:I_TWO*ncre+I_TWO)    = cg%u(iarr_cre, i, j, k)
                  cresp_arguments(2*ncre+3) = 0.e0
                  cresp_arguments(2*ncre+4) = emag(cg%b(xdim,i,j,k), cg%b(ydim,i,j,k), cg%b(zdim,i,j,k))/cg%dvol
                  divvel = cg%q(qna%ind(divv_n))%point([i,j,k])
!                   print *, '(i,j,k) divvel =   (', i,j,k,')  ', divvel

                  call cresp_init_state(dt, cresp_arguments, dt_cre_tmp)

#ifdef VERBOSE
              print *, 'Output of cosmic ray electrons module for grid cell with coordinates i,j,k:', i, j, k
#endif /* VERBOSE */
                  if (dt_cre .ge. dt_cre_tmp) then
                      dt_cre = dt_cre_tmp
                  endif
                  
                  cg%u(iarr_cre, i, j, k) = cresp_arguments(I_ONE:I_TWO*ncre+I_TWO)
                  dt_cre = 0.5
              enddo
           enddo
        enddo
   
   cgl=>cgl%nxt
   enddo

   end subroutine grid_cresp_initialization
   
end module cresp_grid
