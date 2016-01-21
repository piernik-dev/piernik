module cresp_grid
! pulled by COSM_RAY_ELECTRONS
      use initcosmicrays, only: ncre, iarr_cre, iarr_cre_e, iarr_cre_n, crel, cre_eff, f_init, iarr_crn
      use global,         only: dt, t

      
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
!       use fluidtypes,   only: var_numbers
      use fluidindex,     only: flind
      use crhelpers,      only: divv_n
      use named_array_list, only: qna, wna
      

      implicit none

!       class(component_fluid), pointer :: fl
      integer                         :: i, j, k !, icr, ipm, jpm, kpm
!       real                            :: cs_iso, xsn, ysn, zsn, r2, maxv
      type(cg_list_element),  pointer :: cgl
      type(grid_container),   pointer :: cg
      real(kind=8), allocatable, dimension(:)  :: cresp_arguments
      real(kind=8)                             :: dt_cre_tmp
!       type(var_numbers), intent(inout) :: flind
      
      allocate(cresp_arguments(I_ONE:I_TWO*ncre+I_FOUR))
      !       logical, save :: frun = .true.
!       integer       :: cr_id         ! maybe we should make this variable global in the module and do not pass it as an argument?

   i = 0
   j = 0
   k = 0

   dt_cre_tmp = 1.0
   dt_cre = dt_cre_tmp

   cgl => leaves%first
   do while (associated(cgl))
     cg => cgl%cg
     cresp_arguments = 0.0
!      print *, 'plo,pup ', cg%u(flind%cre%plo,34,34,k),cg%u(flind%cre%pup,34,34,k)
        do k = cg%ks, cg%ke
           do j = cg%js, cg%je
              do i = cg%is, cg%ie
                  cresp_arguments(:) = 0.0
                  cresp_arguments(I_ONE:I_TWO*ncre+I_TWO)    = cg%u(iarr_cre, i, j, k)
!                   cresp_arguments(ncre+1:2*ncre) = cg%u(ind_e_beg:ind_e_end, i, j, k)
                  
                  cresp_arguments(2*ncre+3) = emag(cg%b(xdim,i,j,k), cg%b(ydim,i,j,k), cg%b(zdim,i,j,k))/cg%dvol
                  cresp_arguments(2*ncre+4) = cg%q(qna%ind(divv_n))%point([i,j,k])/cg%dvol
                 
#ifdef VERBOSE
              print *, 'Output of cosmic ray electrons module for grid cell with coordinates i,j,k:', i, j, k
#endif /* VERBOSE */
              call cresp_crs_update(2*dt, cresp_arguments, dt_cre_tmp) !cg%u(cr_table(cren)), cg%u(cr_table(cree)), cg%u(cr_table(crepl), &
              cg%u(iarr_cre, i, j, k) = cresp_arguments(I_ONE:I_TWO*ncre+I_TWO)
!              diagnostic:
                if (i.eq.1.and.j.eq.34.and.k.eq.0) then
                      call printer(t)      
!                       print *, cresp_arguments(2*ncre+3)
                      print *, 'cg%u(iarr_cre(e),34,34,:) =', cg%u(iarr_cre_e,34,34,0)
                      print *, 'cg%u(iarr_cre(n),34,34,:) =', cg%u(iarr_cre_n,34,34,0)
                      print *, 'p ', crel%p
!                       print *, 'q ', crel%q
                      print *, 'f ', crel%f
                      print *, 'plo, pup = ', cg%u(flind%cre%plo,i,j,k),cg%u(flind%cre%pup,i,j,k)
                      print *, '-------------------------'
                endif
              
              if (dt_cre .ge. dt_cre_tmp) then
                      dt_cre = dt_cre_tmp
              endif
              
              
           enddo
         enddo
!          print *,'emag = ', cresp_arguments(2*ncre+3)
!          print *,'dvol = ', cg%dvol
       enddo

       cgl=>cgl%nxt

      enddo

   end subroutine grid_cresp_update
   
   subroutine grid_cresp_initialization
   
      use cg_leaves,      only: leaves
      use cg_list,        only: cg_list_element
      use constants,      only: xdim, ydim, zdim, LO, HI, pMAX, I_ONE, I_TWO, I_FOUR, fpi
      use domain,         only: dom, is_multicg
      use func,           only: ekin, emag, operator(.equals.), operator(.notequals.)
      use grid_cont,      only: grid_container
      use cresp_variables, only: ind_p_lo, ind_p_up, cresp_taylor_order, taylor_coeff_2nd, taylor_coeff_3rd, &
                                ind_e_beg, ind_e_end, ind_n_beg, ind_n_end, cnst_c
      use initcosmicrays, only: ncre, iarr_cre, p_lo_init, p_up_init, p_min_fix, q_init
      use cresp_crspectrum, only: cresp_init_state
      use units,          only: clight
      use fluidindex,     only: flind
      implicit none

      integer                         :: i, j, k !, icr, ipm, jpm, kpm
      type(cg_list_element),  pointer :: cgl
      type(grid_container),   pointer :: cg
      real(kind=8), allocatable, dimension(:)  :: cresp_arguments
      real(kind=8)                             :: dt_cre_tmp
      
      allocate(cresp_arguments(I_ONE:I_TWO*ncre+I_FOUR))
      !       logical, save :: frun = .true.
      !       integer       :: cr_id         ! maybe we should make this variable global in the module and do not pass it as an argument?
   i = 0
   j = 0
   k = 0
   dt_cre_tmp = 1.0
   dt_cre = dt_cre_tmp
   
   cgl => leaves%first
   do while (associated(cgl))
     cg => cgl%cg
     cresp_arguments = 0.0
!      print *, 'cg%u(iarr_cre,34,34,:) =',  cg%u(iarr_cre,34,34,:)
         do k = cg%ks, cg%ke
           do j = cg%js, cg%je
              do i = cg%is, cg%ie
!                   f_init = 1/(fpi*clight*(p_lo_init**(4))*((p_lo_init/p_min_fix)**(4-q_init)) + (p_up_init**(q_init))*((p_up_init/p_min_fix)**(4-q_init)))   !!! amplitude and distribution of electron energy density is inherited after those of nucleons
                  f_init = 1/(fpi*clight*(p_lo_init**(I_FOUR))*(((p_up_init/p_lo_init)**(I_FOUR-q_init)) -I_ONE)/(4-q_init))   !!! amplitude and distribution of electron energy density is inherited after those of nucleons, see crspectrum.pdf, eq. 29
!                    f_init = 1.0
                  f_init    = f_init*cg%u(iarr_crn(1),i,j,k)*cre_eff
                  cresp_arguments(:) = 0.0
!                   print *, f_init
!                   cresp_arguments(I_ONE:2*ncre+2)    = cg%u(iarr_cre, i, j, k)
                  cresp_arguments(2*ncre+3) = emag(cg%b(xdim,i,j,k), cg%b(ydim,i,j,k),cg%b(zdim,i,j,k))/cg%dvol
                  call cresp_init_state(dt, cresp_arguments, dt_cre_tmp)
#ifdef VERBOSE
              print *, 'Output of cosmic ray electrons module for grid cell with coordinates i,j,k:', i, j, k
#endif /* VERBOSE */
                  cg%u(iarr_cre, i, j, k) = cresp_arguments(I_ONE:I_TWO*ncre+I_TWO)
              
                if (i.eq.34.and.j.eq.34.and.k.eq.0) then
                      print *, 'cg%u(iarr_cre(e),34,34,:) =', cg%u(iarr_cre_e,i,j,k)
                      print *, 'cg%u(iarr_cre(n),34,34,:) =', cg%u(iarr_cre_n,i,j,k)
                      print *, 'p ', crel%p
                      print *, 'q ', crel%q
                      print *, 'f ', crel%f
                      print *,''
                      print *, 'plo, pup = ', cg%u(flind%cre%plo,i,j,k),cg%u(flind%cre%pup,i,j,k)
                      print *, '-------------------------'
                endif
              

                  if (dt_cre .ge. dt_cre_tmp) then
                      dt_cre = dt_cre_tmp
                  endif
           enddo
         enddo
       enddo
!       print *, 'cg%u(iarr_cre,34,34,:) =',  cg%u(iarr_cre,34,34,:)
      cgl=>cgl%nxt
      print *,'PASS'
      enddo
!       print *,'c = ' ,clight 
   
   end subroutine grid_cresp_initialization
   
end module cresp_grid
