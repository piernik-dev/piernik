module cresp_grid
! pulled by COSM_RAY_ELECTRONS
      use initcosmicrays, only: ncre, iarr_cre, iarr_cre_e, iarr_cre_n
      use global,         only: dt, t

      
      public        dt_cre      
      real(kind=8)                    :: dt_cre
contains


 subroutine grid_cresp_update
      use cg_leaves,      only: leaves
      use cg_list,        only: cg_list_element
!       use cg_level_finest,       only: finest
      use constants,      only: I_ONE, I_TWO, I_FOUR, xdim, ydim, zdim !, LO, HI, pMAX, 
!       use domain,         only: dom!, is_multicg
      use func,           only: ekin, emag, operator(.equals.), operator(.notequals.)
      use grid_cont,      only: grid_container
      use cresp_variables, only: ind_p_lo, ind_p_up, cresp_taylor_order, taylor_coeff_2nd, taylor_coeff_3rd, &
                                ind_e_beg, ind_e_end, ind_n_beg, ind_n_end
      use cresp_crspectrum, only:cresp_crs_update, printer
      use crhelpers,      only: divv_n
      use named_array_list, only: qna
      use units,           only: s_len_u

      implicit none

!       class(component_fluid), pointer :: fl
      integer                         :: i, j, k !, icr, ipm, jpm, kpm
!       real                            :: cs_iso, xsn, ysn, zsn, r2, maxv
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
        do k = cg%ks, cg%ke
           do j = cg%js, cg%je
              do i = cg%is, cg%ie
  
                  cresp_arguments(I_ONE:I_TWO*ncre+I_TWO)    = cg%u(iarr_cre, i, j, k)
!                   cresp_arguments(ncre+1:2*ncre) = cg%u(ind_e_beg:ind_e_end, i, j, k)
!                   cresp_arguments(2*ncre+1) = cg%u(ind_p_lo, i, j, k)
!                   cresp_arguments(2*ncre+2) = cg%u(ind_p_up, i, j, k)
                  
                  cresp_arguments(2*ncre+3) = emag(cg%b(xdim,i,j,k), cg%b(ydim,i,j,k), cg%b(zdim,i,j,k))/cg%dvol  !!! module works properly for small emag. Should emag be > 1e-4, negative values will appear
!                   cresp_arguments(2*ncre+3) = emag(cg%b(xdim,i,j,k), cg%b(ydim,i,j,k), cg%b(zdim,i,j,k))/(4*pi*cg%dvol)
                  cresp_arguments(2*ncre+4) = cg%q(qna%ind(divv_n))%point([i,j,k])/cg%dvol
                  
!                   cresp_arguments(2*ncre+4) = 5.0e-5
#ifdef VERBOSE
              print *, 'Output of cosmic ray electrons module for grid cell with coordinates i,j,k:', i, j, k
#endif /* VERBOSE */
! #ifndef DIFF_TEST
              call cresp_crs_update(2*dt, cresp_arguments, dt_cre_tmp) !cg%u(cr_table(cren)), cg%u(cr_table(cree)), cg%u(cr_table(crepl), &
              cg%u(iarr_cre, i, j, k) = cresp_arguments(I_ONE:I_TWO*ncre+I_TWO)
!              diagnostic:
                if (i.eq.34.and.j.eq.34.and.k.eq.0) then
                      call printer(t)      
                       print *, '   emag = ', emag(cg%b(xdim,i,j,k), cg%b(ydim,i,j,k), cg%b(zdim,i,j,k))/cg%dvol !/(4*pi*cg%dvol)
!                       print *, cresp_arguments(2*ncre+3)
!                       print *, 'cg%u(iarr_cre(e),34,34,:) =', cg%u(iarr_cre_e,34,34,0)
                      print *, 'cg%u(iarr_cre(n),34,34,:) =', cg%u(iarr_cre_n,34,34,0)
!                       print *, 'p ', crel%p
! !                       print *, 'q ', crel%q
!                       print *, 'f ', crel%f
!                       print *, 'plo, pup = ', cg%u(flind%cre%plo,i,j,k),cg%u(flind%cre%pup,i,j,k)
                      print *, '-------------------------'
                endif
!               if(i.eq.1.and.j.eq.1.and.k.eq.0) call printer(t)          
              
              dt_cre = min(dt_cre, dt_cre_tmp)
!               if (dt_cre .ge. dt_cre_tmp) then
!                       dt_cre = dt_cre_tmp
!               endif
           enddo
         enddo
!          print *,'emag = ', cresp_arguments(2*ncre+3)
!          print *,'dvol = ', cg%dvol
       enddo
!        print *, 'cresp_args: = ', cresp_arguments(iarr_cre)
      cgl=>cgl%nxt
!       print *,'  dt_cre =   ', dt_cre !, dt_cre_tmp
!       print *,'dt_cre grid update = ', dt_cre, ' ==0.5!'
      enddo
      
   end subroutine grid_cresp_update
   
   subroutine grid_cresp_initialization
   
      use cg_leaves,      only: leaves
      use cg_list,        only: cg_list_element
      use constants,      only: xdim, ydim, zdim, I_ONE, I_TWO, I_FOUR, fpi !, LO, HI, pMAX,
!       use domain,         only: dom, is_multicg
      use func,           only: ekin, emag, operator(.equals.), operator(.notequals.)
      use grid_cont,      only: grid_container
      use cresp_variables, only: ind_p_lo, ind_p_up, cresp_taylor_order, taylor_coeff_2nd, taylor_coeff_3rd, &
                                ind_e_beg, ind_e_end, ind_n_beg, ind_n_end
      use initcosmicrays, only: ncre, iarr_cre, f_init, p_up_init, p_lo_init, q_init, cre_eff, iarr_crn
      use cresp_crspectrum, only: cresp_init_state
      use units,          only: clight
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
     
         do k = cg%ks, cg%ke
           do j = cg%js, cg%je
              do i = cg%is, cg%ie

                  cresp_arguments(I_ONE:2*ncre+2)    = cg%u(iarr_cre, i, j, k)
                  cresp_arguments(2*ncre+3) = emag(cg%b(xdim,i,j,k), cg%b(ydim,i,j,k),cg%b(zdim,i,j,k))*0.1/cg%dvol
                  
                  f_init = 1/(fpi*clight*(p_lo_init**(I_FOUR))*(((p_up_init/p_lo_init)**(I_FOUR-q_init))-I_ONE)/(I_FOUR-q_init))   !!! amplitude and distribution of electron energy density is inherited after those of nucleons, see crspectrum.pdf, eq. 29
!                    f_init = 1.0
                  f_init    = f_init*cg%u(iarr_crn(1),i,j,k)*cre_eff
                  
                  call cresp_init_state(dt, cresp_arguments, dt_cre_tmp)
#ifdef VERBOSE
              print *, 'Output of cosmic ray electrons module for grid cell with coordinates i,j,k:', i, j, k
#endif /* VERBOSE */
                  cg%u(iarr_cre, i, j, k) = cresp_arguments(I_ONE:I_TWO*ncre+I_TWO)
              
!                 if (i.eq.34.and.j.eq.34.and.k.eq.0) then
!                       print *, 'cg%u(iarr_cre(e),34,34,:) =', cg%u(iarr_cre_e,i,j,k)
!                       print *, 'cg%u(iarr_cre(n),34,34,:) =', cg%u(iarr_cre_n,i,j,k)
!                       print *, 'p ', crel%p
!                       print *, 'q ', crel%q
!                       print *, 'f ', crel%f
!                       print *,''
!                       print *, 'plo, pup = ', cg%u(flind%cre%plo,i,j,k),cg%u(flind%cre%pup,i,j,k)
!                       print *, '-------------------------'
!                 endif

                  if (dt_cre .ge. dt_cre_tmp) then
                      dt_cre = dt_cre_tmp
                  endif
!
           enddo
         enddo
!          print *,'ub = ', cresp_arguments(2*ncre+3)
!          print *, ' dt cre, dt_cre_tmp', dt_cre, dt_cre_tmp
       enddo
      cgl=>cgl%nxt
!       print *,'min_cre_dt = ', dt_cre
!       print *,'dt_cre grid update = ', dt_cre, ' ==0.5!'
!       dt_cre = 0.5
      
      enddo
!       stop
      deallocate(cresp_arguments)
   
   end subroutine grid_cresp_initialization
   
end module cresp_grid