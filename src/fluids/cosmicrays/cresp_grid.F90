module cresp_grid
! pulled by COSM_RAY_ELECTRONS

! This module contains routines necessary to initialize, compute timestep for cre bins and to update spectrum in the whole domain
! as the crspectrum module operates on a single grid cell.
      use initcosmicrays, only: iarr_cre, iarr_cre_e, iarr_cre_n, iarr_cre_pl
      use initcrspectrum, only: ncre
      use global,         only: dt, t

      
      public        dt_cre      
      real(kind=8)                    :: dt_cre
contains


 subroutine cresp_update_grid
      use cg_leaves,      only: leaves
      use cg_list,        only: cg_list_element
      use constants,      only: I_ONE, I_TWO, I_FOUR, xdim, ydim, zdim, LO, HI !, pMAX, 
!       use domain,         only: dom
!       use func,           only: ekin, emag, operator(.equals.), operator(.notequals.)
      use grid_cont,      only: grid_container
      use cresp_crspectrum, only:cresp_update_cell, printer
!       use crhelpers,      only: divv_n
!       use named_array_list, only: qna
      use units,           only: s_len_u
      use initcrspectrum, only: spec_mod_trms, p_lo_nda, p_up_nda
!           use initcrspectrum, only: spec_mod_trms
    use named_array_list, only: qna
    use crhelpers,        only: divv_n
    use func,             only: emag

      implicit none
      integer                         :: i, j, k 
      type(cg_list_element),  pointer :: cgl
      type(grid_container),   pointer :: cg
      real(kind=8), allocatable, dimension(:)  :: cresp_arguments
      real(kind=8)                             :: dt_cre_tmp
      type(spec_mod_trms)  :: sptab
      
      allocate(cresp_arguments(I_ONE:I_TWO*ncre+I_TWO))
   i = 0
   j = 0
   k = 0

   dt_cre_tmp = 1.0
   dt_cre = dt_cre_tmp
   cgl => leaves%first
   do while (associated(cgl))
     cg => cgl%cg
     cresp_arguments = 0.0

          do k = cg%lhn(zdim,LO), cg%lhn(zdim,HI)
             do j = cg%lhn(ydim,LO), cg%lhn(ydim,HI)
               do i = cg%lhn(xdim,LO), cg%lhn(xdim,HI)
                  cresp_arguments(I_ONE:I_TWO*ncre+I_TWO)    = cg%u(iarr_cre, i, j, k)
                  cresp_arguments(2*ncre+1)          = p_lo_nda(i,j,k)
                  cresp_arguments(2*ncre+2)          = p_up_nda(i,j,k)
!                   call  append_dissipative_terms(i,j,k)  ! loads values of magnetic energy and vel divergence
                  sptab%ub = emag(cg%b(xdim,i,j,k), cg%b(ydim,i,j,k), cg%b(zdim,i,j,k))/10.0
                  sptab%ud = cg%q(qna%ind(divv_n))%point([i,j,k])
#ifdef VERBOSE
              print *, 'Output of cosmic ray electrons module for grid cell with coordinates i,j,k:', i, j, k
#endif /* VERBOSE */
! #ifndef DIFF_TEST
              call cresp_update_cell(2*dt, cresp_arguments, sptab) !cg%u(cr_table(cren)), cg%u(cr_table(cree)), cg%u(cr_table(crepl), &
              cg%u(iarr_cre, i, j, k) = cresp_arguments(I_ONE:I_TWO*ncre+I_TWO)
              p_lo_nda(i,j,k)         = cresp_arguments(I_TWO*ncre+I_ONE)
              p_up_nda(i,j,k)         = cresp_arguments(I_TWO*ncre+I_TWO)
!              diagnostic:
                if (i.eq.2.and.j.eq.2.and.k.eq.0) then
                      call printer(t)      
!                        print *, '   emag = ', emag(cg%b(xdim,i,j,k), cg%b(ydim,i,j,k), cg%b(zdim,i,j,k))
!                       print *, cresp_arguments(2*ncre+3)
!                       print *, 'cg%u(iarr_cre(e),34,34,:) =', cg%u(iarr_cre_e,34,34,0)
!                       print *, 'cg%u(iarr_cre(n),34,34,:) =', cg%u(iarr_cre_n,34,34,0)
!                       print *, 'p ', crel%p
! !                       print *, 'q ', crel%q
!                       print *, 'f ', crel%f
!                       print *, 'plo, pup = ', cg%u(flind%cre%plo,i,j,k),cg%u(flind%cre%pup,i,j,k)
!                       print *, '-------------------------'
                      
                endif

           enddo
         enddo
       enddo
       
      cgl=>cgl%nxt
      enddo

   end subroutine cresp_update_grid
   
   
   subroutine cresp_init_grid
   
      use cg_leaves,      only: leaves
      use cg_list,        only: cg_list_element
      use constants,      only: I_ONE, I_TWO, I_FOUR, fpi, LO, HI, xdim, ydim, zdim, ndims !, pMAX,
!       use domain,         only: dom!, is_multicg
!       use func,           only: ekin, emag, operator(.equals.), operator(.notequals.)
      use grid_cont,      only: grid_container
!       use cresp_variables, only: ind_p_lo, ind_p_up, cresp_taylor_order, taylor_coeff_2nd, taylor_coeff_3rd, &
!                                 ind_e_beg, ind_e_end, ind_n_beg, ind_n_end
      use initcrspectrum, only: ncre, f_init, p_up_init, p_lo_init, q_init, cre_eff, spec_mod_trms, p_lo_nda, p_up_nda
      use initcosmicrays, only: iarr_crn, iarr_cre
      use cresp_crspectrum, only: cresp_init_state
      use units,          only: clight
!       use crhelpers,      only: divv_n
!       use named_array_list, only: qna
      
      implicit none

      integer                         :: i, j, k !, icr, ipm, jpm, kpm
!       integer, dimension(ndims)       :: ldm, hdm
      type(cg_list_element),  pointer :: cgl
      type(grid_container),   pointer :: cg
      real(kind=8), allocatable, dimension(:)  :: cresp_arguments
      real(kind=8)                             :: dt_cre_tmp
      type(spec_mod_trms)  :: sptab
      
      allocate(cresp_arguments(I_ONE:I_TWO*ncre+I_TWO))
      !       logical, save :: frun = .true.
      !       integer       :: cr_id         ! maybe we should make this variable global in the module and do not pass it as an argument?
   i = 0
   j = 0
   k = 0
   
   dt_cre_tmp = 1.0
   dt_cre = dt_cre_tmp
   
     p_lo_nda = p_lo_init
     p_up_nda = p_up_init
   
   cgl => leaves%first
   do while (associated(cgl))
     cg => cgl%cg
     cresp_arguments = 0.0

           do k = cg%lhn(zdim,LO), cg%lhn(zdim,HI)
             do j = cg%lhn(ydim,LO), cg%lhn(ydim,HI)
               do i = cg%lhn(xdim,LO), cg%lhn(xdim,HI)

!                   cresp_arguments(I_ONE:2*ncre+2)    = cg%u(iarr_cre, i, j, k)
                  cresp_arguments(2*ncre+1)          = p_lo_nda(i,j,k)
                  cresp_arguments(2*ncre+2)          = p_up_nda(i,j,k)
!                   sptab%ub = emag(cg%b(xdim,i,j,k), cg%b(ydim,i,j,k),cg%b(zdim,i,j,k))
!                   sptab%ud = cg%q(qna%ind(divv_n))%point([i,j,k])
!                  call  append_dissipative_terms(i,j,k)
                  f_init = 1/(fpi*clight*(p_lo_init**(I_FOUR))*(((p_up_init/p_lo_init)**(I_FOUR-q_init))-I_ONE)/(I_FOUR-q_init))   !!! amplitude and distribution of electron energy density is inherited after those of nucleons, see crspectrum.pdf, eq. 29
                  f_init    = f_init*cg%u(iarr_crn(1),i,j,k)*cre_eff
                  
                  call cresp_init_state(cresp_arguments, sptab)
#ifdef VERBOSE
              print *, 'Output of cosmic ray electrons module for grid cell with coordinates i,j,k:', i, j, k
#endif /* VERBOSE */
                  cg%u(iarr_cre, i, j, k) = cresp_arguments(I_ONE:I_TWO*ncre+I_TWO)
                  p_lo_nda(i,j,k)         = cresp_arguments(I_TWO*ncre+I_ONE)
                  p_up_nda(i,j,k)         = cresp_arguments(I_TWO*ncre+I_TWO)
!                 if (i.eq.4.and.j.eq.34.and.k.eq.0) then
!                       print *, 'cg%u(iarr_cre(e),34,34,:) =', cg%u(iarr_cre_e,i,j,k)
!                       print *, 'cg%u(iarr_cre(n),34,34,:) =', cg%u(iarr_cre_n,i,j,k)
!                       print *, 'p ', crel%p
!                       print *, 'q ', crel%q
!                       print *, 'f ', crel%f
!                       print *,''
!                       print *, 'plo, pup = ', cg%u(flind%cre%plo,i,j,k),cg%u(flind%cre%pup,i,j,k)
!                       print *, '-------------------------'
!                 endif
                cresp_arguments = 0.0
           enddo
         enddo

       enddo
      cgl=>cgl%nxt
      
      enddo
      deallocate(cresp_arguments)
   
   end subroutine cresp_init_grid
   
! ------------------------------------------------

   subroutine grid_cresp_timestep
    use cg_leaves,        only: leaves
    use cg_list,          only: cg_list_element
!     use fluidtypes,       only: var_numbers
    use crhelpers,        only: divv_n
    use func,             only: emag !, operator(.equals.), operator(.notequals.)
    use grid_cont,        only: grid_container
    use constants,        only: I_ONE, I_TWO, I_FOUR, xdim, ydim, zdim, pi, LO, HI !, pMAX, 
    use named_array_list, only: qna
    use constants,        only: one
    use initcrspectrum,   only: spec_mod_trms
    use initcosmicrays,   only: iarr_cre_pl, iarr_cre_pu
    
    use timestep_cresp,  only: cresp_timestep_new

     implicit none
     integer                         :: i, j, k 
     type(grid_container), pointer   :: cg
     type(cg_list_element), pointer  :: cgl
!      type(var_numbers)    :: flind
     real   :: p_l, p_u
     real(kind=8)                             :: dt_cre_tmp
     type(spec_mod_trms)   :: sptab
          
     dt_cre = huge(one)
     dt_cre_tmp = huge(one)
     cgl => leaves%first
     do while (associated(cgl))
     cg => cgl%cg
     
           do k = cg%lhn(zdim,LO), cg%lhn(zdim,HI)
             do j = cg%lhn(ydim,LO), cg%lhn(ydim,HI)
               do i = cg%lhn(xdim,LO), cg%lhn(xdim,HI)
               p_l = cg%u(iarr_cre_pl, i, j, k)
               p_u = cg%u(iarr_cre_pu, i, j, k)
!                call append_dissipative_terms(i,j,k)
               sptab%ub = emag(cg%b(xdim,i,j,k), cg%b(ydim,i,j,k), cg%b(zdim,i,j,k))
               sptab%ud = cg%q(qna%ind(divv_n))%point([i,j,k])
               call cresp_timestep_new(dt_cre_tmp, p_l, p_u, sptab)
               dt_cre = min(dt_cre, dt_cre_tmp)
              enddo
           enddo
         enddo
         cgl=>cgl%nxt
     enddo
     dt_cre = 0.5 * dt_cre ! dt comes in to cresp_crspectrum with factor 2
!     print *, '@cresp_timestep_new: ', dt_cre
   
   end subroutine grid_cresp_timestep
   
   subroutine append_dissipative_terms(i,j,k)
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