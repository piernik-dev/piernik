module timestep_cresp
! pulled by COSM_RAY_ELECTRONS

 use initcosmicrays, only: ncre
 use constants,      only: one, zero
!  use global,         only: dt, t
      implicit none
      
      public          dt_cre      
           
      private         p_tmp
 
      real(kind=8) :: dt_cre
      real(kind=8), allocatable, dimension(:) :: p_tmp   ! do not use p variable in this module anywhere outside
      real(kind=8), allocatable, dimension(:) :: dts_new
      real(kind=8) :: u_b

!-------------------      
!      
contains 
!
!-------------------
 subroutine timestep_cresp_update
!      use constansts,     only: I_ZERO, I_ONE, zero, one
     use cg_leaves,      only: leaves
     use cg_list,        only: cg_list_element
     use grid_cont,      only: grid_container
     use func,           only: ekin, emag, operator(.equals.), operator(.notequals.)
     use fluidtypes,     only: var_numbers
     use fluidindex,     only: flind
     use constants,      only: xdim, zdim, ydim, one, zero
     use dataio_pub,     only: msg
     
     use initcosmicrays, only: ncre, p_min_fix, p_max_fix, p_fix_init
     use cresp_arrays_handling, only: allocate_with_index
     implicit none
     
      type(cg_list_element),  pointer         :: cgl
      type(grid_container),   pointer         :: cg
      real(kind=8)                            :: p_l, p_u, c !, p_l_n, p_u_n
      real(kind=8)                            :: dt_tmp 
      integer        :: id,jd,kd, i, k            ! id;kd - domain indices, i,k - sorting indices
  
      call allocate_with_index(p_tmp, 0, ncre)
      call allocate_with_index(dts_new, 1, ncre)

  u_b        = zero
  dt_tmp = huge(one)
  dt_cre = huge(one)

  cgl => leaves%first
  do while (associated(cgl))
     cg => cgl%cg
        do kd = cg%ks, cg%ke
           do jd = cg%js, cg%je
              do id = cg%is, cg%ie
              
               p_l   = cg%u(flind%cre%plo, id, jd, kd)
               p_u   = cg%u(flind%cre%pup, id, jd, kd)
               u_b   = emag(cg%b(xdim,id,jd,kd), cg%b(ydim,id,jd,kd), cg%b(zdim,id,jd,kd))/cg%dvol

               p_tmp = zero                        ! important to have p_tmp = zero after each iteration
               p_tmp = p_fix_init                  ! initial (and constant) values of momentum grid are passed down to momentum array for the current cell
               
               where (p_tmp .gt. p_l)
                  p_tmp          = p_fix_init      ! actual array of p including free edges
               elsewhere
                  p_tmp          = zero
               end where
               
               p_tmp(0) = p_u ! cg%u(flind%cre%plo, id, jd, kd)
               p_tmp(ncre) = p_l  ! cg%u(flind%cre%pup, id, jd, kd)
                                            
               ! Sorting bin edges - arbitrary chosen p_lo and p_up may need to be sorted to appear in growing order
               ! Chunk of code similar to that in cresp_init_state.
               ! To be swapped with subroutine call in the future.
               
               do k = ncre, 1, -1
                    do i = 0, k-1
                      if (p_tmp(i)>p_tmp(i+1)) then 
                       c = p_tmp(i)
                       p_tmp(i) = p_tmp(i+1)
                       p_tmp(i+1) = c
                      endif
                    enddo
               enddo
               
               where (p_tmp .gt. p_u)
                  p_tmp          = zero
               end where

               call cresp_timestep(dt_tmp)  
               
               dt_cre = min(dt_cre, dt_tmp) ! After all iterations gives minimal timestep in whole domain
                  
              enddo
           enddo
        enddo
        cgl=>cgl%nxt
    enddo
        dt_cre = dt_cre*0.5  ! cresp_crspectrum is called with 2*dt; this must be done to compensate it.
!     write(msg,*) '[timestep_cresp:timestep_cresp_update]:                        dt_cre = ', dt_cre

    if (allocated(p_tmp))   deallocate(p_tmp)
    if (allocated(dts_new)) deallocate(dts_new)
!     if (allocated(all_edges)) deallocate(all_edges)
!     if (allocated(p_fix)) deallocate(p_fix)
    
 end subroutine timestep_cresp_update

!----------------------------------------------------------
! remove from cresp_crspectrum when this module starts to work!
 
   subroutine cresp_timestep(dt_calc)
    use initcosmicrays,    only: cfl_cre
    implicit none
    real(kind=8)                  :: dt_calc
    real(kind=8)                  :: dts_min
!     real(kind=8), parameter       :: min_p_diff = 5e-2
      
      dt_calc = huge(one)
      dts_new = huge(one) ! whole dts_new array
      where (abs(b_losses(p_tmp(1:ncre))) .ne. zero)
        dts_new =  (p_tmp(1:ncre)-p_tmp(0:ncre-1))/abs(b_losses(p_tmp(1:ncre)))
        where ((p_tmp(1:ncre)-p_tmp(0:ncre-1)).eq.zero)   !!!
          dts_new = huge(one) 
        end where!!!
!         where  (p_tmp(1:ncre)-p_tmp(0:ncre-1) .le. min_p_diff)
!           dts_new = huge(one)
!         end where
      end where
      
      dts_min = cfl_cre*minval(dts_new)   ! min of array
!        print *, 'p = ', p
!        print *, 'ub  = ', u_b
!        print *,'dts_min = ', dts_min
!        print *,'b_losses(p(1:ncre) = ', b_losses(p(1:ncre))
      
      if (dt_calc.ge.dts_min) then   ! gives minimal timestep in bin space
         dt_calc = dts_min
      endif
      
      dts_new = zero
      
   end subroutine cresp_timestep

!----------------------------------------------------------
! remove from cresp_crspectrum when this module starts to work!
   
  function b_losses(p)
    implicit none
!     real(kind=8), intent(in)                :: u_b
    real(kind=8), dimension(:), intent(in)  :: p
    real(kind=8), dimension(size(p)) :: b_losses
   
    b_losses = u_b*p**2  !!! b_sync_ic = 8.94e-25*(u_b+u_cmb)*gamma_l**2 ! erg/cm

  end function

end module timestep_cresp
