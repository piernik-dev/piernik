module cresp_grid
! pulled by COSMIC_RAY_ELECTRONS

contains


 subroutine grid_cresp_update
      use constants,      only: xdim, ydim, zdim, LO, HI, pMAX
      use dataio_pub,        only: halfstep, warn, printinfo, msg
      use global,            only: dt
      use func,              only: operator(.notequals.), bx, by, bz
      use func,              only: emag
      use grid_cont,         only: grid_container
!       use cresp_variables,   only: cr_table, cren, cree, crepl, crepu, u_d, u_b
      use constants,         only: two
      
      implicit none
      
      integer(kind=4), intent(in)  :: i,j,k ! < loop iterators
!       logical, save :: frun = .true.
!       integer       :: cr_id         ! maybe we should make this variable global in the module and do not pass it as an argument?

      i,j,k    = 0

!       ts =  set_timer(tmr_mgd, .true.)
!       call all_dirty

!       if (diff_explicit .or. (allow_explicit .and. dt/diff_dt_crs_orig<1)) then
! 
!          if (frun) then
!             if (master .and. diff_explicit) call warn("[cresp_multigrid_update:multigrid_solve_diff] Multigrid was initialized but is not used")
!             frun = .false.
!          endif
      
      do k = cg%lhn(zdim,LO), cg%lhn(zdim,HI)
         do j = cg%lhn(ydim,LO), cg%lhn(ydim,HI)
            do i = cg%lhn(xdim,LO), cg%lhn(xdim,HI)
              
!             u_b = emag = pmag(bx(i,j,k), by(i,j,k), bz(i,j,k))
!             u_d = div_v
#ifdef VERBOSE
              print *, 'Output of cosmic ray electrons module for grid cell with coordinates i,j,k:', i, j, k
#endif /* VERBOSE */
!               call cell_cresp_update(two*dt,)
!             call cresp_crsupdate(2*cresp_dt, cg%u(cr_table(cren)), cg%u(cr_table(cree)), cg%u(cr_table(crepl), &
!             cg%u(cr_table(crepu), div_v ,emag(cg%b(i), cg%b(j), cg%b(k)) ! most likely 2 sweeps in one dt are executed in one step, we will test whether it's true or not
          
          enddo
        enddo
      enddo

      
   end subroutine grid_cresp_update
end module cresp_grid