module cresp_grid
! pulled by COSM_RAY_ELECTRONS

contains


 subroutine grid_cresp_update

      use constants,         only: xdim, ydim, zdim, zero, tmr_mgd
!       use crdiffusion,       only: cr_diff
      use dataio_pub,        only: halfstep, warn, printinfo, msg
      use fluidindex,        only: flind
      use global,            only: dt
      use func,              only: operator(.notequals.), bx, by, bz
      use mpisetup,          only: master
      use multigrid_helpers, only: all_dirty
      use multigridvars,     only: ts, tot_ts, stdout
      use timer,             only: set_timer
      use func,              only: emag
      use cresp_variables,   only: cr_table, cren, cree, crepl, crepu, u_d, u_b
      use constants,         only: two
      
      implicit none
      
!       real  :: cresp_dt, u_d1, u_b1
      integer(kind=4), intent(in)  :: i,j,k ! < loop iterators
!       integer(kind=4), intent(in)
      logical, save :: frun = .true.
      integer       :: cr_id         ! maybe we should make this variable global in the module and do not pass it as an argument?
      


      i,j,k    = 0
      
      cresp_dt = dt
      u_d_test = u_d
      u_b_test = u_b
      
!       ts =  set_timer(tmr_mgd, .true.)
!       call all_dirty

!       if (diff_explicit .or. (allow_explicit .and. dt/diff_dt_crs_orig<1)) then
! 
!          if (frun) then
!             if (master .and. diff_explicit) call warn("[cresp_multigrid_update:multigrid_solve_diff] Multigrid was initialized but is not used")
!             frun = .false.
!          endif
      
      
      do i = 1, xdim
        do j = 1, ydim
          do k = 1, zdim
            u_b = emag = pmag(bx(i,j,k), by(i,j,k), bz(i,j,k))
            u_d = div_v
#ifdef VERBOSE
            print *, 'Output of cosmic ray electrons module for grid cell with coordinates i,j,k:', i, j, k
#endif /* VERBOSE */
            call cell_cresp_update(two*dt, cg%u(flind%cre), u_d_test, u_b_test)
!             call cresp_crsupdate(2*cresp_dt, cg%u(cr_table(cren)), cg%u(cr_table(cree)), cg%u(cr_table(crepl), &
!             cg%u(cr_table(crepu), div_v ,emag(cg%b(i), cg%b(j), cg%b(k)) ! most likely 2 sweeps in one dt are executed in one step, we will test whether it's true or not
          
          enddo
        enddo
      enddo

      
   end subroutine grid_cresp_update
end module cresp_grid