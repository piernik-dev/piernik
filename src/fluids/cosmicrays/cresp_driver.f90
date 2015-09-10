program cresp_driver
! pulled by COSM_RAY_ELECTRONS

use types,      only: crel!,tosave
use vars!,       only: nbin, q_init, f_init
use cresp_pectrum, only: crsupdate, b_losses, init_state

implicit none
integer          , parameter      :: ione    = 1
real(kind=8)     , parameter      :: zero    = 0.0d0
real(kind=8)     , parameter      :: half    = 0.5d0
real(kind=8)     , parameter      :: sixth   = 1.6666d-1
real(kind=8)     , parameter      :: one     = 1.d0
!real(kind=8)     , parameter     :: two     = 2.d0
real(kind=8)     , parameter      :: three   = 3.d0
real(kind=8)     , parameter      :: four    = 4.d0
real(kind=8)     , parameter      :: five    = 5.d0
real(kind=8)     , parameter      :: ten     = 10.d0
!real(kind=8)     , parameter     :: twenty  = 20.d0
real(kind=8)     , parameter      :: hundred = 100.d0


logical                           :: first_run
real(kind = 8)                    :: t = zero
real(kind = 8)                    :: dt

integer                           :: i, j


 first_run = .true.
  
 open(10, file="crs.dat")
!   open(20, file="test.dat")
 
 ! ------------ coefficients for order in p_upw_rch and p_upw -- these are declared in vars
 c2nd = (mod(2,order) / 2 + mod(3,order))       ! coefficient which is always equal to 1 when order = 2 or = 3 and 0 if order = 1
 c3rd = (order - 1)*(order - 2) / 2             ! coefficient which is equal to 1 only when order = 3

 
 call init_state(dt)
 
 call timestep_cresp(dt)

 
 do while (t .lt. t_max)
 
   print*,  'Start a new iteration'
   print *, '----------------------------------'

! Compute new dt and update time
   call timestep_cresp(dt)
   t = t+dt
   
   print *, "    Time: ", t
   print *, '----------------------------------'
   
   u_d = u_d0 + div_v*cos(omega_d*t)
   
   call crsupdate(dt)
   
   print *, " "
   print '(A5, 11E18.9)', "t", t
   print '(A5, 11E18.9)', "p_fix", p_fix
   print '(A5, 11E18.9)', "p_act", crel%p
   print '(A5, 11E18.9)', "p_nex", p_next
   print '(A5, 11E18.9)', "p_upw", p_upw

   print *, " "
   print '(A5, 10E18.9)', "    n", n
   print '(A5, 11E18.9)', "nflux", nflux
   print '(A5, 10E18.9)', "  ndt", ndt
!   print '(A5, 10E18.9)', "   n1", n1

   print *, " "
   print '(A5, 10E18.9)', "    e", e
   print '(A5, 11E18.9)', "eflux", eflux
   print '(A5, 10E18.9)', "  edt", edt
!   print '(A5, 10E18.9)', "   e1", e1

   print *, " "
   print '(A5, 10E18.9)', "    r", r
   print '(A5, 10E18.9)', "    q", crel%q
   print '(A5, 11E18.9)', "    f", crel%f

   print *, " "

   call printer
!    call printer_2
   print *, " -------------------------- "

   print *,'coefficients: c2nd = ',c2nd, ', c3rd = ', c3rd
   print *
   print*, 'Accuracy test for adabatic compression/expansion:'
   print*, 'n_tot = ', n_tot, 'n_tot0 = ', n_tot0, 'rel error = ', (n_tot - n_tot0)/n_tot0
   print*, 'e_tot = ', e_tot, 'e_anal = ', e_tot0*exp(-u_d0*t), 'rel error = ', (e_tot-e_tot0*exp(-u_d0*t))/(e_tot0*exp(-u_d0*t))
   print*
   print*,  'End of iteration'
   print *, '=================================='
   print*
   print*
   print*
   print*
   print*,'--------------------'
   print *,''
    
        
    
 enddo
 
 print *,''
 
 close(10)
!   close(20)

!-------------------
!
contains
!
!-------------------
 
   subroutine timestep_cresp(dt)

      implicit none
      real(kind=8), intent(out)  :: dt
      real(kind=8), dimension(nbin) :: dts_new

      !!! beware !!! - reduce the scope to fixed bins
      
      dts_new = dt_ini + one ! Condition which assures that resulting dt will be lowest possible
      
      where (abs(b_losses(crel%p(1:nbin))) .ne. zero)
        dts_new =  (crel%p(1:nbin)-crel%p(0:nbin-1))/abs(b_losses(crel%p(1:nbin)))
      end where

      !print *,'dtsn = ', dts_new
      
      dt = cfl_cr* minval(dts_new)
      if (dt .ge. dt_ini) then 
        dt = dt_ini
      endif
        !dt = dt_ini
      
       print *, 'Timestep: ', dt
   end subroutine timestep
   
!--------------------

   subroutine printer
   implicit none
  
      write(10, '(e16.9, 3(1x,i8), 100(1x,e16.9))') t, nbin, crel%i_lo, crel%i_up, crel%p, crel%f, crel%q
   end subroutine printer
! -------------------- 
!    subroutine printer_2
!    implicit none
!    
!       write(10, '(e16.9, 3(1x,i8), 100(1x,e16.9))') t, nbin, crel
!    end subroutine printer_2
end program cresp_driver