module cresp_crspectrum
! pulled by COSM_RAY_ELECTRONS
 use global,only: dt !!! global timestep calculated by piernik
 use initcosmicrays, only: ncre
 use types, only: crel
 use vars!,  only: nbin, u_b, u_d, c2nd, c3rd, f_init, q_init

 implicit none
 
!   integer          , parameter      :: ione    = 1
!   real(kind=8)     , parameter      :: zero    = 0.0d0
!   real(kind=8)     , parameter      :: half    = 0.5d0
!   real(kind=8)     , parameter      :: sixth   = 1.6666d-1
!   real(kind=8)     , parameter      :: one     = 1.d0
! !real(kind=8)     , parameter     :: two     = 2.d0
!   real(kind=8)     , parameter      :: three   = 3.d0
!   real(kind=8)     , parameter      :: four    = 4.d0
!   real(kind=8)     , parameter      :: five    = 5.d0
!   real(kind=8)     , parameter      :: ten     = 10.d0
! !real(kind=8)     , parameter     :: twenty  = 20.d0
!   real(kind=8)     , parameter      :: hundred = 100.d0
  
  integer, dimension(0:ncre)        :: all_edges, i_act_edges
  integer, dimension(1:ncre)        :: all_bins
  integer                           :: del_i_lo, del_i_up
!   integer                           :: i_lo_act_edge, i_up_act_edge !,i


  
    ! logical arrays / arrays determining use of p/n/e in some lines
  logical, dimension(0:ncre)        :: is_fixed_edge,   is_fixed_edge_next 
  logical, dimension(0:ncre)        :: is_active_edge,  is_active_edge_next 
  logical, dimension(0:ncre)        :: is_cooling_edge, is_cooling_edge_next 
  logical, dimension(0:ncre)        :: is_heating_edge, is_heating_edge_next 
  logical, dimension(1:ncre)        :: is_active_bin,   is_active_bin_next
  
    ! counters
  integer                           :: num_fixed_edges,   num_fixed_edges_next
  integer                           :: num_active_edges,  num_active_edges_next
  integer                           :: num_active_bins,   num_active_bins_next
  integer                           :: num_cooling_edges, num_cooling_edges_next
  integer                           :: num_heating_edges, num_heating_edges_next
  
    ! dynamic arrays
  integer, allocatable              :: fixed_edges(:),   fixed_edges_next(:)
  integer, allocatable              :: active_edges(:),  active_edges_next(:)
  integer, allocatable              :: active_bins(:),   active_bins_next(:)
  integer, allocatable              :: cooling_edges(:), cooling_edges_next(:)
  integer, allocatable              :: heating_edges(:), heating_edges_next(:)
  
 
!-------------------------------------------------------------------------------------------------
!
contains
!
!-------------------------------------------------------------------------------------------------

!----- main subroutine -----

! values saved will be n and e in each cell; f and q have to be calculated

subroutine cresp_update(dt)
  implicit none
!    real(kind = 8)                    :: dt
   
! Update indexes of active bins, fixed edges and active edges at [t]
! Detect heating edges (energy upflow) and cooling edges (energy downflow)
   call cresp_update_bin_index(dt, p_lo, p_up, p_lo_next, p_up_next)    
   
! Compute power indexes for each bin at [t]
   call ne_to_q(n, e, crel%q)
   
! Compute f on left bin faces at [t]
   crel%f = nq_to_f(crel%p(0:ncre-1), crel%p(1:ncre), n(1:ncre), crel%q(1:ncre), active_bins)

! Compute fluxes through fixed edges in time period [t,t+dt]
! using f, q, p_lo and p_up at [t]
! Note that new [t+dt] values of p_lo and p_up 
! in case new fixed edges appear or disappear.
! fill new bins
  call cresp_compute_fluxes(cooling_edges_next,heating_edges_next)
  
! Computing e and n at [t+dt]
   ndt(1:ncre) = n(1:ncre)  - (nflux(1:ncre) - nflux(0:ncre-1))
   edt(1:ncre) = e(1:ncre)  - (eflux(1:ncre) - eflux(0:ncre-1))
   
!   edt(1:nbin) = e(1:nbin) *(one-0.5*dt*r(1:nbin)) - (eflux(1:nbin) - eflux(0:nbin-1))/(one+0.5*dt*r(1:nbin))   !!! oryginalnie u Miniatiego

! Compute coefficients R_i needed to find energy in [t,t+dt]
   call cresp_compute_r(p_next, active_bins_next)                 ! new active bins already received some particles, Ri is needed for those bins too
   
   edt(1:ncre) = edt(1:ncre) *(one-dt*r(1:ncre))

   p_lo = p_lo_next
   p_up = p_up_next

   n = ndt
   e = edt
   
   n_tot = sum(n)
   e_tot = sum(e)

   call cresp_deallocate_active_arrays
   
end subroutine cresp_update


!-------------------------------------------------------------------------------------------------

! all the procedures below are called by crsupdate subroutine or the driver

subroutine update_bin_index(dt,p_lo, p_up, p_lo_next, p_up_next)
   
      implicit none
      real(kind=8)              :: dt
      real(kind=8), intent(in)  :: p_lo, p_up
      real(kind=8), intent(out)  :: p_lo_next, p_up_next
      integer                    :: i_lo_next, i_up_next

! Compute p_lo and p_up at [t+dt]
      call p_update(dt, p_lo, p_lo_next)
      call p_update(dt, p_up, p_up_next)
    
! Locate cut-ofs before and after current timestep
      crel%i_lo = int(floor(dlog10(p_lo/p_fix(1))/w)) + 1
      crel%i_lo = max(0, crel%i_lo)
      crel%i_lo = min(crel%i_lo, ncre - 1)
      
      crel%i_up = int(floor(dlog10(p_up/p_fix(1))/w)) + 2
      crel%i_up = max(1,crel%i_up)
      crel%i_up = min(crel%i_up,ncre)
      
      i_lo_next = int(floor(dlog10(p_lo_next/p_fix(1))/w)) + 1
      i_lo_next = max(0, i_lo_next)
      i_lo_next = min(i_lo_next, ncre - 1)
      
      i_up_next = int(floor(dlog10(p_up_next/p_fix(1))/w)) + 2
      i_up_next = max(1,i_up_next)
      i_up_next = min(i_up_next,ncre)
      
! Detect changes in positions of lower an upper cut-ofs
      del_i_lo = i_lo_next - crel%i_lo
      del_i_up = i_up_next - crel%i_up
      
! Construct index arrays for fixed edges betwen p_lo and p_up, active edges 
! before timestep
      is_fixed_edge = .false.
      is_fixed_edge(crel%i_lo+1:crel%i_up-1) = .true.
      num_fixed_edges = count(is_fixed_edge)
      allocate(fixed_edges(num_fixed_edges))
      fixed_edges = pack(all_edges, is_fixed_edge)
      
      is_fixed_edge_next = .false.
      is_fixed_edge_next(i_lo_next+1:i_up_next-1) = .true.
      num_fixed_edges_next = count(is_fixed_edge_next)
      allocate(fixed_edges_next(num_fixed_edges_next))
      fixed_edges_next = pack(all_edges, is_fixed_edge_next)
      
      is_active_edge = .false.
      is_active_edge(crel%i_lo:crel%i_up) = .true.
      num_active_edges = count(is_active_edge)
      allocate(active_edges(num_active_edges))
      active_edges = pack(all_edges, is_active_edge)
      
      is_active_edge_next = .false.
      is_active_edge_next(i_lo_next:i_up_next) = .true.
      num_active_edges_next = count(is_active_edge_next)
      allocate(active_edges_next(num_active_edges_next))
      active_edges_next = pack(all_edges, is_active_edge_next)
      
! Active bins 
      is_active_bin = .false.
      is_active_bin(crel%i_lo+1:crel%i_up) = .true.
      num_active_bins = count(is_active_bin)
      allocate(active_bins(num_active_bins))
      active_bins = pack(all_bins, is_active_bin)
      
! Active bins after timestep
      is_active_bin_next = .false.
      is_active_bin_next(i_lo_next+1:i_up_next) = .true.
      num_active_bins_next = count(is_active_bin_next)
      allocate(active_bins_next(num_active_bins_next))
      active_bins_next = pack(all_bins, is_active_bin_next)
      
! Compute p for all active edges
      crel%p = zero
      crel%p(fixed_edges) = p_fix(fixed_edges)
      crel%p(crel%i_lo) = p_lo
      crel%p(crel%i_up) = p_up
      
      p_next = zero
      p_next(fixed_edges_next) = p_fix(fixed_edges_next)
      p_next(i_lo_next) = p_lo_next
      p_next(i_up_next) = p_up_next

! Compute upwind momentum p_upw for all fixed edges
      p_upw = zero
      p_upw(1:ncre) = p_fix(1:ncre)*(one+p_upw_rch(dt,p_fix(1:ncre)))
            

      print*, 'Change of  cut index lo,up:', del_i_lo, del_i_up
      
      
! Detect cooling edges and heating edges
      is_cooling_edge_next = .false.
      is_cooling_edge_next(fixed_edges_next)   = (p_upw(fixed_edges_next) &
                                                           > p_fix(fixed_edges_next))
      num_cooling_edges_next = count(is_cooling_edge_next)
      allocate(cooling_edges_next(num_cooling_edges_next))
      cooling_edges_next = pack(all_edges, is_cooling_edge_next)
      
      is_heating_edge_next = .false.
      is_heating_edge_next(fixed_edges_next) = (p_upw(fixed_edges_next) & 
                                                         < p_fix(fixed_edges_next))
      num_heating_edges_next = count(is_heating_edge_next)
      allocate(heating_edges_next(num_heating_edges_next))
      heating_edges_next = pack(all_edges, is_heating_edge_next)
      
      print *, 'In update_bin_index'
      print *, 'active edges: ', is_active_edge_next, active_edges_next
      print *, 'active bins:   ', is_active_bin_next ,'', active_bins_next
      print *, 'fixed  edges: ', is_fixed_edge_next,  fixed_edges_next
      print *, 'cooling edges:', is_cooling_edge_next,  cooling_edges_next
      print *, 'heating edges:', is_heating_edge_next,  heating_edges_next
      print *, 'Exit update_bin_index'

   end subroutine update_bin_index
   
!-------------------------------------------------------------------------------------------------
! 
! update p_range
!
!-------------------------------------------------------------------------------------------------

   subroutine p_update(dt, p_old,  p_new)
      
      implicit none
      real(kind=8), intent(in)  :: dt
      real(kind=8), intent(in)  :: p_old
      real(kind=8), intent(out) :: p_new
      
      p_new = p_old*(one + p_rch(dt, p_old)) ! changed from - to + for the sake of intuitiveness in p_rch subroutine
      
      
   end subroutine p_update

!-------------------------------------------------------------------------------------------------
! 
! arrays initialization
!
!-------------------------------------------------------------------------------------------------
   subroutine cresp_init_state(dt)
      implicit none
      real(kind = 8)    :: dt
      integer           :: i, k
      real(kind=8)      ::  c ! width of bin

      all_edges = (/ (i,i=0,ncre) /)
      all_bins = (/ (i,i=1,ncre) /)

      crel%q = q_init

      w  = (log10(p_max_fix/p_min_fix))/dble(ncre-2)
      
      p_fix(1:ncre-1)  =  p_min_fix*ten**(w*dble(all_edges(1:ncre-1)-1))
      p_fix(0)    = zero     ! initial array of p used to identify fixed edges
      p_fix(ncre) = zero     ! extreme two edges are never fixed edges

      crel%p        = p_fix       ! actual array of p including free edges
      crel%p(0)     = p_lo
      crel%p(ncre)  = p_up
      
! Sorting bin edges - arbitrary chosen p_lo and p_up may need to be sorted to appear in growing order
      do k = ncre, 1, -1
         do i = 0, k-1
            if (crel%p(i)>crel%p(i+1)) then 
               c = crel%p(i)
               crel%p(i) = crel%p(i+1)
               crel%p(i+1) = c
            endif
         enddo
      enddo
      
      crel%f = f_init * (crel%p/p_min_fix)**(-q_init)

      crel%i_lo = 0
      crel%i_up = ncre
      
      dt = 0.0d0
      u_d = 0.0d0
      
!      call timestep(dt) <- will be called in the driver instead

       call cresp_update_bin_index(dt, p_lo, p_up, p_lo_next, p_up_next)

      
      
     print*, 'In init_state'
     print *, num_active_edges, size(active_edges), lbound(active_edges), ubound(active_edges)
     print *, num_active_bins,  size(active_bins), lbound(active_bins), ubound(active_bins)
     print *, 'p      =', active_edges(1:num_active_edges-1), & 
                                                          crel%p(active_edges(1:num_active_edges-1))
     print *, 'p      =', active_edges(2:num_active_edges),   &
                                                          crel%p(active_edges(2:num_active_edges))
     print *, 'init f =', active_edges(1:num_active_edges-1), & 
                                                          crel%f(active_edges(1:num_active_edges-1))
     print *, 'init q =', active_bins(1:num_active_bins),     & 
                                                          crel%q(active_bins(1:num_active_bins))
      
      
       e = fq_to_e(crel%p(0:ncre-1), crel%p(1:ncre), crel%f(0:ncre-1), crel%q(1:ncre), active_bins)
!       
       n = fq_to_n(crel%p(0:ncre-1), crel%p(1:ncre), crel%f(0:ncre-1), crel%q(1:ncre), active_bins)


      print *, 'init n =', n
      print *, 'init e =', e 

       n_tot0 = sum(n)
       e_tot0 = sum(e)

       print *, 'n_tot0 =', n_tot0
       print *, 'e_tot0 =', e_tot0

   call cresp_deallocate_active_arrays
   
!   stop
   
    end subroutine cresp_init_state
  
!-------------------------------------------------------------------------------------------------
! 
! energy integral (eq. 21)
!
!-------------------------------------------------------------------------------------------------
 
   function fq_to_e(p_l, p_r, f_l, q, bins)
      implicit none
      real(kind=8), dimension(:)            :: p_l, p_r, f_l, q
      integer, dimension(:)                 :: bins
      real(kind=8), dimension(size(bins))   :: e_bins
      real(kind=8), dimension(1:ncre)       :: fq_to_e

      fq_to_e = 0.0d0
      e_bins = four*cnst_pi*cnst_c*f_l(bins)*p_l(bins)**4
      where(q(bins) .ne. four) 
         e_bins = e_bins*((p_r(bins)/p_l(bins))**(four-q(bins)) - one)/(four - q(bins))
      elsewhere
         e_bins = e_bins*log(p_r(bins)/p_l(bins))
      end where
      
      fq_to_e(bins) = e_bins
      
   end function
 
!-------------------------------------------------------------------------------------------------
 
 
!-------------------------------------------------------------------------------------------------
! 
! density integral (eq. 9)
!
!-------------------------------------------------------------------------------------------------

   function fq_to_n(p_l, p_r, f_l, q, bins)

      implicit none
      real(kind=8), dimension(:)            :: p_l, p_r, f_l, q
      integer, dimension(:)                 :: bins
      real(kind=8), dimension(size(bins))   :: n_bins
      real(kind=8), dimension(1:ncre)       :: fq_to_n
      
      n_bins = four*cnst_pi*f_l(bins)*p_l(bins)**3
      where(q(bins) .ne. three) 
         n_bins = n_bins*((p_r(bins)/p_l(bins))**(three-q(bins)) - one)/(three - q(bins))
      elsewhere
         n_bins = n_bins*log((p_r(bins)/p_l(bins)))
      end where
      
      fq_to_n = zero
      fq_to_n(bins) = n_bins
      
   end function

!-------------------------------------------------------------------------------------------------
   
   subroutine cresp_deallocate_active_arrays
      implicit none
      if(allocated(fixed_edges)) deallocate(fixed_edges)
      if(allocated(fixed_edges_next)) deallocate(fixed_edges_next)
      if(allocated(active_edges)) deallocate(active_edges)
      if(allocated(active_edges_next)) deallocate(active_edges_next)
      if(allocated(active_bins)) deallocate(active_bins)
      if(allocated(active_bins_next)) deallocate(active_bins_next)
      if(allocated(cooling_edges_next)) deallocate(cooling_edges_next)
      if(allocated(heating_edges_next)) deallocate(heating_edges_next)

   end subroutine cresp_deallocate_active_arrays

!-------------------------------------------------------------------------------------------------
! 
! compute fluxes
!
!-------------------------------------------------------------------------------------------------
   subroutine cresp_compute_fluxes(ce,he)
   
      implicit none
      integer, intent(in), dimension(:) :: ce, he    ! cooling edges, heating edges
      real(kind=8), dimension(1:ncre-1) :: pimh, pimth, fimh,fimth
      
      real(kind=8), dimension(1:ncre-1) :: dn_upw, de_upw, qi,qim1
      
      pimh(1:ncre-1) = crel%p(1:ncre-1)
      pimth(1:ncre-1) = crel%p(0:ncre-2)

      fimh(1:ncre-1) = crel%f(1:ncre-1)
      fimth(1:ncre-1) = crel%f(0:ncre-2)
      
      qi(1:ncre-1)  = crel%q(2:ncre)
      qim1(1:ncre-1) = crel%q(1:ncre-1)
      
      dn_upw = zero
      de_upw = zero
      nflux  = zero
      eflux  = zero
   
      dn_upw(ce) = four*cnst_pi*fimh(ce)*pimh(ce)**3
      where(qi(ce) .ne. three) 
         dn_upw(ce) = dn_upw(ce)*((p_upw(ce)/pimh(ce))**(three-qi(ce)) - one)/(three - qi(ce))
      elsewhere
         dn_upw(ce) = dn_upw(ce)*log((p_upw(ce)/pimh(ce)))
      end where
      nflux(ce) = - dn_upw(ce)
            
      de_upw(ce) = four*cnst_pi*cnst_c*fimh(ce)*pimh(ce)**4
      where(qi(ce) .ne. four) 
         de_upw(ce) = de_upw(ce)*((p_upw(ce)/pimh(ce))**(four-qi(ce)) - one)/(four - qi(ce))
      elsewhere
         de_upw(ce) = de_upw(ce)*log(p_upw(ce)/pimh(ce))
      end where
      eflux(ce) =  - de_upw(ce)
      
      if(del_i_up == -1) then 
         nflux(crel%i_up) = -n(crel%i_up)
         eflux(crel%i_up) = -e(crel%i_up  )
      endif
         
      dn_upw(he) = four*cnst_pi*fimth(he)*p_upw(he)**3*(pimth(he)/p_upw(he))**qim1(he)
      where(qim1(he) .ne. three) 
         dn_upw(he) = dn_upw(he)*((pimh(he)/p_upw(he))**(three-qim1(he)) - one)/(three - qim1(he))
      elsewhere
         dn_upw(he) = dn_upw(he)*log((pimh(he)/p_upw(he)))
      end where
      nflux(he) = dn_upw(he)
            
      de_upw(he) = four*cnst_pi*cnst_c*fimth(he)*p_upw(he)**4*(pimth(he)/p_upw(he))**qim1(he)
      where(qi(he) .ne. four) 
         de_upw(he) = de_upw(he)*((pimh(he)/p_upw(he))**(four-qim1(he)) - one)/(four - qim1(he))
      elsewhere
         de_upw(he) = de_upw(he)*log(pimh(he)/p_upw(he))
      end where
      eflux(he) = de_upw(he)
      
      if(del_i_lo == 1) then 
         nflux(crel%i_lo+1) = n(crel%i_lo+1)
         eflux(crel%i_lo+1) = e(crel%i_lo+1)
      endif
      
   end subroutine cresp_compute_fluxes
      
!-------------------------------------------------------------------------------------------------
! 
! compute R (eq. 25)
!
!-------------------------------------------------------------------------------------------------

   subroutine cresp_compute_r(p, bins)
      implicit none
      
      real(kind=8), dimension(0:ncre), intent(in) :: p
      integer, dimension(:), intent(in)     :: bins
      real(kind=8), dimension(size(bins)) :: r_num, r_den
      
      r = 0.0d0
      where(crel%q(bins) .ne. five) 
         r_num = (p(bins)**(five-crel%q(bins)) - p(bins-1)**(five-crel%q(bins)))/(five - crel%q(bins))
      elsewhere
         r_num = log(p(bins)/p(bins-1))
      end where
      
      where(crel%q(bins) .ne. four) 
         r_den = (p(bins)**(four-crel%q(bins)) - p(bins-1)**(four-crel%q(bins)))/(four - crel%q(bins))
      elsewhere
         r_den = log(p(bins)/p(bins-1))
      end where
      r(bins) = u_d + u_b * r_num/r_den   !!! all cooling effects will come here

   end subroutine cresp_compute_r
   

!-------------------------------------------------------------------------------------------------
! 
! find new q (eq. 29) and new f
!
!-------------------------------------------------------------------------------------------------

subroutine ne_to_q(n, e, q)

   implicit none
   real(kind=8), dimension(ncre), intent(in)  :: n, e
   real(kind=8), dimension(ncre), intent(out) :: q
   integer :: i, j
   real(kind=8)    :: alpha, base, x, dx, delta, q_err
   real(kind=8)    :: dfun1

   dx    = 1.0d-3
   q_err = 1.0d-12
  
   do i = 1, ncre
!      if(is_active_bin_next(i)) then 
!         alpha = e(i)/(n(i)*p_next(i-1)*cnst_c)     ! czy tu ma byc p_next ?
!         base  = p_next(i)/p_next(i-1)                   !

      if(is_active_bin(i).and.n(i).ne.0.and.crel%p(i-1).ne.0.)  then  ! dziala takze z jednym dodatkowym warunkiem, niezaleznie czy jest to n czy p, ktores musi byc wieksze od 0
         alpha = e(i)/(n(i)*crel%p(i-1)*cnst_c)     ! czy tu ma byc p_next ?
         base  = crel%p(i)/crel%p(i-1)                   !
         x     = q_init !q(i)
         do j = 1, 100
            dfun1 = 0.5*(fun1(x+dx, alpha, base)-fun1(x-dx, alpha, base))/dx
     
            if(x .gt. q_big) then 
               q(i) = q_big ! (four*alpha-three)/(alpha-one)
               exit
            endif
            if(x .le. -q_big) then 
                  q(i) = -q_big ! (four*alpha-three)/(alpha-one)
               exit
            endif

            delta = fun1(x, alpha, base)/dfun1
  
            if(abs(delta) .lt. q_err) then
               q(i) = x
!                  q(i) = min(q(i), twenty)
               exit
            else
               x = x - delta
            endif
         enddo
      else
         q(i) = 0.0d0!zero! hundred ! when i-bin is empty 
      endif
!            print *, i, is_active_bin(i), q(i)
   enddo
      
   end subroutine ne_to_q
   
!-------------------------------------------------------------------------------------------------   
   
   function fun1(x, alpha, base)
      implicit none
      real(kind=8) :: x, alpha, base, fun1
      if (x .eq. three) then
         fun1 = -alpha + (-one + base)/log(base) 
      else if (x .eq. four) then
         fun1 = -alpha + base*log(base)/(base - one)
      else
         fun1 = -alpha + ((three-x)/(four-x))*((base**(four-x)-one)/(base**(three-x)-one))
      endif
   end function
    
!-------------------------------------------------------------------------------------------------
! 
! distribution function amplitudes (eq. 9)
!
!-------------------------------------------------------------------------------------------------

   function nq_to_f(p_l, p_r, n, q, bins)

      implicit none
      real(kind=8), dimension(:)            :: p_l, p_r, n, q
      integer, dimension(:)                 :: bins
      real(kind=8), dimension(size(bins))   :: f_bins
      real(kind=8), dimension(0:ncre)       :: nq_to_f

      f_bins = n(bins) / (four*cnst_pi*p_l(bins)**3)
      where(q(bins) .ne. three) 
         f_bins = f_bins*(three - q(bins)) /((p_r(bins)/p_l(bins))**(three-q(bins)) - one)
      elsewhere
         f_bins = f_bins/log((p_r(bins)/p_l(bins)))
      end where

      nq_to_f= zero
      nq_to_f(bins-1) = f_bins


   end function nq_to_f  
   
!-------------------------------------------------------------------------------------------------


  function b_losses(p)
    implicit none
    real(kind=8), dimension(:) :: p
    real(kind=8), dimension(size(p)) :: b_losses
   
    b_losses = u_b*p**2  !!! b_sync_ic = 8.94e-25*(u_b+u_cmb)*gamma_l**2 ! erg/cm

  end function

!-------------------------------------------------------------------------------------------------
! 
! relative change of momentum due to losses (u_b*p*dt) and compression u_d*dt (Taylor expansion up to 3rd order)
!
!-------------------------------------------------------------------------------------------------

  function p_rch(dt, p)
    implicit none
    real(kind=8), intent(in)  :: dt
    real(kind=8), intent(in)  :: p
    real(kind=8)              :: p_rch

   !p_rch =  (u_b*p + u_d) * dt !!! b_sync_ic = 8.94e-25*(u_b+u_cmb)*gamma_l**2 ! erg/cm
    p_rch = (- u_d - p * u_b ) *  dt  + c2nd *( half*u_d**2 + u_b**2 * p**2)*dt**2 + c3rd *(-sixth*u_d**3 - u_b**3 * p**3)*dt**3 ! analitycally correct solution
   

  end function

!-------------------------------------------------------------------------------------------------

  function p_upw_rch(dt, p)
    implicit none
    real(kind=8), intent(in)  :: dt
    real(kind=8), dimension(:)        :: p
    real(kind=8), dimension(size(p))  :: p_upw_rch
   
   !p_upw_rch =  (u_b*p + u_d) * dt      !!! b_sync_ic = 8.94e-25*(u_b+u_cmb)*gamma_l**2 ! erg/cm
    p_upw_rch = (u_d + p * u_b) * dt + c2nd * (half*u_d**2 + u_b**2 * p**2)*dt**2 + c3rd *(sixth*u_d**3 + u_b**3 * p**3)*dt**3 ! analitycally correct
   

  end function
  
end module cresp_crspectrum