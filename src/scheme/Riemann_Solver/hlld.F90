!
! PIERNIK Code Copyright (C) 2006 Michal Hanasz
!
!    This file is part of PIERNIK code.
!
!    PIERNIK is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    PIERNIK is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with PIERNIK.  If not, see <http://www.gnu.org/licenses/>.
!
!    Initial implementation of PIERNIK code was based on TVD split MHD code by
!    Ue-Li Pen
!        see: Pen, Arras & Wong (2003) for algorithm and
!             http://www.cita.utoronto.ca/~pen/MHD
!             for original source code "mhd.f90"
!
!    For full list of developers see $PIERNIK_HOME/license/pdt.txt
!
!  References:
!
!  A multi-state HLLD approximate Riemann solver for ideal magnetohydrodynamics.
!  Takahiro Miyoshi, Kanya Kusano
!  Journal of Computational Physics 208 (2005) 315-344
!
!  ->Solve one dimensional Riemann problem using adiabatic HLLD scheme
!
!  Varadarajan Parthasarathy, CAMK, Warszawa. 2015.
!
!---------------------------------------------------------------------------------------------------------------------------


#include "piernik.def"
!>
!!  \brief This module implements HLLD Riemann solver following the work of Miyoshi & Kusano (2005)
!<
module hlld
! pulled by RIEMANN

  implicit none

  private
  public :: riemann_hlld

contains

  function fluxes(u,b,bb,cs2, cdim) result(f)

    use constants,  only: half, zero, xdim, ydim, zdim, idn, imx, imy, imz, ien
    use fluidindex, only: flind,iarr_mag_swp
    use fluidtypes, only: component_fluid
    use func,       only: ekin
    use dataio_pub, only: die

    implicit none

    real, dimension(:,:),      intent(in)    :: u
    real, dimension(:,:),      intent(in)    :: bb
    real, dimension(:,:),      intent(out)   :: b
    real, dimension(:), pointer, intent(in)  :: cs2
    integer(kind=4),             intent(in)  :: cdim

    real, dimension(size(u,1), size(u,2))    :: f
    real, dimension(size(u,2))               :: vx, vy, vz, p_t, p
    integer :: ip
    class(component_fluid),    pointer       :: fl
    integer(kind=4)                          :: ibx, iby, ibz

    ibx = iarr_mag_swp(cdim,xdim)
    iby = iarr_mag_swp(cdim,ydim)
    ibz = iarr_mag_swp(cdim,zdim)

    do ip = 1, flind%fluids

       fl => flind%all_fluids(ip)%fl

       vx  =  u(fl%imx,:)/u(fl%idn,:)
       vy  =  u(fl%imy,:)/u(fl%idn,:)
       vz  =  u(fl%imz,:)/u(fl%idn,:)

       if(fl%has_energy) then
          p_t = fl%gam_1*(u(fl%ien,:) - ekin(u(fl%imx,:), u(fl%imy,:), u(fl%imz,:), u(fl%idn,:)) - half*(bb(ibx,:)**2 + bb(iby,:)**2 + bb(ibz,:)**2)) + &
                                                                                                           half*(bb(ibx,:)**2 + bb(iby,:)**2 + bb(ibz,:)**2)
       else
          if(associated(cs2)) then
             p = cs2*u(fl%idn,:)
          else
             p = 0.
          endif
       endif

       if(cdim .eq. xdim) then
          f(fl%idn,:)  =  u(fl%imx,:)
          f(fl%imx,:)  =  u(fl%imx,:)*vx(:) + p_t(:) - bb(ibx,:)**2
          f(fl%imy,:)  =  u(fl%imy,:)*vx(:) - bb(ibx,:)*bb(iby,:)
          f(fl%imz,:)  =  u(fl%imz,:)*vx(:) - bb(ibx,:)*bb(ibz,:)
          b(ibx,:)     =  zero
          b(iby,:)     =  bb(iby,:)*vx(:) - bb(ibx,:)*vy(:)
          b(ibz,:)     =  bb(ibz,:)*vx(:) - bb(ibx,:)*vz(:)
          if(fl%has_energy) then
             f(fl%ien,:)  =  (u(fl%ien,:) + p_t(:))*vx(:) - bb(ibx,:)*(bb(ibx,:)*vx(:) + bb(iby,:)*vy(:) + bb(ibz,:)*vz(:))
          endif

       else if(cdim .eq. ydim) then
          f(fl%idn,:)  =  u(fl%imy,:)
          f(fl%imx,:)  =  u(fl%imx,:)*vy(:) - bb(ibx,:)*bb(iby,:)
          f(fl%imy,:)  =  u(fl%imy,:)*vy(:) + p_t(:) - bb(iby,:)**2
          f(fl%imz,:)  =  u(fl%imz,:)*vy(:) - bb(ibz,:)*bb(iby,:)
          b(ibx,:)     =  bb(ibx,:)*vy(:) - bb(iby,:)*vx(:)
          b(iby,:)     =  zero
          b(ibz,:)     =  bb(ibz,:)*vy(:) - bb(iby,:)*vz(:)
          if(fl%has_energy) then
             f(fl%ien,:)  =  (u(fl%ien,:) + p_t(:))*vy(:) -  bb(iby,:)*(bb(ibx,:)*vx(:) + bb(iby,:)*vy(:) + bb(ibz,:)*vz(:))
          endif

       else if(cdim .eq. zdim) then
          f(fl%idn,:)  =  u(fl%imz,:)
          f(fl%imx,:)  =  u(fl%imx,:)*vz(:) - bb(ibx,:)*bb(ibz,:)
          f(fl%imy,:)  =  u(fl%imy,:)*vz(:) - bb(iby,:)*bb(ibz,:)
          f(fl%imz,:)  =  u(fl%imz,:)*vz(:) + p_t(:) - bb(ibz,:)**2
          b(ibx,:)     =  bb(ibx,:)*vz(:) - bb(ibz,:)*vx(:)
          b(iby,:)     =  bb(iby,:)*vz(:) - bb(ibz,:)*vy(:)
          b(ibz,:)     =  zero
          if(fl%has_energy) then
             f(fl%ien,:)  =  (u(fl%ien,:) + p_t(:))*vy(:) - bb(ibz,:)*(bb(ibx,:)*vx(:) + bb(iby,:)*vy(:) + bb(ibz,:)*vz(:))
          endif

       else
          call die("Check the fluxes")
       endif
    enddo

    return

  end function fluxes

  ! Abstraction needed in terms in normal and tangential components that will reduce the if-blocks.

!-------------------------------------------------------------------------------------------------------------------------------------------------

  !subroutine riemann_hlld(n, gamma, uleft, uright, b, cdim, f)
  subroutine riemann_hlld(n, f, b, gamma, cdim)

    ! external procedures
    
    use constants,  only: half, zero, xdim, ydim, zdim, idn, imx, imy, imz, ien
    use fluidindex, only: iarr_mag_swp, flind
    use func,       only: operator(.notequals.)
    !use func,       only: emag, ekin
    !use grid_cont,  only: grid_container
    !use fluxes,     only: all_fluxes, flimiter
    use dataio_pub, only: die

    ! arguments
    
    implicit none

    integer,                       intent(in)    :: n
    real, dimension(:,:),          intent(out)   :: f
    real, dimension(:),            intent(in)    :: b
    real,                          intent(in)    :: gamma
    integer(kind=4),               intent(in)    :: cdim

    ! Local variables

    integer                                      :: i
    real, parameter                              :: four = 4.0
    real, parameter                              :: one  = 1.0   
    integer(kind=4)                              :: ibx, iby, ibz
    
    real                                         :: sm, sm_nr, sm_dr, sl, sr
    real                                         :: alfven_l, alfven_r, c_fastl, c_fastr, gampr_l, gampr_r
    real                                         :: slsm, srsm, slvxl, srvxr, smvxl, smvxr, srmsl, srtsl, dn_l, dn_r
    real                                         :: b_lr, b_lrgam, magprl, magprr, prtl, prtr, prt_star, b_sig
    real                                         :: coeff_1, coeff_2, dn_lsqt, dn_rsqt, add_dnsq, mul_dnsq
    real                                         :: vb_l, vb_starl, vb_r, vb_starr, vb_2star
    
    ! Local arrays

    real, dimension(n, flind%all)                :: ul, ur, fl, fr
    real, dimension(flind%all)                   :: u_starl, u_starr, u_2star
    real, dimension(flind%all)                   :: v_starl, v_starr, b_starl, b_starr, v_2star, b_2star
    
    ibx = iarr_mag_swp(cdim,xdim)
    iby = iarr_mag_swp(cdim,ydim)
    ibz = iarr_mag_swp(cdim,zdim)

    ! SOLVER

    ! Copy normal components of magnetic field to the left and right states

    ul(ibx,:) =  b(:)
    ur(ibx,:) =  b(:)

    do i = 1,n

       gampr_l  =  gamma*ul(ien,i)
       gampr_r  =  gamma*ur(ien,i)

       c_fastl  =   (gampr_l+(ul(ibx,i)**2+ul(iby,i)**2+ul(ibz,i)**2))  &
            + sqrt((gampr_l+(ul(ibx,i)**2+ul(iby,i)**2+ul(ibz,i)**2))**2-(four*gampr_l*ul(ibx,i)**2))

       c_fastl  =  sqrt(half*c_fastl/ul(idn,i))

       c_fastr  =   (gampr_r+(ur(ibx,i)**2+ur(iby,i)**2+ur(ibz,i)**2))  &
            + sqrt((gampr_r+(ur(ibx,i)**2+ur(iby,i)**2+ur(ibz,i)**2))**2-(four*gampr_r*ur(ibx,i)**2))

       c_fastr  =  sqrt(half*c_fastr/ur(idn,i))

       ! Eq. (67)
    
       sl  =  min(ul(imx,i), ur(imx,i)) - max(c_fastl,c_fastr)
       sr  =  max(ur(imx,i), ur(imx,i)) + max(c_fastl,c_fastr)

       ! Speed of contact discontinuity Eq. (38)

       sm_nr  =  (sr*ur(imx,i) - sl*ul(imx,i)) - (fr(imx,i) - fl(imx,i))
       sm_dr  =  (sr*ur(idn,i) - sl*ul(idn,i)) - (fr(idn,i) - fl(idn,i))
       sm     =  sm_nr/sm_dr

       ! Speed differences

       slsm  =  sl - sm
       srsm  =  sr - sm

       slvxl  =  sl - ul(imx,i)
       srvxr  =  sr - ur(imx,i)

       smvxl  =  sm - ul(imx,i)
       smvxr  =  sm - ur(imx,i)

       srmsl  =  sr - sl

       srtsl  =  sr*sl
    
       ! Co-efficients

       dn_l     =  ul(idn,i)*slvxl
       dn_r     =  ur(idn,i)*srvxr
       b_lr     =  ul(ibx,i)*ur(ibx,i)
       b_lrgam  =  b_lr/gamma
       
       ! Magnetic pressure

       magprl  =  half*sum(ul(ibx:ibz,i)*ul(ibx:ibz,i))
       magprr  =  half*sum(ur(ibx:ibz,i)*ur(ibx:ibz,i))

       ! Total pressure

       prtl  =  ul(ien,i) + magprl
       prtr  =  ur(ien,i) + magprr

       ! Pressure of intermediate state Eq. (23)

       prt_star  =  (prtl*dn_l*smvxl) + (prtr*dn_r*smvxr)  !< Check for 0.5.

       

       ! Normal components of velocity and magnetic field

       
       v_starl(imx)  =  sm
       v_starr(imx)  =  sm
    
       b_starl(ibx)  =  ul(ibx,i)
       b_starr(ibx)  =  ur(ibx,i)

       ! HLLD fluxes

       if (sl .ge.  zero) then
          f(:,i)  =  fl(:,i)       !< F_L  

       else if (sr .le.  zero) then
          f(:,i)  =  fr(:,i)       !< F_R

       else

          ! Transversal components of magnetic field for left states (Eq. 45 & 47), taking degeneracy into account
          
          coeff_1  =  dn_l*slsm - b_lr
          
          if ((coeff_1 .notequals. zero) .and. b_lrgam .le. ul(ien,i)) then
             coeff_2  =  (dn_l*slvxl - b_lr)/coeff_1

             b_starl(iby)   =  ul(iby,i)*coeff_2
             b_starl(ibz)  =  ul(ibz,i)*coeff_2

          else

             ! Calculate HLL left states
             
             b_starl(iby)  =  ((sr*ur(iby,i) - sl*ul(iby,i)) - (fr(iby,i) - fl(iby,i)))/srmsl
             b_starl(ibz)  =  ((sr*ur(ibz,i) - sl*ul(ibz,i)) - (fr(ibz,i) - fl(ibz,i)))/srmsl
             
          endif

          ! Transveral components of velocity Eq. 42

          coeff_1  =  ul(ibx,i)/dn_l
          
          v_starl(imy)  =  ul(imy,i) + coeff_1*(ul(iby,i) - u_starl(iby))
          v_starl(imz)  =   ul(imz,i) + coeff_1*(ul(ibz,i) - u_starl(ibz))

          ! Transversal components of magnetic field for right states (Eq. 45 & 47), taking degeneracy into account

          coeff_1  =  dn_r*srsm - b_lr

          if ((coeff_1 .notequals. zero) .and. b_lrgam .le. ur(ien,i)) then
             coeff_2  =  (dn_r*srvxr - b_lr)/coeff_1

             b_starr(iby)  =  ur(iby,i)*coeff_2
             b_starr(ibz)  =  ur(ibz,i)*coeff_2

          else

             ! Calculate HLL right states

             b_starr(iby)  =  ((sr*ur(iby,i) - sl*ul(iby,i)) - (fr(iby,i) - fl(iby,i)))/srmsl
             b_starr(ibz)  =  ((sr*ur(ibz,i) - sl*ul(ibz,i)) - (fr(ibz,i) - fl(ibz,i)))/srmsl

          endif

          ! Transversal components of velocity Eq. 42

          coeff_1  =  ur(ibx,i)/dn_r

          v_starr(imy)  =   ur(imy,i) + coeff_1*(ur(iby,i) - u_starr(iby))
          !v_starr(imz)  =  ur(imy,i) + coeff_1*(ur(iby,i) - u_starr(iby))  ! CHECK: ur(imz,i) + coeff_1*(ur(ibz,i) - u_starr(ibz))
          v_starr(imz)  =  ur(imy,i) + coeff_1*(ur(ibz,i) - u_starr(ibz))

          ! Dot product of velocity and magnetic field

          vb_l  =  sum(ul(imx:imz,i)*ul(ibx:ibz,i))
          vb_r  =  sum(ur(imx:imz,i)*ur(ibx:ibz,i))
          vb_starl  =  sum(v_starl(imx:imz)*b_starl(ibx:ibz))
          vb_starr  =  sum(v_starr(imx:imz)*b_starr(ibx:ibz))
          
          
          ! Left intermediate state conservative form
          
          u_starl(idn)  =  dn_l/slsm   ! Densities for right intermediate state Eq. (43)
          u_starl(imx)  =  u_starl(idn)*v_starl(imx)
          u_starl(imy)  =  u_starl(idn)*v_starl(imy)
          u_starl(imz)  =  u_starl(idn)*v_starl(imz)
          u_starl(ibx)  =  b_starl(ibx)
          u_starl(iby)  =  b_starl(iby)
          u_starl(ibz)  =  b_starl(ibz)

          ! Right intermediate state conservative form

          u_starr(idn)  =  dn_r/srsm ! Density for right intermediate state Eq. (43)
          u_starr(imx)  =  u_starr(idn)*v_starr(imx)
          u_starr(imy)  =  u_starr(idn)*v_starr(imy)
          u_starr(imz)  =  u_starr(idn)*v_starr(imz)
          u_starr(ibx)  =  b_starr(ibx)
          u_starr(iby)  =  b_starr(iby)
          u_starr(ibz)  =  b_starr(ibz)

          

          ! Total energy of left and right intermediate states Eq. (48)

          u_starl(ien)  =  (slvxl*ul(ien,i) - prtl*ul(imx,i) + prt_star*sm + ul(ibx,i)*(vb_l - vb_starl))/slsm
          u_starr(ien)  =  (srvxr*ur(ien,i) - prtr*ur(imx,i) + prt_star*sm + ur(ibx,i)*(vb_r - vb_starr))/srsm

          ! Cases for B_x .ne. and .eq. zero

          if(abs(ul(ibx,i)) > zero) then
                 
             ! Left and right Alfven waves velocity Eq. 51

             dn_lsqt  =  sqrt(u_starl(idn))
             dn_rsqt  =  sqrt(u_starr(idn))

             alfven_l  =  sm - abs(ul(ibx,i))/dn_lsqt
             alfven_r  =  sm + abs(ur(ibx,i))/dn_rsqt

             ! Intermediate discontinuities
             
             if(alfven_l > zero) then

                ! Left intermediate flux Eq. 64

                f(:,i)  =  fl(:,i) + sl*(u_starl(:) - ul(:,i))

             else if (alfven_r < zero) then

                ! Right intermediate flux Eq. 64

                f(:,i)  =  fr(:,i) + sr*(u_starr(:) - ur(:,i))

             else ! alfven_l .le. zero .le. alfven_r

                ! Arragnge for sign of normal component of magnetic field
                
                if(ul(ibx,i) .ge. zero) then

                   b_sig = one

                else

                   b_sig = -one

                endif

                ! Sum and product of density square-root

                add_dnsq  =  dn_lsqt + dn_rsqt
                mul_dnsq  =  dn_lsqt*dn_rsqt

                ! Components of velocity Eq. 39, 59, 60 and magnetic field Eq. 61, 62

                v_2star(imx)  =  sm
                v_2star(imy)  =  ((dn_lsqt*v_starl(imy) + dn_rsqt*v_starr(imy)) + b_sig*(b_starr(iby) - b_starl(iby)))/add_dnsq
                v_2star(imz)  =  ((dn_lsqt*v_starl(imz) + dn_rsqt*v_starr(imz)) + b_sig*(b_starr(ibz) - b_starl(ibz)))/add_dnsq

                b_2star(ibx)  =  ul(ibx,i)
                b_2star(iby)  =  ((dn_lsqt*b_starr(iby) + dn_rsqt*b_starl(iby)) + b_sig*mul_dnsq*(v_starr(imy) - v_starl(imy)))/add_dnsq
                !b_2star(ibz)  =  ((dn_lsqt*b_starr(ibz) + dn_rsqt*b_starl(iby)) + b_sig*mul_dnsq*(v_starr(imz) - v_starl(imy)))/add_dnsq  ! CHECK: dn_rsqt*b_starl(ibz),                                                                                                                                                v_starl(imz)
                b_2star(ibz)  =  ((dn_lsqt*b_starr(ibz) + dn_rsqt*b_starl(ibz)) + b_sig*mul_dnsq*(v_starr(imz) - v_starl(imz)))/add_dnsq
                
                ! Dot product of velocity and magnetic field

                vb_2star  =  sum(v_2star(imx:imz)*b_2star(ibx:ibz))

                ! Choose right Alfven wave according to speed of contact discontinuity

                if(sm .ge. zero) then

                   ! Conservative variables for left Alfven intermediate state
                   
                   u_2star(idn)  =  u_starl(idn)
                   u_2star(imx)  =  u_2star(idn)*v_2star(imx)
                   u_2star(imy)  =  u_2star(idn)*v_2star(imy)
                   u_2star(imz)  =  u_2star(idn)*v_2star(imz)
                   u_2star(ibx)  =  b_2star(ibx)
                   u_2star(iby)  =  b_2star(iby)
                   u_2star(ibz)  =  b_2star(ibz)

                   ! Energy of Alfven intermediate state Eq. 63

                   u_2star(ien)  =  u_starl(ien) - b_sig*dn_lsqt*(vb_starl - vb_2star)

                   ! Left Alfven intermediate flux Eq. 65

                   f(:,i)  =  fl(:,i) + alfven_l*u_2star(:) - (alfven_l - sl)*u_starl(:) - sl*ul(:,i)

                else if (sm .le. zero) then

                   ! Conservative variables for right Alfven intermediate state

                   u_2star(idn)  =  u_starr(idn)
                   u_2star(imx)  =  u_2star(idn)*v_2star(imx)
                   u_2star(imy)  =  u_2star(idn)*v_2star(imy)
                   u_2star(imz)  =  u_2star(idn)*v_2star(imz)
                   u_2star(ibx)  =  b_2star(ibx)
                   u_2star(iby)  =  b_2star(iby)
                   u_2star(ibz)  =  b_2star(ibz)

                   ! Energy of Alfven intermediate state Eq. 63

                   u_2star(ien)  =  u_starr(ien) + b_sig*dn_rsqt*(vb_starr - vb_2star)

                   ! Right Alfven intermediate flux Eq. 65

                   f(:,i)  =  fr(:,i) + alfven_r*u_2star(:) - (alfven_r - sr)*u_starr(:) - sr*ur(:,i)


                else ! sm = 0

                   ! Conservative variables for left Alfven intermediate state
                   
                   u_2star(idn)  =  u_starl(idn)
                   u_2star(imx)  =  u_2star(idn)*v_2star(imx)
                   u_2star(imy)  =  u_2star(idn)*v_2star(imy)
                   u_2star(imz)  =  u_2star(idn)*v_2star(imz)
                   u_2star(ibx)  =  b_2star(ibx)
                   u_2star(iby)  =  b_2star(iby)
                   u_2star(ibz)  =  b_2star(ibz)
                   
                   

                   ! Energy for Alfven intermediate state Eq. 63

                   u_2star(ien)  =  u_starl(ien) - b_sig*dn_lsqt*(vb_starl - vb_2star)

                   ! Left Alfven intermediate flux Eq. 65

                   f(:,i)  =  fl(:,i) + alfven_l*u_2star(:) - (alfven_l - sl)*u_starl(:) - sl*ul(:,i)

                   ! Conservative variables for right Alfven intermediate state

                   u_2star(idn)  =  u_starr(idn)
                   u_2star(imx)  =  u_2star(idn)*v_2star(imx)
                   u_2star(imy)  =  u_2star(idn)*v_2star(imy)
                   u_2star(imz)  =  u_2star(idn)*v_2star(imz)
                   u_2star(ibx)  =  b_2star(ibx)
                   u_2star(iby)  =  b_2star(iby)
                   u_2star(ibz)  =  b_2star(ibz)

                   ! Energy for Alfven intermediate state Eq. 63
                   
                   u_2star(ien)  =  u_starr(ien)  + b_sig*dn_rsqt*(vb_starr - vb_2star)

                   ! Right Alfven intermediate flux Eq. 65

                   f(:,i)  =  half*(f(:,i) + fr(:,i) + alfven_r*u_2star(:) - (alfven_r - sr)*u_starr(:) - sr*ur(:,i))

                endif  ! sm = 0

             endif  ! alfven_l .le. 0 and alfven_r .ge. 0

          else ! B_x = 0

             ! Intermediate state for B_x = 0

             if(sm .ge. zero) then

                ! Left intermediate flux Eq. 64

                f(:,i)  =  fl(:,i) + sl*(u_starl(:) - ul(:,i))

             else if(sm .le. zero) then

                f(:,i)  =  fr(:,i) + sr*(u_starr(:) - ur(:,i))

             else ! sm = 0

                ! Average left and right flux if both sm = 0 = B_x

                f(:,i)  =  half*((fl(:,i) + sl*(u_starl(:) - ul(:,i))) + (fr(:,i) + sr*(u_starr(:) - ur(:,i))))

             endif  ! sm = 0

          endif     ! B_x = 0 

       endif

    end do
    
  end subroutine riemann_hlld

end module hlld
