! Code Copyright (C) 2006 Michal Hanasz
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
!  Dr. Artur Gawryszczak, CAMK, Warszawa.
!---------------------------------------------------------------------------------------------------------------------------


#include "piernik.def"
!>
!!  \brief This module implements HLLD Riemann solver following the work of Miyoshi & Kusano (2005)
!<
module hlld
! pulled by RIEMANN

  implicit none

  private
  public :: riemann_hlld, fluxes

contains

  function fluxes(u,b_cc) result(f)

    use constants,  only: half, xdim, ydim, zdim
    use fluidindex, only: flind
    use fluidtypes, only: component_fluid
    use func,       only: ekin

    implicit none

    real, dimension(:,:),        intent(in)         :: u
    real, dimension(:,:),        intent(inout)      :: b_cc

    real, dimension(size(u,1), size(u,2))           :: f
    real, dimension(size(u,2))                      :: vx, vy, vz, p_t
    integer                                         :: ip
    class(component_fluid),      pointer            :: fl


    do ip = 1, flind%fluids

       fl => flind%all_fluids(ip)%fl

       vx  =  u(fl%imx,:)/u(fl%idn,:)
       vy  =  u(fl%imy,:)/u(fl%idn,:)
       vz  =  u(fl%imz,:)/u(fl%idn,:)

        if (fl%has_energy) then
           p_t = fl%gam_1*(u(fl%ien,:) - ekin(u(fl%imx,:), u(fl%imy,:), u(fl%imz,:), u(fl%idn,:)) - half*sum(b_cc(xdim:zdim,:)**2)) + half*sum(b_cc(xdim:zdim,:)**2)
        endif

       f(fl%idn,:)  =  u(fl%imx,:)
       f(fl%imx,:)  =  u(fl%imx,:)*vx(:) + p_t(:) - b_cc(xdim,:)**2
       f(fl%imy,:)  =  u(fl%imy,:)*vx(:) - b_cc(xdim,:)*b_cc(ydim,:)
       f(fl%imz,:)  =  u(fl%imz,:)*vx(:) - b_cc(xdim,:)*b_cc(zdim,:)
       b_cc(ydim,:) =  b_cc(ydim,:)*vx(:) - b_cc(xdim,:)*vy(:)
       b_cc(zdim,:) =  b_cc(zdim,:)*vx(:) - b_cc(xdim,:)*vz(:)
       if (fl%has_energy) then
          f(fl%ien,:)  =  (u(fl%ien,:) + p_t(:))*vx(:) - b_cc(xdim,:)*(b_cc(xdim,:)*vx(:) + b_cc(ydim,:)*vy(:) + b_cc(zdim,:)*vz(:))
       endif

    enddo

    return

  end function fluxes

 !-------------------------------------------------------------------------------------------------------------------------------------------------


  subroutine riemann_hlld(n,f,ul,ur,b_cc,b_ccl,b_ccr,gamma)

    ! external procedures

    use constants,  only: half, zero, xdim, ydim, zdim, idn, imx, imy, imz, ien
    use fluidindex, only: flind
    use func,       only: operator(.notequals.)

    ! arguments

    implicit none

    integer,                       intent(in)    :: n
    real, dimension(:,:), pointer, intent(inout) :: f
    real, dimension(:,:), pointer, intent(in)    :: ul, ur
    real, dimension(:,:), pointer, intent(in)    :: b_cc
    real, dimension(:,:), pointer, intent(in)    :: b_ccl, b_ccr
    real,                          intent(in)    :: gamma

    ! Local variables

    integer                                      :: i
    real, parameter                              :: four = 4.0
    real, parameter                              :: one  = 1.0
    real                                         :: sm, sm_nr, sm_dr, sl, sr
    real                                         :: alfven_l, alfven_r, c_fastl, c_fastr, gampr_l, gampr_r
    real                                         :: slsm, srsm, slvxl, srvxr, smvxl, smvxr, srmsl, srtsl, dn_l, dn_r
    real                                         :: b_lr, b_lrgam, magprl, magprr, prt_star, b_sig, enl, enr
    real                                         :: coeff_1, coeff_2, dn_lsqt, dn_rsqt, add_dnsq, mul_dnsq
    real                                         :: vb_l, vb_starl, vb_r, vb_starr, vb_2star

    ! Local arrays

    real, dimension(flind%all,n)                 :: fl, fr
    real, dimension(flind%all,n)                 :: prl, prr ! enl, enr
    real, dimension(flind%all)                   :: u_starl, u_starr, u_2star, v_starl, v_starr, v_2star
    real, dimension(xdim:zdim,n)                 :: b_cclf, b_ccrf
    real, dimension(xdim:zdim)                   :: b_starl, b_starr, b_2star

    ! SOLVER

    do i = 1,n

       ! Total left and right pressure

       prl(ien,i) = ul(ien,i)
       prr(ien,i) = ur(ien,i)

       ! Left and right energy Eq. 2

       enl = (prl(ien,i)/(gamma -one)) + half*ul(idn,i)*sum(ul(imx:imz,i)**2) + half*sum(b_ccl(xdim:zdim,i)**2)
       enr = (prr(ien,i)/(gamma -one)) + half*ur(idn,i)*sum(ur(imx:imz,i)**2) + half*sum(b_ccr(xdim:zdim,i)**2)

       gampr_l = gamma*prl(ien,i)
       gampr_r = gamma*prr(ien,i)

       ! Fast magnetosonic waves Eq. 3

       c_fastl  =   (gampr_l+(b_ccl(xdim,i)**2+b_ccl(ydim,i)**2+b_ccl(zdim,i)**2))  &
                                       + sqrt((gampr_l+(b_ccl(xdim,i)**2+b_ccl(ydim,i)**2+b_ccl(zdim,i)**2))**2-(four*gampr_l*b_ccl(xdim,i)**2))

       c_fastl = sqrt(half*c_fastl/ul(idn,i))

       c_fastr  =   (gampr_r+(b_ccr(xdim,i)**2+b_ccr(ydim,i)**2+b_ccr(zdim,i)**2))  &
                                       + sqrt((gampr_r+(b_ccr(xdim,i)**2+b_ccr(ydim,i)**2+b_ccr(zdim,i)**2))**2-(four*gampr_r*b_ccr(xdim,i)**2))

       c_fastr  =  sqrt(half*c_fastr/ur(idn,i))


       ! Eq. (67)

       sl  =  min(ul(imx,i) ,ur(imx,i)) - max(c_fastl,c_fastr)
       sr  =  max(ul(imx,i), ur(imx,i)) + max(c_fastl,c_fastr)

       ! Magnetic pressure

       magprl  =  half*sum(b_ccl(xdim:zdim,i)*b_ccl(xdim:zdim,i))
       magprr  =  half*sum(b_ccr(xdim:zdim,i)*b_ccr(xdim:zdim,i))

       ! Left flux

       fl(idn,i) = ul(idn,i)*ul(imx,i)
       fl(imx,i) = ul(idn,i)*ul(imx,i)**2 + prl(ien,i) - b_ccl(xdim,i)**2
       fl(imy,i) = ul(idn,i)*ul(imy,i)*ul(imx,i) - b_ccl(xdim,i)*b_ccl(ydim,i)
       fl(imz,i) = ul(idn,i)*ul(imz,i)*ul(imx,i) - b_ccl(xdim,i)*b_ccl(zdim,i)
       fl(ien,i) = (enl + prl(ien,i))*ul(imx,i) - b_ccl(xdim,i)*(sum(ul(imx:imz,i)*b_ccl(xdim:zdim,i)))
       b_cclf(ydim,i) = b_ccl(ydim,i)*ul(imx,i) - b_ccl(xdim,i)*ul(imy,i)
       b_cclf(zdim,i) = b_ccl(zdim,i)*ul(imx,i) - b_ccl(xdim,i)*ul(imz,i)

       ! Right flux

       fr(idn,i) = ur(idn,i)*ur(imx,i)
       fr(imx,i) = ur(idn,i)*ur(imx,i)**2 + prr(ien,i) - b_ccr(xdim,i)**2
       fr(imy,i) = ur(idn,i)*ur(imy,i)*ur(imx,i) - b_ccr(xdim,i)*b_ccr(ydim,i)
       fr(imz,i) = ur(idn,i)*ur(imz,i)*ur(imx,i) - b_ccr(xdim,i)*b_ccr(zdim,i)
       fr(ien,i) = (enr + prr(ien,i))*ur(imx,i) - b_ccr(xdim,i)*(sum(ur(imx:imz,i)*b_ccr(xdim:zdim,i)))
       b_ccrf(ydim,i) = b_ccr(ydim,i)*ur(imx,i) - b_ccr(xdim,i)*ur(imy,i)
       b_ccrf(zdim,i) = b_ccr(zdim,i)*ul(imx,i) - b_ccr(xdim,i)*ur(imz,i)

       ! HLLD fluxes

       if (sl .ge.  zero) then
          f(idn,i)  =  fl(idn,i)
          f(imx,i)  =  fl(imx,i)
          f(imy,i)  =  fl(imy,i)
          f(imz,i)  =  fl(imz,i)
          f(ien,i)  =  fl(ien,i)
          b_cc(ydim,i) = b_cclf(ydim,i)
          b_cc(zdim,i) = b_cclf(zdim,i)

       else if (sr .le.  zero) then
          f(idn,i)  =  fr(idn,i)
          f(imx,i)  =  fr(imx,i)
          f(imy,i)  =  fr(imy,i)
          f(imz,i)  =  fr(imz,i)
          f(ien,i)  =  fr(ien,i)
          b_cc(ydim,i) = b_ccrf(ydim,i)
          b_cc(zdim,i) = b_ccrf(zdim,i)

       else

          ! Speed of contact discontinuity Eq. 38

          sm_nr = (sr - ur(imx,i))*ur(idn,i)*ur(imx,i) - (sl - ul(imx,i))*ul(idn,i)*ul(imx,i) - prr(ien,i) + prl(ien,i)
          sm_dr = (sr - ur(imx,i))*ur(idn,i) - (sl - ul(imx,i))*ul(idn,i)
          sm    = sm_nr/sm_dr

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
          b_lr     =  b_ccl(xdim,i)*b_ccr(xdim,i)
          b_lrgam  =  b_lr/gamma

          ! Pressure of intermediate state Eq. (23)

          prt_star  =  half*((prl(ien,i)+dn_l*smvxl) + (prr(ien,i)+dn_r*smvxr))  !< Check for 0.5.

          ! Normal components of velocity and magnetic field

          v_starl(imx)  =  sm
          v_starr(imx)  =  sm

          b_starl(xdim)  =  b_ccl(xdim,i)
          b_starr(xdim)  =  b_ccr(xdim,i)

          ! Transversal components of magnetic field for left states (Eq. 45 & 47), taking degeneracy into account

          coeff_1  =  dn_l*slsm - b_lr

          if ((coeff_1 .notequals. zero) .and. b_lrgam .le. ul(ien,i)) then

             coeff_2  =  (dn_l*slvxl - b_lr)/coeff_1

             b_starl(ydim)   =  b_ccl(ydim,i)*coeff_2
             b_starl(zdim)   =  b_ccl(zdim,i)*coeff_2

          else

             ! Calculate HLL left states

             b_starl(ydim) = ((sr*b_ccr(ydim,i) - sl*b_ccl(ydim,i)) - (b_ccrf(ydim,i) - b_cclf(ydim,i)))/srmsl
             b_starl(zdim) = ((sr*b_ccr(zdim,i) - sl*b_ccl(zdim,i)) - (b_ccrf(zdim,i) - b_cclf(zdim,i)))/srmsl


          endif

          ! Transveral components of velocity Eq. 42

          coeff_1  =  b_ccl(xdim,i)/dn_l

          v_starl(imy) = ul(imy,i)  + coeff_1*(b_ccl(ydim,i) - b_starl(ydim))
          v_starl(imz) = ul(imz,i)  + coeff_1*(b_ccl(zdim,i) - b_starl(zdim))

          ! Transversal components of magnetic field for right states (Eq. 45 & 47), taking degeneracy into account

          coeff_1  =  dn_r*srsm - b_lr

          if ((coeff_1 .notequals. zero) .and. b_lrgam .le. ur(ien,i)) then

             coeff_2  =  (dn_r*srvxr - b_lr)/coeff_1

             b_starr(ydim)  =  b_ccr(ydim,i)*coeff_2
             b_starr(zdim)  =  b_ccr(zdim,i)*coeff_2

          else

             ! Calculate HLL right states

             b_starr(ydim)  =  ((sr*b_ccr(ydim,i) - sl*b_ccl(ydim,i)) - (b_ccrf(ydim,i) - b_cclf(ydim,i)))/srmsl
             b_starr(zdim)  =  ((sr*b_ccr(zdim,i) - sl*b_ccl(zdim,i)) - (b_ccrf(zdim,i) - b_cclf(zdim,i)))/srmsl

          endif

          ! Transversal components of velocity Eq. 42

          coeff_1  =  b_ccr(xdim,i)/dn_r

          v_starr(imy)  =   ur(imy,i) + coeff_1*(b_ccr(ydim,i) - b_starr(ydim))
          v_starr(imz)  =   ur(imz,i) + coeff_1*(b_ccr(zdim,i) - b_starr(zdim))

          ! Dot product of velocity and magnetic field

          vb_l  =  sum(ul(imx:imz,i)*b_ccl(xdim:zdim,i))
          vb_r  =  sum(ur(imx:imz,i)*b_ccr(xdim:zdim,i))
          vb_starl  =  sum(v_starl(imx:imz)*b_starl(xdim:zdim))
          vb_starr  =  sum(v_starr(imx:imz)*b_starr(xdim:zdim))


          ! Left intermediate state conservative form

          u_starl(idn)  =  dn_l/slsm
          u_starl(imx)  =  u_starl(idn)*v_starl(imx)
          u_starl(imy)  =  u_starl(idn)*v_starl(imy)
          u_starl(imz)  =  u_starl(idn)*v_starl(imz)

          ! Right intermediate state conservative form

          u_starr(idn)  =  dn_r/srsm
          u_starr(imx)  =  u_starr(idn)*v_starr(imx)
          u_starr(imy)  =  u_starr(idn)*v_starr(imy)
          u_starr(imz)  =  u_starr(idn)*v_starr(imz)

          ! Total energy of left and right intermediate states Eq. (48)

          u_starl(ien) = (slvxl*enl - prl(ien,i)*ul(imx,i) + prt_star*sm + b_ccl(xdim,i)*(vb_l - vb_starl))/slsm
          u_starr(ien) = (srvxr*enr - prr(ien,i)*ur(imx,i) + prt_star*sm + b_ccr(xdim,i)*(vb_r - vb_starr))/srsm

          ! Cases for B_x .ne. and .eq. zero

          if (abs(b_ccl(xdim,i)) > zero) then

             ! Left and right Alfven waves velocity Eq. 51

             dn_lsqt  =  sqrt(u_starl(idn))
             dn_rsqt  =  sqrt(u_starr(idn))

             alfven_l  =  sm - abs(b_ccl(xdim,i))/dn_lsqt
             alfven_r  =  sm + abs(b_ccr(xdim,i))/dn_rsqt

             ! Intermediate discontinuities

             if (alfven_l > zero) then

                ! Left intermediate flux Eq. 64

                f(idn,i) = fl(idn,i) + sl*(u_starl(idn) - ul(idn,i))
                f(imx,i) = fl(imx,i) + sl*(u_starl(imx) - ul(idn,i)*ul(imx,i))
                f(imy,i) = fl(imy,i) + sl*(u_starl(imy) - ul(idn,i)*ul(imy,i))
                f(imz,i) = fl(imz,i) + sl*(u_starl(imz) - ul(idn,i)*ul(imz,i))
                f(ien,i) = fl(ien,i) + sl*(u_starl(ien) - enl)
                b_cc(ydim,i) = b_cclf(ydim,i) + sl*(b_starl(ydim) - b_ccl(ydim,i))
                b_cc(zdim,i) = b_cclf(zdim,i) + sl*(b_starl(zdim) - b_ccl(zdim,i))

             else if (alfven_r < zero) then

                ! Right intermediate flux Eq. 64

                f(idn,i) = fr(idn,i) + sl*(u_starr(idn) - ur(idn,i))
                f(imx,i) = fr(imx,i) + sl*(u_starr(imx) - ur(idn,i)*ur(imx,i))
                f(imy,i) = fr(imy,i) + sl*(u_starr(imy) - ur(idn,i)*ur(imy,i))
                f(imz,i) = fr(imz,i) + sl*(u_starr(imz) - ur(idn,i)*ur(imz,i))
                f(ien,i) = fr(ien,i) + sl*(u_starr(ien) - enr)
                b_cc(ydim,i) = b_ccrf(ydim,i) + sl*(b_starr(ydim) - b_ccr(ydim,i))
                b_cc(zdim,i) = b_ccrf(zdim,i) + sl*(b_starr(zdim) - b_ccr(zdim,i))

             else ! alfven_l .le. zero .le. alfven_r

                ! Arragnge for sign of normal component of magnetic field

                if (b_ccl(xdim,i) .ge. zero) then

                   b_sig = one

                else

                   b_sig = -one

                endif

                ! Sum and product of density square-root

                add_dnsq  =  dn_lsqt + dn_rsqt
                mul_dnsq  =  dn_lsqt*dn_rsqt

                ! Components of velocity Eq. 39, 59, 60 and magnetic field Eq. 61, 62

                v_2star(imx)  =  sm
                v_2star(imy)  =  ((dn_lsqt*v_starl(imy) + dn_rsqt*v_starr(imy)) + b_sig*(b_starr(ydim) - b_starl(ydim)))/add_dnsq
                v_2star(imz)  =  ((dn_lsqt*v_starl(imz) + dn_rsqt*v_starr(imz)) + b_sig*(b_starr(zdim) - b_starl(zdim)))/add_dnsq

                b_2star(xdim)  =  b_ccl(xdim,i)
                b_2star(ydim)  =  ((dn_lsqt*b_starr(ydim) + dn_rsqt*b_starl(ydim)) + b_sig*mul_dnsq*(v_starr(imy) - v_starl(imy)))/add_dnsq
                b_2star(zdim)  =  ((dn_lsqt*b_starr(zdim) + dn_rsqt*b_starl(zdim)) + b_sig*mul_dnsq*(v_starr(imz) - v_starl(imz)))/add_dnsq

                ! Dot product of velocity and magnetic field

                vb_2star  =  sum(v_2star(imx:imz)*b_2star(xdim:zdim))

                ! Choose right Alfven wave according to speed of contact discontinuity

                if (sm .ge. zero) then

                   ! Conservative variables for left Alfven intermediate state

                   u_2star(idn)  =  u_starl(idn)
                   u_2star(imx)  =  u_2star(idn)*v_2star(imx)
                   u_2star(imy)  =  u_2star(idn)*v_2star(imy)
                   u_2star(imz)  =  u_2star(idn)*v_2star(imz)

                   ! Energy of Alfven intermediate state Eq. 63

                   u_2star(ien)  =  u_starl(ien) - b_sig*dn_lsqt*(vb_starl - vb_2star)

                   ! Left Alfven intermediate flux Eq. 65

                   f(idn,i) = fl(idn,i) + alfven_l*u_2star(idn) - (alfven_l - sl)*u_starl(idn) - sl*ul(idn,i)
                   f(imx,i) = fl(imx,i) + alfven_l*u_2star(imx) - (alfven_l - sl)*u_starl(imx) - sl*ul(idn,i)*ul(imx,i)
                   f(imy,i) = fl(imy,i) + alfven_l*u_2star(imy) - (alfven_l - sl)*u_starl(imy) - sl*ul(idn,i)*ul(imy,i)
                   f(imz,i) = fl(imz,i) + alfven_l*u_2star(imz) - (alfven_l - sl)*u_starl(imz) - sl*ul(idn,i)*ul(imz,i)
                   f(ien,i) = fl(ien,i) + alfven_l*u_2star(ien) - (alfven_l - sl)*u_starl(ien) - sl*enl
                   b_cc(ydim,i) = b_cclf(ydim,i) + alfven_l*b_2star(ydim) - (alfven_l - sl)*b_starl(ydim) - sl*b_ccl(ydim,i)
                   b_cc(zdim,i) = b_cclf(zdim,i) + alfven_l*b_2star(zdim) - (alfven_l - sl)*b_starl(zdim) - sl*b_ccl(zdim,i)


                else if (sm .le. zero) then

                   ! Conservative variables for right Alfven intermediate state

                   u_2star(idn)  =  u_starr(idn)
                   u_2star(imx)  =  u_2star(idn)*v_2star(imx)
                   u_2star(imy)  =  u_2star(idn)*v_2star(imy)
                   u_2star(imz)  =  u_2star(idn)*v_2star(imz)

                   ! Energy of Alfven intermediate state Eq. 63

                   u_2star(ien)  =  u_starr(ien) + b_sig*dn_rsqt*(vb_starr - vb_2star)

                   ! Right Alfven intermediate flux Eq. 65

                   f(idn,i) = fr(idn,i) + alfven_r*u_2star(idn) - (alfven_r - sl)*u_starr(idn) - sr*ur(idn,i)
                   f(imx,i) = fr(imx,i) + alfven_r*u_2star(imx) - (alfven_r - sl)*u_starr(imx) - sr*ur(idn,i)*ur(imx,i)
                   f(imy,i) = fr(imy,i) + alfven_r*u_2star(imy) - (alfven_r - sl)*u_starr(imy) - sr*ur(idn,i)*ur(imy,i)
                   f(imz,i) = fr(imz,i) + alfven_r*u_2star(imz) - (alfven_r - sl)*u_starr(imz) - sr*ur(idn,i)*ur(imz,i)
                   f(ien,i) = fr(ien,i) + alfven_r*u_2star(ien) - (alfven_r - sl)*u_starr(ien) - sr*enr
                   b_cc(ydim,i) = b_ccrf(ydim,i) + alfven_r*b_2star(ydim) - (alfven_r - sr)*b_starr(ydim) - sr*b_ccr(ydim,i)
                   b_cc(zdim,i) = b_ccrf(zdim,i) + alfven_r*b_2star(zdim) - (alfven_r - sr)*b_starr(zdim) - sr*b_ccr(zdim,i)


                else ! sm = 0

                   ! Conservative variables for left Alfven intermediate state

                   u_2star(idn)  =  u_starl(idn)
                   u_2star(imx)  =  u_2star(idn)*v_2star(imx)
                   u_2star(imy)  =  u_2star(idn)*v_2star(imy)
                   u_2star(imz)  =  u_2star(idn)*v_2star(imz)


                   ! Energy for Alfven intermediate state Eq. 63

                   u_2star(ien)  =  u_starl(ien) - b_sig*dn_lsqt*(vb_starl - vb_2star)

                   ! Left Alfven intermediate flux Eq. 65

                   f(idn,i) = fl(idn,i) + alfven_l*u_2star(idn) - (alfven_l - sl)*u_starl(idn) - sl*ul(idn,i)
                   f(imx,i) = fl(imx,i) + alfven_l*u_2star(imx) - (alfven_l - sl)*u_starl(imx) - sl*ul(idn,i)*ul(imx,i)
                   f(imy,i) = fl(imy,i) + alfven_l*u_2star(imy) - (alfven_l - sl)*u_starl(imy) - sl*ul(idn,i)*ul(imy,i)
                   f(imz,i) = fl(imz,i) + alfven_l*u_2star(imz) - (alfven_l - sl)*u_starl(imz) - sl*ul(idn,i)*ul(imz,i)
                   f(ien,i) = fl(ien,i) + alfven_l*u_2star(ien) - (alfven_l - sl)*u_starl(ien) - sl*enl
                   b_cc(ydim,i) = b_cclf(ydim,i) + alfven_l*b_2star(ydim) - (alfven_l - sl)*b_2star(ydim) - sl*b_ccl(ydim,i)
                   b_cc(zdim,i) = b_cclf(zdim,i) + alfven_l*b_2star(zdim) - (alfven_l - sl)*b_2star(zdim) - sl*b_ccl(zdim,i)

                   ! Conservative variables for right Alfven intermediate state

                   u_2star(idn)  =  u_starr(idn)
                   u_2star(imx)  =  u_2star(idn)*v_2star(imx)
                   u_2star(imy)  =  u_2star(idn)*v_2star(imy)
                   u_2star(imz)  =  u_2star(idn)*v_2star(imz)

                   ! Energy for Alfven intermediate state Eq. 63

                   u_2star(ien)  =  u_starr(ien)  + b_sig*dn_rsqt*(vb_starr - vb_2star)

                   ! Right Alfven intermediate flux Eq. 65

                   f(idn,i)  =  half*(f(idn,i) + (fr(idn,i) + alfven_r*u_2star(idn) - (alfven_r - sr)*u_starr(idn) - sr*ur(idn,i)))
                   f(imx,i)  =  half*(f(imx,i) + (fr(imx,i) + alfven_r*u_2star(imx) - (alfven_r - sr)*u_starr(imx) - sr*ur(idn,i)*ur(imx,i)))
                   f(imy,i)  =  half*(f(imy,i) + (fr(imy,i) + alfven_r*u_2star(imy) - (alfven_r - sr)*u_starr(imy) - sr*ur(idn,i)*ur(imy,i)))
                   f(imz,i)  =  half*(f(imz,i) + (fr(imz,i) + alfven_r*u_2star(imz) - (alfven_r - sr)*u_starr(imz) - sr*ur(idn,i)*ur(imz,i)))
                   f(ien,i)  =  half*(f(ien,i) + (fr(ien,i) + alfven_r*u_2star(ien) - (alfven_r - sr)*u_starr(ien) - sr*enr))
                   b_cc(ydim,i) = half*(b_cc(ydim,i) + (b_ccrf(ydim,i) + alfven_r*b_2star(ydim) - (alfven_r - sr)*b_starr(ydim) - sr*b_ccr(ydim,i)))
                   b_cc(zdim,i) = half*(b_cc(zdim,i) + (b_ccrf(zdim,i) + alfven_r*b_2star(zdim) - (alfven_r - sr)*b_starr(zdim) - sr*b_ccr(zdim,i)))

                endif  ! sm = 0

             endif  ! alfven_l .le. 0 and alfven_r .ge. 0

          else ! B_x = 0

             ! Intermediate state for B_x = 0

             if (sm .ge. zero) then

                ! Left intermediate flux Eq. 64

                f(idn,i)  =  fl(idn,i) + sl*(u_starl(idn) - ul(idn,i))
                f(imx,i)  =  fl(imx,i) + sl*(u_starl(imx) - ul(idn,i)*ul(imx,i))
                f(imy,i)  =  fl(imy,i) + sl*(u_starl(imy) - ul(idn,i)*ul(imy,i))
                f(imz,i)  =  fl(imz,i) + sl*(u_starl(imz) - ul(idn,i)*ul(imz,i))
                f(ien,i)  =  fl(ien,i) + sl*(u_starl(ien) - enl)
                b_cc(ydim,i) = b_cclf(ydim,i) + sl*(b_starl(ydim) - b_ccl(ydim,i))
                b_cc(zdim,i) = b_cclf(zdim,i) + sl*(b_starl(zdim) - b_ccl(zdim,i))

             else if (sm .le. zero) then

                f(idn,i)  =  fr(idn,i) + sr*(u_starr(idn) - ur(idn,i))
                f(imx,i)  =  fr(imx,i) + sr*(u_starr(imx) - ur(idn,i)*ur(imx,i))
                f(imy,i)  =  fr(imy,i) + sr*(u_starr(imy) - ur(idn,i)*ur(imy,i))
                f(imz,i)  =  fr(imz,i) + sr*(u_starr(imz) - ur(idn,i)*ur(imz,i))
                f(ien,i)  =  fr(ien,i) + sr*(u_starr(ien) - enr)
                b_cc(ydim,i) = b_ccrf(ydim,i) + sr*(b_starr(ydim) - b_ccr(ydim,i))
                b_cc(zdim,i) = b_ccrf(zdim,i) + sr*(b_starr(zdim) - b_ccr(zdim,i))

             else ! sm = 0

                ! Average left and right flux if both sm = 0 = B_x

                f(idn,i)  =  half*((fl(idn,i) + sl*(u_starl(idn) - ul(idn,i))) + (fr(idn,i) + sr*(u_starr(idn) - ur(idn,i))))
                f(imx,i)  =  half*((fl(imx,i) + sl*(u_starl(imx) - ul(idn,i)*ul(imx,i))) + (fr(imx,i) + sr*(u_starr(imx) - ur(idn,i)*ur(imx,i))))
                f(imy,i)  =  half*((fl(imy,i) + sl*(u_starl(imy) - ul(idn,i)*ul(imy,i))) + (fr(imy,i) + sr*(u_starr(imy) - ur(idn,i)*ur(imy,i))))
                f(imz,i)  =  half*((fl(imz,i) + sl*(u_starl(imz) - ul(idn,i)*ul(imz,i))) + (fr(imz,i) + sr*(u_starr(imz) - ur(idn,i)*ur(imz,i))))
                f(ien,i)  =  half*((fl(ien,i) + sl*(u_starl(ien) - enl)) + (fr(ien,i) + sr*(u_starr(ien) - enr)))
                b_cc(ydim,i) = half*((b_cclf(ydim,i) + sl*(b_starl(ydim) - b_ccl(ydim,i))) + (b_ccrf(ydim,i) + sr*(b_starr(ydim) - b_ccr(ydim,i))))
                b_cc(zdim,i) = half*((b_cclf(zdim,i) + sl*(b_starl(zdim) - b_ccl(zdim,i))) + (b_ccrf(zdim,i) + sr*(b_starr(zdim) - b_ccr(zdim,i))))

             endif  ! sm = 0

          endif     ! B_x = 0

       endif

    enddo

  end subroutine riemann_hlld




end module hlld
