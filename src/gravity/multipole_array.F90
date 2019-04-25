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

#include "piernik.h"

!> \brief This module implements an array containing radial distribution of multipole moments

module multipole_array
! pulled by MULTIGRID && SELF_GRAV

   implicit none

   private
   public :: mpole_container
   public :: interp_pt2mom, interp_mom2pot  ! initialized in multigrid_gravity

   ! namelist parameters for MULTIGRID_GRAVITY
   logical :: interp_pt2mom  !< Distribute contribution from a cell between two adjacent radial bins (linear interpolation in radius)
   logical :: interp_mom2pot !< Compute the potential from moments from two adjacent radial bins (linear interpolation in radius)

   type :: mpole_container

      real, dimension(:,:,:), allocatable          :: Q       !< The whole moment array with dependence on radius
      integer(kind=4),                     private :: lmax    !< Maximum l-order of multipole moments
      integer(kind=4),                     private :: mmax    !< Maximum m-order of multipole moments. Equal to lmax by default.

      ! auxiliary factors
      real, dimension(:,:,:), allocatable, private :: k12     !< array of Legendre recurrence factors
      real, dimension(:),     allocatable, private :: ofact   !< arrays of Legendre normalization factor (compressed)
      ! radial discretization
      integer,                             private :: irmin   !< minimum Q(:, :, r) index in use
      integer,                             private :: irmax   !< maximum Q(:, :, r) index in use
      integer,                             private :: rqbin   !< number of radial samples of multipoles
      real,                                private :: drq     !< radial resolution of multipoles
      real,                                private :: rscale  !< scaling factor that limits risk of floating poin
      real, dimension(:),     allocatable, private :: rn      !< 0..l_max sized array of positive powers of r
      real, dimension(:),     allocatable, private :: irn     !< 0..l_max sized array of negative powers of r

   contains

      procedure          :: init_once
      procedure          :: refresh
      procedure          :: reset
      procedure          :: red_int_norm
      procedure          :: cleanup
      procedure          :: point2moments
      procedure          :: moments2pot
      procedure, private :: lm
      procedure, private :: geomfac4moments

   end type mpole_container

   enum, bind(C)
      enumerator :: INSIDE = 1, OUTSIDE  !< indices for interior and exterior multipole expansion, respectively
   end enum

contains

   subroutine init_once(this, lmax, mmax)

      use dataio_pub, only: die

      implicit none

      class(mpole_container), intent(inout) :: this  !< object invoking type-bound procedure
      integer(kind=4),        intent(in)    :: lmax  !< Maximum l-order of multipole moments
      integer(kind=4),        intent(in)    :: mmax  !< Maximum m-order of multipole moments. Equal to lmax by default.

      integer :: l, m

      if (lmax < 0) call die("[multipole_array:init_once] lmax < 0")
      if (mmax < 0) call die("[multipole_array:init_once] mmax < 0")
      if (mmax > lmax) call die("[multipole_array:init_once] mmax > lmax")

      this%lmax = lmax
      this%mmax = mmax

      if (allocated(this%k12) .or. allocated(this%ofact)) call die("[multipole_array:init_once] k12 or ofact already allocated")
      allocate(this%k12(2, 1:lmax, 0:mmax), this%ofact(0:this%lm(int(lmax), int(2*mmax))))

      if (allocated(this%rn) .or. allocated(this%irn)) call die("[multipole_array:init_once] rn or irn already allocated")
      allocate(this%rn(0:lmax), this%irn(0:lmax))

      this%ofact(:) = 0. ! prevent spurious FP exceptions in multipole:img_mass2moments
      do l = 1, this%lmax
         do m = 0, min(l, int(this%mmax))
            if (m == 0) then
               this%ofact(l) = 2. ! this%lm(l,0)
            else                     ! this%ofact(l, 2*m) = ((-1)^m (2m-1)!!)**2 (l-m)! / (l+m)! ; Should work up to m=512 and even beyond
               this%ofact(this%lm(l, 2*m))   = this%ofact(this%lm(l, 2*(m-1))) * real(1 - 2*m)**2 / real((l+m) * (l-m+1))
               this%ofact(this%lm(l, 2*m-1)) = this%ofact(this%lm(l, 2*m))
            endif
            if (l>m) then
               this%k12(1, l, m) = real(2 * l - 1) / real(l - m)
               this%k12(2, l, m) = real(l + m - 1) / real(l - m)
            endif
         enddo
      enddo
      this%ofact(0:this%lmax) = 1. ! this%lm(0:this%lmax,0)

   end subroutine init_once

   subroutine refresh(this)

      use cg_level_finest, only: finest
      use constants,       only: xdim, zdim, GEO_XYZ, GEO_RPZ, HI, pMIN
      use dataio_pub,      only: die
      use domain,          only: dom
      use mpisetup,        only: piernik_MPI_Allreduce

      implicit none

      class(mpole_container), intent(inout) :: this  !< object invoking type-bound procedure

      if (associated(finest%level%first)) then
         select case (dom%geometry_type)
            case (GEO_XYZ)
               this%drq = minval(finest%level%first%cg%dl(:), mask=dom%has_dir(:)) / 2.
            case (GEO_RPZ)
               this%drq = min(finest%level%first%cg%dx, dom%C_(xdim)*finest%level%first%cg%dy, finest%level%first%cg%dz) / 2.
            case default
               call die("[multipole_array:refresh] Unsupported geometry.")
         end select
      else
         this%drq = maxval(dom%L_(:))
      endif
      call piernik_MPI_Allreduce(this%drq, pMIN)

      select case (dom%geometry_type)
         case (GEO_XYZ)
            this%rqbin = int(sqrt(sum(dom%L_(:)**2))/this%drq) + 1
            ! arithmetic average of the closest and farthest points of computational domain with respect to its center
            !>
            !!\todo check what happens if there are points that are really close to the domain center (maybe we should use a harmonic average?)
            !! Issue a warning or error if it is known that given lmax leads to FP overflows in rn(:) and irn(:)
            !<
            this%rscale = ( minval(dom%L_(:)) + sqrt(sum(dom%L_(:)**2)) )/4.
         case (GEO_RPZ)
            this%rqbin = int(sqrt((2.*dom%edge(xdim, HI))**2 + dom%L_(zdim)**2)/this%drq) + 1
            this%rscale = ( min(2.*dom%edge(xdim, HI), dom%L_(zdim)) + sqrt((2.*dom%edge(xdim, HI))**2 + dom%L_(zdim)**2) )/4.
         case default
            call die("[multipole_array:refresh] Unsupported geometry.")
      end select

      ! The particles don't affect this%rqbin by design.
      ! The particles are either close enough to contribute to the multipole radial distribution or
      ! they're far enough to be assigned to last radial bin.
      ! Please note that the flag cg%pset%p(i)%outside does not always imply contributing to the farthest bin.
      !
      ! Recovering radial profile of the multipoles far from the domain is possible but may cost way too much memory.
      !
      ! This simplification holds as long as the mass of particles that lie further than this%rqbin * this%drq
      ! is small compared to the total mass of gas and all particles.
      !
      ! The simulation with significant mass outside that radiuss (or even outside the domain) seems to be
      ! a bit pathological anyway.
      !
      ! ToDo: on each call integrate mass outside the domain and call warn() or even die()
      ! when their weight gets too big.
      !
      ! Another solution is to expand base level when too many particles go outside but still remain bound.

      if (allocated(this%Q)) deallocate(this%Q)
      allocate(this%Q(0:this%lm(int(this%lmax), int(2*this%mmax)), INSIDE:OUTSIDE, 0:this%rqbin))

   end subroutine refresh

!> \brief reset the multipole data

   subroutine reset(this)

      implicit none

      class(mpole_container), intent(inout) :: this  !< object invoking type-bound procedure

      this%Q = 0.
      this%irmax = 0
      this%irmin = this%rqbin

   end subroutine reset

!> \brief Integrate radially and apply normalization factor (the (4 \pi)/(2 l  + 1) terms cancel out)

   subroutine red_int_norm(this)

      use constants, only: pSUM, pMIN, pMAX
      use mpisetup,  only: piernik_MPI_Allreduce

      implicit none

      class(mpole_container), intent(inout) :: this  !< object invoking type-bound procedure

      integer :: r, rr

      call piernik_MPI_Allreduce(this%irmin, pMIN)
      call piernik_MPI_Allreduce(this%irmax, pMAX)

      ! integrate radially and apply normalization factor (the (4 \pi)/(2 l  + 1) terms cancel out)
      rr = max(1, this%irmin)
      this%Q(:, INSIDE, rr-1) = this%Q(:, INSIDE, rr-1) * this%ofact(:)
      do r = rr, this%irmax
         this%Q(:, INSIDE, r) = this%Q(:, INSIDE, r) * this%ofact(:) + this%Q(:, INSIDE, r-1)
      enddo

      this%Q(:, OUTSIDE, this%irmax+1) = this%Q(:, OUTSIDE, this%irmax+1) * this%ofact(:)
      do r = this%irmax, rr, -1
         this%Q(:, OUTSIDE, r) = this%Q(:, OUTSIDE, r) * this%ofact(:) + this%Q(:, OUTSIDE, r+1)
      enddo

      call piernik_MPI_Allreduce(this%Q(:, :, this%irmin:this%irmax), pSUM)

   end subroutine red_int_norm

   subroutine cleanup(this)

      implicit none

      class(mpole_container), intent(inout) :: this  !< object invoking type-bound procedure

      if (allocated(this%rn))    deallocate(this%rn)
      if (allocated(this%irn))   deallocate(this%irn)
      if (allocated(this%Q))     deallocate(this%Q)
      if (allocated(this%k12))   deallocate(this%k12)
      if (allocated(this%ofact)) deallocate(this%ofact)

   end subroutine cleanup


!>
!! \brief This function returns array index for "compressed" this%Q(:, inout, r) array replacing "plain" this%Q(l, m, inout, r) array
!! \details Originally there were separate indices for l and m multipole numbers in Q and ofact arrays. Since there are no Y_{l,m} harmonics for l<|m|, approximately half of the Q array
!! was left unused. For each m there is only lmax-|m|+1 valid Y_{l,m} harmonic contributions (for l = |m| .. lmax). The function lm(l, m) converts each valid (l,m) pair into an
!! unique index leaving no unused entries in the Q array.
!! Let store harmonics in the following way
!!    Y_{l,m}  in this%Q(this%lm(l, 2*|m|  ), :, :) for even (cosine, positive-m) harmonic,
!!    Y_{l,-m} in this%Q(this%lm(l, 2*|m|-1), :, :) for odd  (sine,   negative-m) harmonic.
!! Then for best performance consecutive entries in compressed this%Q(:, inout, r) array should be related to Y_{l,m} harmonics in the same order as they're referenced in point2moments and moments2pot.
!! So the ordering goes from index 0 as follows
!! Y_{0,0}, Y_{1,0}, ..., Y_{lmax, 0}, Y_{1,-1}, Y_{1,1}, ..., Y_{lmax,-1}, Y_{lmax,1}, Y_{2,-2}, Y{2,2}, ..., Y(lmax-1,-(mmax-1)), Y(lmax-1,mmax-1), Y(lmax,-(mmax-1)), Y(lmax,mmax-1), Y(lmax,-mmax), Y(lmax,mmax)
!! Does it looks a bit cryptic? I agree.
!!
!! Since this ordering does not depend on mmax we can safely truncate highest azimuthal moments and use mmax < lmax.
!<

   pure integer function lm(this, l, m)

      implicit none

      class(mpole_container), intent(in) :: this  !< object invoking type-bound procedure
      integer,                intent(in) :: l, m  !< harmonic indices

      integer :: mu

      mu = this%lmax - int((m+1)/2)

      lm = merge(this%lmax**2 - 1 + 2*l + mod(m+1, 2) - mu*(mu + 1), l, m /= 0)

   end function lm

!>
!! \brief Compute multipole moments for a single point
!!
!! \todo try to improve accuracy with linear interpolation over radius
!<

   subroutine point2moments(this, mass, x, y, z)

      use constants, only: zero
      use func,      only: operator(.notequals.)

      implicit none

      class(mpole_container), intent(inout) :: this  !< object invoking type-bound procedure
      real,                   intent(in)    :: mass  !< mass of the contributing point
      real,                   intent(in)    :: x     !< x coordinate of the contributing point
      real,                   intent(in)    :: y     !< y coordinate of the contributing point
      real,                   intent(in)    :: z     !< z coordinate of the contributing point

      real    :: sin_th, cos_th, sin_ph, cos_ph, del, cfac, sfac, tmpfac
      real    :: Ql, Ql1, Ql2
      integer :: l, m, ir, m2c

      call this%geomfac4moments(mass, x, y, z, sin_th, cos_th, sin_ph, cos_ph, ir, del)

      if (.not. interp_pt2mom) del = 0.

      ! monopole, the (0,0) moment; P_0 = 1.
      ! this%lm(0, 0) == 0
      this%Q(0, INSIDE,  ir) = this%Q(0, INSIDE,  ir) +  this%rn(0) * (1.-del)
      this%Q(0, OUTSIDE, ir) = this%Q(0, OUTSIDE, ir) + this%irn(0) * (1.-del)
      if (del .notequals. zero) then
         this%Q(0, INSIDE,  ir+1) = this%Q(0, INSIDE,  ir+1) +  this%rn(0) * del
         this%Q(0, OUTSIDE, ir+1) = this%Q(0, OUTSIDE, ir+1) + this%irn(0) * del
      endif

      ! axisymmetric (l,0) moments
      ! Legendre polynomial recurrence: l P_l = x (2l-1) P_{l-1} - (l-1) P_{l-2}, x \eqiv \cos(\theta)
      Ql2 = 0.
      Ql1 = 1.
      do l = 1, this%lmax
         ! this%lm(l, 0) == l
         Ql = cos_th * this%k12(1, l, 0) * Ql1 - this%k12(2, l, 0) * Ql2
         this%Q(l, INSIDE,  ir) = this%Q(l, INSIDE,  ir) +  this%rn(l) * Ql * (1.-del)
         this%Q(l, OUTSIDE, ir) = this%Q(l, OUTSIDE, ir) + this%irn(l) * Ql * (1.-del)
         if (del .notequals. zero) then
            this%Q(l, INSIDE,  ir+1) = this%Q(l, INSIDE,  ir+1) +  this%rn(l) * Ql * del
            this%Q(l, OUTSIDE, ir+1) = this%Q(l, OUTSIDE, ir+1) + this%irn(l) * Ql * del
         endif
         Ql2 = Ql1
         Ql1 = Ql
      enddo

      cfac = 1.
      sfac = 0.
      ! ph = atan2(y, x); cfac(m) = cos(m * ph); sfac(m) = sin(m * ph)

      ! non-axisymmetric (l,m) moments for 1 <= m <= this%mmax, m <= l <= this%lmax.
      do m = 1, this%mmax

         ! The (m,m) moment
         ! Associated Legendre polynomial: P_m^m = (-1)^m (2m-1)!! (1-x^2)^{m/2}
         ! The (2m-1)!! factor is integrated in ofact(:) array, where it mostly cancels out, note that (2m-1)!! \simeq m! exp(m/sqrt(2)) so it grows pretty fast with m
         Ql1 = sin_th ** m

         tmpfac = cfac
         cfac = cos_ph * cfac - sin_ph * sfac
         sfac = cos_ph * sfac + sin_ph * tmpfac

         m2c = this%lm(m, 2*m)  ! "cosine" coefficient
         ! The "sine" coefficient has index this%lm(m, 2*m-1) in our cmapping convention.
         ! Here we also exploit the fact that it equals m2c -1
         this%Q(m2c-1:m2c, INSIDE,  ir) = this%Q(m2c-1:m2c, INSIDE,  ir) +  this%rn(m) * Ql1 * (1.-del) * [ sfac, cfac ]
         this%Q(m2c-1:m2c, OUTSIDE, ir) = this%Q(m2c-1:m2c, OUTSIDE, ir) + this%irn(m) * Ql1 * (1.-del) * [ sfac, cfac ]
         if (del .notequals. zero) then
            this%Q(m2c-1:m2c, INSIDE,  ir+1) = this%Q(m2c-1:m2c, INSIDE,  ir+1) +  this%rn(m) * Ql1 * del * [ sfac, cfac ]
            this%Q(m2c-1:m2c, OUTSIDE, ir+1) = this%Q(m2c-1:m2c, OUTSIDE, ir+1) + this%irn(m) * Ql1 * del * [ sfac, cfac ]
         endif

         !>
         !! \deprecated BEWARE: most of computational cost of multipoles is here
         !! from (m+1,m) to (this%lmax,m)
         !! Associated Legendre polynomial: P_{m+1}^m = x (2m+1) P_m^m
         !! Associated Legendre polynomial recurrence: (l-m) P_l^m = x (2l-1) P_{l-1}^m - (l+m-1) P_{l-2}^m
         !<
         Ql2 = 0.
         do l = m + 1, this%lmax
            Ql = cos_th * this%k12(1, l, m) * Ql1 - this%k12(2, l, m) * Ql2
            m2c = m2c + 2  ! = this%lm(l, 2*m)
            ! Here we exploit that this%lm(l, 2*m) + 2 == this%lm(l+1, 2*m) for current implementation
            this%Q(m2c-1:m2c, INSIDE,  ir) = this%Q(m2c-1:m2c, INSIDE,  ir) +  this%rn(l) * Ql * (1.-del) * [ sfac, cfac ]
            this%Q(m2c-1:m2c, OUTSIDE, ir) = this%Q(m2c-1:m2c, OUTSIDE, ir) + this%irn(l) * Ql * (1.-del) * [ sfac, cfac ]
            if (del .notequals. zero) then
               this%Q(m2c-1:m2c, INSIDE,  ir+1) = this%Q(m2c-1:m2c, INSIDE,  ir+1) +  this%rn(l) * Ql * del * [ sfac, cfac ]
               this%Q(m2c-1:m2c, OUTSIDE, ir+1) = this%Q(m2c-1:m2c, OUTSIDE, ir+1) + this%irn(l) * Ql * del * [ sfac, cfac ]
            endif
            Ql2 = Ql1
            Ql1 = Ql
         enddo

      enddo

   end subroutine point2moments

!>
!! \brief Compute potential from multipole moments at a single point
!!
!! \todo improve accuracy with linear interpolation over radius
!<

   real function moments2pot(this, x, y, z) result(potential)

      use constants, only: zero
      use func,      only: operator(.notequals.)
      use units,     only: newtong

      implicit none

      class(mpole_container), intent(inout) :: this  !< object invoking type-bound procedure
      real,                   intent(in)    :: x     !< x coordinate of the contributing point
      real,                   intent(in)    :: y     !< y coordinate of the contributing point
      real,                   intent(in)    :: z     !< z coordinate of the contributing point

      real :: sin_th, cos_th, sin_ph, cos_ph, del, cfac, sfac, tmpfac
      real :: Ql, Ql1, Ql2
      integer :: l, m, ir, m2c

      call this%geomfac4moments(-newtong, x, y, z, sin_th, cos_th, sin_ph, cos_ph, ir, del)

      if (.not. interp_mom2pot) del = 0.

      ! monopole, the (0,0) moment; P_0 = 1.
      ! this%lm(0, 0) == 0
      potential = (1.-del) * ( &
           this%Q(0, INSIDE,  ir)   * this%irn(0) + &
           this%Q(0, OUTSIDE, ir+1) *  this%rn(0) )
      if (del .notequals. zero) potential = potential + del * ( &
           this%Q(0, INSIDE,  ir-1) * this%irn(0) + &
           this%Q(0, OUTSIDE, ir)   *  this%rn(0) )
      ! ir+1 to prevent duplicate accounting contributions from ir bin; alternatively one can modify radial integration

      ! axisymmetric (l,0) moments
      Ql2 = 0.
      Ql1 = 1.
      do l = 1, this%lmax
         ! this%lm(l, 0) == l
         Ql = cos_th * this%k12(1, l, 0) * Ql1 - this%k12(2, l, 0) * Ql2
         potential = potential + Ql * (1.-del) * ( &
              &      this%Q(l, INSIDE,  ir)   * this%irn(l) + &
              &      this%Q(l, OUTSIDE, ir+1) *  this%rn(l) )
         if (del .notequals. zero) potential = potential + Ql  * del * ( &
              &      this%Q(l, INSIDE,  ir-1) * this%irn(l) + &
              &      this%Q(l, OUTSIDE, ir)   *  this%rn(l) )
         Ql2 = Ql1
         Ql1 = Ql
      enddo

      cfac = 1.
      sfac = 0.

      ! non-axisymmetric (l,m) moments for 1 <= m <= this%mmax, m <= l <= this%lmax.
      do m = 1, this%mmax
         ! The (m,m) moment
         Ql1 = sin_th ** m
         tmpfac = cfac
         cfac = cos_ph * cfac - sin_ph * sfac
         sfac = cos_ph * sfac + sin_ph * tmpfac
         m2c = this%lm(m, 2*m)  ! see comments in point2moments
         potential = potential + Ql1 * (1.-del) * ( &
              &       this%irn(m) * (this%Q(m2c-1, INSIDE,  ir)   * sfac  + &
              &                      this%Q(m2c,   INSIDE,  ir)   * cfac) + &
              &        this%rn(m) * (this%Q(m2c-1, OUTSIDE, ir+1) * sfac  + &
              &                      this%Q(m2c,   OUTSIDE, ir+1) * cfac) )
         if (del .notequals. zero) potential = potential + Ql1 * del * ( &
              &       this%irn(m) * (this%Q(m2c-1, INSIDE,  ir-1) * sfac  + &
              &                      this%Q(m2c,   INSIDE,  ir-1) * cfac) + &
              &        this%rn(m) * (this%Q(m2c-1, OUTSIDE, ir)   * sfac  + &
              &                      this%Q(m2c,   OUTSIDE, ir)   * cfac) )

         !> \deprecated BEWARE: lots of computational cost of multipoles is here
         ! from (m+1,m) to (this%lmax,m)
         Ql2 = 0.
         do l = m+1, this%lmax
            m2c = m2c + 2  ! see comments in point2moments
            Ql = cos_th * this%k12(1, l, m) * Ql1 - this%k12(2, l, m) * Ql2
            potential = potential + Ql * (1.-del) * ( &
                 &       this%irn(l) * (this%Q(m2c-1, INSIDE,  ir)   * sfac  + &
                 &                      this%Q(m2c,   INSIDE,  ir)   * cfac) + &
                 &        this%rn(l) * (this%Q(m2c-1, OUTSIDE, ir+1) * sfac  + &
                 &                      this%Q(m2c,   OUTSIDE, ir+1) * cfac) )
            if (del .notequals. zero) potential = potential + Ql * del * ( &
                 &       this%irn(l) * (this%Q(m2c-1, INSIDE,  ir-1) * sfac  + &
                 &                      this%Q(m2c,   INSIDE,  ir-1) * cfac) + &
                 &        this%rn(l) * (this%Q(m2c-1, OUTSIDE, ir)   * sfac  + &
                 &                      this%Q(m2c,   OUTSIDE, ir)   * cfac) )
            Ql2 = Ql1
            Ql1 = Ql
         enddo
      enddo

   end function moments2pot

!>
!! \brief This routine calculates various geometrical numbers required for multipole evaluation
!!
!! \details It modifies the this%rn(:) and this%irn(:) arrays. Scalars are passed through argument list.
!<

   subroutine geomfac4moments(this, factor, x, y, z, sin_th, cos_th, sin_ph, cos_ph, ir, delta)

      use constants,  only: GEO_XYZ, GEO_RPZ, zero
      use dataio_pub, only: die, msg
      use domain,     only: dom
      use func,       only: operator(.notequals.)

      implicit none

      class(mpole_container), intent(inout) :: this  !< object invoking type-bound procedure
      real,                   intent(in)  :: factor  !< scaling factor (e.g. mass) of the contributing point
      real,                   intent(in)  :: x       !< x coordinate of the contributing point
      real,                   intent(in)  :: y       !< x coordinate of the contributing point
      real,                   intent(in)  :: z       !< x coordinate of the contributing point
      real,                   intent(out) :: sin_th  !< sine of the vertical angle
      real,                   intent(out) :: cos_th  !< cosine of the vertical angle
      real,                   intent(out) :: sin_ph  !< sine of the azimuthal angle
      real,                   intent(out) :: cos_ph  !< cosine of the azimuthal angle
      integer,                intent(out) :: ir      !< radial index for the this%Q(:, :, r) array
      real,                   intent(out) :: delta   !< fraction of the radial cell for interpolation between ir and ir+1

      real    :: rxy, r, rinv
      integer :: l

      ! radius and its projection onto XY plane
      select case (dom%geometry_type)
         case (GEO_XYZ)
            rxy = x**2 + y**2
         case (GEO_RPZ)
            rxy = x**2
         case default
            call die("[multipole_array:geomfac4moments] Unsupported geometry.")
            rxy = 0.
      end select
      r    = sqrt(rxy + z**2)
      rxy  = sqrt(rxy)
      if (r .notequals. zero) then
         rinv = 1. / r
      else
         rinv = 0.
      endif

      !radial index for the this%Q(:, :, r) array
      ir = int(r / this%drq)
      if (ir < this%rqbin) then
         delta = r/this%drq - ir
      else
         delta = 0
      endif
      if (ir > this%rqbin .or. ir < 0) then
         write(msg,'(2(a,i7),a)')"[multipole_array:geomfac4moments] radial index = ", ir, " outside this%Q(:, :, ", this%rqbin, ") range"
         call die(msg)
      endif
      this%irmax = max(this%irmax, ir)
      this%irmin = min(this%irmin, ir)

      ! azimuthal angle sine and cosine tables
      ! ph = atan2(y, x); this%cfac(m) = cos(m * ph); this%sfac(m) = sin(m * ph)
      ! this%cfac(0) and this%sfac(0) are set in init_multigrid
      if (rxy .notequals. zero) then
         select case (dom%geometry_type)
            case (GEO_XYZ)
               cos_ph = x / rxy
               sin_ph = y / rxy
            case (GEO_RPZ)
               cos_ph = cos(y)
               sin_ph = sin(y)
            case default
               call die("[multipole_array:geomfac4moments] Unsupported geometry.")
               cos_ph = 0. ; sin_ph = 0.
         end select
      else
         cos_ph = 1.
         sin_ph = 0.
      endif
!> \todo Possible optimization: number of computed elements can be doubled on each loop iteration (should give better pipelining)

      ! vertical angle
      cos_th = z   * rinv
      sin_th = rxy * rinv

!> \todo check how much it would degrade solution and improve performance to move the multiplications by this%rn(:) and this%irn(:) to img_mass2moments (before and after integration)
      ! this%rn(l) = factor * r ** l; this%irn(l) = factor * r ** -(l+1)
      this%rn(0)  = factor
      this%irn(0) = factor * rinv
      ! scale r to reduce risk of occurring a Floating Point Overflow for high l_max and large r values. Leave the original value of r only in this%irn(0)
      r = r / this%rscale
      rinv = rinv * this%rscale
!> \todo Possible optimization: number of computed elements can be doubled on each loop iteration (should give better pipelining)
      do l = 1, this%lmax
         this%rn(l)  =  this%rn(l-1) * r
         this%irn(l) = this%irn(l-1) * rinv
      enddo

   end subroutine geomfac4moments

end module multipole_array
