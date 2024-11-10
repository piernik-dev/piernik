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

!>
!! \brief Unified refinement criteria for filters that provide scalar indicator
!! of the need for refinement based on selected fields.
!<

module unified_ref_crit_var

   use constants,               only: INVALID, cbuff_len
   use unified_ref_crit_filter, only: urc_filter

   implicit none

   private
   public :: urc_var, decode_urcv

!>
!! \brief Things that should be common for all refinement criteria based on filters
!! that look for shockwaves or do other checks based on selected fields.
!<

   type, extends(urc_filter) :: urc_var
      private
      !unified_ref_crit_list:create_plotfields needs these components
      character(len=cbuff_len), public      :: rname         !< name of the refinement routine
      integer, public                       :: iv = INVALID  !< field index in cg%q or cg%w array
      integer, public                       :: ic = INVALID  !< component index (cg%w(iv)%arr(ic,:,:,:)) or INVALID for 3D arrays
      character(len=cbuff_len)              :: rvar          !< name of the refinement variable
      real                                  :: aux           !< auxiliary parameter (can be smoother or filter strength)
      procedure(refine_crit), pass, pointer :: refine        !< refinement routine
   contains
      procedure :: mark => mark_var
   end type urc_var

   interface urc_var
      procedure :: init
   end interface urc_var

   interface

      subroutine refine_crit(this, cg, p3d)

         use grid_cont, only: grid_container
         import urc_var

         implicit none

         class(urc_var),                  intent(in)    :: this  !<  an object invoking the type-bound procedure
         type(grid_container), pointer,   intent(inout) :: cg    !< current grid piece
         real, dimension(:,:,:), pointer, intent(in)    :: p3d   !< pointer to array to be examined for (de)refinement

      end subroutine refine_crit

   end interface

contains

!>
!! \brief Decode field-based refinement criteria from problem.par.
!! Return a pointer or chain of pointers if necessary
!<

   function decode_urcv(rf) result(this)

      use constants,  only: INVALID
      use dataio_pub, only: warn
      use mpisetup,   only: master
      use refinement, only: ref_auto_param

      implicit none

      type(ref_auto_param), intent(in) :: rf  !< the data read from problem.par

      type(urc_var), pointer :: this  !< a pointer to first of the constructed objects

      type(urc_var), pointer :: link
      integer(kind=4) :: iv
      integer(kind=4), dimension(:), allocatable :: ic
      integer :: i

      this => null()
      call identify_field(rf%rvar, iv, ic)

      if (iv == INVALID) then
         if (master) call warn("[URC_var:decode_urcv] ignoring '" // trim(rf%rvar) // "'")
         return
      endif

      if (allocated(ic)) then
         link => null()
         do i = lbound(ic, dim=1), ubound(ic, dim=1)
            allocate(this)
            this = init(rf, iv, ic(i))
            if (associated(link)) this%next => link
            link => this
         enddo
      else
         allocate(this)
         this = init(rf, iv, INVALID)
      endif

   end function decode_urcv

!> \brief Identify field name and return indices to cg%q or cg%w arrays

   subroutine identify_field(vname, iv, ic)

      use constants,        only: INVALID, cbuff_len
      use dataio_pub,       only: msg, warn
      use fluidindex,       only: iarr_all_dn, iarr_all_mx, iarr_all_my, iarr_all_mz, iarr_all_en
      use mpisetup,         only: master
      use named_array_list, only: qna, wna
      use refinement,       only: inactive_name

      implicit none

      character(len=cbuff_len),                   intent(in)  :: vname !< string specifying the field on
      integer(kind=4),                            intent(out) :: iv    !< field index in cg%q or cg%w array
      integer(kind=4), dimension(:), allocatable, intent(out) :: ic    !< component index array (cg%w(iv)%arr(ic,:,:,:)) or INVALID for 3D arrays

      iv = INVALID

      if (trim(vname) == trim(inactive_name)) return ! ignore this

      if (qna%exists(trim(vname))) then
         iv = qna%ind(trim(vname))
         return ! this is a 3d array name
      endif

      if (trim(vname) == "dens") then
         call alloc_ic(iarr_all_dn)
         return
      else if (trim(vname) == "velx") then
         call alloc_ic(iarr_all_mx)
         return
      else if (trim(vname) == "vely") then
         call alloc_ic(iarr_all_my)
         return
      else if (trim(vname) == "velz") then
         call alloc_ic(iarr_all_mz)
         return
      else if (trim(vname) == "ener") then
         call alloc_ic(iarr_all_en)
         return
      endif
      !> \todo identify here all {den,vl[xyz],ene}{d,n,i}
      !> \todo introduce possibility to operate on pressure or other indirect fields

      write(msg,'(3a)')"[URC_var:identify_field] Unidentified refinement variable: '",trim(vname),"'"
      if (master) call warn(msg)

   contains

      subroutine alloc_ic(tab)

         implicit none

         integer(kind=4), dimension(:), intent(in) :: tab

         iv = wna%fi

         allocate(ic(size(tab)))
         ic = tab

      end subroutine alloc_ic

   end subroutine identify_field

!> \brief A constructor for single scalar fields

   function init(rf, iv, ic) result(this)

      use constants,        only: V_VERBOSE, INVALID
      use dataio_pub,       only: printinfo, msg, die, warn
      use func,             only: operator(.notequals.)
      use mpisetup,         only: master
      use named_array_list, only: qna, wna
      use refinement,       only: ref_auto_param, inactive_name

      implicit none

      type(ref_auto_param), intent(in) :: rf  !< the data read from problem.par
      integer(kind=4),      intent(in) :: iv  !< index in qna or wna
      integer(kind=4),      intent(in) :: ic  !< sub index in wna or INVALID if qna

      type(urc_var) :: this  !< an object to be constructed

      if (master) then
         write(msg, '(5a,g13.5)')"[URC var]   Initializing refinement on variable '", trim(rf%rvar), "', method: '", trim(rf%rname), "', threshold = ", rf%ref_thr
         if (rf%aux .notequals. 0.) write(msg(len_trim(msg)+1:), '(a,g13.5)') ", with parameter = ", rf%aux
         if (rf%plotfield)  write(msg(len_trim(msg)+1:), '(a)') ", with plotfield"
         if (ic /= INVALID) then
            write(msg(len_trim(msg)+1:), '(a, i3,a,i3,a)') ", wna index: ", iv,"(", ic, ")"
         else
            write(msg(len_trim(msg)+1:), '(a, i3)') ", qna index: ", iv
         endif
         call printinfo(msg, V_VERBOSE)
         if (ic /= INVALID) then
            if (.not. wna%lst(iv)%vital) call warn("[URC_var:init] 4D field '" // trim(wna%lst(iv)%name) // "' is not vital. Please make sure that the guardcells are properly updater for refinement update.")
         else
            if (.not. qna%lst(iv)%vital) call warn("[URC_var:init] 3D field '" // trim(qna%lst(iv)%name) // "' is not vital. Please make sure that the guardcells are properly updater for refinement update.")
         endif
      endif

      ! urc_filter
      this%ref_thr   = rf%ref_thr
      this%plotfield = rf%plotfield

      ! own components
      this%iv = iv
      if (ic /= INVALID) this%ic = ic
      this%rvar  = rf%rvar
      this%rname = rf%rname
      this%aux   = rf%aux

!> \todo Implement Richardson extrapolation method, as described in M. Berger papers

      select case (trim(this%rname))
         case ("grad")
            this%refine => refine_on_gradient
         case ("relgrad")
            this%refine => refine_on_relative_gradient
         case ("Loechner", "second_order", "d2")
            this%refine => refine_on_second_derivative
         case (trim(inactive_name)) ! do nothing
         case default
            call die("[URC_var:init] unknown refinement detection routine")
      end select

   end function init

!> \brief Mark regions for refinement

   subroutine mark_var(this, cg)

      use constants, only: INVALID
      use grid_cont, only: grid_container

      implicit none

      class(urc_var),                intent(inout) :: this  !< an object invoking the type-bound procedure
      type(grid_container), pointer, intent(inout) :: cg    !< current grid piece

      real, dimension(:,:,:), pointer :: p3d

      if (this%ic == INVALID) then
         p3d => cg%q(this%iv)%arr
      else
         associate (a => cg%w(this%iv)%arr)
            p3d(lbound(a, dim=2):, lbound(a, dim=3):, lbound(a, dim=4):) => cg%w(this%iv)%arr(this%ic, :, :, :)
         end associate
      endif
      call this%refine(cg, p3d)

   end subroutine mark_var

!>
!! \brief R. Loechner criterion
!! Original paper: https://www.researchgate.net/publication/222452974_An_adaptive_finite_element_scheme_for_transient_problems_in_CFD
!! Cartesian grid implementation:
!! http://flash.uchicago.edu/~jbgallag/2012/flash4_ug/node14.html#SECTION05163100000000000000
!! (note that some indices in the left part of denominator seem to be slightly messed up)
!!
!! this%aux is the noise filter (epsilon in Loechner's paper)
!!
!! Known limitations:
!! * Extremely small and big cell sizes easily cause FP over- and underflows. Don't go much beyond 1e±70.
!!
!! OPT: this routine seems to take way too much CPU! Implement 3D variant without functions?
!!
!<

   subroutine refine_on_second_derivative(this, cg, p3d)

      use constants,  only: xdim, ydim, zdim, GEO_XYZ, INVALID, I_ONE, PPP_AMR, PPP_CG
      use dataio_pub, only: die, msg
      use domain,     only: dom
      use grid_cont,  only: grid_container
      use ppp,        only: ppp_main

      implicit none

      class(urc_var),                  intent(in)    :: this !< this contains refinement parameters
      type(grid_container), pointer,   intent(inout) :: cg   !< current grid piece
      real, dimension(:,:,:), pointer, intent(in)    :: p3d  !< pointer to array to be examined for (de)refinement needs (should contain at least two layers of up-to-date guardcells)

      integer(kind=4) :: i, j, k
      real :: sn, sd, r
      integer(kind=4), parameter :: how_far = 2
      character(len=*), parameter :: L_label = "Loechner_mark"

      call ppp_main%start(L_label, PPP_AMR + PPP_CG)
      if (dom%geometry_type /= GEO_XYZ) call die("[URC_var:refine_on_second_derivative] noncartesian geometry not supported yet")
      if (dom%nb <= how_far+I_ONE) then
         write(msg, '(a,i1,a)')"[URC_var:refine_on_second_derivative] at least ", how_far+I_ONE+I_ONE, " guardcells are required"
         ! dux(i+I_ONE, j, k) is accessing p3d(i+I_ONE+I_ONE, j, k), so for i = cg%ie + how_far*dom%D_x means that dom%nb has to be at least 4
         ! ToDo: check if it is safe to reduce how_far and avoid refinement flickering
         call die(msg)
      endif

      ! Perhaps it will be a bit faster with arrays for storing first-order differences
      ! but let's see if it works first and then how expensive it is.

      do k = cg%ks - how_far*dom%D_z, cg%ke + how_far*dom%D_z
         do j = cg%js - how_far*dom%D_y, cg%je + how_far*dom%D_y
            do i = cg%is - how_far*dom%D_x, cg%ie + how_far*dom%D_x

               sn = 0.
               sd = 0.

               if (dom%has_dir(xdim)) then
                  ! d/dx d/dx
                  sn = sn + ((dux(i+I_ONE, j, k) - dux(i-I_ONE, j, k)) * cg%idx)**2
                  sd = sd + ((abs(dux(i+I_ONE, j, k)) + abs(dux(i-I_ONE, j, k)) + this%aux * (daux(i+I_ONE, j, k) + daux(i-I_ONE, j, k))) * cg%idx)**2
                  if (dom%has_dir(ydim)) then
                     ! d/dy d/dx
                     sn = sn + 2 * ((dux(i, j+I_ONE, k) - dux(i, j-I_ONE, k)) * cg%idy)**2  ! == d/dx d/dy
                     sd = sd + ((abs(dux(i, j+I_ONE, k)) + abs(dux(i, j-I_ONE, k)) + this%aux * (daux(i, j+I_ONE, k) + daux(i, j-I_ONE, k))) * cg%idy)**2 + &
                          &    ((abs(duy(i+I_ONE, j, k)) + abs(duy(i-I_ONE, j, k)) + this%aux * (dauy(i+I_ONE, j, k) + dauy(i-I_ONE, j, k))) * cg%idx)**2     ! d/dx d/dy
                     ! to be exploited: (dauy(i+I_ONE, j, k) + dauy(i-I_ONE, j, k)) * cg%idx == (daux(i, j+I_ONE, k) + daux(i, j-I_ONE, k)) * cg%idy
                  endif
                  if (dom%has_dir(zdim)) then
                     ! d/dz d/dx
                     sn = sn + 2 * ((dux(i, j, k+I_ONE) - dux(i, j, k-I_ONE)) * cg%idz)**2  ! == d/dx d/dz
                     sd = sd + ((abs(dux(i, j, k+I_ONE)) + abs(dux(i, j, k-I_ONE)) + this%aux * (daux(i, j, k+I_ONE) + daux(i, j, k-I_ONE))) * cg%idz)**2 + &
                          &    ((abs(duz(i+I_ONE, j, k)) + abs(duz(i-I_ONE, j, k)) + this%aux * (dauz(i+I_ONE, j, k) + dauz(i-I_ONE, j, k))) * cg%idx)**2     ! d/dx d/dz
                     ! to be exploited: dauz(i+I_ONE, j, k) + dauz(i-I_ONE, j, k)) * cg%idx == (daux(i, j, k+I_ONE) + daux(i, j, k-I_ONE)) * cg%idz
                  endif
               endif

               if (dom%has_dir(ydim)) then
                  ! d/dy d/dy
                  sn = sn + ((duy(i, j+I_ONE, k) - duy(i, j-I_ONE, k)) * cg%idy)**2
                  sd = sd + ((abs(duy(i, j+I_ONE, k)) + abs(duy(i, j-I_ONE, k)) + this%aux * (dauy(i, j+I_ONE, k) + dauy(i, j-I_ONE, k))) * cg%idy)**2
                  if (dom%has_dir(zdim)) then
                     ! d/dz d/dy
                     sn = sn + 2 * ((duy(i, j, k+I_ONE) - duy(i, j, k-I_ONE)) * cg%idz)**2  ! == d/dy d/dz
                     sd = sd + ((abs(duy(i, j, k+I_ONE)) + abs(duy(i, j, k-I_ONE)) + this%aux * (dauy(i, j, k+I_ONE) + dauy(i, j, k-I_ONE))) * cg%idz)**2 + &
                          &    ((abs(duz(i, j+I_ONE, k)) + abs(duz(i, j-I_ONE, k)) + this%aux * (dauz(i, j+I_ONE, k) + dauz(i, j-I_ONE, k))) * cg%idy)**2     ! d/dy d/dz
                     ! to be exploited: (dauz(i, j+I_ONE, k) + dauz(i, j-I_ONE, k)) * cg%idy == (dauy(i, j, k+I_ONE) + dauy(i, j, k-I_ONE)) * cg%idz
                  endif
               endif

               if (dom%has_dir(zdim)) then
                  ! d/dz d/dz
                  sn = sn + ((duz(i, j, k+I_ONE) - duz(i, j, k-I_ONE)) * cg%idz)**2
                  sd = sd + ((abs(duz(i, j, k+I_ONE)) + abs(duz(i, j, k-I_ONE)) + this%aux * (dauz(i, j, k+I_ONE) + dauz(i, j, k-I_ONE))) * cg%idz)**2
               endif

               if (sd > 0.) then
                  r = sn / sd
               else  ! sd == 0 because it should never be < 0.
                  if (sn > 0.) then
                     r = 1.  ! strange case
                  else
                     r = 0.  ! most likely constant field == 0.
                  endif
               endif
               if (this%iplot /= INVALID) cg%q(this%iplot)%arr(i, j, k) = r

               if (r >= this%ref_thr) call cg%flag%set(i, j, k)

            enddo
         enddo
      enddo
      call ppp_main%stop(L_label, PPP_AMR + PPP_CG)

   contains

      elemental real function dux(i, j, k)
         use constants, only: I_ONE
         implicit none
         integer(kind=4), intent(in) :: i, j, k !< (x, y, z)-indices
         dux = (p3d(i+I_ONE, j, k) - p3d(i-I_ONE, j, k)) * cg%idx
      end function dux

      elemental real function daux(i, j, k)
         use constants, only: I_ONE
         implicit none
         integer(kind=4), intent(in) :: i, j, k !< (x, y, z)-indices
         daux = (abs(p3d(i+I_ONE, j, k)) + abs(p3d(i-I_ONE, j, k))) * cg%idx
      end function daux

      elemental real function duy(i, j, k)
         use constants, only: I_ONE
         implicit none
         integer(kind=4), intent(in) :: i, j, k !< (x, y, z)-indices
         duy = (p3d(i, j+I_ONE, k) - p3d(i, j-I_ONE, k)) * cg%idy
      end function duy

      elemental real function dauy(i, j, k)
         use constants, only: I_ONE
         implicit none
         integer(kind=4), intent(in) :: i, j, k !< (x, y, z)-indices
         dauy = (abs(p3d(i, j+I_ONE, k)) + abs(p3d(i, j-I_ONE, k))) * cg%idy
      end function dauy

      elemental real function duz(i, j, k)
         use constants, only: I_ONE
         implicit none
         integer(kind=4), intent(in) :: i, j, k !< (x, y, z)-indices
         duz = (p3d(i, j, k+I_ONE) - p3d(i, j, k-I_ONE)) * cg%idz
      end function duz

      elemental real function dauz(i, j, k)
         use constants, only: I_ONE
         implicit none
         integer(kind=4), intent(in) :: i, j, k !< (x, y, z)-indices
         dauz = (abs(p3d(i, j, k+I_ONE)) + abs(p3d(i, j, k-I_ONE))) * cg%idz
      end function dauz

   end subroutine refine_on_second_derivative

!>
!! \brief Refinement based on ||grad u||
!! This is sensitive to gradients, but the thresholds must be rescaled, when you change units of the problem.
!<

   subroutine refine_on_gradient(this, cg, p3d)

      use constants, only: INVALID
      use grid_cont, only: grid_container

      implicit none

      class(urc_var),                  intent(in)    :: this !< this contains refinement parameters
      type(grid_container), pointer,   intent(inout) :: cg   !< current grid piece
      real, dimension(:,:,:), pointer, intent(in)    :: p3d  !< pointer to array to be examined for (de)refinement needs (should contain at least one layer of updated guardcells)

      integer(kind=4) :: i, j, k
      real :: r

      !> \todo implement how far we should look for refinements

      do k = cg%ks, cg%ke
         do j = cg%js, cg%je
            do i = cg%is, cg%ie
               r = grad2(i, j, k)
               if (this%iplot /= INVALID) cg%q(this%iplot)%arr(i, j, k) = r
               if (r >= this%ref_thr**2) call cg%flag%set(i, j, k)
               ! we can avoid calculating square root here
            enddo
         enddo
      enddo

   contains

      !> \brief scan given area of indices for maximum value of specified function

      elemental real function maxgradoverarea(i1, i2, j1, j2, k1, k2)

         implicit none

         integer(kind=4), intent(in) :: i1, i2, j1, j2, k1, k2

         integer(kind=4) :: i, j, k

         maxgradoverarea = -huge(1.)

         do k = k1, k2
            do j = j1, j2
               do i = i1, i2
                  maxgradoverarea = max(maxgradoverarea, grad2(i, j, k))
               enddo
            enddo
         enddo

      end function maxgradoverarea

      !>
      !! \brief Square of a gradient without taking into account cell sizes
      !!
      !! \details Perhaps it is not gradient but rather a norm of 3D difference
      !<

      elemental real function grad2(i, j, k)

         use domain, only: dom

         implicit none

         integer(kind=4), intent(in) :: i !< x-index
         integer(kind=4), intent(in) :: j !< y-index
         integer(kind=4), intent(in) :: k !< z-index

         grad2 = (p3d(i+dom%D_x, j, k) - p3d(i-dom%D_x, j, k))**2 + &
              &  (p3d(i, j+dom%D_y, k) - p3d(i, j-dom%D_y, k))**2 + &
              &  (p3d(i, j, k+dom%D_z) - p3d(i, j, k-dom%D_z))**2

      end function grad2

   end subroutine refine_on_gradient

!>
!! \brief Refinement based on ||grad u||/||u||
!! This is sensitive to changes of sign or strong gradients (such as fast approaching 0.)
!<

   subroutine refine_on_relative_gradient(this, cg, p3d)

      use constants, only: INVALID
      use grid_cont, only: grid_container

      implicit none

      class(urc_var),                  intent(in)    :: this !< this contains refinement parameters
      type(grid_container), pointer,   intent(inout) :: cg   !< current grid piece
      real, dimension(:,:,:), pointer, intent(in)    :: p3d  !< pointer to array to be examined for (de)refinement needs (should contain at least one layer of updated guardcells)

      integer(kind=4) :: i, j, k
      real :: r

      !> \todo implement how far we should look for (de)refinements

      do k = cg%ks, cg%ke
         do j = cg%js, cg%je
            do i = cg%is, cg%ie
               r = rel_grad2(i, j, k)
               if (this%iplot /= INVALID) cg%q(this%iplot)%arr(i, j, k) = r
               if (r >= this%ref_thr**2) call cg%flag%set(i, j, k)
               ! we can avoid calculating square root here
            enddo
         enddo
      enddo

   contains

      !> \brief scan given area of indices for maximum value of specified function

      elemental real function maxrelgradoverarea(i1, i2, j1, j2, k1, k2)

         implicit none

         integer(kind=4), intent(in) :: i1, i2, j1, j2, k1, k2

         integer(kind=4) :: i, j, k

         maxrelgradoverarea = -huge(1.)

         do k = k1, k2
            do j = j1, j2
               do i = i1, i2
                  maxrelgradoverarea = max(maxrelgradoverarea, rel_grad2(i, j, k))
               enddo
            enddo
         enddo

      end function maxrelgradoverarea

      !>
      !! \brief Square of a 'relative gradient' without taking into account cell sizes
      !!
      !! \details This routine should return value between 0 and 1.
      !<

      elemental real function rel_grad2(i, j, k)

         use constants, only: xdim, ydim, zdim
         use domain,    only: dom

         implicit none

         integer(kind=4), intent(in) :: i !< x-index
         integer(kind=4), intent(in) :: j !< y-index
         integer(kind=4), intent(in) :: k !< z-index

         rel_grad2 = 0.

         if (dom%has_dir(xdim)) rel_grad2 = sumsq01( &
              rel_grad_1pair(p3d(i, j, k), p3d(i+dom%D_x, j, k)), &
              rel_grad_1pair(p3d(i, j, k), p3d(i-dom%D_x, j, k)) )

         if (dom%has_dir(ydim)) rel_grad2 = sumsq01( rel_grad2, sumsq01( &
              rel_grad_1pair(p3d(i, j, k), p3d(i, j+dom%D_y, k)), &
              rel_grad_1pair(p3d(i, j, k), p3d(i, j-dom%D_y, k)) ) )

         if (dom%has_dir(zdim)) rel_grad2 = sumsq01( rel_grad2, sumsq01( &
              rel_grad_1pair(p3d(i, j, k), p3d(i, j, k+dom%D_z)), &
              rel_grad_1pair(p3d(i, j, k), p3d(i, j, k-dom%D_z)) ) )

      end function rel_grad2

      !>
      !! \brief Return 'relative gradient' defined as |a-b|/(|a| + |b|)
      !!
      !! \details This routine should return value between 0 and 1.
      !! rel_grad_1pair == 0 when a == b, also when a == b == 0
      !! rel_grad_1pair == 1 when one argument == 0 and the other /= 0
      !! rel_grad_1pair == 1 when a*b < 0
      !<

      elemental real function rel_grad_1pair(a, b)

         use func, only: operator(.notequals.)

         implicit none

         real, intent(in) :: a, b

         if ((abs(a) .notequals. 0.) .and. (abs(b) .notequals. 0.)) then
            rel_grad_1pair = abs(a-b) / (abs(a) + abs(b))
         else
            rel_grad_1pair = 0.
         endif

      end function rel_grad_1pair

      !>
      !! \brief square of an 'addition' of two numbers from [0, 1] range.
      !!
      !! \details Let us define a+b so it obeys the following rules:
      !! a (+) b = b (+) a
      !! a (+) 0 = a
      !! a (+) 1 = 1
      !! if a << 1 and b << 1, a (+) b \simeq sqrt(a**2+ b**2)
      !!
      !! One of possible formulas is:
      !! a (+) b = sqrt(a**2 - a**2 * b**2 + b**2)
      !! Note that: (a (+) b) (+) c = a (+) (b (+) c) = sqrt( (a (+) b)**2 * (1 - c**2) + c**2)
      !! That's why this routine returns (a (+) b)**2 for best performance when 'adding' multiple numbers.
      !<

      elemental real function sumsq01(a, b)

         implicit none

         real, intent(in) :: a, b

         sumsq01 = a**2 * (1-b**2) + b**2

      end function sumsq01

   end subroutine refine_on_relative_gradient

end module unified_ref_crit_var
