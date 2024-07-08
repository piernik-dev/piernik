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

!> \brief This module contains variables and initialization routines related to refinement

module refinement

   use constants, only: ndims, LO, HI, cbuff_len

   implicit none

   private
   public :: n_updAMR, oop_thr, ref_point, refine_points, ref_auto_param, refine_vars, level_min, level_max, inactive_name, bsize, &
        &    ref_box, refine_boxes, refine_zcyls, init_refinement, emergency_fix, set_n_updAMR, prefer_n_bruteforce, jeans_ref, jeans_plot, &
        &    nbody_ref, updAMR_after

   integer(kind=4), protected :: n_updAMR            !< How often to update the refinement structure, 0 to disable it even at IC
   real,            protected :: oop_thr             !< Maximum allowed ratio of Out-of-Place grid pieces (according to current ordering scheme)
   logical,         protected :: prefer_n_bruteforce !< If .false. then try SFC algorithms for neighbor searches
   integer(kind=4), protected :: level_min           !< Minimum allowed refinement, base level by default.
   integer(kind=4), protected :: level_max           !< Maximum allowed refinement, if set > 0 then AMR is requested (don't need to be reached if not necessary)
   integer(kind=4), dimension(ndims), protected :: bsize  !< The size of cg for multiblock decomposition, if not set explicitly, then a heuristic values will be used, if possible.
   logical,         protected :: auto_bsize          !< Allows to compute bsize automatically when bsize is not explicitly provided

   ! some refinement primitives
   integer, parameter :: nshapes = 10 !< number of shapes of each kind allowed to be predefined by user in problem.par

   integer, parameter :: n_upd_steps = 10 !< number of entries in updAMR_after
   integer(kind=4), dimension(n_upd_steps), protected :: updAMR_after

   !> \brief Refinement point
   type :: ref_point
      integer(kind=4)        :: level  !< desired level of refinement
      real, dimension(ndims) :: coords !< coordinates, where to refine
   end type ref_point
   type(ref_point), dimension(nshapes), protected :: refine_points !< Points of refinement to be used from problem.par: level and (x-y-z)-coordinates

   !> \brief Refinement box
   type :: ref_box
      integer(kind=4)               :: level  !< desired level of refinement
      real, dimension(ndims, LO:HI) :: coords !< coordinates, where to refine
   end type ref_box
   type(ref_box), dimension(nshapes), protected :: refine_boxes !< Areas (boxes) of refinement to be used from problem.par: level, (x-y-z)-coordinates of lower left corner and (x-y-z)-coordinates of upper right corner
   type(ref_box), dimension(nshapes), protected :: refine_zcyls !< z-axis cylinders that fit the specified boxes

   !> \brief Parameters of automagic refinement
   type :: ref_auto_param
      character(len=cbuff_len) :: rvar  !< name of the refinement variable
      character(len=cbuff_len) :: rname !< name of the refinement routine
      real :: ref_thr                   !< refinement threshold
      real :: aux                       !< auxiliary parameter (can be smoother or filter strength)
      logical :: plotfield              !< create a 3D array to keep the value of refinement criterion when set to .true.
   end type ref_auto_param
   integer, parameter :: n_ref_auto_param = 10                                 !< number of automatic refinement criteria available to user
   type(ref_auto_param), dimension(n_ref_auto_param), protected :: refine_vars !< Definitions of user-supplied automatic refinement criteria: refinement variable, refinement algorithm, refinement threshold, derefinement threshold, auxiliary parameter

   ! \brief Parameters of Jeans length based refinement
   real    :: jeans_ref   !< minimum resolution in cells per Jeans wavelengths
   logical :: jeans_plot  !<create a 3D array to keep the value of Jeans resolution

   ! \brief Parameter for particle-based refinement
   integer(kind=4) :: nbody_ref  !< maximum allowed number of particles per cg

   character(len=cbuff_len), parameter :: inactive_name = "none"               !< placeholder for inactive refinement criterion

   logical :: emergency_fix                                                    !< set to .true. if you want to call update_refinement ASAP

   namelist /AMR/ level_min, level_max, bsize, auto_bsize, n_updAMR, updAMR_after, prefer_n_bruteforce, oop_thr, &
        &         refine_points, refine_boxes, refine_zcyls, refine_vars, jeans_ref, jeans_plot, &
        &         nbody_ref

contains

!>
!! \brief Initialization of parameters of refinement mechanics
!!
!! @b AMR
!! \n \n
!! <table border="+1">
!!   <tr><td> bsize(3)            </td><td> 0       </td><td> integer          </td><td> \copydoc refinement::bsize               </td></tr>
!!   <tr><td> auto_bsize          </td><td> .false. </td><td> logical          </td><td> \copydoc refinement::auto_bsize          </td></tr>
!!   <tr><td> level_min           </td><td> 0       </td><td> integer          </td><td> \copydoc refinement::level_min           </td></tr>
!!   <tr><td> level_max           </td><td> 0       </td><td> integer          </td><td> \copydoc refinement::level_max           </td></tr>
!!   <tr><td> n_updAMR            </td><td> HUGE    </td><td> integer          </td><td> \copydoc refinement::n_updAMR            </td></tr>
!!   <tr><td> updAMR_after        </td><td> 0       </td><td> integer(10)      </td><td> \copydoc refinement::updAMR_after        </td></tr>
!!   <tr><td> oop_thr             </td><td> 0.1     </td><td> real             </td><td> \copydoc refinement::oop_thr             </td></tr>
!!   <tr><td> refine_points(10)   </td><td> none    </td><td> integer, 3*real  </td><td> \copydoc refinement::refine_points       </td></tr>
!!   <tr><td> refine_boxes(10)    </td><td> none    </td><td> integer, 6*real  </td><td> \copydoc refinement::refine_boxes        </td></tr>
!!   <tr><td> refine_zcyls(10)    </td><td> none    </td><td> integer, 6*real  </td><td> \copydoc refinement::refine_zcyls        </td></tr>
!!   <tr><td> refine_vars(10)     </td><td> none    </td><td> 2*string, 3*real </td><td> \copydoc refinement::refine_vars         </td></tr>
!!   <tr><td> jeans_ref           </td><td> 0.      </td><td> real             </td><td> \copydoc refinement::jeans_ref           </td></tr>
!!   <tr><td> jeans_plot          </td><td> .false. </td><td> logical          </td><td> \copydoc refinement::jeans_plot          </td></tr>
!!   <tr><td> nbody_ref           </td><td> INVALID </td><td> integer          </td><td> \copydoc refinement::nbody_ref           </td></tr>
!!   <tr><td> prefer_n_bruteforce </td><td> .false. </td><td> logical          </td><td> \copydoc refinement::prefer_n_bruteforce </td></tr>
!! </table>
!! \n \n
!<
   subroutine init_refinement

      use constants,  only: base_level_id, PIERNIK_INIT_DOMAIN, xdim, ydim, zdim, I_ZERO, I_ONE, LO, HI, cbuff_len, refinement_factor, INVALID
      use dataio_pub, only: die, code_progress, warn, msg, printinfo, nh
      use domain,     only: dom
      use mpisetup,   only: cbuff, ibuff, lbuff, rbuff, master, slave, piernik_MPI_Bcast
      use user_hooks, only: problem_domain_update

      implicit none

      integer :: d, level_crit
      logical :: do_refine
      integer(kind=4), save :: level_insane = 64  ! That's absolute limit for 1D and current implementation of ordering::Morton_id based on 8-bit integers

      if (code_progress < PIERNIK_INIT_DOMAIN) call die("[refinement:init_refinement] Domain not initialized.")

      level_min = base_level_id
      level_max = level_min
      bsize(:)  = I_ZERO
      auto_bsize = associated(problem_domain_update)  ! allow by default for problems with domain expansion configured
      n_updAMR  = huge(I_ONE)
      prefer_n_bruteforce = .false.
      oop_thr = 0.1
      refine_points(:) = ref_point(base_level_id-1, [ 0., 0., 0.] )
      refine_boxes (:) = ref_box  (base_level_id-1, reshape([ 0., 0., 0., 0., 0., 0.], [ndims, HI-LO+I_ONE] ) )
      refine_zcyls(:)  = refine_boxes (:)
      refine_vars  (:) = ref_auto_param (inactive_name, inactive_name, 0., 0., .false.)
      jeans_ref = 0.       !< inactive by default, 4. is the absolute minimum for reasonable use
      jeans_plot = .false.
      nbody_ref = INVALID
      updAMR_after = I_ZERO

      if (2*n_ref_auto_param              > ubound(cbuff, dim=1)) call die("[refinement:init_refinement] increase cbuff size")
      if (10+3*nshapes                    > ubound(ibuff, dim=1)) call die("[refinement:init_refinement] increase ibuff size")
      if (2+15*nshapes+2*n_ref_auto_param > ubound(rbuff, dim=1)) call die("[refinement:init_refinement] increase rbuff size")

      if (master) then

         if (.not.nh%initialized) call nh%init()
         open(newunit=nh%lun, file=nh%tmp1, status="unknown")
         write(nh%lun,nml=AMR)
         close(nh%lun)
         open(newunit=nh%lun, file=nh%par_file)
         nh%errstr=""
         read(unit=nh%lun, nml=AMR, iostat=nh%ierrh, iomsg=nh%errstr)
         close(nh%lun)
         call nh%namelist_errh(nh%ierrh, "AMR")
         read(nh%cmdl_nml,nml=AMR, iostat=nh%ierrh)
         call nh%namelist_errh(nh%ierrh, "AMR", .true.)
         open(newunit=nh%lun, file=nh%tmp2, status="unknown")
         write(nh%lun,nml=AMR)
         close(nh%lun)
         call nh%compare_namelist()

         if (any(bsize(:) > 0 .and. bsize(:) < dom%nb .and. dom%has_dir(:))) call die("[refinement:init_refinement] bsize(:) is too small.")

         ! minimal sanitizing
         level_min = max(level_min, base_level_id)
         level_max = max(level_max, level_min)

         cbuff(1                 :  n_ref_auto_param) = refine_vars(:)%rvar
         cbuff(1+n_ref_auto_param:2*n_ref_auto_param) = refine_vars(:)%rname

         ibuff(1) = level_min
         ibuff(2) = level_max
         ibuff(3) = n_updAMR
         ibuff(4) = nbody_ref
         ibuff(5:4+ndims) = bsize
         ! it is safe because 4+ndims < 11
         ibuff(11          :10+  nshapes) = refine_points(:)%level
         ibuff(11+  nshapes:10+2*nshapes) = refine_boxes (:)%level
         ibuff(11+2*nshapes:10+3*nshapes) = refine_zcyls (:)%level
         ibuff(11+3*nshapes:10+3*nshapes+n_upd_steps) = updAMR_after

         lbuff(1) = jeans_plot
         lbuff(2) = prefer_n_bruteforce
         lbuff(3) = auto_bsize
         lbuff(4:3+n_ref_auto_param) = refine_vars(:)%plotfield

         rbuff(1) = oop_thr
         rbuff(2) = jeans_ref
         rbuff(3          :2+  nshapes) = refine_points(:)%coords(xdim)
         rbuff(3+  nshapes:2+2*nshapes) = refine_points(:)%coords(ydim)
         rbuff(3+2*nshapes:2+3*nshapes) = refine_points(:)%coords(zdim)

         rbuff(3+3*nshapes:2+4*nshapes) = refine_boxes (:)%coords(xdim, LO)
         rbuff(3+4*nshapes:2+5*nshapes) = refine_boxes (:)%coords(xdim, HI)
         rbuff(3+5*nshapes:2+6*nshapes) = refine_boxes (:)%coords(ydim, LO)
         rbuff(3+6*nshapes:2+7*nshapes) = refine_boxes (:)%coords(ydim, HI)
         rbuff(3+7*nshapes:2+8*nshapes) = refine_boxes (:)%coords(zdim, LO)
         rbuff(3+8*nshapes:2+9*nshapes) = refine_boxes (:)%coords(zdim, HI)

         rbuff(3+ 9*nshapes:2+10*nshapes) = refine_zcyls (:)%coords(xdim, LO)
         rbuff(3+10*nshapes:2+11*nshapes) = refine_zcyls (:)%coords(xdim, HI)
         rbuff(3+11*nshapes:2+12*nshapes) = refine_zcyls (:)%coords(ydim, LO)
         rbuff(3+12*nshapes:2+13*nshapes) = refine_zcyls (:)%coords(ydim, HI)
         rbuff(3+13*nshapes:2+14*nshapes) = refine_zcyls (:)%coords(zdim, LO)
         rbuff(3+14*nshapes:2+15*nshapes) = refine_zcyls (:)%coords(zdim, HI)

         rbuff(3+15*nshapes                   :2+15*nshapes+  n_ref_auto_param) = refine_vars(:)%ref_thr
         rbuff(3+15*nshapes+  n_ref_auto_param:2+15*nshapes+2*n_ref_auto_param) = refine_vars(:)%aux

      endif

      call piernik_MPI_Bcast(cbuff, cbuff_len)
      call piernik_MPI_Bcast(ibuff)
      call piernik_MPI_Bcast(lbuff)
      call piernik_MPI_Bcast(rbuff)

      if (slave) then

         refine_vars(:)%rvar  = cbuff(1                 :  n_ref_auto_param)
         refine_vars(:)%rname = cbuff(1+n_ref_auto_param:2*n_ref_auto_param)

         level_min = ibuff(1)
         level_max = ibuff(2)
         n_updAMR  = ibuff(3)
         nbody_ref = ibuff(4)
         bsize     = ibuff(5:4+ndims)
         refine_points(:)%level = ibuff(11          :10+  nshapes)
         refine_boxes (:)%level = ibuff(11+  nshapes:10+2*nshapes)
         refine_zcyls (:)%level = ibuff(11+2*nshapes:10+3*nshapes)
         updAMR_after           = ibuff(11+3*nshapes:10+3*nshapes+n_upd_steps)

         jeans_plot               = lbuff(1)
         prefer_n_bruteforce      = lbuff(2)
         auto_bsize               = lbuff(3)
         refine_vars(:)%plotfield = lbuff(4:3+n_ref_auto_param)

         oop_thr   = rbuff(1)
         jeans_ref = rbuff(2)
         refine_points(:)%coords(xdim)     = rbuff(3          :2+  nshapes)
         refine_points(:)%coords(ydim)     = rbuff(3+  nshapes:2+2*nshapes)
         refine_points(:)%coords(zdim)     = rbuff(3+2*nshapes:2+3*nshapes)

         refine_boxes (:)%coords(xdim, LO) = rbuff(3+3*nshapes:2+4*nshapes)
         refine_boxes (:)%coords(xdim, HI) = rbuff(3+4*nshapes:2+5*nshapes)
         refine_boxes (:)%coords(ydim, LO) = rbuff(3+5*nshapes:2+6*nshapes)
         refine_boxes (:)%coords(ydim, HI) = rbuff(3+6*nshapes:2+7*nshapes)
         refine_boxes (:)%coords(zdim, LO) = rbuff(3+7*nshapes:2+8*nshapes)
         refine_boxes (:)%coords(zdim, HI) = rbuff(3+8*nshapes:2+9*nshapes)

         refine_zcyls (:)%coords(xdim, LO) = rbuff(3+ 9*nshapes:2+10*nshapes)
         refine_zcyls (:)%coords(xdim, HI) = rbuff(3+10*nshapes:2+11*nshapes)
         refine_zcyls (:)%coords(ydim, LO) = rbuff(3+11*nshapes:2+12*nshapes)
         refine_zcyls (:)%coords(ydim, HI) = rbuff(3+12*nshapes:2+13*nshapes)
         refine_zcyls (:)%coords(zdim, LO) = rbuff(3+13*nshapes:2+14*nshapes)
         refine_zcyls (:)%coords(zdim, HI) = rbuff(3+14*nshapes:2+15*nshapes)

         refine_vars  (:)%ref_thr          = rbuff(3+15*nshapes                   :2+15*nshapes+  n_ref_auto_param)
         refine_vars  (:)%aux              = rbuff(3+15*nshapes+  n_ref_auto_param:2+15*nshapes+2*n_ref_auto_param)

      endif

      emergency_fix = .false.

      do_refine = (level_max > base_level_id) .or. all((bsize /= I_ZERO) .or. .not. dom%has_dir)

      if ((do_refine .or. auto_bsize) .and. all(bsize == I_ZERO)) call automagic_bsize

      where (.not. dom%has_dir) bsize = I_ONE

      ! If bsize was set then check if it is sane and fail if it is wrong
      do d = xdim, zdim
         if (dom%has_dir(d)) then
            if (bsize(d) < dom%nb) then
               if (do_refine .and. master) then
                  if (bsize(d) > 1) then
                     call die("[refinement:init_refinement] Refinements disabled bsize smaller than nb")
                  else
                     call warn("[refinement:init_refinement] any(bsize == 1) disables refinement")
                  endif
               endif
               do_refine = .false.
            else
               if (mod(dom%n_d(d), bsize(d)) /= 0) then
                  if (do_refine .and. master) call die("[refinement:init_refinement] Refinements disabled (domain not divisible by bsize)")
                  do_refine = .false.
               endif
            endif
            if (mod(dom%n_d(d), refinement_factor) /= I_ZERO) then
               if (do_refine .and. master) call warn("[refinement:init_refinement] Refinements disabled (domain not divisible by refinement factor)")
               do_refine = .false.
            endif
         endif
      enddo

      if (any(dom%has_dir .and. modulo(bsize, refinement_factor) /= 0)) then
         write(msg, '(a,3i5,a,i2)')"[refinement:init_refinement] bsize = [", bsize, "] not divisible by ", refinement_factor
         call die(msg)
         ! A block with odd bsize cannot be refined to blocks with odd bsizes because we prefer to allow for incomplete coverage.
         ! Then partially refined block would have to handle fluxes in partially refined cells in some places.
         ! A workaround for this would involve handling blocks with different bsize which is not worth the added complexity.
         do_refine = .false.  ! just in case someone wants to experiment with odd bsize anyway
      endif

      if (dom%eff_dim == 0) then
         level_crit = huge(1)
      else
         level_crit = (63 - int(log(maxval(dom%n_d)-1.)/log(2.)+1)*dom%eff_dim)/dom%eff_dim  ! these are the limits of ordering::Morton_id
      endif

      ! Such large refinements may require additional work in I/O routines, visualization, computing MPI tags and so on.
      if (level_max >  level_crit) then
         write(msg, '(a,i2,a)')"[refinement:init_refinement] BEWARE: Refinement levels above ", level_crit, " may lead to integer overflows in this run under certain conditions."
         if (master) call warn(msg)
      endif
      if (level_max > level_insane) then
         write(msg, '(a,i8,a)')"[refinement:init_refinement] level_max = ", level_max, " is way too much. Reducing. Expect failure somewhat earlier anyway."
         level_max = level_insane
         if (master) call warn(msg)
      endif

      if (.not. do_refine .and. .not. auto_bsize) bsize = I_ZERO

      if (do_refine) then
         write(msg, '(a)')"[refinement]"
         if (level_min /= base_level_id) write(msg(len_trim(msg)+1:), '(a,i2,a)')" minimum level = ", level_min, ","
         write(msg(len_trim(msg)+1:), '(a,i2)')" maximum allowed level = ", level_max
         write(msg(len_trim(msg)+1:), '(a,3i5,a)')", block size = [", bsize, "]"
      else
         if (level_max > base_level_id .and. master) call warn("[refinement] refinements were requested but are disabled by sanity checks")
         msg = "[refinement] No static or adaptive refinement"
         if (any(bsize(:) > 0)) write(msg(len_trim(msg)+1:), '(a,3i5,a)')", block size = [", bsize, "]"
         ! switch off refinement options
         level_min = base_level_id
         level_max = base_level_id
         n_updAMR  = huge(I_ONE)
      endif
      if (master) call printinfo(msg)

      if ((n_updAMR == 0) .and. (level_max > base_level_id)) then
         if (master) call warn("[refinement:init_refinement] refinements disabled (n_updAMR = 0)")
      endif

   end subroutine init_refinement

!> \brief Change the protected parameter n_updAMR

   subroutine set_n_updAMR(n)

      use dataio_pub, only: die
      use mpisetup,   only: piernik_MPI_Bcast

      implicit none

      integer(kind=4), intent(in) :: n

      integer(kind=4) :: nn

      nn = n
      call piernik_MPI_Bcast(nn)
      if (nn /= n) call die("[refinement:set_n_updAMR] n_updAMR goes out of sync")
      n_updAMR = nn

   end subroutine set_n_updAMR

!> \brief Guess some safe bsize

   subroutine automagic_bsize

      use constants,  only: ndims, xdim, ydim, zdim, refinement_factor, I_ONE, I_TWO, INVALID
      use dataio_pub, only: msg, die
      use domain,     only: dom
      use mpisetup,   only: nproc

      implicit none

      integer(kind=4), parameter :: not_too_small = 16 ! bsize below that tends to be inefficient due to huge memory and computation overhead
      integer(kind=4) :: d, i
      integer(kind=4), dimension(ndims) :: b1, b2
      integer(kind=4) :: sq

      if (dom%eff_dim == 0) then
         bsize(:) = 1
         return
      endif

      b1 = INVALID
      b2 = b1

      ! start with size that results with roughly one block per process on highest full level, but don't go below 16 cells per dimension
      ! also divide each dimension to at least 4 pieces even for low thread count
      sq = max(not_too_small, int(((product(dom%n_d, mask=dom%has_dir) * refinement_factor**(dom%eff_dim * level_min))/max(nproc, int((2*refinement_factor)**dom%eff_dim, kind=4)))**(1./dom%eff_dim), kind=4))
      if (mod(sq, I_TWO) == I_ONE) sq = sq + I_ONE
      ! find divisible values in each dim in the range [nb .. dom%n_d] starting from sq in both directions
      if (all(mod(dom%n_d, sq) == 0 .or. .not. dom%has_dir)) then
         bsize = sq
      else
         do d = xdim, zdim
            if (dom%has_dir(d)) then
               if (sq > dom%n_d(d)) then
                  b1(d) = dom%n_d(d)
                  b2(d) = dom%n_d(d)
               else
                  do i = sq, dom%nb, -2
                     if (mod(dom%n_d(d), i) == 0) then
                        b2(d) = i
                        exit
                     endif
                  enddo
                  do i = sq, dom%n_d(d), 2
                     if (i == b2(d)) cycle
                     if (mod(dom%n_d(d), i) == 0) then
                        b1(d) = i
                        exit
                     endif
                  enddo
               endif
            else
               b1(d) = I_ONE
               b2(d) = I_ONE
            endif
         enddo

         ! pick the best one somehow: a bit longer block in x-direction usually doesn't hurt, sometimes gives good performance
         where (b1 == INVALID) b1 = b2
         where (b2 == INVALID) b2 = b1
         bsize = [ b1(xdim), b2(ydim:zdim) ]
      endif

      if (any(bsize == INVALID .and. dom%has_dir)) then
         write(msg, '(a,3i5,a)')"[refinement:automagic_bsize] somewhat invalid block size = [", bsize, "]"
         call die(msg)
      endif

   end subroutine automagic_bsize

end module refinement
