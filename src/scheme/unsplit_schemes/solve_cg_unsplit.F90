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
!! \brief HLLD Riemann solver for ideal magnetohydrodynamics
!!
!! Varadarajan Parthasarathy, CAMK, Warszawa. 2015.
!! Dr. Artur Gawryszczak, CAMK, Warszawa.
!!
!! RK(N) with N .GE. 3 could be helpful for WENO3 ( this statement to be tested )
!!
!! Reference:Relativistic Hydrodynamics, L. Rezzolla, O. Zanotti
!! ---------------------------------------------------------------------------
!! L (or dtodx)--> discretization of spatial differential operator (Eq. 9.135)
!! ---------------------------------------------------------------------------
!! RK2 (Eq. 9.140)
!! u^(1)   = u^(n) + \Delta t L(u^(n))
!! u^(n+1) = 1/2 ( u^(n) + u^(1) + \Delta t L(u^(1)  )
!! ---------------------------------------------------------------------------
!! RK3 (Eq. 9.141)
!! u^(1)   = u(n) + \Delta t L(u^(n))
!! u^(2)   = 1/4 ( 3 u^(n) + u^(1) + \Delta t L(u^(1) ) )
!! u^(n+1) = 1/3 u^(n) + 2/3 u^(2) + 2/3 \Delta t (u^(2))
!! ---------------------------------------------------------------------------
!!
!! Energy fix up routines for CT and its related comments are not used in the current version.
!! The algorithm is simply present for experimental purposes.
!<

module solvecg_unsplit

! pulled by ANY

   implicit none

   private
   public  :: solve_cg_unsplit

contains

! This routine has to conform to the interface defined in sweeps::sweep

   subroutine solve_cg_unsplit(cg,istep)

      use constants,        only: mag_n, GEO_XYZ
      use dataio_pub,       only: die
      use domain,           only: dom
      use fluidindex,       only: flind
      use grid_cont,        only: grid_container
      use global,           only: use_fargo
      use named_array_list, only: wna
      use sources,          only: prepare_sources
       
      implicit none

      type(grid_container), pointer, intent(in) :: cg
      integer,                       intent(in) :: istep     ! stage in the time integration scheme
      integer :: nmag, i

      if (dom%geometry_type /= GEO_XYZ) call die("[solve_cg_unsplit:solve_cg_unsplit] Non-cartesian geometry is not implemented yet in this Unsplit solver.")

      call prepare_sources(cg)

      if (wna%exists(mag_n)) then
         call die("[solve_cg_unsplit:solve_cg_unsplit] Magnetic field is still unsafe for the flux named arrays")
         nmag = 0
         do i = 1, flind%fluids
            if (flind%all_fluids(i)%fl%is_magnetized) nmag = nmag + 1
         enddo
         if (nmag > 1) call die("[solve_cg_riemann:solve_cg_riemann] At most one magnetized fluid is implemented")
         !call solve_cg_ub(cg, ddim, istep)
      else
         call solve_cg_u(cg,istep)
      endif

      cg%processed = .true.

   end subroutine solve_cg_unsplit

   subroutine solve_cg_u(cg,istep)
      use grid_cont,        only: grid_container
      use named_array_list, only: wna, qna
      use constants,        only: pdims, ORTHO1, ORTHO2, I_ONE, LO, HI, uh_n, rk_coef, cs_i2_n, first_stage, last_stage, xdim, ydim, zdim
      use global,           only: dt, integration_order, nstep
      use domain,           only: dom
      use fluidindex,       only: flind, iarr_all_dn, iarr_all_mx, iarr_all_swp
      use fluxtypes,        only: ext_fluxes
      use unsplit_source,   only: apply_source
      implicit none

      type(grid_container), pointer, intent(in) :: cg
      integer,                       intent(in) :: istep    

      integer                                    :: i1, i2, ddim
      integer(kind=4)                            :: uhi
      real, dimension(:,:),allocatable           :: u
      real, dimension(:,:), pointer              :: pu
      real, dimension(:,:), pointer              :: pflux
      real, dimension(:),   pointer              :: cs2
      real, dimension(:,:),allocatable           :: flux
      real, dimension(:,:),allocatable           :: tflux
      type(ext_fluxes)                           :: eflx
      integer                                    :: i_cs_iso2
      integer :: dummy
      uhi = wna%ind(uh_n)
      if (qna%exists(cs_i2_n)) then
         i_cs_iso2 = qna%ind(cs_i2_n)
      else
         i_cs_iso2 = -1
      endif
      cs2 => null()
      do ddim=xdim,zdim
         if (.not. dom%has_dir(ddim)) cycle
         if (.not. allocated(u)) then
            allocate(u(cg%n_(ddim), size(cg%u,1)))
         else
            deallocate(u)
            allocate(u(cg%n_(ddim), size(cg%u,1)))
         endif

         if (.not. allocated(flux)) then
            allocate(flux(size(u, 1)-1,size(u, 2)))
         else
            deallocate(flux)
            allocate(flux(size(u, 1)-1,size(u, 2)))
         endif    
         if (.not. allocated(tflux)) then
            allocate(tflux(size(u, 2),size(u, 1)))
         else
            deallocate(tflux)
            allocate(tflux(size(u, 2),size(u, 1)))
         endif  
         do i2 = cg%ijkse(pdims(ddim, ORTHO2), LO), cg%ijkse(pdims(ddim, ORTHO2), HI)
            do i1 = cg%ijkse(pdims(ddim, ORTHO1), LO), cg%ijkse(pdims(ddim, ORTHO1), HI)  
               if (ddim==xdim) then
                  pflux => cg%w(wna%xflx)%get_sweep(xdim,i1,i2)
               else if (ddim==ydim) then
                  pflux => cg%w(wna%yflx)%get_sweep(ydim,i1,i2)
               else if (ddim==zdim) then 
                  pflux => cg%w(wna%zflx)%get_sweep(zdim,i1,i2)
               endif


               pu => cg%w(uhi)%get_sweep(ddim,i1,i2)
               if (istep == first_stage(integration_order) .or. integration_order < 2 ) pu => cg%w(wna%fi)%get_sweep(ddim,i1,i2)
            

               u(:, iarr_all_swp(ddim,:)) = transpose(pu(:,:))


               
               if (i_cs_iso2 > 0) cs2 => cg%q(i_cs_iso2)%get_sweep(ddim,i1,i2)

               
               call cg%set_fluxpointers(ddim, i1, i2, eflx)
               call solve_u(u,cs2, eflx, flux)


               call cg%save_outfluxes(ddim, i1, i2, eflx)
               tflux(:,2:) = transpose(flux(:, iarr_all_swp(ddim,:)))
               tflux(:,1) = 0.0
               pflux(:,:) = tflux


               !pflux(:,1)  = pflux(:,2)      ! simple outflow BC on the first interface
               !pflux(:,size(pflux,2)) = pflux(:,size(pflux,2)-1)
            end do
         end do
      end do
      call apply_flux(cg,istep)
      call apply_source(cg,istep)
      deallocate(u,flux,tflux)
      nullify(cs2)

   end subroutine solve_cg_u

   subroutine solve_u(ui, cs2, eflx, flx)

      use fluxtypes,      only: ext_fluxes
      use hlld,           only: riemann_wrap_u
      use interpolations, only: interpol
      use dataio_pub,     only: die

      implicit none

      !real, dimension(:,:),        intent(in)    :: u0     !< cell-centered initial fluid states
      real, dimension(:,:),        intent(in) :: ui       !< cell-centered intermediate fluid states
      real, dimension(:), pointer, intent(in)    :: cs2     !< square of local isothermal sound speed
      !real,                        intent(in)    :: dtodx   !< timestep advance: RK-factor * timestep / cell length
      type(ext_fluxes),            intent(inout) :: eflx    !< external fluxes
      real, dimension(:,:),        intent(inout) :: flx     !< Output flux of a 1D chain of a domain at a fixed ortho location of that dimension

      ! left and right states at interfaces 1 .. n-1
      real, dimension(size(ui, 1)-1, size(ui, 2)), target :: ql, qr
      integer, parameter :: in = 1  ! index for cells

      ! fluxes through interfaces 1 .. n-1
      !real, dimension(size(u0, 1)-1, size(u0, 2)), target :: flx

      ! updates required for higher order of integration will likely have shorter length
      if (size(flx,dim=1) /= size(ui, 1)-1 .or. size(flx,dim=2) /= size(ui, 2)  ) then
         call die("[solvecg_unsplit:solve_u] flux array dimension does not match the expected dimensions")
      endif


      call interpol(ui, ql, qr)
      call riemann_wrap_u(ql, qr, cs2, flx) ! Now we advance the left and right states by a timestep.

      if (associated(eflx%li)) flx(eflx%li%index, :) = eflx%li%uflx
      if (associated(eflx%ri)) flx(eflx%ri%index, :) = eflx%ri%uflx
      if (associated(eflx%lo)) eflx%lo%uflx = flx(eflx%lo%index, :)
      if (associated(eflx%ro)) eflx%ro%uflx = flx(eflx%ro%index, :)

      !associate (nx => size(u0, in))
      !   u1(2:nx-1, :) = u0(2:nx-1, :) + dtodx * (flx(:nx-2, :) - flx(2:, :))
      !end associate
   end subroutine solve_u

   subroutine apply_flux(cg, istep)
      use domain,             only: dom
      use grid_cont,          only: grid_container
      use constants,          only: xdim, ydim, zdim, first_stage, last_stage, rk_coef, uh_n, ORTHO1, LO, HI, I_ONE, I_TWO, I_THREE, pdims
      use global,             only: integration_order, dt, nstep
      use named_array_list,   only: wna

      implicit none

      type(grid_container), pointer, intent(in)    :: cg
      integer,                       intent(in)    :: istep

      integer                                      :: afdim, uhi, iul, iuh, jul, juh, kul, kuh
      !real, dimension(:,:,:,:), pointer                :: fld
      !nb = dom%nb

      iul = lbound(cg%w(wna%fi)%arr , 2) ; iuh = ubound(cg%w(wna%fi)%arr , 2) 
      jul = lbound(cg%w(wna%fi)%arr , 3) ; juh = ubound(cg%w(wna%fi)%arr , 3)
      kul = lbound(cg%w(wna%fi)%arr , 4) ; kuh = ubound(cg%w(wna%fi)%arr , 4)
      
      if (dom%has_dir(xdim))  then
         iul = iul + I_ONE
         iuh = iuh - I_ONE
      endif
      if (dom%has_dir(ydim)) then 
         jul = jul + I_ONE
         juh = juh - I_ONE
      endif
      if (dom%has_dir(zdim)) then
         kul = kul + I_ONE
         kuh = kuh - I_ONE
      endif

      !write(111,*) igli,ighi,jgli,jghi,kgli,kghi,iglo,igho,jglo,jgho,kglo,kgho,iul,iuh,jul,juh,kul,kuh

      uhi = wna%ind(uh_n)
      ! if (istep==first_stage(integration_order) .and. integration_order>I_ONE) then
      !    cg%w(uhi)%arr(:,:,:,:) = cg%w(wna%fi)%arr(:,:,:,:)
      !    fld => cg%w(uhi)%arr
      ! else if (istep==last_stage(integration_order) .or. integration_order<2) then
      !    fld => cg%w(wna%fi)%arr


      if (istep==first_stage(integration_order) .and. integration_order>I_ONE) then
         cg%w(uhi)%arr(:,:,:,:) = cg%w(wna%fi)%arr(:,:,:,:)
         do afdim=xdim,zdim
               if (afdim==xdim .and. dom%has_dir(afdim)) then
                  if (dom%has_dir(ydim)) then
                     jul = jul + I_THREE ; juh = juh - I_THREE
                  endif
                  if (dom%has_dir(zdim)) then
                     kul = kul + I_THREE ; kuh = kuh - I_THREE
                  endif

                  cg%w(uhi)%arr(:,iul:iuh,jul:juh,kul:kuh) = cg%w(uhi)%arr(:,iul:iuh,jul:juh,kul:kuh) + &
                                                         dt/cg%dl(xdim) * rk_coef(istep) * (cg%fx(:,iul:iuh,jul:juh,kul:kuh) &
                                                                                          - cg%fx(:,iul+I_ONE:iuh+I_ONE,jul:juh,kul:kuh) )
                  if (dom%has_dir(ydim)) then
                     jul = jul - I_THREE ; juh = juh + I_THREE
                  endif
                  if (dom%has_dir(zdim)) then
                     kul = kul - I_THREE ; kuh = kuh + I_THREE
                  endif 

               else if (afdim==ydim .and. dom%has_dir(afdim)) then 

                  if (dom%has_dir(xdim)) then
                     iul = iul + I_THREE ; iuh = iuh - I_THREE
                  endif
                  if (dom%has_dir(zdim)) then
                     kul = kul + I_THREE ; kuh = kuh - I_THREE
                  endif           

                  cg%w(uhi)%arr(:,iul:iuh,jul:juh,kul:kuh) = cg%w(uhi)%arr(:,iul:iuh,jul:juh,kul:kuh) + &
                                                         dt/cg%dl(ydim) * rk_coef(istep) * (cg%gy(:,iul:iuh,jul:juh,kul:kuh) &
                                                                                          - cg%gy(:,iul:iuh,jul+I_ONE:juh+I_ONE,kul:kuh) )  
                  if (dom%has_dir(xdim)) then
                     iul = iul - I_THREE ; iuh = iuh + I_THREE
                  endif
                  if (dom%has_dir(zdim)) then
                     kul = kul - I_THREE ; kuh = kuh + I_THREE
                  endif   

               else if (afdim==zdim .and. dom%has_dir(afdim)) then

                  if (dom%has_dir(xdim)) then
                     iul = iul + I_THREE ; iuh = iuh - I_THREE
                  endif
                  if (dom%has_dir(ydim)) then
                     jul = jul + I_THREE ; juh = juh - I_THREE
                  endif 

                  cg%w(uhi)%arr(:,iul:iuh,jul:juh,kul:kuh) = cg%w(uhi)%arr(:,iul:iuh,jul:juh,kul:kuh) + &
                                                         dt/cg%dl(zdim) * rk_coef(istep) * (cg%hz(:,iul:iuh,jul:juh,kul:kuh) &
                                                                                          - cg%hz(:,iul:iuh,jul:juh,kul+I_ONE:kuh+I_ONE) )

                  if (dom%has_dir(xdim)) then
                     iul = iul - I_THREE ; iuh = iuh + I_THREE
                  endif
                  if (dom%has_dir(ydim)) then
                     jul = jul - I_THREE ; juh = juh + I_THREE
                  endif 

               end if
         end do
         
      else if (istep==last_stage(integration_order) .or. integration_order<2) then
         do afdim=xdim,zdim
               if (afdim==xdim .and. dom%has_dir(afdim)) then
                  if (dom%has_dir(ydim)) then
                     jul = jul + I_THREE ; juh = juh - I_THREE
                  endif
                  if (dom%has_dir(zdim)) then
                     kul = kul + I_THREE ; kuh = kuh - I_THREE
                  endif

                  cg%w(wna%fi)%arr(:,iul:iuh,jul:juh,kul:kuh) = cg%w(wna%fi)%arr(:,iul:iuh,jul:juh,kul:kuh) + &
                                                         dt/cg%dl(xdim) * rk_coef(istep) * (cg%fx(:,iul:iuh,jul:juh,kul:kuh) &
                                                                                          - cg%fx(:,iul+I_ONE:iuh+I_ONE,jul:juh,kul:kuh) )
                  if (dom%has_dir(ydim)) then
                     jul = jul - I_THREE ; juh = juh + I_THREE
                  endif
                  if (dom%has_dir(zdim)) then
                     kul = kul - I_THREE ; kuh = kuh + I_THREE
                  endif 

               else if (afdim==ydim .and. dom%has_dir(afdim)) then 
                  if (dom%has_dir(xdim)) then
                     iul = iul + I_THREE ; iuh = iuh - I_THREE
                  endif
                  if (dom%has_dir(zdim)) then
                     kul = kul + I_THREE ; kuh = kuh - I_THREE
                  endif           

                  cg%w(wna%fi)%arr(:,iul:iuh,jul:juh,kul:kuh) = cg%w(wna%fi)%arr(:,iul:iuh,jul:juh,kul:kuh) + &
                                                         dt/cg%dl(ydim) * rk_coef(istep) * (cg%gy(:,iul:iuh,jul:juh,kul:kuh) &
                                                                                          - cg%gy(:,iul:iuh,jul+I_ONE:juh+I_ONE,kul:kuh) )  
                  if (dom%has_dir(xdim)) then
                     iul = iul - I_THREE ; iuh = iuh + I_THREE
                  endif
                  if (dom%has_dir(zdim)) then
                     kul = kul - I_THREE ; kuh = kuh + I_THREE
                  endif   

               else if (afdim==zdim .and. dom%has_dir(afdim)) then
                  if (dom%has_dir(xdim)) then
                     iul = iul + I_THREE ; iuh = iuh - I_THREE
                  endif
                  if (dom%has_dir(ydim)) then
                     jul = jul + I_THREE ; juh = juh - I_THREE
                  endif 

                  cg%w(wna%fi)%arr(:,iul:iuh,jul:juh,kul:kuh) = cg%w(wna%fi)%arr(:,iul:iuh,jul:juh,kul:kuh) + &
                                                         dt/cg%dl(zdim) * rk_coef(istep) * (cg%hz(:,iul:iuh,jul:juh,kul:kuh) &
                                                                                          - cg%hz(:,iul:iuh,jul:juh,kul+I_ONE:kuh+I_ONE) )

                  if (dom%has_dir(xdim)) then
                     iul = iul - I_THREE ; iuh = iuh + I_THREE
                  endif
                  if (dom%has_dir(ydim)) then
                     jul = jul - I_THREE ; juh = juh + I_THREE
                  endif 

               end if
         end do
      end if
   end subroutine apply_flux
end module solvecg_unsplit
