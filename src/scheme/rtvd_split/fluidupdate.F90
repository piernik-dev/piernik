! $Id$
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
!    Initial implemetation of PIERNIK code was based on TVD split MHD code by
!    Ue-Li Pen 
!        see: Pen, Arras & Wong (2003) for algorithm and
!             http://www.cita.utoronto.ca/~pen/MHD 
!             for original source code "mhd.f90" 
!   
!    For full list of developers see $PIERNIK_HOME/license/pdt.txt
!
#include "piernik.def"
module fluidupdate   ! SPLIT

  implicit none

  contains

subroutine fluid_update
  use dataio, only : nlog,ntsl,write_log,write_timeslice, dt_log, dt_tsl, nstep
  use timestep,   only : time_step
  use sweeps, only : sweepx,sweepy,sweepz
  use mpisetup, only : proc,dt,t
  use grid, only : nxd,nyd,nzd

#ifdef DEBUG
  use dataio, only : nhdf,write_hdf
#endif /* DEBUG */

#ifdef COSM_RAYS
  use crdiffusion, only : cr_diff_x, cr_diff_y, cr_diff_z
#endif /* COSM_RAYS */

#ifdef SHEAR
  use shear, only : yshift
  use fluidboundaries, only : bnd_u
#endif /* SHEAR */

#ifdef SN_SRC
  use snsources
#endif /* SN_SRC */

#ifdef SNE_DISTR
  use sndistr
#endif /* SNE_DISTR */

#ifdef SELF_GRAV
  use poissonsolver, only : poisson
#endif /* SELF_GRAV */

  implicit none
  logical, save :: first_run = .true.
#ifdef DEBUG
  integer system, syslog
#endif /* DEBUG */

  call time_step

#ifdef RESISTIVE
   if(first_run) then 
      dt = 0.0
      first_run = .false.
   endif 
#endif /* RESISTIVE */

  if(dt_log .gt. 0.0) then
    if(nlog .lt. (int(t / dt_log) + 1)) then
      call write_log
      nlog = nlog + 1
    endif
  endif

  if(dt_tsl .gt. 0.0) then
    if(ntsl .lt. (int(t / dt_tsl) + 1)) then
      call write_timeslice
      ntsl = ntsl + 1
    endif
  endif

  if(proc.eq.0) write(*,900) nstep,dt,t
900      format('   nstep = ',i7,'   dt = ',e22.16,'   t = ',e22.16)

      t=t+dt

#ifdef SHEAR
      call yshift(t)
      if(nxd /= 1) call bnd_u('xdim')
      if(nyd /= 1) call bnd_u('ydim')
#endif /* SHEAR */
#ifdef SELF_GRAV
      call poisson
#endif /* SELF_GRAV */

!------------------- X->Y->Z ---------------------
      if(nxd /=1) then
         call sweepx
      
#ifdef MAGNETIC      
         call magfieldbyzx
#endif /* MAGNETIC */      
      
#ifdef COSM_RAYS
         call cr_diff_x
#endif /* COSM_RAYS */

#ifdef DEBUG
         syslog = system('echo -n sweep x')
         call write_hdf
         nhdf = nhdf + 1     
#endif /* DEBUG */
      endif

      if(nyd /= 1) then
         call sweepy

#ifdef MAGNETIC            
         call magfieldbzxy
#endif /* MAGNETIC */      

#ifdef COSM_RAYS
         call cr_diff_y
#endif /* COSM_RAYS */

#ifdef DEBUG
         syslog = system('echo -n sweep y')
         call write_hdf
         nhdf = nhdf + 1
#endif /* DEBUG */
      endif

      if(nzd /= 1) then
         call sweepz
      
#ifdef MAGNETIC            
         call magfieldbxyz
#endif /* MAGNETIC */      
      
#ifdef COSM_RAYS
         call cr_diff_z
#endif  /* COSM_RAYS */

#ifdef DEBUG
         syslog = system('echo -n sweep z')
         call write_hdf
         nhdf = nhdf + 1
#endif /* DEBUG */
      endif

! Sources ----------------------------------------

#ifdef SN_SRC
#ifndef SNE_DISTR
      call random_sn
#endif /* SNE_DISTR */
#endif /* SN_SRC */

      t=t+dt
      
#ifdef SHEAR
      call yshift(t)
      if(nxd /= 1) call bnd_u('xdim')
      if(nyd /= 1) call bnd_u('ydim')
#endif /* SHEAR */

#ifdef SELF_GRAV
      call poisson
#endif /* SELF_GRAV */
!-------------------------------------------------



!------------------- Z->Y->X ---------------------
      if(nzd /= 1) then
#ifdef COSM_RAYS
         call cr_diff_z
#endif /* COSM_RAYS */

#ifdef MAGNETIC      
         call magfieldbxyz
#endif /* MAGNETIC */      
      
         call sweepz
      
#ifdef DEBUG
         syslog = system('echo -n sweep z')
         call write_hdf
         nhdf = nhdf + 1
#endif /* DEBUG */
      endif
      
      if(nyd /= 1) then
#ifdef COSM_RAYS
         call cr_diff_y
#endif /* COSM_RAYS */

#ifdef MAGNETIC
         call magfieldbzxy
#endif /* MAGNETIC */      
      
         call sweepy
      
#ifdef DEBUG
         syslog = system('echo -n sweep y')
         call write_hdf
         nhdf = nhdf + 1
#endif /* DEBUG */
      endif

      if(nxd /= 1) then
#ifdef COSM_RAYS
         call cr_diff_x
#endif /* COSM_RAYS */

#ifdef MAGNETIC
         call magfieldbyzx
#endif /* MAGNETIC */
      
         call sweepx
      
#ifdef DEBUG
         syslog = system('echo -n sweep x')
         call write_hdf
         nhdf = nhdf + 1
#endif /* DEBUG */
      endif

#ifdef SNE_DISTR
      call supernovae_distribution
#ifdef DEBUG
      syslog = system('echo -n sne_dis')
      call write_hdf
      nhdf = nhdf + 1
#endif /* DEBUG */
#endif /* SNE_DISTR */

end subroutine fluid_update

!------------------------------------------------------------------------------------------

#ifdef MAGNETIC
  subroutine magfieldbyzx
    use fluidindex, only : ibx,iby,ibz
    use arrays,  only : b
    use grid, only : xdim,ydim,zdim,nyd,nzd
    use advects, only : advectby_x,advectbz_x

#ifdef RESISTIVE
    use resistivity, only : diffuseby_x,diffusebz_x
#endif /* RESISTIVE */

      call advectby_x
      
#ifdef RESISTIVE
      call diffuseby_x
#endif /* RESISTIVE */

      call mag_add(iby,xdim,ibx,ydim)     

      call advectbz_x

#ifdef RESISTIVE
      call diffusebz_x
#endif /* RESISTIVE */

      call mag_add(ibz,xdim,ibx,zdim)

  end subroutine magfieldbyzx

!------------------------------------------------------------------------------------------

  subroutine magfieldbzxy
    use fluidindex, only : ibx,iby,ibz
    use arrays,  only : b
    use grid, only : xdim,ydim,zdim,nzd,nxd
    use advects, only : advectbx_y,advectbz_y

#ifdef RESISTIVE
    use resistivity, only : diffusebx_y,diffusebz_y
#endif /* RESISTIVE */

      call advectbz_y
      
#ifdef RESISTIVE
      call diffusebz_y
#endif /* RESISTIVE */

      call mag_add(ibz,ydim,iby,zdim)

      call advectbx_y
      
#ifdef RESISTIVE
      call diffusebx_y      
#endif /* RESISTIVE */

      call mag_add(ibx,ydim,iby,xdim)

  end subroutine magfieldbzxy

!------------------------------------------------------------------------------------------

  subroutine magfieldbxyz
    use fluidindex, only : ibx,iby,ibz
    use arrays,  only : b
    use grid, only : xdim,ydim,zdim,nxd,nyd
    use advects, only : advectbx_z,advectby_z

#ifdef RESISTIVE
    use resistivity, only : diffusebx_z,diffuseby_z
#endif /* RESISTIVE */

      call advectbx_z
#ifdef RESISTIVE
      call diffusebx_z
#endif /* RESISTIVE */
      call mag_add(ibx,zdim,ibz,xdim)

      call advectby_z
#ifdef RESISTIVE
      call diffuseby_z
#endif /* RESISTIVE */
    
      call mag_add(iby,zdim,ibz,ydim)
  end subroutine magfieldbxyz

!------------------------------------------------------------------------------------------

  subroutine mag_add(ib1,dim1,ib2,dim2)
    use func,   only : pshift, mshift
    use arrays, only : b,wa,wcu
    use grid, only   : dl
    use magboundaries, only : all_mag_boundaries


    implicit none
    integer             :: ib1,ib2,dim1,dim2

#ifdef RESISTIVE
! DIFFUSION FULL STEP

    b(ib1,:,:,:) = b(ib1,:,:,:) - wcu/dl(dim1)
!   wcu = cshift(wcu,shift= 1,dim=dim1)
    wcu = pshift(wcu,dim1)
    b(ib1,:,:,:) = b(ib1,:,:,:) + wcu/dl(dim1)
!   wcu = cshift(wcu,shift=-1,dim=dim1)
    wcu = mshift(wcu,dim1)
    b(ib2,:,:,:) = b(ib2,:,:,:) + wcu/dl(dim2)
!   wcu = cshift(wcu,shift= 1,dim=dim2)
    wcu = pshift(wcu,dim2)
    b(ib2,:,:,:) = b(ib2,:,:,:) - wcu/dl(dim2)
#endif /* RESISTIVE */
! ADVECTION FULL STEP

    b(ib1,:,:,:) = b(ib1,:,:,:) - wa/dl(dim1)
!   wa = cshift(wa,shift=-1,dim=dim1)
    wa = mshift(wa,dim1)
    b(ib1,:,:,:) = b(ib1,:,:,:) + wa/dl(dim1)
    b(ib2,:,:,:) = b(ib2,:,:,:) - wa/dl(dim2)
!   wa = cshift(wa,shift=1,dim=dim2)
    wa = pshift(wa,dim2)
    b(ib2,:,:,:) = b(ib2,:,:,:) + wa/dl(dim2)


    call all_mag_boundaries

  end subroutine mag_add
#endif /* MAGNETIC */
!------------------------------------------------------------------------------------------

end module fluidupdate
