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
module shear
  real    :: ts, dely, eps, omega,qshear
  integer :: delj

  contains

  subroutine init_shear
    use mpisetup
    use errh, only : namelist_errh
    implicit none
    integer :: ierrh
    character(LEN=100) :: par_file, tmp_log_file

    namelist /SHEARING/ omega, qshear
    
    par_file = trim(cwd)//'/problem.par'
    tmp_log_file = trim(cwd)//'/tmp.log'

    omega  = 0.0
    qshear = 0.0
    if(proc .eq. 0) then
       open(1,file=par_file)
          read(unit=1,nml=SHEARING,iostat=ierrh)
          call namelist_errh(ierrh,'SHEARING')
       close(1)
       open(3, file=tmp_log_file, position='append')
          write(unit=3,nml=SHEARING)
       close(3)
  
       rbuff(1) = omega
       rbuff(2) = qshear

       call MPI_BCAST(rbuff,    buffer_dim, MPI_DOUBLE_PRECISION, 0, comm, ierr)

    else

       call MPI_BCAST(rbuff,    buffer_dim, MPI_DOUBLE_PRECISION, 0, comm, ierr)

       omega  = rbuff(1)
       qshear = rbuff(2) 

    endif
  end subroutine init_shear

  subroutine yshift(ts)
    use grid, only : dy,xmin,xmax,nyd
    implicit none
    real, intent(in) :: ts

    dely  = ts*qshear*omega*(xmax-xmin)
    delj  = mod(int(dely/dy),nyd)
    eps   = mod(dely,dy)/dy
  end subroutine yshift
#ifdef FFTW
  subroutine unshear_fft_b(qty,xx,lb,rb,inv)
    use mpisetup, only  : smalld
    use grid, only   : dy,nb,xmax,xmin,nyd,ymax,ymin,nxd
    use constants, only : dpi
    implicit none
    logical, optional               :: inv
    real, dimension(:,:,:)          :: qty
    real, dimension(:), intent(in)  :: xx
    real ,dimension(size(xx) - 2*nb):: x
    real, dimension(nb,size(qty,2),size(qty,3)),intent(out)  :: lb,rb
    real, dimension(size(qty,1),size(qty,2),size(qty,3))     :: unshear_fft
    integer :: nx,ny,nz,p,np
    real    :: St

    integer(kind=8) :: planf,planb

    complex*16      , dimension(size(qty,2)/2 +1) :: ctmp
    real(kind=8)    , dimension(size(qty,2))      :: rtmp
    real(kind=8)    , dimension(size(qty,2)/2+1)  :: ky

    integer     , parameter :: FFTW_ESTIMATE = 64

    x(:) = xx(nb+1:nxd+nb)

    St = dely / (xmax - xmin)
!    St = dely

    nx = size(qty,1)
    ny = size(qty,2)
    nz = size(qty,3)

    np = ny / 2 + 1


!    if(.not.allocated(ctmp)) allocate(ctmp(np))
!    if(.not.allocated(rtmp)) allocate(rtmp(ny))
!    if(.not.allocated(ky))   allocate(ky(np))

    ky(1) = 0.0
    do p = 2, np
      ky(p) = dpi * (p-1) / ny
    enddo
    
    St = -St * dpi / (ymax - ymin) * (nyd / dpi)
    if (.not.present(inv)) St = -St 


    call dfftw_plan_dft_r2c_1d(planf, size(rtmp), rtmp, ctmp, FFTW_ESTIMATE)
    call dfftw_plan_dft_c2r_1d(planb, size(rtmp), ctmp, rtmp, FFTW_ESTIMATE)

    do p = 1,nx
      rtmp(:)  = qty(p,:,1)
      call dfftw_execute(planf)
      ctmp(:)  = ctmp(:)*exp(cmplx( 0.0, St*ky(:)*x(p) ) )
      call dfftw_execute(planb)
      unshear_fft(p,:,1)  = rtmp(:) / ny
    enddo

    if (present(inv)) then 
    do p = 1,nb
      rtmp(:)  = qty(p,:,1)
      call dfftw_execute(planf)
      ctmp(:)  = ctmp(:)*exp(cmplx( 0.0, St*ky(:)*x(p) ) )
      call dfftw_execute(planb)
      rb(p,:,1)  = rtmp(:) / ny

      rtmp(:)  = qty(nx-nb+p,:,1)
      call dfftw_execute(planf)
      ctmp(:)  = ctmp(:)*exp(cmplx( 0.0, St*ky(:)*x(nx-nb+p) ) )
      call dfftw_execute(planb)
      lb(p,:,1)  = rtmp(:) / ny
    enddo
    else
    do p = 1,nb
      rtmp(:)  = qty(p,:,1)
      call dfftw_execute(planf)
      ctmp(:)  = ctmp(:)*exp(cmplx( 0.0, St*ky(:)*xx(nxd+nb+p) ) )
      call dfftw_execute(planb)
      rb(p,:,1)  = rtmp(:) / ny

      rtmp(:)  = qty(nx-nb+p,:,1)
      call dfftw_execute(planf)
      ctmp(:)  = ctmp(:)*exp(cmplx( 0.0, St*ky(:)*xx(p) ) )
      call dfftw_execute(planb)
      lb(p,:,1)  = rtmp(:) / ny
    enddo
    endif


    call dfftw_destroy_plan(planf)
    call dfftw_destroy_plan(planb)

!    if (allocated(rtmp))   deallocate(rtmp)
!    if (allocated(ctmp))   deallocate(ctmp)
!    if (allocated(ky))     deallocate(ky)
    return

  end subroutine unshear_fft_b
  function unshear_fft(qty,x,inv)
    use mpisetup, only  : smalld 
    use grid, only   : dy,nb,xmax,xmin,nyd,ymax,ymin
    use constants, only : dpi
    implicit none
    logical, optional               :: inv
    real, dimension(:,:,:)          :: qty
    real, dimension(:), intent(in)  :: x
    real, dimension(size(qty,1),size(qty,2),size(qty,3)) :: unshear_fft
    integer :: nx,ny,nz,p,np
    real    :: St

    integer(kind=8) :: planf,planb

    complex*16      , dimension(:)    , allocatable :: ctmp
    real(kind=8)    , dimension(:)    , allocatable :: rtmp
    real(kind=8)    , dimension(:)    , allocatable :: ky

    integer     , parameter :: FFTW_ESTIMATE = 64

    St = dely / (xmax - xmin)
!    St = dely

    nx = size(qty,1)
    ny = size(qty,2)
    nz = size(qty,3)

    np = ny / 2 + 1

    if(.not.allocated(ctmp)) allocate(ctmp(np))
    if(.not.allocated(rtmp)) allocate(rtmp(ny))
    if(.not.allocated(ky))   allocate(ky(np))

    ky(1) = 0.0
    do p = 2, np
      ky(p) = dpi * (p-1) / ny
    enddo
    
    St = -St * dpi / (ymax - ymin) * (nyd / dpi)
    if (.not.present(inv)) St = -St 


    call dfftw_plan_dft_r2c_1d(planf, size(rtmp), rtmp, ctmp, FFTW_ESTIMATE)
    call dfftw_plan_dft_c2r_1d(planb, size(rtmp), ctmp, rtmp, FFTW_ESTIMATE)

    do p = 1,nx
      rtmp(:)  = qty(p,:,1)
      call dfftw_execute(planf)
      ctmp(:)  = ctmp(:)*exp(cmplx( 0.0, St*ky(:)*x(p) ) )
      call dfftw_execute(planb)
      unshear_fft(p,:,1)  = rtmp(:) / ny
    enddo

    call dfftw_destroy_plan(planf)
    call dfftw_destroy_plan(planb)

    if (allocated(rtmp))   deallocate(rtmp)
    if (allocated(ctmp))   deallocate(ctmp)
    if (allocated(ky))     deallocate(ky)
    return

  end function unshear_fft
#endif /* FFTW */

  function unshear(qty,x,inv)
    use grid,  only  : nb,xmax,xmin,nyd,dy
    use mpisetup, only  : smalld
    
    logical, optional               :: inv
    real, dimension(:,:,:)          :: qty
    real, dimension(:), intent(in)  :: x
    real, dimension(size(qty,1),size(qty,2),size(qty,3)) :: unshear
    real, dimension(:,:), allocatable:: temp
    integer :: i,sg,my,nx,ny,nz
    real    :: fx

    nx = size(qty,1)
    ny = size(qty,2)
    nz = size(qty,3)

    my = 3*nyd+2*nb

    fx = dely / (xmax - xmin)
    sg = -1

    if(.not.allocated(temp)) allocate(temp(my,nz))

    unshear = 0.0

    if(present(inv)) then
       fx = - fx
    endif
    do i = 1, nx
      dl  = fx * x(i)
      ndl = mod(int(dl/dy),nyd)
      ddl = mod(dl,dy)/dy

      temp(         1:  nyd+nb,:)   = qty(i,   1:nyd+nb ,:)
      temp(  nyd+nb+1:2*nyd+nb,:)   = qty(i,nb+1:nyd+nb,:)
      temp(2*nyd+nb+1:3*nyd+2*nb,:) = qty(i,nb+1:ny    ,:)

      temp = cshift(temp,dim=1,shift=ndl)

!      temp(1:nb,:) = temp(nyb+1:nyb+nb,:)          ! not needed
!      temp(nb+nyb+1:nyb+2*nb,:) = temp(nb+1:2*nb,:)

      temp(:,:) = (1.0+ddl)*(1.0-ddl) * temp(:,:) &
            - 0.5*(ddl)*(1.0-ddl) * cshift(temp(:,:),shift= sg,dim=1) &
            + 0.5*(ddl)*(1.0+ddl) * cshift(temp(:,:),shift=-sg,dim=1) 

      unshear(i,nb+1:nb+nyd,:) = temp(nb+nyd+1:nb+2*nyd,:)

      unshear(i,1:nb,:)          = unshear(i,nyd+1:nyd+nb,:)
      unshear(i,nyd+nb+1:ny,:)   = unshear(i,nb+1 :2*nb,:)

!      unshear(i,:,:) = max(unshear(i,:,:), smalld)
    enddo
    if (allocated(temp)) deallocate(temp)
    return
  end function unshear
end module shear

