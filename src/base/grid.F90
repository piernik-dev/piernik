! $Id$
#include "piernik.def"
module grid

! Written by: M. Hanasz, January/February 2006

   implicit none
   real    :: dx, dy, dz, dxmn, dvol
   integer :: nxd, nyd, nzd, nb
   integer :: nx, ny, nz
   integer :: nxb, nyb, nzb
   integer :: nxt, nyt, nzt
   integer :: is, ie, js, je, ks, ke
   integer :: maxxyz

   real    :: xmin, xmax, ymin, ymax, zmin, zmax
   real    :: xminb, xmaxb, yminb, ymaxb, zminb, zmaxb
   real    :: Lx, Ly, Lz
   integer,parameter  :: xdim=1, ydim=2, zdim=3

   real, allocatable :: dl(:)
   real, allocatable, dimension(:)  :: x, xl, xr
   real, allocatable, dimension(:)  :: y, yl, yr
   real, allocatable, dimension(:)  :: z, zl, zr


contains

   subroutine init_grid
      use mpi_setup
      implicit none
      character(LEN=100) :: par_file, tmp_log_file
      
      namelist /DOMAIN_SIZES/ nxd, nyd, nzd, nb
      namelist /DOMAIN_LIMITS/ xmin, xmax, ymin, ymax, zmin, zmax

      nxd  = 1
      nyd  = 1
      nzd  = 1
      nb   = 4

      if(proc == 0) then
         par_file = trim(cwd)//'/problem.par'
         tmp_log_file = trim(cwd)//'/tmp.log'
         open(1,file=par_file)
            read(unit=1,nml=DOMAIN_SIZES)
         close(1)
         open(1,file=par_file)
            read(unit=1,nml=DOMAIN_LIMITS)
         close(1)
         open(3, file='tmp.log', position='append')
           write(3,nml=DOMAIN_SIZES)
           write(3,nml=DOMAIN_LIMITS)
           write(3,*)
         close(3)
      endif


      if(proc == 0) then

         ibuff(1)   = nxd
         ibuff(2)   = nyd
         ibuff(3)   = nzd
         ibuff(4)   = nb

         rbuff(1)   = xmin
         rbuff(2)   = xmax
         rbuff(3)   = ymin
         rbuff(4)   = ymax
         rbuff(5)   = zmin
         rbuff(6)   = zmax

         call MPI_BCAST(ibuff,    buffer_dim, MPI_INTEGER,          0, comm, ierr)
         call MPI_BCAST(rbuff,    buffer_dim, MPI_DOUBLE_PRECISION, 0, comm, ierr)

      else

         call MPI_BCAST(ibuff,    buffer_dim, MPI_INTEGER,          0, comm, ierr)
         call MPI_BCAST(rbuff,    buffer_dim, MPI_DOUBLE_PRECISION, 0, comm, ierr)

         nxd  = ibuff(1)
         nyd  = ibuff(2)
         nzd  = ibuff(3)
         nb   = ibuff(4)
         
         xmin = rbuff(1)
         xmax = rbuff(2)
         ymin = rbuff(3)
         ymax = rbuff(4)
         zmin = rbuff(5)
         zmax = rbuff(6)

      endif

      if((mod(nxd, pxsize) .ne. 0) .or. &
         (mod(nyd, pysize) .ne. 0) .or. &
         (mod(nzd, pzsize) .ne. 0)) then
         call mpistop
         if (proc .eq. 0) then
            write(*,*) 'One of: (mod(nxd,pxsize) .or. mod(nyd,pysize) .or. mod(nzd,pzsize)) .ne. 0'
         endif
         stop
      endif

      nxb = nxd/pxsize     !
      nyb = nyd/pysize     ! Block 'physical' grid sizes
      nzb = nzd/pzsize     !

      nx=nxb+2*nb          !
      ny=nyb+2*nb          ! Block total grid sizes
      nz=nzb+2*nb          !

      nxt=nxd+2*nb         !
      nyt=nyd+2*nb         ! Domain total grid sizes
      nzt=nzd+2*nb         !

      if(nxd == 1) then
         nx     = 1
         nxb    = 1
         nxt    = 1
         pxsize = 1
         is     = 1
         ie     = 1
      else
         is = nb+1
         ie = nb+nxb
      endif

      if(nyd == 1) then
         ny     = 1
         nyb    = 1
         nyt    = 1
         pysize = 1
         js     = 1
         je     = 1
      else
         js = nb+1
         je = nb+nyb
      endif

      if(nzd == 1) then
         nz     = 1
         nzb    = 1
         nzt    = 1
         pzsize = 1
         ks     = 1
         ke     = 1
      else 
         ks = nb+1
         ke = nb+nzb
      endif
      allocate(dl(3))
      allocate(x(nx), xl(nx), xr(nx))
      allocate(y(ny), yl(ny), yr(ny))
      allocate(z(nz), zl(nz), zr(nz))

   end subroutine init_grid

  subroutine grid_xyz
    use mpi_setup

    implicit none
    integer i,j,k

    maxxyz = max(size(x),size(y))
    maxxyz = max(size(z),maxxyz)

    xminb = xmin + real(pcoords(1)  )*(xmax-xmin)/real(psize(1))
    xmaxb = xmin + real(pcoords(1)+1)*(xmax-xmin)/real(psize(1))
    yminb = ymin + real(pcoords(2)  )*(ymax-ymin)/real(psize(2))
    ymaxb = ymin + real(pcoords(2)+1)*(ymax-ymin)/real(psize(2))
    zminb = zmin + real(pcoords(3)  )*(zmax-zmin)/real(psize(3))
    zmaxb = zmin + real(pcoords(3)+1)*(zmax-zmin)/real(psize(3))


    if(nxd /= 1) then
       dx = (xmaxb-xminb)/nxb
    else
       dx = 1.0
    endif
    if(nyd /= 1) then
       dy = (ymaxb-yminb)/nyb
    else
       dy = 1.0
    endif
    if(nzd /= 1) then
       dz = (zmaxb-zminb)/nzb
    else 
       dz = 1.0
    endif

    dl(xdim) = dx
    dl(ydim) = dy
    dl(zdim) = dz

    dvol = dx*dy*dz
    dxmn = min(dx,dy,dz)


!    write(*,*) 'proc=',proc, zminb, zmaxb, dl(zdim)


!--- Asignments -----------------------------------------------------------
    ! left zone boundaries:  xl, yl, zl
    ! zone centers:          x,  y,  z
    ! right zone boundaries: xr, yr, zr

!--- x-grids --------------------------------------------------------------

    if(nxd /= 1) then
       do i= 1, nx
          x(i)  = xminb + 0.5*dx + (i-nb-1)*dx
          xl(i) = x(i)  - 0.5*dx
          xr(i) = x(i)  + 0.5*dx
       enddo
    else
       x  =  0.0
       xl = -0.5*dx
       xr =  0.5*dx
    endif

!--- y-grids --------------------------------------------------------------

    if(nyd /= 1) then
       do j= 1, ny
          y(j)  = yminb + 0.5*dy + (j-nb-1)*dy
          yl(j) = y(j)  - 0.5*dy
          yr(j) = y(j)  + 0.5*dy
       enddo
    else
       y  =  0.0
       yl = -0.5*dy
       yr =  0.5*dy
    endif

!--- z-grids --------------------------------------------------------------

    if(nzd /= 1) then
       do k= 1, nz
          z(k)  = zminb + 0.5*dz + (k-nb-1) * dz
          zl(k) = z(k)  - 0.5*dz
          zr(k) = z(k)  + 0.5*dz
       enddo
    else 
       z  =  0.0  
       zl = -0.5*dz
       zr =  0.5*dz
    endif
!--------------------------------------------------------------------------

    Lx = xmax - xmin
    Ly = ymax - ymin
    Lz = zmax - zmin

  end subroutine grid_xyz

   subroutine cleanup_grid
      implicit none
      
      deallocate(dl)
      deallocate(x, xl, xr)
      deallocate(y, yl, yr)
      deallocate(z, zl, zr)


   end subroutine cleanup_grid

end module grid
