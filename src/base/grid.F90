#include "mhd.def"
module grid

! Written by: M. Hanasz, January/February 2006
      
  implicit none
  real dx,dy,dz,dxmn,dvol
  
  real xminb,xmaxb, yminb,ymaxb, zminb,zmaxb
  real Lx, Ly, Lz

contains

  subroutine grid_xyz 
    use mpi_setup
    use start, only  : xmin,ymin,zmin,xmax,ymax,zmax, dimensions, nb
    use arrays, only : dl,xdim,ydim,zdim,xl,yl,zl,x,y,z,xr,yr,zr, &
       nxb,nyb,nzb,nx,ny,nz
                    
    integer i,j,k 

    xminb = xmin + real(pcoords(1))*(xmax-xmin)/real(psize(1))
    xmaxb = xmin + real(pcoords(1)+1)*(xmax-xmin)/real(psize(1))
    yminb = ymin + real(pcoords(2))*(ymax-ymin)/real(psize(2))
    ymaxb = ymin + real(pcoords(2)+1)*(ymax-ymin)/real(psize(2))
    zminb = zmin + real(pcoords(3))*(zmax-zmin)/real(psize(3))
    zmaxb = zmin + real(pcoords(3)+1)*(zmax-zmin)/real(psize(3))


    dx = (xmaxb-xminb)/nxb
    dy = (ymaxb-yminb)/nyb
    dl(xdim) = dx
    dl(ydim) = dy
    if(dimensions .eq. '3d') then
      dz = (zmaxb-zminb)/nzb
      dl(zdim) = dz
      dxmn = min(dx,dy,dz)
    else if (dimensions .eq. '2dxy') then
      dz = 1.0
      dl(zdim) = dz
      dxmn = min(dx,dy)
    endif
      
!    write(*,*) 'proc=',proc, zminb, zmaxb, dl(zdim)


!--- Asignments -----------------------------------------------------------
    ! left zone boundaries:  xl, yl, zl
    ! zone centers:          x,  y,  z 
    ! right zone boundaries: xr, yr, zr

!--- x-grids --------------------------------------------------------------
    
    do i= 1, nx
      x(i)  = xminb + 0.5*dx + (i-nb-1)*dx
      xl(i) = x(i)  - 0.5*dx
      xr(i) = x(i)  + 0.5*dx
    enddo
       
!--- y-grids --------------------------------------------------------------


    do j= 1, ny

      y(j)  = yminb + 0.5*dy + (j-nb-1)*dy
      yl(j) = y(j)  - 0.5*dy
      yr(j) = y(j)  + 0.5*dy
    enddo
       

!--- z-grids --------------------------------------------------------------

    if(dimensions .eq. '3d') then

      do k= 1, nz
        z(k)  = zminb + 0.5*dz + (k-nb-1) * dz
        zl(k) = z(k)  - 0.5*dz
        zr(k) = z(k)  + 0.5*dz
      enddo      
       
!--------------------------------------------------------------------------

    else if(dimensions .eq. '2dxy') then
   
      z = 0.0  ! dz =1
      zl = -dz/2.
      zr =  dz/2.
     
    endif

    dvol = dx*dy*dz

    Lx = xmax - xmin
    Ly = ymax - ymin
    Lz = zmax - zmin

  end subroutine grid_xyz

end module grid
