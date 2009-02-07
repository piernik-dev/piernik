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

module snsources
   use initproblem, only : amp_ecr_sn,ethu,f_sn, h_sn, r_sn

#ifdef SHEAR
   use shear
#endif /* SHEAR */

   implicit none
   real epsi, epso
   real ysna, ysni, ysno
   integer, save :: nsn, nsn_last
   real, save    :: dt_sn_prev, ecr_supl, decr_supl
   real,    save :: gset
   integer, save :: irand, iset


   contains

   subroutine random_sn
! Written by: M. Hanasz
      use mpisetup, only  : t
      use constants, only : small

      implicit none
      real :: dt_sn
!      real, dimension(2) :: orient
      real, dimension(3) :: snpos
      integer :: isn, nsn_per_timestep

      dt_sn = 1./(f_sn+small)

      nsn = t/dt_sn
      nsn_per_timestep = nsn - nsn_last
      nsn_last = nsn

      do isn = 0, nsn_per_timestep

         call rand_coords(snpos)

#ifdef COSM_RAYS
         call cr_sn(snpos)
#endif /* COSM_RAYS */

#ifdef DIPOLS
!        call magn_multipole_sn(rand_angles(),snpos,)
#endif /* DIPOLS */
      enddo ! isn
      return
   end subroutine random_sn

!--------------------------------------------------------------------------

   subroutine cr_sn(pos)
! Written by: M. Hanasz
      use grid,   only : nx,ny,nz,x,y,z,Lx,Ly
      use arrays, only : u
      use initcosmicrays, only : iecr
      implicit none
      real, dimension(3), intent(in) :: pos
      integer i,j,k, ipm, jpm
      real decr, xsn,ysn,zsn
      xsn = pos(1)
      ysn = pos(2)
      zsn = pos(3)

      do k=1,nz
         do j=1,ny
            do i=1,nx

               do ipm=-1,1

                  if(ipm .eq. -1) ysna = ysno
                  if(ipm .eq.  0) ysna = ysn
                  if(ipm .eq.  1) ysna = ysni

                  do jpm=-1,1

                     decr = amp_ecr_sn * ethu  &
                           * EXP(-((x(i)-xsn+real(ipm)*Lx)**2  &
                           + (y(j)-ysna+real(jpm)*Ly)**2  &
                           + (z(k)-zsn)**2)/r_sn**2)

                     u(iecr,i,j,k) = u(iecr,i,j,k) + decr

                  enddo ! jpm
               enddo ! ipm

            enddo ! i
         enddo ! j
      enddo ! k

      return

   end subroutine cr_sn

!--------------------------------------------------------------------------

   subroutine rand_coords(pos)
! Written by M. Hanasz
      use grid,   only : Lx,Ly,xmin,ymin,nzd 
#ifdef SHEAR
      use grid,   only : dy,nyd,nzd,y,js,je

      integer :: jsn,jremap
      real :: dysn
#endif /* SHEAR */

      real, dimension(3), intent(out) :: pos
      real, dimension(4) :: rand
      real :: xsn,ysn,zsn,znorm,znorm2

      call random_number(rand)
      xsn = xmin+ Lx*rand(1)
      ysn = ymin+ Ly*rand(2)

      if(nzd /= 1) then
         irand = irand+4
         znorm = gasdev(rand(3),rand(4))
         zsn = h_sn*znorm
      else
         zsn = 0.0
      endif

#ifdef SHEAR
      jsn  = js+int((ysn-ymin)/dy)
      dysn  = dmod(ysn,dy)

      epsi   = eps*dy
      epso   = -epsi

!  outer boundary
      jremap = jsn - delj
      jremap = mod(mod(jremap, nyd)+nyd,nyd)
      if (jremap .le. (js-1)) jremap = jremap + nyd

      ysno = y(jremap) + epso + dysn

!  inner boundary
      jremap = jsn + delj
      jremap = mod(jremap, nyd)+nyd
      if (jremap .ge. (je+1)) jremap = jremap - nyd

      ysni = y(jremap) + epsi + dysn
#else /* SHEAR */
      ysno = ysn
      ysni = ysn
#endif /* SHEAR */

      pos(1) = xsn
      pos(2) = ysn
      pos(3) = zsn
      return

   end subroutine rand_coords

!-----------------------------------------------------------------------

!>
!! \brief Function generate point on the surface of unit
!! sphere with uniform distribution and returns its latidude and longitude
!<
   function rand_angles()

      use constants, only : pi
      implicit none
      real :: rand(2)
      real :: rnx,rny,rnz               !> Point's position in cartesian coordinates
      real, dimension(2) :: rand_angles !> Latidue and longitude

      call random_number(rand)
      rnz = (1.0-2.0*rand(1))
      rnx = sqrt(1.0-rnz**2)*cos(2.*pi*rand(2))
      rny = sqrt(1.0-rnz**2)*sin(2.*pi*rand(2))
      rand_angles(1) = 2.0*pi*rand(2)
      rand_angles(2) = acos(rnz/sqrt(rnz**2+rny**2+rnx**2) )

   end function rand_angles

    !=======================================================================
    !
    !      \\\\\\\         B E G I N   S U B R O U T I N E S        ///////
    !      ///////    F R O M   N U M E R I C A L   R E C I P E S   \\\\\\\
    !
    !=======================================================================
    
    
       function gasdev(x,y)
    
          implicit none
          real x, y, x1, y1,  r
          real gasdev, rand(2)
          real fac
          real, save :: gset
          integer, save :: iset, irand
    
          if (iset.eq.0) then
1            x1=2.*x-1.
             y1=2.*y-1.
             r=x1**2+y1**2
             if(r.ge.1.) then
                call random_number(rand)
                x = rand(1)
                y = rand(2)
                irand = irand+2
                go to 1
             endif
             fac=sqrt(-2.*log(r)/r)
             gset=x1*fac
             gasdev=y1*fac
             iset=1
          else
             gasdev=gset
             iset=0
          endif
          return
       end function gasdev
    
    !=======================================================================
    !
    !      \\\\\\\          E N D   S U B R O U T I N E S           ///////
    !      ///////    F R O M   N U M E R I C A L   R E C I P E S   \\\\\\\
    !
    !=======================================================================

end module snsources
