! $Id$
#include "piernik.def"

module gravity_user

   character*7 gravpart
#ifdef GALACTIC_DISK
   real, allocatable :: gpotdisk(:,:,:), gpothalo(:,:,:), gpotbulge(:,:,:)
#endif /* GALACTIC_DISK */

   contains

   subroutine grav_pot_user(gpot,sweep,i1,i2,xsw,n,status,temp_log)
      use arrays, only : x,y,z
      use constants, only : gmu, kpc, newtong
      implicit none
      character, intent(in)    :: sweep*6
      character, intent(inout) :: status*9
      integer, intent(in)      :: i1, i2, n
      logical, optional     :: temp_log
      real, dimension(n)    :: xsw
      real, dimension(n),intent(out) :: gpot
      integer i
      real, allocatable :: gpdisk(:), gphalo(:), gpblg(:)
      real rgcs, seccoord, zet
      real, dimension(n)    :: rgcv
      real Mhalo, Mbulge, Mdisk, ahalo, bbulge, adisk, bdisk
      real, dimension(n)    :: rgsv

! galactic case as in vollmer'01 (gravitational potential of Allen & Santillan '01)
         allocate(gpdisk(n),gphalo(n),gpblg(n))
         Mhalo = 4615.0*gmu !g_z*Msun !2.43e11*Msun   !Milky Way
         Mbulge= 0.0 ! 606.0*gmu !ptmass*Msun !0.8e10*Msun  !Milky Way (estimation)
         Mdisk = 3690.0*gmu !dg_dz*Msun !3.7e10*Msun  !Milky Way
         ahalo =   12.0*kpc !n_gravr2*kpc !35.*kpc
         bbulge= 3.0*kpc !0.3873*kpc !ptm_x*pc !2100.*pc
         adisk = 5.3178*kpc !ptm_y*kpc !4.9*kpc
         bdisk = 0.2500*kpc !ptm_z*pc !150.*pc
         if (sweep == 'zsweep') then
            rgcs = sqrt(x(i1)**2+y(i2)**2)
            rgsv = sqrt(x(i1)**2+y(i2)**2+xsw**2)
            gpdisk=-Mdisk/sqrt(rgcs**2+(adisk+sqrt(xsw**2+bdisk**2))**2)
            gphalo=-1./rgsv*Mhalo*(rgsv/ahalo)**2.02/(1.+(rgsv/ahalo)**1.02)
            gphalo=gphalo-(Mhalo/1.02/ahalo)*(-1.02/(1.+(100.0*kpc/ahalo)**1.02)+log(1.0+(100.0*kpc/ahalo)**1.02))
            gphalo=gphalo+(Mhalo/1.02/ahalo)*(-1.02/(1.+(rgsv/ahalo)**1.02)+log(1.0+(rgsv/ahalo)**1.02))
            gpblg=-Mbulge/sqrt(rgsv**2+bbulge**2)
            gpot=(gpdisk+gpblg+gphalo)*newtong
            if(.not.present(temp_log)) then
               gpotdisk(i1,i2,:)=gpdisk*newtong
               gpothalo(i1,i2,:)=gphalo*newtong
               gpotbulge(i1,i2,:)=gpblg*newtong
            endif
         else ! y or x
            if (sweep == 'xsweep') then
               seccoord = y(i1)
               zet = z(i2)
            elseif (sweep == 'ysweep') then
               seccoord = x(i2)
               zet = z(i1)
            endif
            do i=1,n
               rgcv(i) = sqrt(xsw(i)**2+seccoord**2)
               rgsv(i) = sqrt(xsw(i)**2+seccoord**2+zet**2)
               gpot(i)=-Mdisk/sqrt(rgcv(i)**2+(adisk+sqrt(zet**2+bdisk**2))**2)
               gpot(i)=gpot(i)-1./rgsv(i)*Mhalo*(rgsv(i)/ahalo)**2.02/(1.+(rgsv(i)/ahalo))
               gpot(i)=gpot(i)-(Mhalo/1.02/ahalo)*(-1.02/(1.+(100.0*kpc/ahalo)**1.02)+log(1.0+(100.0*kpc/ahalo)**1.02))
               gpot(i)=gpot(i)+(Mhalo/1.02/ahalo)*(-1.02/(1.+(rgsv(i)/ahalo)**1.02)+log(1.0+(rgsv(i)/ahalo)**1.02))
               gpot(i)=gpot(i)-Mbulge/sqrt(rgsv(i)**2+bbulge**2)
               gpot(i)=gpot(i)*newtong
            enddo
         endif
         deallocate(gpdisk,gphalo,gpblg)


   end subroutine grav_pot_user

end module gravity_user