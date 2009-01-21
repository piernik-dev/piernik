! $Id: cr_diffusion.F90 594 2009-01-21 12:33:44Z xarth $
#include "piernik.def"

! Written by M. Hanasz, 2003
! Adapted for this code M. Hanasz, October 2007

module crdiffusion


  use arrays
  use constants
  use grid
  use start


contains


!
      subroutine cr_diff_x
!
!  PURPOSE:  Diffusive transport of ecr in 1-direction
!
      implicit none
      integer i,j,k, n
      real b1b, b2b, b3b, n1b, n2b, n3b
      real, dimension(COSM_RAYS) :: decr1, decr2, decr3, fcrdif1
      real, dimension(COSM_RAYS) :: dqp2, dqm2, dqp3, dqm3

!=======================================================================

do n=1,COSM_RAYS
      do 30 k=ks,ke
        do 20 j=js,je
          do 10 i=is,ie+1

             b1b =  b(ibx,i,  j,  k)
             b2b = (b(iby,i,  j,  k) + b(iby,i-1,j,  k)       &
                  + b(iby,i-1,j+1,k) + b(iby,i,  j+1,k))/4.
        if(dimensions .eq. '3d') then
             b3b = (b(ibz,i,  j,  k) + b(ibz,i-1,j,  k)       &
                  + b(ibz,i-1,j,k+1) + b(ibz,i,  j,k+1))/4.
        endif

             n1b = b1b/SQRT(b1b**2 + b2b**2 + b3b**2+small)
             n2b = b2b/SQRT(b1b**2 + b2b**2 + b3b**2+small)
             n3b = b3b/SQRT(b1b**2 + b2b**2 + b3b**2+small)

             decr1(n) =  (u(iecr(n),i,  j,  k) - u(iecr(n),i-1,j,  k))/dx

             dqm2(n)  = 0.5*((u(iecr(n),i-1,j  ,k ) + u(iecr(n),i ,j  ,k ))    &
                         -(u(iecr(n),i-1,j-1,k ) + u(iecr(n),i ,j-1,k )))/dy
             dqp2(n)  = 0.5*((u(iecr(n),i-1,j+1,k ) + u(iecr(n),i ,j+1,k ))    &
                         -(u(iecr(n),i-1,j  ,k ) + u(iecr(n),i ,j  ,k )))/dy

             decr2(n) = (dqp2(n)+dqm2(n))* (1.0 + sign(1.0, dqm2(n)*dqp2(n)))/4.

        if(dimensions .eq. '3d') then
             dqm3(n)  = 0.5*((u(iecr(n),i-1,j ,k  ) + u(iecr(n),i ,j ,k  ))    &
                         -(u(iecr(n),i-1,j ,k-1) + u(iecr(n),i ,j ,k-1)))/dz
             dqp3(n)  = 0.5*((u(iecr(n),i-1,j ,k+1) + u(iecr(n),i ,j ,k+1))    &
                         -(u(iecr(n),i-1,j ,k  ) + u(iecr(n),i ,j ,k  )))/dz

              decr3(n) = (dqp3(n)+dqm3(n))* (1.0 + sign(1.0, dqm3(n)*dqp3(n)))/4.
        endif

              fcrdif1(n) = K_cr_paral * n1b *   &
                   (n1b*decr1(n) + n2b*decr2(n) + n3b*decr3(n))  &
                   + K_cr_perp * decr1(n)


             wa(i,j,k) = - fcrdif1(n) * dt

10          continue
20        continue
30      continue

      do 60 k=ks,ke
        do 50 j=js,je
          do 40 i=is,ie
            u(iecr(n),i,j,k) = u(iecr(n),i,j,k)  &
                        -(wa(i+1,j,k)-wa(i,j,k))/dx
40        continue
50      continue
60    continue
enddo
      return
      end subroutine cr_diff_x


!****************************************************************************

!
      subroutine cr_diff_y
!
!  PURPOSE:   Diffusive transport of ecr in 2-direction
!
!-----------------------------------------------------------------------
      implicit none
      integer i,j,k,n
      real b1b, b2b, b3b, n1b, n2b, n3b
      real, dimension(COSM_RAYS) :: decr1, decr2, decr3, fcrdif2
      real, dimension(COSM_RAYS) :: dqp2, dqm1, dqm2, dqp3, dqm3, dqp1

!=======================================================================
!

do n=1,COSM_RAYS
      do 30 k=ks,ke
        do 20 j=js,je+1
          do 10 i=is,ie

             b1b = (b(ibx,i,j,k) + b(ibx,i,j-1,k)   &
                  + b(ibx,i+1,j-1,k) + b(ibx,i+1,j,k))/4.

             b2b =  b(iby,i,j,k)
             if(dimensions .eq. '3d') then
                b3b = (b(ibz,i,j,k) + b(ibz,i,j-1,k)   &
                     + b(ibz,i,j-1,k+1) + b(ibz,i,j,k+1))/4.
             endif

             n1b = b1b/SQRT(b1b**2 + b2b**2 + b3b**2+small)
             n2b = b2b/SQRT(b1b**2 + b2b**2 + b3b**2+small)
             n3b = b3b/SQRT(b1b**2 + b2b**2 + b3b**2+small)

             dqm1(n)  = 0.5*((u(iecr(n),i  ,j-1,k ) + u(iecr(n),i  ,j  ,k ))   &
                         -(u(iecr(n),i-1,j-1,k ) + u(iecr(n),i-1,j  ,k )))/dx
             dqp1(n)  = 0.5*((u(iecr(n),i+1,j-1,k ) + u(iecr(n),i+1,j  ,k ))   &
                         -(u(iecr(n),i  ,j-1,k ) + u(iecr(n),i  ,j  ,k )))/dx

             decr1(n) = (dqp1(n)+dqm1(n))* (1.0 + sign(1.0, dqm1(n)*dqp1(n)))/4.

             decr2(n) = (u(iecr(n),i,j,k) - u(iecr(n),i,j-1,k))/dy

             if(dimensions .eq. '3d') then
                dqm2(n)  = 0.5*((u(iecr(n),i ,j-1,k  ) + u(iecr(n),i ,j  ,k  ))   &
                            -(u(iecr(n),i ,j-1,k-1) + u(iecr(n),i ,j  ,k-1)))/dz
                dqp2(n)  = 0.5*((u(iecr(n),i ,j-1,k+1) + u(iecr(n),i ,j  ,k+1))   &
                            -(u(iecr(n),i ,j-1,k  ) + u(iecr(n),i ,j  ,k  )))/dz

                decr3(n) = (dqp2(n)+dqm2(n))* (1.0 + sign(1.0, dqm2(n)*dqp2(n)))/4.
             endif


             fcrdif2(n) = K_cr_paral * n2b * &
                      (n1b*decr1(n) + n2b*decr2(n) + n3b*decr3(n)) &
                      + K_cr_perp * decr2(n)

             wa(i,j,k) = - fcrdif2(n) * dt

10          continue
20        continue
30      continue

      do 60 k=ks,ke
        do 50 j=js,je
          do 40 i=is,ie
            u(iecr(n),i,j,k) = u(iecr(n),i,j,k)    &
                      - (wa(i,j+1,k)-wa(i,j,k))/dy
40        continue
50      continue
60    continue
enddo

      return
      end subroutine cr_diff_y

!************************************************************************

      subroutine cr_diff_z
!
!  PURPOSE:   Diffusive transport of ecr in 3-direction
!
      implicit none
      integer i,j,k,n
      real b1b, b2b, b3b, n1b, n2b, n3b
      real, dimension(COSM_RAYS) :: decr1, decr2, decr3, fcrdif3
      real, dimension(COSM_RAYS) :: dqp2, dqm2, dqp3, dqm3, dqm1, dqp1

!=======================================================================
!
do n=1,COSM_RAYS
   if(dimensions .eq. '3d') then
      do 100 k=ks,ke+1
        do 100 j=js,je
          do 100 i=is,ie

             b1b = (b(ibx,i,  j,k  ) + b(ibx,i,  j,k-1)  &
                  + b(ibx,i+1,j,k-1) + b(ibx,i+1,j,k  ))/4.
             b2b = (b(iby,i,j,  k  ) + b(iby,i,j,  k-1)  &
                  + b(iby,i,j+1,k-1) + b(iby,i,j+1,k  ))/4.
             b3b =  b(ibz,i,j,k)

             n1b = b1b/SQRT(b1b**2 + b2b**2 + b3b**2+small)
             n2b = b2b/SQRT(b1b**2 + b2b**2 + b3b**2+small)
             n3b = b3b/SQRT(b1b**2 + b2b**2 + b3b**2+small)

             dqm1(n)  = 0.5*((u(iecr(n),i  ,j ,k-1) + u(iecr(n),i  ,j  ,k ))   &
                         -(u(iecr(n),i-1,j ,k-1) + u(iecr(n),i-1,j  ,k )))/dx
             dqp1(n)  = 0.5*((u(iecr(n),i+1,j ,k-1) + u(iecr(n),i+1,j  ,k ))   &
                         -(u(iecr(n),i  ,j ,k-1) + u(iecr(n),i  ,j  ,k )))/dx

             decr1(n) = (dqp1(n)+dqm1(n))* (1.0 + sign(1.0, dqm1(n)*dqp1(n)))/4.

             dqm2(n)  = 0.5*((u(iecr(n),i ,j  ,k-1) + u(iecr(n),i  ,j  ,k ))   &
                         -(u(iecr(n),i ,j-1,k-1) + u(iecr(n),i  ,j-1,k )))/dy
             dqp2(n)  = 0.5*((u(iecr(n),i ,j+1,k-1) + u(iecr(n),i  ,j+1,k ))   &
                         -(u(iecr(n),i ,j  ,k-1) + u(iecr(n),i  ,j  ,k )))/dy

             decr2(n) = (dqp2(n)+dqm2(n))* (1.0 + sign(1.0, dqm2(n)*dqp2(n)))/4.

             decr3(n) =  (u(iecr(n),i,j,k    ) - u(iecr(n),i,j,  k-1))/dz

             fcrdif3(n) = K_cr_paral * n3b * &
                  (n1b*decr1(n) + n2b*decr2(n) + n3b*decr3(n)) &
                   + K_cr_perp * decr3(n)

             wa(i,j,k) = - fcrdif3(n) * dt

100   continue

      do 60 k=ks,ke
        do 50 j=js,je
          do 40 i=is,ie
            u(iecr(n),i,j,k) = u(iecr(n),i,j,k)  &
                        -(wa(i,j,k+1)-wa(i,j,k))/dz
40        continue
50      continue
60    continue

      endif
enddo

      return
      end subroutine cr_diff_z

  subroutine div_v
    implicit none
    real, dimension(nx)  :: dvx,dvy,dvz, tmp
    integer i,j,k

    wa(:,:,:) = 0.0
    if(dimensions .eq. '2dxy') then
      k=1
      do j=2,ny-1
        tmp=u(imxi(1),:,j,k)/u(idni(1),:,j,k)
        wa(:,j,k) = (eoshift(tmp,1) - eoshift(tmp,-1))/(2.*dx)
        tmp=(u(imyi(1),:,j+1,k)/u(idni(1),:,j+1,k)-u(imyi(1),:,j-1,k)/u(idni(1),:,j-1,k))/(2.*dy)
        wa(:,j,k) = wa(:,j,k) + tmp
      enddo
    else if (dimensions .eq. '3d') then
      do k=2,nz-1
        do j=2,ny-1
          tmp=u(imxi(1),:,j,k)/u(idni(1),:,j,k)
          wa(:,j,k) = (eoshift(tmp,1) - eoshift(tmp,-1))/(2.*dx)
          tmp=(u(imyi(1),:,j+1,k)/u(idni(1),:,j+1,k)-u(imyi(1),:,j-1,k)/u(idni(1),:,j-1,k))/(2.*dy)
          wa(:,j,k) = wa(:,j,k) + tmp
          tmp=(u(imzi(1),:,j,k+1)/u(idni(1),:,j,k+1)-u(imzi(1),:,j,k-1)/u(idni(1),:,j,k-1))/(2.*dz)
          wa(:,j,k) = wa(:,j,k) + tmp
        enddo
      enddo
    endif
  end subroutine div_v


end module crdiffusion
