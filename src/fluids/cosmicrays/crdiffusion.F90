! $Id: cr_diffusion.F90 594 2009-01-21 12:33:44Z xarth $
#include "piernik.def"

! Written by M. Hanasz, 2003
! Adapted for this code M. Hanasz, October 2007

module crdiffusion

  use initcosmicrays, only : iecr
  use initcosmicrays, only : K_cr_paral,K_cr_perp 
  
  use fluidindex,     only : ibx,iby,ibz
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
      real :: decr1, decr2, decr3, fcrdif1
      real :: dqp2, dqm2, dqp3, dqm3

!=======================================================================

      do 30 k=ks,ke
        do 20 j=js,je
          do 10 i=is,ie+1

             b1b =  b(ibx,i,  j,  k)
             b2b = (b(iby,i,  j,  k) + b(iby,i-1,j,  k)       &
                  + b(iby,i-1,j+1,k) + b(iby,i,  j+1,k))/4.
        if(nzd .gt. 1) then
             b3b = (b(ibz,i,  j,  k) + b(ibz,i-1,j,  k)       &
                  + b(ibz,i-1,j,k+1) + b(ibz,i,  j,k+1))/4.
        endif

             n1b = b1b/SQRT(b1b**2 + b2b**2 + b3b**2+small)
             n2b = b2b/SQRT(b1b**2 + b2b**2 + b3b**2+small)
             n3b = b3b/SQRT(b1b**2 + b2b**2 + b3b**2+small)

             decr1 =  (u(iecr,i,  j,  k) - u(iecr,i-1,j,  k))/dx

             dqm2  = 0.5*((u(iecr,i-1,j  ,k ) + u(iecr,i ,j  ,k ))    &
                         -(u(iecr,i-1,j-1,k ) + u(iecr,i ,j-1,k )))/dy
             dqp2  = 0.5*((u(iecr,i-1,j+1,k ) + u(iecr,i ,j+1,k ))    &
                         -(u(iecr,i-1,j  ,k ) + u(iecr,i ,j  ,k )))/dy

             decr2 = (dqp2+dqm2)* (1.0 + sign(1.0, dqm2*dqp2))/4.

        if(nzd .gt. 1) then
             dqm3  = 0.5*((u(iecr,i-1,j ,k  ) + u(iecr,i ,j ,k  ))    &
                         -(u(iecr,i-1,j ,k-1) + u(iecr,i ,j ,k-1)))/dz
             dqp3  = 0.5*((u(iecr,i-1,j ,k+1) + u(iecr,i ,j ,k+1))    &
                         -(u(iecr,i-1,j ,k  ) + u(iecr,i ,j ,k  )))/dz

              decr3 = (dqp3+dqm3)* (1.0 + sign(1.0, dqm3*dqp3))/4.
        endif

              fcrdif1 = K_cr_paral * n1b *   &
                   (n1b*decr1 + n2b*decr2 + n3b*decr3)  &
                   + K_cr_perp * decr1


             wa(i,j,k) = - fcrdif1 * dt

10          continue
20        continue
30      continue

      do 60 k=ks,ke
        do 50 j=js,je
          do 40 i=is,ie
            u(iecr,i,j,k) = u(iecr,i,j,k)  &
                        -(wa(i+1,j,k)-wa(i,j,k))/dx
40        continue
50      continue
60    continue

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
      real :: decr1, decr2, decr3, fcrdif2
      real :: dqp2, dqm1, dqm2, dqp3, dqm3, dqp1

!=======================================================================
!

      do 30 k=ks,ke
        do 20 j=js,je+1
          do 10 i=is,ie

             b1b = (b(ibx,i,j,k) + b(ibx,i,j-1,k)   &
                  + b(ibx,i+1,j-1,k) + b(ibx,i+1,j,k))/4.

             b2b =  b(iby,i,j,k)
             if(nzd .gt. 1) then
                b3b = (b(ibz,i,j,k) + b(ibz,i,j-1,k)   &
                     + b(ibz,i,j-1,k+1) + b(ibz,i,j,k+1))/4.
             endif

             n1b = b1b/SQRT(b1b**2 + b2b**2 + b3b**2+small)
             n2b = b2b/SQRT(b1b**2 + b2b**2 + b3b**2+small)
             n3b = b3b/SQRT(b1b**2 + b2b**2 + b3b**2+small)

             dqm1  = 0.5*((u(iecr,i  ,j-1,k ) + u(iecr,i  ,j  ,k ))   &
                         -(u(iecr,i-1,j-1,k ) + u(iecr,i-1,j  ,k )))/dx
             dqp1  = 0.5*((u(iecr,i+1,j-1,k ) + u(iecr,i+1,j  ,k ))   &
                         -(u(iecr,i  ,j-1,k ) + u(iecr,i  ,j  ,k )))/dx

             decr1 = (dqp1+dqm1)* (1.0 + sign(1.0, dqm1*dqp1))/4.

             decr2 = (u(iecr,i,j,k) - u(iecr,i,j-1,k))/dy

             if(nzd .gt. 1) then
                dqm2  = 0.5*((u(iecr,i ,j-1,k  ) + u(iecr,i ,j  ,k  ))   &
                            -(u(iecr,i ,j-1,k-1) + u(iecr,i ,j  ,k-1)))/dz
                dqp2  = 0.5*((u(iecr,i ,j-1,k+1) + u(iecr,i ,j  ,k+1))   &
                            -(u(iecr,i ,j-1,k  ) + u(iecr,i ,j  ,k  )))/dz

                decr3 = (dqp2+dqm2)* (1.0 + sign(1.0, dqm2*dqp2))/4.
             endif


             fcrdif2 = K_cr_paral * n2b * &
                      (n1b*decr1 + n2b*decr2 + n3b*decr3) &
                      + K_cr_perp * decr2

             wa(i,j,k) = - fcrdif2 * dt

10          continue
20        continue
30      continue

      do 60 k=ks,ke
        do 50 j=js,je
          do 40 i=is,ie
            u(iecr,i,j,k) = u(iecr,i,j,k)    &
                      - (wa(i,j+1,k)-wa(i,j,k))/dy
40        continue
50      continue
60    continue

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
      real :: decr1, decr2, decr3, fcrdif3
      real :: dqp2, dqm2, dqp3, dqm3, dqm1, dqp1

!=======================================================================
!
   if(nzd .gt. 1) then
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

             dqm1  = 0.5*((u(iecr,i  ,j ,k-1) + u(iecr,i  ,j  ,k ))   &
                         -(u(iecr,i-1,j ,k-1) + u(iecr,i-1,j  ,k )))/dx
             dqp1  = 0.5*((u(iecr,i+1,j ,k-1) + u(iecr,i+1,j  ,k ))   &
                         -(u(iecr,i  ,j ,k-1) + u(iecr,i  ,j  ,k )))/dx

             decr1 = (dqp1+dqm1)* (1.0 + sign(1.0, dqm1*dqp1))/4.

             dqm2  = 0.5*((u(iecr,i ,j  ,k-1) + u(iecr,i  ,j  ,k ))   &
                         -(u(iecr,i ,j-1,k-1) + u(iecr,i  ,j-1,k )))/dy
             dqp2  = 0.5*((u(iecr,i ,j+1,k-1) + u(iecr,i  ,j+1,k ))   &
                         -(u(iecr,i ,j  ,k-1) + u(iecr,i  ,j  ,k )))/dy

             decr2 = (dqp2+dqm2)* (1.0 + sign(1.0, dqm2*dqp2))/4.

             decr3 =  (u(iecr,i,j,k    ) - u(iecr,i,j,  k-1))/dz

             fcrdif3 = K_cr_paral * n3b * &
                  (n1b*decr1 + n2b*decr2 + n3b*decr3) &
                   + K_cr_perp * decr3

             wa(i,j,k) = - fcrdif3 * dt

100   continue

      do 60 k=ks,ke
        do 50 j=js,je
          do 40 i=is,ie
            u(iecr,i,j,k) = u(iecr,i,j,k)  &
                        -(wa(i,j,k+1)-wa(i,j,k))/dz
40        continue
50      continue
60    continue

      endif

      return
      end subroutine cr_diff_z

  subroutine div_v
  
    use initionized, only : idni,imxi,imyi,imzi
  
    implicit none
    real, dimension(nx)  :: dvx,dvy,dvz, tmp
    integer i,j,k

    wa(:,:,:) = 0.0
    if(nzd .eq. 1) then
      k=1
      do j=2,ny-1
        tmp=u(imxi,:,j,k)/u(idni,:,j,k)
        wa(:,j,k) = (eoshift(tmp,1) - eoshift(tmp,-1))/(2.*dx)
        tmp=(u(imyi,:,j+1,k)/u(idni,:,j+1,k)-u(imyi,:,j-1,k)/u(idni,:,j-1,k))/(2.*dy)
        wa(:,j,k) = wa(:,j,k) + tmp
      enddo
    else if (nzd .gt. 1) then
      do k=2,nz-1
        do j=2,ny-1
          tmp=u(imxi,:,j,k)/u(idni,:,j,k)
          wa(:,j,k) = (eoshift(tmp,1) - eoshift(tmp,-1))/(2.*dx)
          tmp=(u(imyi,:,j+1,k)/u(idni,:,j+1,k)-u(imyi,:,j-1,k)/u(idni,:,j-1,k))/(2.*dy)
          wa(:,j,k) = wa(:,j,k) + tmp
          tmp=(u(imzi,:,j,k+1)/u(idni,:,j,k+1)-u(imzi,:,j,k-1)/u(idni,:,j,k-1))/(2.*dz)
          wa(:,j,k) = wa(:,j,k) + tmp
        enddo
      enddo
    endif
  end subroutine div_v


end module crdiffusion
