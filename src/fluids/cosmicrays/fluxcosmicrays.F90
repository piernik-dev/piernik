! $Id: 
#include "piernik.def"
#define RNG 2:n-1


module fluxcosmicrays
  implicit none

  contains 
!==========================================================================================

  subroutine flux_crs(fluxc,vion,uuc,n)
  
    use constants,       only : small
    use fluidindex,      only : nvar_crs,icr   
  
    implicit none
    integer n

! locals
    real, dimension(nvar_crs,n):: fluxc,uuc
    real, dimension(n) :: vion 
    
    fluxc   = 0.0

    fluxc(icr,RNG)= uuc(icr,RNG)*vion(RNG)

  end subroutine flux_crs


end module fluxcosmicrays
