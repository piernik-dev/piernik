! $Id$
#include "piernik.def"

module ionizeds

contains

subroutine flux_ionized(iflux,icfr,iuu,bb,n)
use arrays, only : nui
implicit none
integer n
real,dimension(nui,n) :: iflux,icfr,iuu
real,dimension(3,n)   :: bb


end subroutine flux_ionized

subroutine src_ionized
implicit none
end subroutine src_ionized

end module ionizeds
