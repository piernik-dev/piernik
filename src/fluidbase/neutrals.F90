! $Id$
#include "piernik.def"

module neutrals

contains

subroutine flux_neutral(nflux,ncfr,nuu,n)
use arrays, only : nun
implicit none
integer n
real,dimension(nun,n) :: nflux,ncfr,nuu

end subroutine flux_neutral

subroutine src_neutral
implicit none
end subroutine src_neutral

end module neutrals
