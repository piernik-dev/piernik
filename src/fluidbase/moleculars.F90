! $Id$
#include "piernik.def"

module moleculars

contains

subroutine flux_molecular(mflux,mcfr,muu,n)
use arrays, only : num
implicit none
integer n
real,dimension(num,n) :: mflux,mcfr,muu

end subroutine flux_molecular

subroutine src_molecular
implicit none
end subroutine src_molecular

end module moleculars
