! $Id$
#include "piernik.def"

module dusts

contains

subroutine flux_dust(dflux,dcfr,duu,n)
use arrays, only : nud
implicit none
integer n
real,dimension(nud,n) :: dflux,dcfr,duu

end subroutine flux_dust

subroutine src_dust
implicit none
end subroutine src_dust

end module dusts
