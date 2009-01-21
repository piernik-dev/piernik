! $Id: 
#include "piernik.def"

module timestepcosmicrays

  real :: dt_crs

 contains 

  subroutine timestep_crs

    use grid,           only : dxmn
    use initcosmicrays, only : cfl_cr,K_cr_paral,K_cr_perp 

    implicit none

    real dt_crs_proc, dt_crs_all
    
      dt_crs = cfl_cr * 0.5*dxmn**2/(K_cr_paral+K_cr_perp+small)

  end subroutine timestep_crs

!-------------------------------------------------------------------------------
end module timestepcosmicrays

