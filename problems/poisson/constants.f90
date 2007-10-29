module constants

  real, parameter :: pi = 3.141592653589793238
  real, parameter :: dpi = 2.*pi
  real, parameter :: fpi = 4.*pi
  real, parameter :: e  = 2.718281828459045235
  real, parameter :: cmps2 = 3.23e8
  real, parameter :: kpc = 1000.0
  real, parameter :: r_gc_sun = 8.5 ! kpc 
  
  real, parameter :: G_cgs = 6.6725985E-8  ! (cgs)
  real, parameter :: G_one = 1.0e-4  ! (scaled units)


  real, parameter :: fpiG  = fpi*G_one
  real, parameter :: newtong = G_one
end module constants
