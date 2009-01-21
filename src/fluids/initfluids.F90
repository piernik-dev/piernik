! $Id$
#include "piernik.def"

module initfluids
  
#ifdef IONIZED    
  use initionized, only : init_ionized
#endif /* IONIZED */  

#ifdef NEUTRAL    
  use initneutral, only : init_neutral
#endif /* NEUTRAL */  

#ifdef DUST    
  use initdust, only : init_dust
#endif /* DUST */  

#ifdef COSM_RAYS    
  use initcosmicrays, only : init_cosmicrays
#endif /* COSM_RAYS */  

  contains

  subroutine init_fluids
  

#ifdef IONIZED
  call init_ionized
#endif /* IONIZED */  

#ifdef NEUTRAL
  call init_neutral
#endif /* NEUTRAL */  

#ifdef DUST
  call init_dust
#endif /* DUST */  

#ifdef COSM_RAYS
  call init_cosmicrays
#endif /* COSM_RAYS */  


  end subroutine init_fluids

end module initfluids
