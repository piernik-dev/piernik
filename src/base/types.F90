!$Id$
#include "piernik.def"
module types
   type :: hdf
      integer :: nhdf, ntsl, nres, nlog, step_hdf, step_res, log_lun, &
         nstep
      real    :: last_hdf_time
      character(len=128) :: log_file
   end type hdf

   type :: value
      real    :: val
      integer, dimension(3) :: loc
      integer :: proc
   end type value

   type :: tsl_container
#ifdef NEUTRAL 
      real :: denn_min, denn_max, vxn_max, vyn_max, vzn_max, &
              pren_min, pren_max, temn_min, temn_max, csn_max
#endif /* NEUTRAL */

#ifdef DUST
      real :: dend_min, dend_max, vxd_max, vyd_max, vzd_max
#endif /* DUST */

#ifdef COSM_RAYS
      real :: encr_min, encr_max
#endif /* COSM_RAYS */

#ifdef RESISTIVE
      real :: etamax
#endif /* RESISTIVE */

#ifdef MAGNETIC
    real :: b_min, b_max, divb_max
#endif /* MAGNETIC */ 

#ifdef IONIZED 
    real :: deni_min, deni_max, vxi_max, vyi_max, vzi_max, &
            prei_min, prei_max, temi_min, temi_max, vai_max, csi_max      
#endif /* IONIZED */


   end type tsl_container

end module types
