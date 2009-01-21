module types
   type :: hdf
      integer :: nhdf, ntsl, nres, nlog, step_hdf, step_res, log_lun, &
         nstep
      real    :: last_hdf_time
      character(len=128) :: log_file
   end type hdf

end module types
