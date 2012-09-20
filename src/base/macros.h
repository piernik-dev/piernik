#define _MACROS_H_VERSION_STR '$Id$'
/*

  The diff_nml macro is a workaround (DIRTY HACK)
  We would like to pass a namelist name to compare_namelist subroutine, but Fortran currently does not permit it.


  Requires:
    use dataio_pub, only: nh ! QA_WARN required for diff_nml


  It does not work with cpp -traditional-cpp (the default way gfortran calls cpp).

 */
#define diff_nml(namelist)\
  if (.not.nh%initialized) call nh%init();\
  open(newunit=nh%lun, file="temp1.dat", status="unknown");\
  write(nh%lun,nml=namelist);\
  close(nh%lun);\
  open(newunit=nh%lun, file=nh%par_file);\
  nh%errstr="";\
  read(unit=nh%lun, nml=namelist, iostat=nh%ierrh, iomsg=nh%errstr);\
  close(nh%lun);\
  call nh%namelist_errh(nh%ierrh, #namelist);\
  read(nh%cmdl_nml,nml=namelist, iostat=nh%ierrh);\
  call nh%namelist_errh(nh%ierrh, #namelist, .true.);\
  open(newunit=nh%lun, file="temp2.dat", status="unknown");\
  write(nh%lun,nml=namelist);\
  close(nh%lun);\
  call nh%compare_namelist("temp1.dat", "temp2.dat")
