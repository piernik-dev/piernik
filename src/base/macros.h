#define _MACROS_H_VERSION_STR "$Id$"
/*

  The diff_nml macro is a workaround (DIRTY HACK)
  We would like to pass a namelist name to compare_namelist subroutine, but Fortran currently does not permit it.


  Requires:
    use dataio_pub, only: par_file, ierrh, namelist_errh, compare_namelist, cmdl_nml, lun, get_lun  ! QA_WARN required for diff_nml


  It does not work with cpp -traditional-cpp (the default way gfortran calls cpp).

 */
#define diff_nml(namelist)\
  lun = getlun();\
  open(lun, file="temp1.dat", status="unknown");\
  write(lun,nml=namelist);\
  close(lun);\
  lun = getlun();\
  open(lun, file=par_file);\
  read(unit=lun, nml=namelist, iostat=ierrh);\
  close(lun);\
  call namelist_errh(ierrh, #namelist);\
  read(cmdl_nml,nml=namelist, iostat=ierrh);\
  call namelist_errh(ierrh, #namelist, .true.);\
  lun = getlun();\
  open(lun, file="temp2.dat", status="unknown");\
  write(lun,nml=namelist);\
  close(lun);\
  call compare_namelist("temp1.dat", "temp2.dat")
