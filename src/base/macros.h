/* $Id$ */
/*

  The diff_nml macro is a workaround (DIRTY HACK)
  We would like to pass a namelist name to compare_namelist subroutine, but Fortran currently does not permit it.


  Requires:
    use dataio_pub, only: par_file, ierrh, namelist_errh, compare_namelist  ! QA_WARN


  It does not work with cpp -traditional-cpp (the default way gfortran calls cpp).

 */
#define diff_nml(namelist)\
  open(501, file="temp1.dat", status="unknown");\
  write(501,nml=namelist);\
  close(501);\
  open(1, file=par_file);\
  read(unit=1, nml=namelist, iostat=ierrh);\
  call namelist_errh(ierrh, #namelist);\
  close(1);\
  open(502, file="temp2.dat", status="unknown");\
  write(502,nml=namelist);\
  close(502);\
  call compare_namelist("temp1.dat", "temp2.dat")
