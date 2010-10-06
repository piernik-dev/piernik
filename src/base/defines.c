/*

  The diff_nml macro is a workaround (DIRTY HACK)
  We would like to pass a namelist name to compare_namelist subroutine, but Fortran currently does not permit it.


  Requires:
    use mpisetup, only : cwd
    use func,     only : compare_namelist
    integer :: ierrh
    par_file = trim(cwd)//'/problem.par'


  Potential problem: for gnu cpp I'd write:

    call namelist_errh(ierrh, #namelist);

  but it looks like gfortran calls cpp -traditional-cpp and does not perform stringifying

 */
#define diff_nml(namelist)\
  open(501, file="temp1.dat", status="unknown");\
  write(501,nml=namelist);\
  close(501);\
  open(1, file=par_file);\
  read(unit=1, nml=namelist, iostat=ierrh);\
  call namelist_errh(ierrh, "\
namelist\
");\
  close(1);\
  open(502, file="temp2.dat", status="unknown");\
  write(502,nml=namelist);\
  close(502);\
  call compare_namelist("temp1.dat", "temp2.dat", cwd)

/*
   Here we check for consistency of the piernik.def file.

   This is a .c file to make sure that #error directive stops the make process.
*/

/*
   Flux limiters

   Exclusive: VANLEER, MONCEN, MINMOD, SUPERBEE

   Used by: src/fluids/fluxes.F90
*/

#include "piernik.def"

#undef FLIM
#undef FLIM2

#ifdef VANLEER
#define FLIM
#endif

#ifdef MONCEN
#ifdef FLIM
#define FLIM2
#else
#define FLIM
#endif
#endif

#ifdef MINMOD
#ifdef FLIM
#define FLIM2
#else
#define FLIM
#endif
#endif

#ifdef SUPERBEE
#ifdef FLIM
#define FLIM2
#else
#define FLIM
#endif
#endif

#ifdef FLIM2
#error Too many flux limiters have been defined. Only one is allowed.
#endif

#ifndef FLIM
#error No flux limiters have been defined.
#endif

/*
  Freezing speed

  Exclusive: LOCAL_FR_SPEED, GLOBAL_FR_SPEED

  Used by: src/fluids/ionized/fluxionized.F90, src/fluids/neutral/fluxneutral.F90, src/fluids/dust/fluxdust.F90
*/

#undef FR_SPEED
#undef FR_SPEED2

#ifdef LOCAL_FR_SPEED
#define FR_SPEED
#endif

#ifdef GLOBAL_FR_SPEED
#ifdef FR_SPEED
#define FR_SPEED2
#else
#define FR_SPEED
#endif
#endif

#ifdef FR_SPEED2
#error Both freezing speeds defined
#endif

#ifndef FR_SPEED
#error No freezing speed defined.
#endif

/*
  Units

  Exclusive: PGM, SSY, SI, CGS, WT4, PSM, PLN, KSG, KSM

  Code generated automatically by:
  ./setup -units | grep "\!" | awk '{print $2}' |\
  awk 'BEGIN {\
    def="UNIT";\
    n=0;\
  } {\
    print "#ifdef", $1;\
    if (n>0) print "#ifdef",def,"\n#define",def"2\n#else";\
    print "#define",def;\
    if (n>0) print "#endif";\
    print "#endif\n";\
    a[n++]=$1;\
  } END {\
    printf "%s %s", "#ifdef",def"2\n#error Only one of { ";\
    for (b in a) printf "%s ", a[b];\
    printf "%s %s %s","} is allowed.\n#endif\n\n#ifndef",def,"\n#error None of { ";\
    for (b in a) printf "%s ", a[b];\
    printf "%s\n","} were defined.\n#endif";\
  }'
*/

#ifdef PSM
#define UNIT
#endif

#ifdef PLN
#ifdef UNIT
#define UNIT2
#else
#define UNIT
#endif
#endif

#ifdef KSG
#ifdef UNIT
#define UNIT2
#else
#define UNIT
#endif
#endif

#ifdef KSM
#ifdef UNIT
#define UNIT2
#else
#define UNIT
#endif
#endif

#ifdef PGM
#ifdef UNIT
#define UNIT2
#else
#define UNIT
#endif
#endif

#ifdef SSY
#ifdef UNIT
#define UNIT2
#else
#define UNIT
#endif
#endif

#ifdef SI
#ifdef UNIT
#define UNIT2
#else
#define UNIT
#endif
#endif

#ifdef CGS
#ifdef UNIT
#define UNIT2
#else
#define UNIT
#endif
#endif

#ifdef WT4
#ifdef UNIT
#define UNIT2
#else
#define UNIT
#endif
#endif

#ifdef UNIT2
#error Only one of { PGM SSY SI CGS WT4 PSM PLN KSG KSM } is allowed
#endif

#ifndef UNIT
#warning None of { PGM SSY SI CGS WT4 PSM PLN KSG KSM } were defined.
#warning Assuming SCALED
#endif

/* basic sanity check for isothermal fluid */

#ifdef ISO_LOCAL
#  ifndef ISO
#     error ISO must be defined with ISO_LOCAL
#  endif
#  ifndef IONIZED
#     error ISO_LOCAL currently works only with ionized fluid
#  endif
#endif

/* at least one of { ionized, neutral, dust } must be defined */

#undef FLUID

#ifdef IONIZED
#define FLUID
#endif

#ifdef DUST
#define FLUID
#endif

#ifdef NEUTRAL
#define FLUID
#endif

#ifndef FLUID
#error None of { IONIZED DUST NEUTRAL } were defined.
#endif

#ifdef GRAV
/*
  Constant gravity field

  Exclusive: GRAV_NULL, GRAV_UNIFORM, GRAV_LINEAR, GRAV_PTMASS, GRAV_PTMASSSTIFF, GRAV_PTMASSPURE, GRAV_PTFLAT, GRAV_USER

  Used by: src/gravity/gravity.F90
*/

#ifdef GRAV_NULL
#define GRVC
#endif

#ifdef GRAV_UNIFORM
#ifdef GRVC
#define GRVC2
#else
#define GRVC
#endif
#endif

#ifdef GRAV_LINEAR
#ifdef GRVC
#define GRVC2
#else
#define GRVC
#endif
#endif

#ifdef GRAV_PTMASS
#ifdef GRVC
#define GRVC2
#else
#define GRVC
#endif
#endif

#ifdef GRAV_PTMASSSTIFF
#ifdef GRVC
#define GRVC2
#else
#define GRVC
#endif
#endif

#ifdef GRAV_PTMASSPURE
#ifdef GRVC
#define GRVC2
#else
#define GRVC
#endif
#endif

#ifdef GRAV_PTFLAT
#ifdef GRVC
#define GRVC2
#else
#define GRVC
#endif
#endif

#ifdef GRAV_USER
#ifdef GRVC
#define GRVC2
#else
#define GRVC
#endif
#endif

#ifdef GRVC2
#error Only one of { GRAV_PTMASSSTIFF GRAV_PTMASSPURE GRAV_PTFLAT GRAV_USER GRAV_NULL GRAV_UNIFORM GRAV_LINEAR GRAV_PTMASS } is allowed.
#endif

#ifndef GRVC
#warning None of { GRAV_PTMASSSTIFF GRAV_PTMASSPURE GRAV_PTFLAT GRAV_USER GRAV_NULL GRAV_UNIFORM GRAV_LINEAR GRAV_PTMASS } were defined. Relying on grav_accel.
/*
  ToDo: add GRAV_ACC symbol then change the #warning above to #error
*/
#endif

#endif /* GRAV */

/*
  Gravity solvers

  Exclusive: MULTIGRID, POISSON_FFT
*/

#if defined(MULTIGRID) && defined(POISSON_FFT)
#error MULTIGRID and POISSON_FFT are not meant to work together
#endif
