!>
!! \page precompiler Precompiler directives
!!
!! The code uses precompiler directives. All of them are supposed to be set in the file \b piernik.def in the problem directory.
!!
!! Precompiler directives:<ul>
!! <li>are used to select parts of PIERNIK code, which are relevant for given problem</li>
!! <li>instruct setup on which routines should be selected for compilation in ./obj directory</li>
!! <li>instruct precompiler on which parts of selected routines should be retained for compilation</li>
!! </ul>
!!
!! @n @n
!!
!! \b FLUID \b COMPONENTS:
!!
!! \b \#define \b "IONIZED" - to include parts of the code related to ionized fluid
!!
!! \b \#define \b "NEUTRAL" - to include parts of the code related to neutral fluid
!!
!! \b \#define \b "DUST" - to include parts of the code related to dust fluid
!!
!! \b \#define \b "COSM_RAYS" - to include parts of the code related to cosmic ray component
!!
!! \b \#define \b "COSM_RAYS_SOURCES" - to take into account several different species of cosmic rays particles
!!
!! \b \#define \b "TRACER" - to include parts of the code related to tracer
!! @n @n
!!
!! \b PHYSICS \b INCLUDED \b (i.e. \b FIELDS, \b INTERACTIONS):
!!
!! \anchor isothermal
!! \b \#define \b "ISO" - the choice between adiabatic and isothermal fluids (if defined then isothermal fluid is set, otherwise adiabatic)
!!
!! \b \#define \b "ISO_LOCAL" - to choose isothermal fluids with locally specified soundspeed and temperature (required ISO; currently works only with IONIZED)
!!
!! \b \#define \b "MAGNETIC" - to include parts of the code related to magnetic field (otherwise magnetic field is switched off even if ionized fluid is simulated)
!!
!! \b \#define \b "RESISTIVE" - to include resistive dissipation of magnetic field
!!
!! \b \#define \b "FLUID_INTERACTIONS" - to take into account interactions between different fluids (eg. drag force)
!!
!! \b \#define \b "BALSARA"
!!
!! \b \#define \b "GRAV" - to include gravitational forces. The type of gravity is govern by external_gp value. Then chosen gravitational potential is computed.
!! If no %gravity type is defined by external_gp character value then grav_accel defined by user is used and from given gravitational acceleration potential is computed.
!!
!! \b \#define \b "SELF_GRAV" - switch on self-gravity
!!
!! \b \#define \b "VARIABLE_USER_GP" - user flag to switch on gravity maxima dumps to log and tsl files (for SELF_GRAV that is as default)
!!
!! \b \#define \b "CORIOLIS" - switch on Coriolis force (RTVD only)
!!
!! \b \#define \b "SN_SRC" - switch on supernovae insert (RTVD only; hooking to problem_customize_solution is strongly recommended)
!! @n @n
!!
!! \b NUMERICAL \b SCHEMES \b AND \b SOLUTIONS:
!!
!! \b \#define \b "RTVD"  - to choose Relaxing TVD scheme (then HLLC must be undefined)
!!
!! \b \#define \b "HLLC" - to choose HLLC scheme (without magnetic fields; RTVD must be undefined)
!!
!! \b \#define \b "LOCAL_FR_SPEED" - choose locally computed freezing speed (i.e. in each cell)
!!
!! \b \#define \b "GLOBAL_FR_SPEED" - choose globally computed freezing speed (i.e. constant for the whole domain)
!!
!! \b \#define \b "MULTIGRID" - to include self–gravity (multigrid solver, uses FFT when possible, recommended)
!!
!! \b \#define \b "POISSON_FFT" - to include self–gravity (pure FFT solver)
!!
!! \b \#define \b "FFTW"
!!
!! \b \#define \b "SHEAR" - to include Coriolis and tidal forces in gas equation of motion and use  <a href="http://cdsads.u-strasbg.fr/abs/1995ApJ...440..742H">Hawley, Gammie and Balbus (1995)</a> type approach for shearing BCS
!!
!! \b \#define \b "SHEAR_BND"
!!
!! \b \#define \b "NEW_HYDROSTATIC" - change hydrostatic_zeq scheme and find midplane corresponding to zero gravity (practically the least gravity absolute value)
!!
!! \b \#define \b "ZERO_BND_EMF" - switch to putting 0 on emf boundaries for outflow type boundary conditions instead of constant emf gradient as default
!! @n @n
!!
!! \b OTHERS:
!!
!! \b \#define \b "DEBUG" - enable modules: piernikdebug (src/base/debug.F90) and piernikiodebug (src/hdf5/io_debug.F90)
!!
!! \b \#define \b "VERBOSE" - print additional diagnostic information on stdout
!!
!! \b \#define \b "PGPLOT"
!!
!! \b \#define \b "JEANS_PROBLEM" - switch on jeans problem-specific quirks
!!
!! \b \#define \b "INDEPENDENT_ATOUTB" - switch to independent (instead of collective) write to files (restart files) for arrays with outer boundary area type
!!
!! \b \#define \b "PERFMON"
!!
!! \b \#define \b "PIERNIK_OPENCL"
!!
!! \b \#define \b "USER_RULES" - include instructions from user_rules.h file from problem directory at the beginning of compilation
!!
!! \b \#define \b "__INTEL_COMPILER" - automatic flag placed by Intel Compiler
!!
!! \b \#define \b "__PGI" - flag for Portland Compilers
!!
!<


