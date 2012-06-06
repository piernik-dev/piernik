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
!! \b \#define \b "ISO_LOCAL" - to choose isothermal fluids with locally specified soundspeed and temperature
!!
!! \b \#define \b "MAGNETIC" - to include parts of the code related to magnetic field (otherwise magnetic field is switched off even if ionized fluid is simulated)
!!
!! \b \#define \b "RESISTIVE" - to include resistive dissipation of magnetic field
!!
!! \b \#define \b "FLUID_INTERACTIONS" - to take into account interactions between different fluids (eg. drag force)
!!
!! \b \#define \b "FLUID_INTERACTIONS_DW" - older type of fluids interaction scheme
!!
!! \b \#define \b "BALSARA"
!!
!! \b \#define \b "GRAV" - to include gravitational forces
!!
!! If no %gravity type is defined then vertical component of galactic %gravity is used:
!!
!! local Galactic %gravity only in z-direction (see <a href="http://cdsads.u-strasbg.fr/abs/1998ApJ...497..759F">Ferriere K., 1998, Astrophys. Journal, 497, 759</a>)\n
!! \f[
!! F_z = 3.23 \cdot 10^8 \cdot \left[\left(-4.4 \cdot 10^{-9} \cdot exp\left(-\frac{(r_{gc}-r_{gc_{}Sun})}{(4.9kpc)}\right) \cdot \frac{z}{\sqrt{(z^2+(0.2kpc)^2)}}\right)-\left( 1.7 \cdot 10^{-9} \cdot \frac{(r_{gc_{}Sun}^2 + (2.2kpc)^2)}{(r_{gc}^2 + (2.2kpc)^2)} \cdot \frac{z}{1kpc}\right) \right]
!! \f]
!! where \f$r_{gc}\f$ is galactocentric radius and \f$r_{gcSun}\f$ is the galactocentric radius of Sun.
!! </td></tr></table>
!! \n \n
!!
!! \b \#define \b "SELF_GRAV"
!!
!! \b \#define \b "VARIABLE_GP"
!!
!! \b \#define \b "CORIOLIS"
!!
!! \b \#define \b "SN_SRC"
!! @n @n
!!
!! \b NUMERICAL \b SCHEMES \b AND \b SOLUTIONS:
!!
!! \b \#define \b "RTVD"  - to choose Relaxing TVD scheme (then HLLC must be undefined)
!!
!! \b \#define \b "HLLC" - to choose HLLC scheme (without magnetic fields; RTVD must be undefined)
!!
!! \b \#define \b "LOCAL_FR_SPEED"
!!
!! \b \#define \b "GLOBAL_FR_SPEED"
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
!! \b \#define \b "NEW_HYDROSTATIC"
!!
!! \b \#define \b "ZERO_BND_EMF"
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
!! \b \#define \b "INDEPENDENT_ATOUTB"
!!
!! \b \#define \b "PERFMON"
!!
!! \b \#define \b "PIERNIK_OPENCL"
!!
!! \b \#define \b "USER_RULES"
!!
!! \b \#define \b "__INTEL_COMPILER" - automatic flag placed by Intel Compiler
!!
!! \b \#define \b "__PGI" - flag for Portland Compilers
!!
!<


