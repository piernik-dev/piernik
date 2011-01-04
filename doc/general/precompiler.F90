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
!! \anchor isothermal
!! \b \#define \b ISO - the choice between adiabatic and isothermal fluids (if defined then isothermal fluid is set, otherwise adiabatic)
!!
!! \b \#define \b IONIZED - to include parts of the code related to ionized fluid
!!
!! \b \#define \b NEUTRAL - to include parts of the code related to neutral fluid
!!
!! \b \#define \b DUST - to include parts of the code related to dust fluid
!!
!! \b \#define \b COSM_RAYS - to include parts of the code related to cosmic ray component
!!
!! \b \#define \b MAGNETIC - to include parts of the code related to magnetic field (otherwise magnetic field is switched off even if ionized fluid is simulated)
!!
!! \b \#define \b FLUID_INTERACTIONS - to take into account interactions between different fluids (eg. drag force)
!!
!! \b \#define \b RESIST - to include resistive dissipation of magnetic field
!!
!! \b \#define \b GRAV - to include gravitational forces
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
!! \b \#define \b MULTIGRID - to include self–gravity (multigrid solver, uses FFT when possible, recommended)
!! \b \#define \b POISSON_FFT - to include self–gravity (pure FFT solver)
!!
!! \b \#define \b SHEAR - to include Coriolis and tidal forces in gas equation of motion and use  <a href="http://cdsads.u-strasbg.fr/abs/1995ApJ...440..742H">Hawley, Gammie and Balbus (1995)</a> type approach for shearing BCS
!! @n @n
!!
!! There is a list of precompiler directives to choose unit system used in simulation:
!!
!! <table border="0">
!! <tr><td width="40pt"><b></b></td><td width="760pt">
!! \copydetails constants
!! </td></tr></table>
!! \n \n
!<
