
#undef DEBUG
#undef GNU
#undef VERBOSE

/*======= re-setup needed ==========*/
/* choosing scheme */
#undef UNSPLIT_1D      /*scheme*/
#undef  UNSPLIT_3D      /*scheme*/
#define  SPLIT           /*scheme*/

/* choosing solver */
#define  ORIG            /*solver*/
#undef  SSP            /*solver*/

/* choosing unitsystem */
#define STANDARD		/*genuine unitsystem*/
#undef PSM		/*pc,Msun,Myr unitsystem*/
#undef PGM		/*pc,G=1,Myr unitsystem*/
#undef DMY		/*cm^16,Msun,year unitsystem*/
#undef CGS		/*cm,gram,sek unitsystem*/
#undef SI		/*metr,kg,sek unitsystem*/
/*===================================*/

/* adiabatic or isothermal */
#define ISO

/* choosing freezing speed */
#define LOCAL_FR_SPEED
#undef  GLOBAL_FR_SPEED
#define OLDLOCALCFR

/* choosing flux limiter */
#define VANLEER
#undef MONCEN
#undef MINMOD
#undef SUPERBEE

/* physics including */
#define RESIST 		
#undef VISC			
#undef HEAT_COND	
#undef COOL_HEAT
#undef COSM_RAYS

/* gravity type choosing */
#define GRAV 
#define GRAV_NULL
#undef SELF_GRAV

/* structure features including */
#undef GALAXY
#undef GALACTIC_DISK
#undef ARMS_POTENTIAL
#undef SNE_DISTR
#undef MASS_COMPENS
#undef SN_SRC

/* numerical solutions and limitations */
#undef SHEAR
#undef SHEAR_MPI
#undef FLX_BND
#undef VZ_LIMITS

/* problem specific directives */
