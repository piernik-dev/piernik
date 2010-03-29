/*---------------------------------------- gslwr.c --*/
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

static gsl_rng* r;            /*  r is (almost) never used explicitly and does
                                  not need to be seen by the Fortran program.  */

void rng_initialise_(int* s){
   r = gsl_rng_alloc(gsl_rng_taus);    /*  constant gsl_rng_taus is unknown to Fortran  */
   gsl_rng_set(r, (unsigned long int)(*s));  /*  s is cast from int to (unsigned long int)  */
}

void rng_sample_gaussian_(double* x, int* n, double* sigma){
   int i;
   for(i=0; i<*n; i++)
      x[i] = gsl_ran_gaussian(r,*sigma);
}


