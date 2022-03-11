

/* Maximum energy [eV] used in the integrals */
#define EMAX 1.0e4

/* Secondary electron energy table in log(E/eV) */
#define NESECTAB 1000
#define ESECMIN  1.0
#define ESECMAX  4.0

/* Electron fraction table in log xe = log(ne/ne_max). Used for
   secondary ionisations by fast electrons */
#define NXTAB 100
#define XMIN -4.0
#define XMAX  0.0

double *weight, *absc;
double gJH0,gJHe0,gJHep,epsH0,epsHe0,epsHep;
double redshift,Ha;

/* Define physical constants, cgs units
https://physics.nist.gov/cuu/Constants/index.html  */
#define GRAVITY       6.67430e-8
#define BOLTZMANN     1.380649e-16
#define C             2.99792458e10
#define PLANCK        6.62607015e-27
#define PROTONMASS    1.67262192369e-24
#define ELECTRONVOLT  1.602176634e-12

/* Derived from IAU definition of astronomical unit:
   https://www.iau.org/public/themes/measuring/
   https://www.iau.org/static/resolutions/IAU2012_English.pdf */
#define MPC 3.085677581e24

#define GAMMA_MINUS1  (2.0/3.0)

/* Ionisation thresholds in eV 
   https://physics.nist.gov/PhysRefData/ASD/ionEnergy.html */
#define EION_H1  13.598434599702
#define EION_HE1 24.587389011 
#define EION_HE2 54.417765486
