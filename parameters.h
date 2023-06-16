

/****************************************************************/
/* >>> STANDARD PARAMETERS <<< */

/* The UV background file (should be placed in ./uvb_models/ directory).
   Used as the default, unless the PLAW_UVB flag is raised in the
   Makefile. A selection of different models are available.  */
#define UVBFILE "puchwein2019.txt"

/* The output file (appears in ./outputs/ directory) */
//#define OUTFILE "Plaw_J0.1_zH6.5_zHe3.0_a3.5_d1.0.dat"
#define OUTFILE "puchwein2019_m8.0e-14_e5.0e-15.dat"

#define ZSTART     20.0  /* initial redshift */
#define ZEND       0.0   /* final redshift */
#define LOGOVERDEN 0.0   /* the gas overdensity, log(rho/<rho>)*/

/* Total number of time steps, equally spaced in ln(a) */
#define NT 10000 

/* Cosmological parameters */
#define OMEGAM 0.308
#define OMEGAL 0.692
#define OMEGAB 0.0482
#define HUBBLE 0.678

/* Hydrogen mass fraction */
#define XH     0.76


/****************************************************************/
/* >>> DARK_PHOTON PARAMETERS <<< */

/* Extra parameters if optionally running with the DARK_PHOTON flag
   raised in the Makefile.  Calculates the heating from the resonant
   conversion of dark photons.   See e.g.

   Caputo et al. 2020, PhRvD, 102, 103533                                       
   Caputo et al. 2020, PRL, 125, 221303   */


#define MASS_DPHOT 8.0e-14  /* dark photon mass [eV/c^2] */
#define EPS_DPHOT  5.0e-15 /* kinetic mixing parameter */


/****************************************************************/
/* >>> PLAW_UVB PARAMETERS <<< */

/* Extra parameters if optionally running with the PLAW_UVB flag
   raised in the Makefile.  This calculates photoionisation and
   heating rates directly from the single power-law input spectrum
   specified below using Gauss-Legendre quadrature */

#define J22       0.1   /* Specific intensity at the HI ionisation
			   edge, [10^-22 erg s-1 cm^-2 sr^-1 Hz^-1] */
#define ALPHA_UV  1.17  /* Power-law spectral index */
#define ZRH1      6.5   /* Redshift of HI (and HeI) reionisation */
#define ZRHE2     3.0   /* Redshift of HeII reionisation */

/* Photon energy integrals */
#define NWEIGHTS 1000   /* Number of Gauss-Legendre weights used */

/****************************************************************/
