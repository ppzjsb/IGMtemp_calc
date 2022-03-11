#define fnID 2
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "parameters.h"
#include "proto.h"
#include "global_vars.h"


#include <cvode/cvode.h>
#include <nvector/nvector_serial.h>
#include <cvode/cvode_dense.h>

#define Ith(v,i) NV_Ith_S(v,i-1)

static int f_solve_rates(realtype dt, N_Vector y, N_Vector ydot, void *user_data);


double nHp,nHep,nHepp,nH0,nHe0,ne,necgs;
double aHp,aHep,aHepp,ad,geH0,geHe0,geHep;
double bH0,bHe0,bHep,bff;
double logT,yhelium,nHcgs,ratefact;

static double deltaT;
static double *AlphaHp, *AlphaHep, *AlphaHepp, *Alphad;
static double *BetaH0, *BetaHe0, *BetaHep, *Betaff;
static double *GammaeH0, *GammaeHe0, *GammaeHep;
static double logTmin, logTmax;

#define SMALLNUM 1.0e-60
#define TMIN 1.0
#define TMAX 1.0e9
#define NCOOLTAB 2000



/* External UVB table */
#define TABLESIZE 2000

static float inlogz[TABLESIZE];
static float gH0[TABLESIZE], gHe[TABLESIZE], gHep[TABLESIZE];
static float eH0[TABLESIZE], eHe[TABLESIZE], eHep[TABLESIZE];
static float log_gH0[TABLESIZE], log_gHe[TABLESIZE], log_gHep[TABLESIZE];
static float log_eH0[TABLESIZE], log_eHe[TABLESIZE], log_eHep[TABLESIZE];
static int nheattab;		/* length of table */



/******************************************************************************************/

/* Solve the non-equilibrium network using SUNDIALS */

double DoCooling(double u,double rho,double dt,double *nHp_guess,double *nHep_guess,double *nHepp_guess)
{
  
  if (dt > 0)
    u = solve_rates(u, rho, dt, nHp_guess, nHep_guess, nHepp_guess);
  
  return u;
}

/******************************************************************************************/

double solve_rates(double u, double rho, double dt, double *nHp_guess, double *nHep_guess, double *nHepp_guess)
{
  
  double t_init = 0.0;
  
  void *Cvode_mem;
  N_Vector y = NULL, abstol = NULL;
  realtype dt_out, reltol;
  
  int NEQ = 4;

 
  yhelium  = (1.0-XH)/(4.0*XH);
  nHcgs    = XH * rho / PROTONMASS;	
  
    
  /* Initial guess for implicit solver */
  nHp   = dmin(1.0, dmax(*nHp_guess, 0.0));
  nHep  = dmin(yhelium, dmax(*nHep_guess, 0.0));
  nHepp = dmin(yhelium, dmax(*nHepp_guess, 0.0));
  
  if(nHep + nHepp > yhelium)
    {
      nHep  *= yhelium/(nHep + nHepp);
      nHepp *= yhelium/(nHep + nHepp);
    }
  
  
  /* Vector of initial values */
  y = N_VNew_Serial(NEQ);
  
  Ith(y,1) = u;
  Ith(y,2) = nHp;
  Ith(y,3) = nHep;
  Ith(y,4) = nHepp;


  /* Set tolerances */
  reltol = RCONST(1.0e-5);
  abstol = N_VNew_Serial(NEQ);
  
  Ith(abstol,1) = 1.0e-4;
  Ith(abstol,2) = 1.0e-7 * dmax(nHp,1.0e-7); 
  Ith(abstol,3) = 1.0e-7 * dmax(nHep,1.0e-7);
  Ith(abstol,4) = 1.0e-7 * dmax(nHepp,1.0e-7);
  
  
  /* Instantiates a cvode solver object and specifies the solution method */
  Cvode_mem = CVodeCreate(CV_BDF, CV_NEWTON);

  /* Provide required problem and solution specications, allocate
     internal memory, and initialize cvode */
  CVodeInit(Cvode_mem,f_solve_rates, t_init, y);

  /* Specifies scalar relative tolerance and vector absolute tolerances */
  CVodeSVtolerances(Cvode_mem, reltol, abstol);

  //  long int mxsteps = 5000;
  //int flag_mxsteps;
  //flag_mxsteps = CVodeSetMaxNumSteps(Cvode_mem, mxsteps);
  
  /* Select linear solver */
  CVDense(Cvode_mem, NEQ);
  
  /* Integrate rate equations */
  CVode(Cvode_mem, dt, y, &dt_out, CV_NORMAL);
  
  u     = Ith(y,1);
  nHp   = dmin(1.0, dmax(Ith(y,2), 0.0));
  nHep  = dmin(yhelium, dmax(Ith(y,3), 0.0));
  nHepp = dmin(yhelium, dmax(Ith(y,4), 0.0));
  
  if(nHep + nHepp > yhelium)
    {
      nHep  *= yhelium/(nHep + nHepp);
      nHepp *= yhelium/(nHep + nHepp);
    }
  
  /* Deallocate solution vector and solver memory */
  N_VDestroy_Serial(y);
  N_VDestroy_Serial(abstol);
  CVodeFree(&Cvode_mem);
 
  
  *nHp_guess   = nHp;
  *nHep_guess  = nHep;
  *nHepp_guess = nHepp;
  
  return u;

}


/******************************************************************************************/

static int f_solve_rates(realtype dt, N_Vector y, N_Vector ydot, void *user_data)
{
  
  double u;
  double dnHp_dt,dnHepp_dt ,dnHep_dt, du_dt;
  double mu, T, ztime2;
  double Lambda, Heat;
  double LambdaExc, LambdaIon, LambdaRec, LambdaFF, LambdaCmptn;
  double LambdaExcH0, LambdaExcHe0, LambdaExcHep, LambdaIonH0, LambdaIonHe0, LambdaIonHep;
  double LambdaRecHp, LambdaRecHep, LambdaRecHepp, LambdaRecHepd;
  double t, fhi, flow, Tlow, Thi;
  int j;

  
  u     = (double)Ith(y,1);
  nHp   = (double)Ith(y,2);
  nHep  = (double)Ith(y,3);
  nHepp = (double)Ith(y,4);
  nHe0  = yhelium - nHep - nHepp;
  nH0   = 1.0 - nHp;
  ne    = nHp + nHep + 2.0*nHepp;
  
  
  /* Get temperature */
  mu = (1.0 + 4.0*yhelium)/(1.0 + yhelium + ne);
  T  = GAMMA_MINUS1 / BOLTZMANN * u * PROTONMASS * mu;

  
  /* Check table limits */
  if(T < TMIN) T = TMIN;
  if(T > TMAX) T = TMAX;
  
   
  /* Obtain interpolations from the cooling table*/
  t = (log10(T) - logTmin) / deltaT;
  j = (int) t;
  Tlow = logTmin + deltaT * j;
  Thi = Tlow + deltaT;
  fhi = t - j;
  flow = 1 - fhi;


 /* get photo ionization and heating rates */
  
#ifndef PLAW_UVB
  
  int i, ilow;
  double logz, dzlow, dzhi;
    
  logz = log10(1.0 + redshift);
  ilow = 0;
  
  for(i = 0; i < nheattab; i++)
    {
      if(inlogz[i] < logz)
  	ilow = i;
      else
  	break;
    }

  dzlow = logz - inlogz[ilow];
  dzhi = inlogz[ilow + 1] - logz;
  
  if(logz > inlogz[nheattab - 1] || gH0[ilow] == 0 || gH0[ilow + 1] == 0 || nheattab == 0)
    {
      gJHe0 = gJHep = gJH0 = 0;
      epsHe0 = epsHep = epsH0 = 0;
    }
  else
    {
      gJH0   = pow(10., (dzhi * log_gH0[ilow] + dzlow  * log_gH0[ilow + 1])  / (dzlow + dzhi));
      gJHe0  = pow(10., (dzhi * log_gHe[ilow] + dzlow  * log_gHe[ilow + 1])  / (dzlow + dzhi));
      gJHep  = pow(10., (dzhi * log_gHep[ilow] + dzlow * log_gHep[ilow + 1]) / (dzlow + dzhi));
      
      epsH0  = pow(10., (dzhi * log_eH0[ilow]  + dzlow * log_eH0[ilow + 1])  / (dzlow + dzhi));
      epsHe0 = pow(10., (dzhi * log_eHe[ilow]  + dzlow * log_eHe[ilow + 1])  / (dzlow + dzhi));
      epsHep = pow(10., (dzhi * log_eHep[ilow] + dzlow * log_eHep[ilow + 1]) / (dzlow + dzhi));
    }
  
#endif

#ifdef NO_HE2_HEAT
  epsHep = 0.0;
#endif
  

  
  /* Recombination rates [cm^3 s^-1] */
  aHp   = flow * AlphaHp[j]   + fhi * AlphaHp[j + 1];
  aHep  = flow * AlphaHep[j]  + fhi * AlphaHep[j + 1];
  aHepp = flow * AlphaHepp[j] + fhi * AlphaHepp[j + 1];
  ad    = flow * Alphad[j]    + fhi * Alphad[j + 1];
  
  /* Collisional ionisation rates  [cm^3 s^-1] */
  geH0  = flow * GammaeH0[j]  + fhi * GammaeH0[j + 1];
  geHe0 = flow * GammaeHe0[j] + fhi * GammaeHe0[j + 1];
  geHep = flow * GammaeHep[j] + fhi * GammaeHep[j + 1];
  

  /* Collisional excitation cooling [erg s^-1 cm^3]*/
  bH0  = flow * BetaH0[j]  + fhi * BetaH0[j + 1];
  bHe0 = flow * BetaHe0[j] + fhi * BetaHe0[j + 1];
  bHep = flow * BetaHep[j] + fhi * BetaHep[j + 1];
  
  /* Free-free cooling  [erg s^-1 cm^3] */
  bff = flow * Betaff[j] + fhi * Betaff[j + 1];


  
  /* Recombination cooling [erg s^-1 cm^3] */
  LambdaRecHp   = 1.035e-16 * T * aHp  * ne * nHp;
  LambdaRecHep  = 1.035e-16 * T * aHep  * ne * nHep;
  LambdaRecHepp = 1.035e-16 * T * aHepp * ne * nHepp;
  LambdaRecHepd = 6.539e-11 * ad * ne * nHep;
  LambdaRec     = LambdaRecHp + LambdaRecHep + LambdaRecHepp + LambdaRecHepd; 

  
  /* Collisional excitation cooling [erg s^-1 cm^3]*/
  LambdaExcH0  = bH0 * ne * nH0; 
  LambdaExcHe0 = bHe0 * ne * nHe0; 
  LambdaExcHep = bHep * ne * nHep;
  LambdaExc    = LambdaExcH0 +  LambdaExcHe0 + LambdaExcHep; 
      
  /* Collisional ionisation cooling [erg s^-1 cm^3]*/
  LambdaIonH0  = 2.179e-11 * geH0 * ne * nH0;
  LambdaIonHe0 = 3.939e-11 * geHe0 * ne * nHe0;
  LambdaIonHep = 8.720e-11 * geHep * ne * nHep;
  LambdaIon    = LambdaIonH0 + LambdaIonHe0 + LambdaIonHep;	
      
  /* Free-free cooling [erg s^-1 cm^3] */
  LambdaFF = bff * (nHp + nHep + 4.0 * nHepp) * ne; 
  
  /* Compton scattering off CMB photons [erg s^-1 cm^3], 
     Weymann 1966, ApJ, 145, 560*/
  ztime2 = (1.0+redshift)*(1.0+redshift);
  LambdaCmptn = 5.653e-36*ne*(T-2.73*(1.0+redshift))*ztime2*ztime2/nHcgs;
  
  /* Total cooling rate /n_H^2  [erg s^-1 cm^3] */
  Lambda = LambdaExc + LambdaIon + LambdaRec + LambdaFF + LambdaCmptn;
  
  /* Total heating rate / n_H^2  [erg s^-1 cm^3]*/
  Heat = (nH0 * epsH0 + nHe0 * epsHe0 + nHep * epsHep) / nHcgs;

  /* Species rate of change / n_H [s^-1]*/
  dnHp_dt   = (nH0*(gJH0 + geH0*ne*nHcgs) - aHp*ne*nHcgs*nHp);
  dnHepp_dt = (nHep*(gJHep + geHep*ne*nHcgs) - aHepp*ne*nHcgs*nHepp);
  dnHep_dt  = (nHe0*(gJHe0 + geHe0*ne*nHcgs) - (aHep + ad)*ne*nHcgs*nHep) - dnHepp_dt;
  
  /* Rate of change of internal energy u [erg s^-1 g^-1] */
  du_dt = (Heat - Lambda) * (nHcgs*XH/PROTONMASS) - 2.0*Ha*u;

  Ith(ydot,1) = du_dt;
  Ith(ydot,2) = dnHp_dt;
  Ith(ydot,3) = dnHep_dt;
  Ith(ydot,4) = dnHepp_dt;
  
  return 0;
}


/********************************************************************************/

/* Initialise cooling routine */

void InitCool(double temp, double logxe)
{
  char fname[400];
  
  sprintf(fname,"./uvb_models/%s",UVBFILE);
  
  InitCoolMemory();
  MakeCoolingTable();

#ifdef PLAW_UVB
  InitUVB();
#else  
  ReadIonizeParams(fname);
#endif

}

/******************************************************************************/

/* Allocate the arrays for the cooling function table */

void InitCoolMemory(void)
{
  BetaH0         = (double *) calloc((NCOOLTAB + 1), sizeof(double));
  if(NULL==BetaH0){free(BetaH0);printf("Memory allocation failed for BetaH0.\n");
    endrun(fnID);}

  BetaHe0        = (double *) calloc((NCOOLTAB + 1), sizeof(double));
  if(NULL==BetaHe0){free(BetaHe0);printf("Memory allocation failed for BetaHe0.\n");
    endrun(fnID);}
  
  BetaHep        = (double *) calloc((NCOOLTAB + 1), sizeof(double));
  if(NULL==BetaHep){free(BetaHep);printf("Memory allocation failed for BetaHep.\n");
    endrun(fnID);}
  
  AlphaHp       = (double *) calloc((NCOOLTAB + 1), sizeof(double));
  if(NULL==AlphaHp){free(AlphaHp);printf("Memory allocation failed for AlphaHp.\n");
    endrun(fnID);}
  
  AlphaHep      = (double *) calloc((NCOOLTAB + 1), sizeof(double));
  if(NULL==AlphaHep){free(AlphaHep);printf("Memory allocation failed for AlphaHep.\n");
    endrun(fnID);}
  
  AlphaHepp     = (double *) calloc((NCOOLTAB + 1), sizeof(double));
  if(NULL==AlphaHepp){free(AlphaHepp);printf("Memory allocation failed for AlphaHepp.\n");
    endrun(fnID);}
  
  Alphad     = (double *) calloc((NCOOLTAB + 1), sizeof(double));
  if(NULL==Alphad){free(Alphad);printf("Memory allocation failed for Alphad.\n");
    endrun(fnID);}
  
  GammaeH0       = (double *) calloc((NCOOLTAB + 1), sizeof(double));
  if(NULL==GammaeH0){free(GammaeH0);printf("Memory allocation failed for GammaeH0.\n");
    endrun(fnID);}
  
  GammaeHe0      = (double *) calloc((NCOOLTAB + 1), sizeof(double));
  if(NULL==GammaeHe0){free(GammaeHe0);printf("Memory allocation failed for GammaeHe0.\n");
    endrun(fnID);}
  
  GammaeHep      = (double *) calloc((NCOOLTAB + 1), sizeof(double));
  if(NULL==GammaeHep){free(GammaeHep);printf("Memory allocation failed for GammaeHep.\n");
    endrun(fnID);}
  
  Betaff         = (double *) calloc((NCOOLTAB + 1), sizeof(double));
  if(NULL==Betaff){free(Betaff);printf("Memory allocation failed for Betaff.\n");
    endrun(fnID);}
}
/******************************************************************************************/

/* Make the look-up table for the cooling function */

void MakeCoolingTable(void)  
{
  int i;
  double T, T_eV;
  
  /* Verner & Ferland case-A recombination rate fit coefficients */
  double a_vf96[3] = {7.982e-11, 9.356e-10, 1.891e-10};
  double b_vf96[3] = {0.7480   , 0.7892   , 0.7524};
  double T0_vf96[3]= {3.148e0  , 4.266e-2 , 9.370e0};
  double T1_vf96[3]= {7.036e5  , 4.677e6  , 2.7674e6};
  
  /* Voronov collisional ionisation rate fit coefficients */
  double dE_vor97[3] = {13.6    , 24.6    , 54.4};
  double P_vor97[3]  = {0       , 0       , 1};
  double A_vor97[3]  = {0.291e-7, 0.175e-7, 0.205e-8};
  double X_vor97[3]  = {0.232   , 0.180   , 0.265};
  double K_vor97[3]  = {0.39    , 0.35    , 0.25};
  
  
  logTmin = log10(TMIN);
  logTmax = log10(TMAX);
  deltaT  = (logTmax - logTmin) / NCOOLTAB;
 
  for(i = 0; i <= NCOOLTAB; i++)
    {
      
      T      = pow(10.0, log10(TMIN) + deltaT * i);
      T_eV   = T * BOLTZMANN/ELECTRONVOLT;
      
      /* Case-A recombination rates, [cm^3 s^-1]
	 Verner & Ferland, 1996, ApJS, 103, 467 */
      AlphaHp[i]  = a_vf96[0]/(sqrt(T/T0_vf96[0])*pow(1.0+sqrt(T/T0_vf96[0]),1.0-b_vf96[0])
			       *pow(1.0+sqrt(T/T1_vf96[0]),1.0+b_vf96[0]));
      
      AlphaHep[i] = a_vf96[1]/(sqrt(T/T0_vf96[1])*pow(1.0+sqrt(T/T0_vf96[1]),1.0-b_vf96[1])
			       *pow(1.0+sqrt(T/T1_vf96[1]),1.0+b_vf96[1]));
	
      AlphaHepp[i]= a_vf96[2]/(sqrt(T/T0_vf96[2])*pow(1.0+sqrt(T/T0_vf96[2]),1.0-b_vf96[2])
			       *pow(1.0+sqrt(T/T1_vf96[2]),1.0+b_vf96[2]));
      
      
      /* He+ dielectronic recombination rate [cm^3 s^-1]
	 Aldrovandi & Pequignot, 1973, A&A, 25, 137 */
      if(4.7e5/T < 70) 
        Alphad[i] = 1.9e-3*(1.0 + 0.3*exp(-9.4e4/T))*exp(-4.7e5/T)*pow(T,-1.5);

      
      /* Collisional ionisation rates [cm^3 s^-1]
	 Voronov 1997, ADNDT, 65, 1 */
      if(dE_vor97[0]/T_eV < 70)
	GammaeH0[i] = A_vor97[0]*pow(dE_vor97[0]/T_eV, K_vor97[0])*exp(-dE_vor97[0]/T_eV) 
	  *(1.0+P_vor97[0]*sqrt(dE_vor97[0]/T_eV)) / (X_vor97[0] + dE_vor97[0]/T_eV);
      
      if(dE_vor97[1]/T_eV < 70)
        GammaeHe0[i] = A_vor97[1]*pow(dE_vor97[1]/T_eV, K_vor97[1])*exp(-dE_vor97[1]/T_eV) 
	  *(1.0+P_vor97[1]*sqrt(dE_vor97[1]/T_eV)) / (X_vor97[1] + dE_vor97[1]/T_eV);
      
      if(dE_vor97[2]/T_eV < 70)
        GammaeHep[i] = A_vor97[2]*pow(dE_vor97[2]/T_eV, K_vor97[2])*exp(-dE_vor97[2]/T_eV) 
	  *(1.0+P_vor97[2]*sqrt(dE_vor97[2]/T_eV)) / (X_vor97[2] + dE_vor97[2]/T_eV);
      
      
      /* Collisional excitation cooling [erg s^-1 cm^3]
	 Cen 1992, ApJS, 78, 341 */
      if(118348/T < 70)
        BetaH0[i] = 7.5e-19*exp(-118348/T)/(1.0+sqrt(T/1.0e5)); 

      if(13179/T < 70)
        BetaHe0[i] = 9.1e-27*pow(T,-0.1687)*exp(-13179/T)/(1.0+sqrt(T/1.0e5));
      
      if(473638/T < 70) 
        BetaHep[i] = 5.54e-17*pow(T,-0.397)*exp(-473638/T)/(1.0+sqrt(T/1.0e5)); 
      

      /* Bremsstrahlung (free-free) [erg cm^3 s^-1]
	 Theuns et al. 1998, MNRAS, 301, 478 */
      Betaff[i] = 1.42e-27*sqrt(T)*(1.1+0.34*exp(-(5.5-log10(T))*(5.5-log10(T))/3.0));
    }
  
}

/******************************************************************************************/

/* Reads in the information from TREECOOL */


void ReadIonizeParams(char *fname)
{
  int i;
  FILE *fdcool;

  printf("\nUsing %s:\n",UVBFILE);

  
  if(!(fdcool = fopen(fname, "r")))
    {
      printf(" Cannot read ionization table in file `%s'\n", fname);
      endrun(fnID);
    }
  
  for(i = 0; i < TABLESIZE; i++)
    if(fscanf(fdcool, "%g %g %g %g %g %g %g", &inlogz[i], &gH0[i], &gHe[i], &gHep[i], &eH0[i], &eHe[i], &eHep[i]) == EOF)
      break;
  fclose(fdcool);
  
  /*  nheattab is the number of entries in the table */
  for(i = 0, nheattab = 0; i < TABLESIZE; i++)
    if(gH0[i] != 0.0)
      {
	nheattab++;
	log_gH0[i]  = log10(gH0[i]);
	log_gHe[i]  = log10(gHe[i]);
	log_gHep[i] = log10(gHep[i]);
	log_eH0[i]  = log10(eH0[i]);
	log_eHe[i]  = log10(eHe[i]);
	log_eHep[i] = log10(eHep[i]);
      }
    else
      break;
  
}

/******************************************************************************************/

