
/**************************************************************************/

/**! \file uvb.c
 *
 * \brief If using power-law UVB spectrum, use Gauss-Legendre
 * quadrature to solve the integrals for the H0, He0 and He+
 * photo-ionisation and photo-heating rates.  Optionally includes
 * secondary ionisations by fast photo-electrons.
 *
 */

/***************************************************************************/

#ifdef PLAW_UVB

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "parameters.h"
#include "proto.h"
#include "global_vars.h"


static double *E_tab, *JE_tab;
static double *sH1, *sHe1, *sHe2;
static double dlogE_inv, logE_min, logE_max; 

static double **fheat_tab, **nionH1_tab, **nionHe1_tab, **nionHe2_tab;
static double delta_logxe_inv,delta_logEsec_inv;


static double logEion_H1, logEion_He1, logEion_He2;

#define NETAB 100000


/** \brief Allocate the memory and calculate the spectrum and
 * cross-sections look-up tables at the start of the run
 *
 */

void InitUVB(void)
{
  InitUVBMemory();
  MakeUVBTable();
}



/** \brief Evaluate the integrals for the photo-ionisation and
 * photo-heating rates using Gauss-Legendre quadrature.
 *
 * \param logxe log10(nelec / nH)
 */

void IonizationRates(double logxe)
{

  double sigma_H1, sigma_He1, sigma_He2;
  double JE, E, logE;
  double lim_plus,lim_minus;
  double gfac,efac;
  double eval_H1,eval_He1,eval_He2;
  double gJH0_H0,gJH0_He0,gJH0_Hep;
  double gJHe0_H0,gJHe0_He0,gJHe0_Hep;
  double gJHep_H0,gJHep_He0,gJHep_Hep;
  double flow,fhi,t;
  
  int i, j;
 
#ifdef SECONDARY
  double evalH1_H1, evalH1_He1, evalH1_He2;
  double evalHe1_H1, evalHe1_He1, evalHe1_He2;
  double evalHe2_H1, evalHe2_He1, evalHe2_He2;
  double flow1,fhi1,t1,flow2,fhi2,t2;
  double fheat,nionH1,nionHe1,nionHe2;
  double logEsec;
  
  int j1,j2;
#endif

   
  gJH0 = gJHe0 = gJHep = epsH0 = epsHe0 = epsHep = 0.0;
 
  gJH0_H0  = gJH0_He0  = gJH0_Hep  = 0.0;
  gJHe0_H0 = gJHe0_He0 = gJHe0_Hep = 0.0;
  gJHep_H0 = gJHep_He0 = gJHep_Hep = 0.0;
  


  if(redshift <= ZRH1)
    {
      
      /* HI rates.  Integrate from EionHI to Emax in steps of log10(E) */
      lim_plus  = 0.5*(logE_max + logEion_H1);
      lim_minus = 0.5*(logE_max - logEion_H1);
      
      for(i=0; i<NWEIGHTS; i++)
	{
	  
	  logE = lim_plus + lim_minus*absc[i]; 
	  E = pow(10.0,logE);
	  
	  t    = (logE - logE_min) * dlogE_inv;
	  j    = (int) t;
	  fhi  = t - j;
	  flow = 1 - fhi;
	  
	  /* Set spectrum to zero at E>54.4 eV if z>ZRHE2 */
	  JE = logE >= logEion_He2 && redshift > ZRHE2 ? 0.0 : flow * JE_tab[j] + fhi * JE_tab[j + 1];
	  
	  sigma_H1 = flow * sH1[j]  + fhi * sH1[j + 1];
	  
	  eval_H1 = lim_minus * weight[i]  * JE * sigma_H1;
	  gJH0   += eval_H1;
	  	  
#ifdef SECONDARY
	  
	  logEsec = log10(E - EION_H1);
	  
	  if(logEsec > ESECMIN && logxe < XMAX)
	    {
	      logxe = dmax(XMIN,logxe);
	      t1    = (logxe - XMIN) * delta_logxe_inv;
	      j1    = (int)t1;
	      fhi1  = t1 - j1;
	      flow1 = 1.0 - fhi1;
	      
	      logEsec = dmin(ESECMAX,logEsec);
	      t2      = (logEsec - ESECMIN) * delta_logEsec_inv;
	      j2      = (int)t2;
	      fhi2    = t2 - j2;
	      flow2   = 1.0 - fhi2;
	      
	      fheat = flow2*(flow1*fheat_tab[j1][j2] + fhi1*fheat_tab[j1+1][j2]) 
		+ fhi2*(flow1*fheat_tab[j1][j2+1] + fhi1*fheat_tab[j1+1][j2+1]);
	      
	      nionH1 = flow2*(flow1*nionH1_tab[j1][j2] + fhi1*nionH1_tab[j1+1][j2]) 
		+ fhi2*(flow1*nionH1_tab[j1][j2+1] + fhi1*nionH1_tab[j1+1][j2+1]);
	      
	      nionHe1 = flow2*(flow1*nionHe1_tab[j1][j2] + fhi1*nionHe1_tab[j1+1][j2]) 
		+ fhi2*(flow1*nionHe1_tab[j1][j2+1] + fhi1*nionHe1_tab[j1+1][j2+1]);
	      
	      nionHe2 = flow2*(flow1*nionHe2_tab[j1][j2]  + fhi1*nionHe2_tab[j1+1][j2]) 
		+ fhi2*(flow1*nionHe2_tab[j1][j2+1]  + fhi1*nionHe2_tab[j1+1][j2+1]);
	      
	      /* HI ionisations from HI photo-electrons */
	      evalH1_H1 = eval_H1 * nionH1;
	      gJH0_H0  += evalH1_H1;
	      
	      /* HeI ionisations from HI photo-electrons */
	      evalHe1_H1 = eval_H1 * nionHe1;
	      gJHe0_H0  += evalHe1_H1;
	      
	      /* HeII ionisations from HI photo-electrons */
	      evalHe2_H1 = eval_H1 * nionHe2;
	      gJHep_H0  += evalHe2_H1;
	      
	    }
	  else
	    fheat  = 1.0;
	  
	  epsH0 += eval_H1 * (E - EION_H1) * fheat;
#else
	  epsH0 += eval_H1 * (E - EION_H1); 
#endif
	  
	}
 
      
      
      /* HeI rates.  Integrate from EionHeI to Emax */
      lim_plus  = 0.5*(logE_max + logEion_He1);
      lim_minus = 0.5*(logE_max - logEion_He1);
      
      for(i=0; i<NWEIGHTS; i++)
	{
	  logE = lim_plus + lim_minus*absc[i]; 
	  E = pow(10.0,logE);
	  
	  t    = (logE - logE_min) * dlogE_inv;
	  j    = (int) t;
	  fhi  = t - j;
	  flow = 1 - fhi;

	  JE = logE >= logEion_He2 && redshift > ZRHE2 ? 0.0 : flow * JE_tab[j] + fhi * JE_tab[j + 1];
	
	  sigma_He1 = flow * sHe1[j] + fhi * sHe1[j + 1];
	  
	  eval_He1 = lim_minus * weight[i]  * JE * sigma_He1;
	  gJHe0   += eval_He1;
	  
#ifdef SECONDARY
	  
	  logEsec = log10(E - EION_HE1);
	  
	  if(logEsec > ESECMIN && logxe < XMAX)
	    {
	      logxe = dmax(XMIN,logxe);
	      t1    = (logxe - XMIN) * delta_logxe_inv;
	      j1    = (int)t1;
	      fhi1  = t1 - j1;
	      flow1 = 1.0 - fhi1;
	      	      
	      logEsec  = dmin(ESECMAX,logEsec);
	      t2    = (logEsec - ESECMIN) * delta_logEsec_inv;
	      j2    = (int)t2;
	      fhi2  = t2 - j2;
	      flow2 = 1.0 - fhi2;
	      
	      fheat = flow2*(flow1*fheat_tab[j1][j2] + fhi1*fheat_tab[j1+1][j2]) 
		+ fhi2*(flow1*fheat_tab[j1][j2+1] + fhi1*fheat_tab[j1+1][j2+1]);
	      
	      nionH1 = flow2*(flow1*nionH1_tab[j1][j2] + fhi1*nionH1_tab[j1+1][j2]) 
		+ fhi2*(flow1*nionH1_tab[j1][j2+1] + fhi1*nionH1_tab[j1+1][j2+1]);
	      
	      nionHe1 = flow2*(flow1*nionHe1_tab[j1][j2] + fhi1*nionHe1_tab[j1+1][j2]) 
		+ fhi2*(flow1*nionHe1_tab[j1][j2+1] + fhi1*nionHe1_tab[j1+1][j2+1]);
	      
	      nionHe2 = flow2*(flow1*nionHe2_tab[j1][j2]  + fhi1*nionHe2_tab[j1+1][j2]) 
		+ fhi2*(flow1*nionHe2_tab[j1][j2+1]  + fhi1*nionHe2_tab[j1+1][j2+1]);
	      
	      /* HI ionisations from HeI photo-electrons */
	      evalH1_He1 = eval_He1 * nionH1;
	      gJH0_He0  += evalH1_He1;
	  
	      /* HeI ionisations from HeI photo-electrons */
	      evalHe1_He1 = eval_He1 * nionHe1;
	      gJHe0_He0  += evalHe1_He1;
	      
	      /* HeII ionisations from HeI photo-electrons */
	      evalHe2_He1 = eval_He1 * nionHe2;
	      gJHep_He0  += evalHe2_He1;
	      
	    }
	  else
	    fheat  = 1.0;
	  
	  
	  epsHe0 += eval_He1 * (E - EION_HE1) * fheat;
#else
	  epsHe0 += eval_He1 * (E - EION_HE1);
#endif
	}

    }

  
  
  if(redshift <= ZRHE2)
    {
      
      /* HeII rates. Integrate from EionHeII to Emax */
      lim_plus  = 0.5*(logE_max + logEion_He2);
      lim_minus = 0.5*(logE_max - logEion_He2);
      
      for(i=0;i<NWEIGHTS;i++)
	{
	  logE = lim_plus + lim_minus*absc[i]; 
	  E = pow(10.0, logE);
      
	  t    = (logE - logE_min) * dlogE_inv;
	  j    = (int) t;
	  fhi  = t - j;
	  flow = 1 - fhi;
	  
	  JE        = flow * JE_tab[j] + fhi * JE_tab[j + 1];

	  sigma_He2 = flow * sHe2[j] + fhi * sHe2[j + 1];
	  
	  eval_He2 = lim_minus * weight[i]  * JE * sigma_He2;
	  gJHep   += eval_He2;
	  

#ifdef SECONDARY
	  
	  logEsec = log10(E - EION_HE2);
	  
	  if(logEsec > ESECMIN && logxe < XMAX)
	    {
	      logxe = dmax(XMIN,logxe);
	      t1    = (logxe - XMIN) * delta_logxe_inv;
	      j1    = (int)t1;
	      fhi1  = t1 - j1;
	      flow1 = 1.0 - fhi1;
	      
	      logEsec  = dmin(ESECMAX,logEsec);
	      t2    = (logEsec - ESECMIN) * delta_logEsec_inv;
	      j2    = (int)t2;
	      fhi2  = t2 - j2;
	      flow2 = 1.0 - fhi2;
	      
	      fheat = flow2*(flow1*fheat_tab[j1][j2] + fhi1*fheat_tab[j1+1][j2]) 
		+ fhi2*(flow1*fheat_tab[j1][j2+1] + fhi1*fheat_tab[j1+1][j2+1]);
	      
	      nionH1 = flow2*(flow1*nionH1_tab[j1][j2] + fhi1*nionH1_tab[j1+1][j2]) 
		+ fhi2*(flow1*nionH1_tab[j1][j2+1] + fhi1*nionH1_tab[j1+1][j2+1]);
	      
	      nionHe1 = flow2*(flow1*nionHe1_tab[j1][j2] + fhi1*nionHe1_tab[j1+1][j2]) 
		+ fhi2*(flow1*nionHe1_tab[j1][j2+1] + fhi1*nionHe1_tab[j1+1][j2+1]);
	      
	      nionHe2 = flow2*(flow1*nionHe2_tab[j1][j2]  + fhi1*nionHe2_tab[j1+1][j2]) 
		+ fhi2*(flow1*nionHe2_tab[j1][j2+1]  + fhi1*nionHe2_tab[j1+1][j2+1]);
   
	      /* HI ionisations from HeII photo-electrons */
	      evalH1_He2 = eval_He2 * nionH1;
	      gJH0_Hep  += evalH1_He2;
	      
	      /* HeI ionisations from HeII photo-electrons */
	      evalHe1_He2 = eval_He2 * nionHe1;
	      gJHe0_Hep  += evalHe1_He2;
	      
	      /* HeII ionisations from HeII photo-electrons */
	      evalHe2_He2 = eval_He2 * nionHe2;
	      gJHep_Hep  += evalHe2_He2;
	      
	    }
	  else
	    fheat  = 1.0;
	  
	  
	  epsHep += eval_He2 * (E - EION_HE2) * fheat;
#else
	  epsHep += eval_He2 * (E - EION_HE2);
#endif
	}
    }

  gfac  = 4.0 * M_PI * log(10.0) / PLANCK;
  efac  = gfac * ELECTRONVOLT;
  
  gJH0  = gfac * (gJH0  + gJH0_H0  + gJH0_He0  + gJH0_Hep);
  gJHe0 = gfac * (gJHe0 + gJHe0_H0 + gJHe0_He0 + gJHe0_Hep);
  gJHep = gfac * (gJHep + gJHep_H0 + gJHep_He0 + gJHep_Hep);

  epsH0  *= efac;
  epsHe0 *= efac;
  epsHep *= efac;
  
}



/** \brief Calculate the look-up tables for the photo-ionisation
 * cross-sections and the power-law UVB spectrum.  The cross-sections
 * for H0, He0 and He+ are from Verner et al. 1996, ApJ, 465, 487
 *
 */

void MakeUVBTable()
{
  
  double E, dlogE;
  
  double E0_v96[3] = {4.298e-1, 1.361e1,  1.720e0};
  double s0_v96[3] = {5.475e4 , 9.492e2,  1.369e4};
  double ya_v96[3] = {3.288e1 , 1.469e0,  3.288e1};
  double P_v96[3]  = {2.963e0 , 3.188e0,  2.963e0};
  double yw_v96[3] = {0.0     , 2.039e0,  0.0    };
  double y0_v96[3] = {0.0     , 4.434e-1, 0.0    };
  double y1_v96[3] = {0.0     , 2.136e0,  0.0    };

  double x_v96,y_v96,Mb_to_cm2=1.0e-18;
  

  int i;
  
  logEion_H1  = log10(EION_H1);
  logEion_He1 = log10(EION_HE1);
  logEion_He2 = log10(EION_HE2);
  
  logE_min  = log10(EION_H1);
  logE_max  = log10(EMAX); 
  
  dlogE     = (logE_max-logE_min) / NETAB;
  dlogE_inv = 1.0 / dlogE;
  
  
  for(i = 0; i <= NETAB; i++)
    {
      E_tab[i] = logE_min  + dlogE * i;
      E        = pow(10.0, E_tab[i]);
      
      
      /* Specific intensity, erg s^-1 cm^-2 sr^-1 Hz^-1 */
      JE_tab[i] = J22*1.0e-22*pow(E/EION_H1, -ALPHA_UV);
      
      /* Photo-ionisation cross-sections, cm^-2. */
      if(E >= EION_H1)
	{
	  x_v96 = E/E0_v96[0] - y0_v96[0];
	  y_v96 = sqrt(x_v96*x_v96 + y1_v96[0]*y1_v96[0]);
	  
	  sH1[i] = s0_v96[0] * Mb_to_cm2 * ((x_v96 - 1.0) * (x_v96 - 1.0) + yw_v96[0] * yw_v96[0])
	    * pow(y_v96, 0.5 * P_v96[0]  - 5.5) * pow(1.0 + sqrt(y_v96/ya_v96[0]), -P_v96[0]);
	}
      
      if(E >= EION_HE1)
	{
	  x_v96 = E/E0_v96[1] - y0_v96[1];
	  y_v96 = sqrt(x_v96*x_v96 + y1_v96[1]*y1_v96[1]);

	  sHe1[i] = s0_v96[1] * Mb_to_cm2 * ((x_v96 - 1.0) * (x_v96 - 1.0) + yw_v96[1]  * yw_v96[1])
	    * pow(y_v96, 0.5 * P_v96[1] - 5.5) * pow(1.0 + sqrt(y_v96/ya_v96[1]), -P_v96[1]);
	}
      
      if(E >= EION_HE2)
	{
	  x_v96 = E/E0_v96[2] - y0_v96[2];
	  y_v96 = sqrt(x_v96*x_v96 + y1_v96[2]*y1_v96[2]);
	  
	  sHe2[i] = s0_v96[2] * Mb_to_cm2 * ((x_v96 - 1.0) * (x_v96 - 1.0) + yw_v96[2] * yw_v96[2])
	    * pow(y_v96, 0.5 * P_v96[2] - 5.5) * pow(1.0 + sqrt(y_v96/ya_v96[2]), -P_v96[2]);
	}

    
    }
  
}


/** \brief Allocate the memory for the UVB and cross-section look-up tables.
 *
 */

void InitUVBMemory(void)
{
  
  E_tab = (double *) calloc((NETAB + 1), sizeof(double));
  if(NULL==E_tab){free(E_tab); printf("Memory allocation failed for E_tab.\n"); exit(0);}
  
  JE_tab = (double *) calloc((NETAB + 1), sizeof(double));
  if(NULL==JE_tab){free(JE_tab); printf("Memory allocation failed for JE_tab.\n"); exit(0);}

  sH1  = (double *) calloc((NETAB + 1), sizeof(double));
  if(NULL==sH1){free(sH1); printf("Memory allocation failed for sH1.\n"); exit(0);}
  
  sHe1 = (double *) calloc((NETAB + 1), sizeof(double));
  if(NULL==sHe1){free(sHe1); printf("Memory allocation failed for sHe1.\n"); exit(0);}
  
  sHe2 = (double *) calloc((NETAB + 1), sizeof(double));
  if(NULL==sHe2){free(sHe2); printf("Memory allocation failed for sHe2.\n"); exit(0);}
}


/** \brief Allocate the memory for the secondary ionisation look-up
 * table.  This is a 2D table that depends on the photo-electron
 * energy and the free electron fraction.
 *
 */

void InitSecondaryMemory(void)
{
  
  fheat_tab   = (double **) calloc((NXTAB+1) , sizeof(double *));
  if(NULL==fheat_tab){free(fheat_tab); printf("Memory allocation failed for fheat_tab.\n"); exit(0);}
  
  nionH1_tab  = (double **) calloc((NXTAB+1) , sizeof(double *));
  if(NULL==nionH1_tab){free(nionH1_tab); printf("Memory allocation failed for nionH1_tab.\n"); exit(0);}
  
  nionHe1_tab = (double **) calloc((NXTAB+1) ,sizeof(double *));
  if(NULL==nionHe1_tab){free(nionHe1_tab); printf("Memory allocation failed for nionHe1_tab.\n"); exit(0);}
  
  nionHe2_tab = (double **) calloc((NXTAB+1),sizeof(double *));
  if(NULL==nionHe2_tab){free(nionHe2_tab); printf("Memory allocation failed for nionHe2_tab.\n"); exit(0);}

  int i; 
  
  for(i = 0; i <= NXTAB; i++)
    {
      fheat_tab[i]   = (double *) calloc((NESECTAB+1) , sizeof(double));
      if(NULL == fheat_tab[i]){free(fheat_tab[i]); printf("Memory allocation failed for fheat_tab[%d][].\n",i); exit(0);}
      
      nionH1_tab[i]  = (double *) calloc((NESECTAB+1) , sizeof(double));
      if(NULL == nionH1_tab[i]){free(nionH1_tab[i]); printf("Memory allocation failed for nionH1_tab[%d][].\n",i); exit(0);}
      
      nionHe1_tab[i] = (double *) calloc((NESECTAB+1) , sizeof(double));
      if(NULL == nionHe1_tab[i]){free(nionHe1_tab[i]); printf("Memory allocation failed for nionHe1_tab[%d][].\n",i); exit(0);}
      
      nionHe2_tab[i] = (double *) calloc((NESECTAB+1) , sizeof(double));
      if(NULL == nionHe2_tab[i]){free(nionHe2_tab[i]); printf("Memory allocation failed for nionHe2_tab[%d][].\n",i); exit(0);}
    }
}


/** \brief Read in the pre-computed secondary ionisation tables. These
 *   are from Furlanetto & Johnson-Stoever 2010, MNRAS, 404, 1869
 */

void ReadSecondaryTable(void)
{
  double logxe,delta_logxe,delta_logEsec;
  int i;
  
  char filename[200];
  FILE *input;
  
  delta_logxe = (XMAX - XMIN)/ (double)(NXTAB);
  delta_logxe_inv = 1.0/delta_logxe;
  
  delta_logEsec     = (ESECMAX - ESECMIN)/ (double)(NESECTAB);
  delta_logEsec_inv = 1.0/delta_logEsec;
  
  for(i = 0; i <= NXTAB; i++)
    {
      logxe = XMIN + delta_logxe * i;
      
      sprintf(filename,"./secondaries/logxe_%1.3f.dat",logxe);
      if(!(input = fopen(filename, "rb")))
	{
	  printf("Cannot read file:\n");
	  printf("%s\n", filename);
	  exit(0);
	}
      fread(fheat_tab[i],sizeof(double),NESECTAB + 1,input);
      fread(nionH1_tab[i],sizeof(double),NESECTAB + 1,input);
      fread(nionHe1_tab[i],sizeof(double),NESECTAB + 1,input);
      fread(nionHe2_tab[i],sizeof(double),NESECTAB + 1,input);
      fclose(input);
    }
}


#endif
