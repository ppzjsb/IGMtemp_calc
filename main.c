/********************************************************************************/

/**! \file main.c 
 *
 * \brief Main routine that drives the calculation.
 *
 */

/********************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "parameters.h"
#include "proto.h"
#include "global_vars.h"

/*struct output_data
{
  double z;
  double t;
  double temp;
  double  nH0, nHep, nHepp, ne, heat;
  double gJH0, gJHe0, gJHep, epsH0, epsHe0, epsHep;
} *IO;
*/

double *z_out, *t_out;
double  *nH0_out, *nHep_out, *nHepp_out, *ne_out, *temp_out, *heat_out;
double *gJH0_out, *gJHe0_out, *gJHep_out, *epsH0_out, *epsHe0_out, *epsHep_out;



/** \brief Initialises the calculation, drives the time-stepping, writes the output.
 *
 * \return (0)                                       
 */

int main()
{
  
  double start = clock();
    
  redshift = ZSTART;
  
  double H       = 1.0e7 / MPC;
  double yhelium = (1.0 - XH) / (4.0 * XH);
  double rhocb   = 3.0 * pow(H*HUBBLE, 2.0) * OMEGAB * pow(10.0,LOGOVERDEN) / (8.0 * M_PI * GRAVITY);
  
#ifdef DARK_PHOTON
  double rho_old=0.0;
#endif
  
  double atime     = 1.0 / (1.0 + ZSTART);
  double atime_end = 1.0 / (1.0 + ZEND);
  double dlna      = (log(atime_end) - log(atime)) / ((double)NT - 1.0);
  
  double nHp   = 0.0; 
  double nHep  = 0.0;
  double nHepp = 0.0;
  double ne    = nHp + nHep + 2.0 * nHepp;;
 
  
  /* Thermal decoupling redshift for CMB, appropriate for z<zdec: see
     Eq (66) in Furlanetto et al. 2006, PhR, 433, 181.  Use RECFAST if
     higher accuracy desired. */
  double zdec  = 150.0 * pow((OMEGAB * HUBBLE * HUBBLE) / 0.023, 0.4);
  double temp  = 2.726 / (atime * atime * (1.0 + zdec));
  

#ifdef PLAW_UVB
  double logxe = ne > 0.0 ? log10(ne/(1.0+2.0*yhelium)) : -10.0;

  gausslegendre();
#endif
  
  /* Load the external Furlanetto & Johnson-Stoever (2010, MNRAS, 404,
     1869) secondaries tables */
#ifdef SECONDARY
  MakeSecondaryTable();
  InitSecondaryMemory();
  ReadSecondaryTable();
#endif
 
  
  InitCool();

  InitDumpMemory();
  
  int i;
  double pcdone = 0.0;
  double t = 0.0;

  for(i=0; i<NT; i++)
    {
     
      Ha = H*HUBBLE*sqrt(OMEGAM/(atime*atime*atime) + OMEGAL);

      double dtdlna   = 1.0 / Ha;
      double dt       = dtdlna * dlna;

      redshift = 1.0/atime - 1.0;
    
      
#ifdef PLAW_UVB
      IonizationRates(logxe);
#endif  
      
      double rho   = rhocb / (atime*atime*atime);
      double nHcgs = XH * rho / PROTONMASS;
      
      z_out[i]      = redshift;
      t_out[i]      = t;
      temp_out[i]   = temp;
      nH0_out[i]    = 1.0-nHp;
      nHep_out[i]   = nHep;
      nHepp_out[i]  = nHepp;
      ne_out[i]     = ne;
      gJH0_out[i]   = gJH0;
      gJHe0_out[i]  = gJHe0;
      gJHep_out[i]  = gJHep;
      epsH0_out[i]  = epsH0;
      epsHe0_out[i] = epsHe0;
      epsHep_out[i] = epsHep;

     
      /*       IO[i].z      = redshift;
       IO[i].t      = t;
       IO[i].temp   = temp;
       IO[i].nH0    = 1.0-nHp;
       IO[i].nHep   = nHep;
       IO[i].nHepp = nHepp;
       IO[i].ne     = ne;
       IO[i].gJH0   = gJH0;
       IO[i].gJHe0  = gJHe0;
       IO[i].gJHep  = gJHep;
       IO[i].epsH0  = epsH0;
       IO[i].epsHe0 = epsHe0;
       IO[i].epsHep = epsHep;
      */
    
#ifdef DARK_PHOTON
      heat_out[i] = (Heat_dark + nHcgs*((1.0-nHp)*epsH0 + (yhelium-nHep-nHepp)*epsHe0 + nHep*epsHep)) / rho; /* erg s^-1 g^-1 */
      // IO[i].heat = (Heat_dark + nHcgs*((1.0-nHp)*epsH0 + (yhelium-nHep-nHepp)*epsHe0 + nHep*epsHep)) / rho; /* erg s^-1 g^-1 */
#else
      heat_out[i] = nHcgs*((1.0-nHp)*epsH0 + (yhelium-nHep-nHepp)*epsHe0 + nHep*epsHep) / rho;
      //IO[i].heat = nHcgs*((1.0-nHp)*epsH0 + (yhelium-nHep-nHepp)*epsHe0 + nHep*epsHep) / rho;
#endif
      
      double mu = (1.0 + 4.0 * yhelium) / (1.0 + yhelium + ne);
      double u  = (BOLTZMANN * temp)/(GAMMA_MINUS1 * mu * PROTONMASS); /*  erg g^-1 */

      /* Solves for the specific internal energy [erg g^-1] and abundances */
#ifdef DARK_PHOTON
      double u_new = DoCooling(u, rho, rho_old, dt, &nHp, &nHep, &nHepp);
      rho_old = rho;
#else
      double u_new = DoCooling(u, rho, dt, &nHp, &nHep, &nHepp);
#endif
      
      ne = nHp + nHep + 2.0*nHepp;

#ifdef PLAW_UVB
      logxe = ne > 0.0 ? log10(ne/(1.0+2.0*yhelium)) : XMIN;
#endif

      mu    = (1.0 + 4.0 * yhelium) / (1.0 + yhelium + ne);
      temp  = GAMMA_MINUS1 / BOLTZMANN * u_new * PROTONMASS * mu;
      
      if(atime/atime_end > pcdone || i == NT-1)
	{
	  printf("z=%1.3f, T=%1.1fK, fH0=%1.3e, fHep=%1.3e\n",
		 redshift,temp_out[i],nH0_out[i],nHep_out[i]/yhelium);
	  pcdone += 0.05;
	}
      
      
      atime = exp(log(atime) + dlna);
      t += dt; /* Note that t=0 corresponds to ZSTART, not z --> infty */
    }
  
  write_data();
  
  double end = clock();
  printf("\nTotal CPU time used %lf s\n",(end - start) / CLOCKS_PER_SEC);
  
  return(0);
}



/** \brief Save the data to a binary file                                   
 */

void write_data(void)
{
 
  FILE *output;  
  char fname[400];
  
  sprintf(fname,"./output/%s",OUTFILE);
  if(!(output = fopen(fname, "wb")))
    {
      printf("Cannot read file:\n");
      printf("%s\n", fname);
      exit(0);
    }
  
  int nbins = NT;
  
  double density_out = LOGOVERDEN;
  double om_out      = OMEGAM;
  double ol_out      = OMEGAL;
  double ob_out      = OMEGAB;
  double h_out       = HUBBLE;
  double xh_out      = XH;


  /* Header */
  fwrite(&om_out,sizeof(double),1,output);
  fwrite(&ol_out,sizeof(double),1,output);
  fwrite(&ob_out,sizeof(double),1,output);
  fwrite(&h_out,sizeof(double),1,output);
  fwrite(&xh_out,sizeof(double),1,output);
  fwrite(&density_out,sizeof(double),1,output);
  fwrite(&nbins,sizeof(int),1,output);
  
  /* Data */
  fwrite(z_out,sizeof(double),NT,output);
  fwrite(t_out,sizeof(double),NT,output);
  fwrite(temp_out,sizeof(double),NT,output);
  fwrite(nH0_out,sizeof(double),NT,output);
  fwrite(nHep_out,sizeof(double),NT,output);
  fwrite(nHepp_out,sizeof(double),NT,output);
  fwrite(ne_out,sizeof(double),NT,output);
  fwrite(heat_out,sizeof(double),NT,output);
  fwrite(gJH0_out,sizeof(double),NT,output);
  fwrite(gJHe0_out,sizeof(double),NT,output);
  fwrite(gJHep_out,sizeof(double),NT,output);
  fwrite(epsH0_out,sizeof(double),NT,output);
  fwrite(epsHe0_out,sizeof(double),NT,output);
  fwrite(epsHep_out,sizeof(double),NT,output);

  // fwrite(IO, sizeof(struct output_data), NT, output);

  fclose(output);

 
  
}



/** \brief Allocate memory for output file.  Would be better written as a structure.                                   
 */

void InitDumpMemory(void)
{

  /*  if(!(IO=calloc(NT,sizeof(struct output_data))))
    {
      printf("Memory allocation failed in read_snapshot.c\n\n");
      exit(0);
    }
  */

  z_out = (double *)calloc(NT, sizeof(double));
  if(NULL==z_out){free(z_out); printf("Memory allocation failed for z_out.\n"); exit(0);}
  
  t_out = (double *)calloc(NT, sizeof(double));
  if(NULL==t_out){free(t_out); printf("Memory allocation failed for t_out.\n"); exit(0);}
  
  temp_out = (double *)calloc(NT, sizeof(double));
  if(NULL==temp_out){free(temp_out); printf("Memory allocation failed for temp_out.\n"); exit(0);}
  
  nH0_out = (double *)calloc(NT, sizeof(double));
  if(NULL==nH0_out){free(nH0_out); printf("Memory allocation failed for nH0_out.\n"); exit(0);}

  nHep_out  = (double *)calloc(NT, sizeof(double));
  if(NULL==nHep_out){free(nHep_out);printf("Memory allocation failed for nHep_out.\n"); exit(0);}

  nHepp_out = (double *)calloc(NT, sizeof(double));
  if(NULL==nHepp_out){free(nHepp_out);printf("Memory allocation failed for nHepp_out.\n"); exit(0);}
  
  ne_out    = (double *)calloc(NT, sizeof(double));
  if(NULL==ne_out){free(ne_out);printf("Memory allocation failed for ne_out.\n"); exit(0);}
  
  heat_out    = (double *)calloc(NT, sizeof(double));
  if(NULL==heat_out){free(heat_out);printf("Memory allocation failed for heat_out.\n"); exit(0);} 

  gJH0_out    = (double *)calloc(NT, sizeof(double));
  if(NULL==gJH0_out){free(gJH0_out);printf("Memory allocation failed for gJH0_out.\n"); exit(0);} 
  
  gJHe0_out    = (double *)calloc(NT, sizeof(double));
  if(NULL==gJHe0_out){free(gJHe0_out);printf("Memory allocation failed for gJHe0_out.\n"); exit(0);} 
  
  gJHep_out    = (double *)calloc(NT, sizeof(double));
  if(NULL==gJHep_out){free(gJHep_out);printf("Memory allocation failed for gJHep_out.\n"); exit(0);} 
  
  epsH0_out    = (double *)calloc(NT, sizeof(double));
  if(NULL==epsH0_out){free(epsH0_out);printf("Memory allocation failed for epsH0_out.\n"); exit(0);} 
  
  epsHe0_out    = (double *)calloc(NT, sizeof(double));
  if(NULL==epsHe0_out){free(epsHe0_out);printf("Memory allocation failed for epsHe0_out.\n"); exit(0);} 
  
  epsHep_out    = (double *)calloc(NT, sizeof(double));
  if(NULL==epsHep_out){free(epsHep_out);printf("Memory allocation failed for epsHep_out.\n"); exit(0);} 
}
