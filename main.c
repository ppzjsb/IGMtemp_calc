#define fnID 1
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>


#include "parameters.h"
#include "proto.h"
#include "global_vars.h"


double *z_out,*temp_out,*nH0_out,*nHep_out,*nHepp_out,*ne_out,*heat_out,*t_out;
double *gJH0_out,*gJHe0_out,*gJHep_out;
double *epsH0_out,*epsHe0_out,*epsHep_out;

/**********************************************************************************/

int main()
{
  
  double dlna,dtdlna,atime,atime_end,dt,t=0.0;
  double temp,rhocb,H,yhelium,rho,u,pcdone=0.0;
  double nHp,nHep,nHepp,ne,mu,nHcgs,logxe,zdec;
  double start,end,tot_cpu_used;
  int i;
  
  start = clock();
 
  
  H       = 1.0e7/MPC;
  yhelium = (1.0 - XH) / (4.0 * XH);
  rhocb   = 3.0*(H*HUBBLE)*(H*HUBBLE)*OMEGAB*pow(10.0,LOGOVERDEN)/(8.0*M_PI*GRAVITY);
     
  /* Initialise */
  redshift  = ZSTART;
  atime     = 1.0/(1.0+ZSTART);
  atime_end = 1.0/(1.0+ZEND);
  dlna      = (log(atime_end)-log(atime))/((double)NT-1.0);
  
  nHp   = 0.0; 
  nHep  = 0.0;
  nHepp = 0.0;
  ne    = nHp + nHep + 2.0*nHepp;;
 
  
  /* Thermal decoupling redshift for CMB, appropriate for z<zdec: see
     Eq (66) in Furlanetto et al. 2006, PhR, 433, 181.  Use RECFAST if
     higher accuracy desired. */
  zdec  = 150.0 * pow((OMEGAB*HUBBLE*HUBBLE)/0.023, 0.4);
  temp  = 2.726/(atime*atime*(1.0+zdec));

  logxe = ne > 0.0 ? log10(ne/(1.0+2.0*yhelium)) : -10.0;

  
  /* Load the external Furlanetto & Johnson-Stoever (2010, MNRAS, 404,
     1869) secondaries tables */
#ifdef SECONDARY
  MakeSecondaryTable();
  InitSecondaryMemory();
  ReadSecondaryTable();
#endif
  
  gausslegendre();
  
  InitCool(temp, logxe);

  InitDumpMemory();
  

  for(i=0;i<NT;i++)
    {
     
      Ha       = H*HUBBLE*sqrt(OMEGAM/(atime*atime*atime) + OMEGAL);
      dtdlna   = 1.0 / Ha;
      dt       = dtdlna * dlna;
      redshift = 1.0/atime - 1.0;
    
      
#ifdef PLAW_UVB
      IonizationRates(temp, logxe);
#endif  
      
      rho   = rhocb / (atime*atime*atime);
      nHcgs = XH * rho / PROTONMASS;
      
      
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
    
      
      /* Specific heating rate, erg s^-1 g^-1 */
      heat_out[i] = nHcgs*((1.0-nHp)*epsH0 + (yhelium-nHep-nHepp)*epsHe0 + nHep*epsHep) / rho;
      
      
      /* Specific internal energy erg g^-1 or (cm/s)^2, since energy is g cm^2 s^-2 */
      mu = (1.0 + 4.0 * yhelium) / (1.0 + yhelium + ne);
      u  = (BOLTZMANN * temp)/(GAMMA_MINUS1 * mu * PROTONMASS);

      /* Solve for specific internal energy [erg g^-1] and abundances */
      u = DoCooling(u, rho, dt, &nHp, &nHep, &nHepp);
      
      ne    = nHp + nHep + 2.0*nHepp;
      logxe = ne > 0.0 ? log10(ne/(1.0+2.0*yhelium)) : XMIN;
      
      mu    = (1.0 + 4.0 * yhelium) / (1.0 + yhelium + ne);
      temp  = GAMMA_MINUS1 / BOLTZMANN * u * PROTONMASS * mu;
      
      if(atime/atime_end > pcdone || i == NT-1)
	{
	  printf("z=%1.3f, T=%1.1fK, fH0=%1.3e, fHep=%1.3e, t=%1.2f Gyr\n",
		 redshift,temp_out[i],nH0_out[i],nHep_out[i]/yhelium,t_out[i]/(1.0e9*3600*24*365.0));
	  // printf("%e %e\n",gJH0_out[i],epsH0_out[i]);
	  // printf("%e %e\n",gJHep_out[i],epsHep_out[i]);
	  pcdone += 0.05;
	}
      
      
      atime = exp(log(atime) + dlna);
      t += dt; /* Note that t=0 corresponds to ZSTART, not z --> infty */
    }
  
  write_data();
  
  
  end = clock();
  tot_cpu_used = (end - start) / CLOCKS_PER_SEC;
  printf("\nTotal CPU time used %lf s\n",tot_cpu_used);
  
  return 0;
}


/*****************************************************************************/

void write_data(void)
{
  
  double om_out,ol_out,ob_out,h_out,xh_out,density_out;
  int nbins;
  FILE *output;  
  char fname[400];
  
  
  //sprintf(fname,"../plaw_data/%s",OUTFILE);

  sprintf(fname,"./output/%s",OUTFILE);
  if(!(output = fopen(fname, "wb")))
    {
      printf("Cannot read file:\n");
      printf("%s\n", fname);
      exit(0);
    }
  
  nbins       = NT;
  density_out = LOGOVERDEN;
  om_out      = OMEGAM;
  ol_out      = OMEGAL;
  ob_out      = OMEGAB;
  h_out       = HUBBLE;
  xh_out      = XH;


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
  fclose(output);

}


/*****************************************************************************/

void InitDumpMemory(void)
{
  z_out     = (double *)calloc(NT, sizeof(double));
  if(NULL==z_out){free(z_out);printf("Memory allocation failed for z_out.\n");
    exit(0);}
  
  t_out     = (double *)calloc(NT, sizeof(double));
  if(NULL==t_out){free(t_out);printf("Memory allocation failed for t_out.\n");
    exit(0);}
  
  temp_out  = (double *)calloc(NT, sizeof(double));
  if(NULL==temp_out){free(temp_out);printf("Memory allocation failed for temp_out.\n");
    exit(0);}
  
  nH0_out   = (double *)calloc(NT, sizeof(double));
  if(NULL==nH0_out){free(nH0_out);printf("Memory allocation failed for nH0_out.\n");
    exit(0);}

  nHep_out  = (double *)calloc(NT, sizeof(double));
  if(NULL==nHep_out){free(nHep_out);printf("Memory allocation failed for nHep_out.\n");
    exit(0);}

  nHepp_out = (double *)calloc(NT, sizeof(double));
  if(NULL==nHepp_out){free(nHepp_out);printf("Memory allocation failed for nHepp_out.\n");
    exit(0);}
  
  ne_out    = (double *)calloc(NT, sizeof(double));
  if(NULL==ne_out){free(ne_out);printf("Memory allocation failed for ne_out.\n");
    exit(0);}
  
  heat_out    = (double *)calloc(NT, sizeof(double));
  if(NULL==heat_out){free(heat_out);printf("Memory allocation failed for heat_out.\n");
    exit(0);} 

  gJH0_out    = (double *)calloc(NT, sizeof(double));
  if(NULL==gJH0_out){free(gJH0_out);printf("Memory allocation failed for gJH0_out.\n");
    exit(0);} 
  
  gJHe0_out    = (double *)calloc(NT, sizeof(double));
  if(NULL==gJHe0_out){free(gJHe0_out);printf("Memory allocation failed for gJHe0_out.\n");
    exit(0);} 
  
  gJHep_out    = (double *)calloc(NT, sizeof(double));
  if(NULL==gJHep_out){free(gJHep_out);printf("Memory allocation failed for gJHep_out.\n");
    exit(0);} 
  
  epsH0_out    = (double *)calloc(NT, sizeof(double));
  if(NULL==epsH0_out){free(epsH0_out);printf("Memory allocation failed for epsH0_out.\n");
    exit(0);} 
  
  epsHe0_out    = (double *)calloc(NT, sizeof(double));
  if(NULL==epsHe0_out){free(epsHe0_out);printf("Memory allocation failed for epsHe0_out.\n");
    exit(0);} 
  
  epsHep_out    = (double *)calloc(NT, sizeof(double));
  if(NULL==epsHep_out){free(epsHep_out);printf("Memory allocation failed for epsHep_out.\n");
    exit(0);} 
}

/*****************************************************************************/
