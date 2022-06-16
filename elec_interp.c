#define fnID 6
#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "elec_interp.h"

/********************************************************************************/

/* This routine is based on the publicly available code provided by
 Furlanetto & Johnson-Stoever, 2010, MNRAS, 404, 1869.  It has been
 modified to output a more finely spaced set of look-up tables that are
 faster to access.  The modifications to the public code are confined
 to MakeSecondaryTable. */

#include "parameters.h"
#include "global_vars.h"
#include "proto.h"

void MakeSecondaryTable(void)
{
  double xe,logxe,delta_logxe,E,logE,delta_logE;
  double fheat[NESECTAB+1],nionH1[NESECTAB+1],nionHe1[NESECTAB+1],nionHe2[NESECTAB+1];
  double tab_old[6],tab_new[6];
  
  int i,j,itab,maketable;
  
  FILE *output,*tab;
  
  char filename_out[200],filename_tab[200];
  
  
  tab_new[0] = (double)NXTAB;   
  tab_new[1] = (double)XMIN; 
  tab_new[2] = (double)XMAX;
  
  tab_new[3] = (double)NESECTAB; 
  tab_new[4] = (double)ESECMIN;        
  tab_new[5] = (double)ESECMAX;
  
  
  sprintf(filename_tab,"./secondaries/tablesize.dat");
  
  /* First check whether we need to change the table size */
  if(!(tab = fopen(filename_tab, "r")))
    {
      printf("Cannot read file:\n");
      printf("%s\n", filename_tab);
      maketable = 1;
    }
  else
    {
      fscanf(tab, "%lf %lf %lf\n", &tab_old[0],&tab_old[1],&tab_old[2]);
      fscanf(tab, "%lf %lf %lf\n", &tab_old[3],&tab_old[4],&tab_old[5]);
      fclose(tab);
      
      /* Check the existing table boundaries against those in global_vars.h */
      for(itab=0;itab<6;itab++)
	{
	  if(tab_old[itab]-tab_new[itab] != 0)
	    {
	      maketable = 1; 
	      break;
	    }
	  else 
	    maketable = 0;
	}
    }
  
  switch(maketable)
    {
      /* Use existing tables */
    case 0:
      break;
      /* Recompute the tables */
    case 1:
      printf("Recomputing secondaries table...");
      initialize_interp_arrays();
      
      delta_logxe = (XMAX - XMIN) / (double)(NXTAB);
      delta_logE  = (ESECMAX - ESECMIN)/ (double)(NESECTAB);
      
      /* Loop over electron fraction */
      for(i = 0; i <= NXTAB; i++)
	{
	  logxe = XMIN + delta_logxe * i;
	  xe = pow(10.0, logxe);
	  
	  sprintf(filename_out,"./secondaries/logxe_%1.3f.dat",logxe);
	  
	  /* Loop over photo-electron energy */
	  for(j = 0; j <= NESECTAB; j++)
	    {
	      logE = ESECMIN + delta_logE * j;

	      E = pow(10.0,logE);
	      
	      fheat[j]   = interp_fheat(E,xe);
	      nionH1[j]  = interp_nion_HI(E,xe);
	      nionHe1[j] = interp_nion_HeI(E,xe);
	      nionHe2[j] = interp_nion_HeII(E,xe);
	    }
	  
	  /* Write the look-up table */
	  if(!(output = fopen(filename_out, "wb")))
	    {
	      printf("Cannot write file:\n");
	      printf("%s\n", filename_out);
	      endrun(fnID);
	    }
	  fwrite(fheat,sizeof(double),NESECTAB+1,output);
	  fwrite(nionH1,sizeof(double),NESECTAB+1,output);
	  fwrite(nionHe1,sizeof(double),NESECTAB+1,output);
	  fwrite(nionHe2,sizeof(double),NESECTAB+1,output);
	  fclose(output);
	}
      printf("done\n\n");
      
      /* Save the table boundaries */
      if(!(tab = fopen(filename_tab, "w")))
	{
	  printf("Cannot write file:\n");
	  printf("%s\n", filename_tab);
	  endrun(fnID);
	}
      fprintf(tab, "%lf %lf %lf\n",tab_new[0],tab_new[1],tab_new[2]);
      fprintf(tab, "%lf %lf %lf\n",tab_new[3],tab_new[4],tab_new[5]);
      fclose(tab);
      break;
    }
  
  printf("Secondaries table:\n");
  printf("log(xe) table: N=%g min=%g max=%g\n",tab_new[0],tab_new[1],tab_new[2]);
  printf("log(E/eV) table: N=%g min=%g max=%g\n\n",tab_new[3],tab_new[4],tab_new[5]);
}

/****************************************************************************************/



// Functions to interpolate the energy deposition fractions of high-energy secondary electrons
// in the IGM.
// A call to initialize_interp_arrays() is necessary to read in and set up the tables.
// Then the functions
//     interp_fheat = heating fraction
//     interp_n_Lya = number of HI Lya photons produced
//     interp_nion_HI = number of HI ionizations produced
//     interp_nion_HeI = number of HeI ionizations produced
//     interp_nion_HeII = number of HeII ionizations produced
// Note all energies are in eV.
// So for example the fraction of energy in Lya photons is just (n_Lya*10.2)/(electron energy)

// The interpolation arrays.  Note all are defined in global scope, but x_int_ prefix should
// ensure no conflicts with other code
float x_int_Energy[x_int_NENERGY];
float x_int_XHII[x_int_NXHII];
float x_int_fheat[x_int_NXHII][x_int_NENERGY];
float x_int_n_Lya[x_int_NXHII][x_int_NENERGY];
float x_int_nion_HI[x_int_NXHII][x_int_NENERGY];
float x_int_nion_HeI[x_int_NXHII][x_int_NENERGY];
float x_int_nion_HeII[x_int_NXHII][x_int_NENERGY];

// Call once to read in data files and set up arrays for interpolation.
// All data files should be in a local subdirectory "x_int_tables/"; if moved, change input_base
// below to new location.
void initialize_interp_arrays()
{
  FILE *input_file;
  char input_file_name[100];
  char input_base[100] = "./secondaries/x_int_tables/";
  char input_tail[100] = ".dat";
  char mode[10] = "r";

  float xHI,xHeI,xHeII,z,T;
  float trash;
  //char label[] = " ";
  char label[100]; /* JSB 16/06/22: fix for fscanf error that occurs for some compiler options */
  
  int i;
  int n_ion;

  // Initialize array of ionized fractions
  x_int_XHII[0] = 1.0e-4;
  x_int_XHII[1] = 2.318e-4;
  x_int_XHII[2] = 4.677e-4;
  x_int_XHII[3] = 1.0e-3;
  x_int_XHII[4] = 2.318e-3;
  x_int_XHII[5] = 4.677e-3;
  x_int_XHII[6] = 1.0e-2;
  x_int_XHII[7] = 2.318e-2;
  x_int_XHII[8] = 4.677e-2;
  x_int_XHII[9] = 1.0e-1;
  x_int_XHII[10] = 0.5;
  x_int_XHII[11] = 0.9;
  x_int_XHII[12] = 0.99;
  x_int_XHII[13] = 0.999;

  for (n_ion=0;n_ion<x_int_NXHII;n_ion++) {

    // Construct filename
    if (x_int_XHII[n_ion] < 0.3) {
      sprintf(input_file_name,"%slog_xi_%1.1f%s",input_base,log10(x_int_XHII[n_ion]),input_tail);
    } else {
      sprintf(input_file_name,"%sxi_%1.3f%s",input_base,x_int_XHII[n_ion],input_tail);      
    }

    //    printf("%s\n",input_file_name);
    
    input_file = fopen(input_file_name, mode);

    if (input_file == NULL) {
      fprintf(stderr, "Can't open input file %s!\n",input_file_name);
      exit(1);
    }

    // Read in first line
    for (i=1;i<=5;i++) {
      fscanf(input_file,"%s", label);
      //      printf("%s\n",label);
    }
    
    // Read in second line (ionized fractions info)
    fscanf(input_file,"%g %g %g %g %g", &xHI, &xHeI, &xHeII, &z, &T);
    
    // Read in column headings
    for (i=1;i<=11;i++) {
      fscanf(input_file,"%s", label);
      //      printf("%s\n",label);
    }
    
    // Read in data table
    for (i=0;i<x_int_NENERGY;i++) {
      fscanf(input_file,"%g %g %g %g %g %g %g %g %g",
	     &x_int_Energy[i],
	     &trash,
	     &x_int_fheat[n_ion][i],
	     &trash,
	     &x_int_n_Lya[n_ion][i],
	     &x_int_nion_HI[n_ion][i],
	     &x_int_nion_HeI[n_ion][i],
	     &x_int_nion_HeII[n_ion][i],
	     &trash);
      //      printf("%g\t%g\t%g\t%g\t%g\t%g\n", x_int_Energy[i], x_int_fheat[n_ion][i], 
      //      	     x_int_n_Lya[n_ion][i], x_int_nion_HI[n_ion][i], x_int_nion_HeI[n_ion][i], 
      //      	     x_int_nion_HeII[n_ion][i]);
    }

    fclose(input_file);

  }
  return;
}

// Function to compute fheat for an interacting electron, given energy En and IGM ionized fraction
// xHII_call; note we assume xHeI=xHI and all the rest of the helium is in HeII (though that makes
// almost no difference in practice).
// Function returns the fraction of the this electron energy that is returned to heat, thus
// it does NOT include the primary photo-ionization energy loss
// Note that if En>highest element of array, it just uses that highest value.  Similarly, if 
// xHII_call is less than or smaller than the limits of the ionized fraction array (10^-4 and 
// 0.999), it just uses those values.
float interp_fheat(float En, float xHII_call)
{
  int n_low,n_high;
  int m_xHII_low,m_xHII_high;

  float elow_result,ehigh_result,final_result;

  // Check if En is inside interpolation boundaries
  if (En > 0.999*x_int_Energy[x_int_NENERGY-1]) {
    // If it is above the upper limit, we just assume that it is near the upper limit, which
    // has anyway reached the asymptotic limit
    En = x_int_Energy[x_int_NENERGY-1]*0.999;
  } else if (En < x_int_Energy[0]) {
    return 1.0;
  }

  // Check if ionized fraction is within boundaries; if not, adjust to be within
  if (xHII_call > x_int_XHII[x_int_NXHII-1]*0.999) {
    xHII_call = x_int_XHII[x_int_NXHII-1]*0.999;
  } else if (xHII_call < x_int_XHII[0]) {
    xHII_call = 1.001*x_int_XHII[0];
  }

  n_low = locate_energy_index(En);
  n_high = n_low + 1;
  //  printf("nlow=%d %f\n",n_low,En);

  m_xHII_low = locate_xHII_index(xHII_call);
  m_xHII_high = m_xHII_low + 1;
  //    printf("%f %f %f\n",xHII_call,XHII[m_xHII_low],XHII[m_xHII_high]);
  
  // First linear interpolation in energy
  elow_result = ((x_int_fheat[m_xHII_low][n_high]-x_int_fheat[m_xHII_low][n_low])/
		 (x_int_Energy[n_high]-x_int_Energy[n_low]));
  elow_result *= (En - x_int_Energy[n_low]);
  elow_result += x_int_fheat[m_xHII_low][n_low];
  //  printf("elow= %f\n",elow_result);
  
  // Second linear interpolation in energy
  ehigh_result = ((x_int_fheat[m_xHII_high][n_high]-x_int_fheat[m_xHII_high][n_low])/
		  (x_int_Energy[n_high]-x_int_Energy[n_low]));
  ehigh_result *= (En - x_int_Energy[n_low]);
  ehigh_result += x_int_fheat[m_xHII_high][n_low];
  //  printf("ehigh= %f\n",ehigh_result);
  
  // Final interpolation over the ionized fraction
  final_result = (ehigh_result - elow_result)/(x_int_XHII[m_xHII_high] - x_int_XHII[m_xHII_low]);
  final_result *= (xHII_call - x_int_XHII[m_xHII_low]);
  final_result += elow_result;
  //  printf("final= %f\n",final_result);
  
  return final_result;
}

// Function to compute nLya for an interacting electron, given energy En and IGM ionized fraction
// xHII_call; note we assume xHeI=xHI and all the rest of the helium is in HeII (though that makes
// almost no difference in practice).
// Function returns the number of Lyman-alpha photons produced as a function of the secondary
// electron energy, thus it does NOT include the primary photo-ionization energy loss
// Note that if En>highest element of array, it just uses that high value.  Similarly, if xHII_call
// is less than or smaller than the limits of the ionized fraction array (10^-4 and 0.999),
// it just uses those values.
float interp_n_Lya(float En, float xHII_call)
{
  int n_low,n_high;
  int m_xHII_low,m_xHII_high;

  float elow_result,ehigh_result,final_result;

  // Check if En is inside interpolation boundaries
  if (En > 0.999*x_int_Energy[x_int_NENERGY-1]) {
    // If it is above the upper limit, we just assume that it is near the upper limit, which
    // has anyway reached the asymptotic limit
    En = x_int_Energy[x_int_NENERGY-1]*0.999;
  } else if (En < x_int_Energy[0]) {
    return 0.0;
  }

  // Check if ionized fraction is within boundaries; if not, adjust to be within
  if (xHII_call > x_int_XHII[x_int_NXHII-1]*0.999) {
    xHII_call = x_int_XHII[x_int_NXHII-1]*0.999;
  } else if (xHII_call < x_int_XHII[0]) {
    xHII_call = 1.001*x_int_XHII[0];
  }

  n_low = locate_energy_index(En);
  n_high = n_low + 1;
  
  m_xHII_low = locate_xHII_index(xHII_call);
  m_xHII_high = m_xHII_low + 1;
  //    printf("%f %f %f\n",xHII_call,XHII[m_xHII_low],XHII[m_xHII_high]);
  
  // First linear interpolation in energy
  elow_result = ((x_int_n_Lya[m_xHII_low][n_high]-x_int_n_Lya[m_xHII_low][n_low])/
		 (x_int_Energy[n_high]-x_int_Energy[n_low]));
  elow_result *= (En - x_int_Energy[n_low]);
  elow_result += x_int_n_Lya[m_xHII_low][n_low];
  
  // Second linear interpolation in energy
  ehigh_result = ((x_int_n_Lya[m_xHII_high][n_high]-x_int_n_Lya[m_xHII_high][n_low])/
		  (x_int_Energy[n_high]-x_int_Energy[n_low]));
  ehigh_result *= (En - x_int_Energy[n_low]);
  ehigh_result += x_int_n_Lya[m_xHII_high][n_low];
  
  // Final interpolation over the ionized fraction
  final_result = (ehigh_result - elow_result)/(x_int_XHII[m_xHII_high] - x_int_XHII[m_xHII_low]);
  final_result *= (xHII_call - x_int_XHII[m_xHII_low]);
  final_result += elow_result;
  
  return final_result;
}

// Function to compute nHI for an interacting electron, given energy En and IGM ionized fraction
// xHII_call; note we assume xHeI=xHI and all the rest of the helium is in HeII (though that makes
// almost no difference in practice).
// Function returns the number of hydrogen ionizations produced as a function of the secondary
// electron energy, thus it does NOT include the primary photo-ionization energy loss
// Note that if En>highest element of array, it just uses that high value.  Similarly, if xHII_call
// is less than or smaller than the limits of the ionized fraction array (10^-4 and 0.999),
// it just uses those values.
float interp_nion_HI(float En, float xHII_call)
{
  int n_low,n_high;
  int m_xHII_low,m_xHII_high;

  float elow_result,ehigh_result,final_result;

  // Check if En is inside interpolation boundaries
  if (En > 0.999*x_int_Energy[x_int_NENERGY-1]) {
    // If it is above the upper limit, we just assume that it is near the upper limit, which
    // has anyway reached the asymptotic limit
    En = x_int_Energy[x_int_NENERGY-1]*0.999;
  } else if (En < x_int_Energy[0]) {
    return 0.0;
  }

  // Check if ionized fraction is within boundaries; if not, adjust to be within
  if (xHII_call > x_int_XHII[x_int_NXHII-1]*0.999) {
    xHII_call = x_int_XHII[x_int_NXHII-1]*0.999;
  } else if (xHII_call < x_int_XHII[0]) {
    xHII_call = 1.001*x_int_XHII[0];
  }

  n_low = locate_energy_index(En);
  n_high = n_low + 1;
  
  m_xHII_low = locate_xHII_index(xHII_call);
  m_xHII_high = m_xHII_low + 1;
  //    printf("%f %f %f\n",xHII_call,XHII[m_xHII_low],XHII[m_xHII_high]);
  
  // First linear interpolation in energy
  elow_result = ((x_int_nion_HI[m_xHII_low][n_high]-x_int_nion_HI[m_xHII_low][n_low])/
		 (x_int_Energy[n_high]-x_int_Energy[n_low]));
  elow_result *= (En - x_int_Energy[n_low]);
  elow_result += x_int_nion_HI[m_xHII_low][n_low];
  
  // Second linear interpolation in energy
  ehigh_result = ((x_int_nion_HI[m_xHII_high][n_high]-x_int_nion_HI[m_xHII_high][n_low])/
		  (x_int_Energy[n_high]-x_int_Energy[n_low]));
  ehigh_result *= (En - x_int_Energy[n_low]);
  ehigh_result += x_int_nion_HI[m_xHII_high][n_low];
  
  // Final interpolation over the ionized fraction
  final_result = (ehigh_result - elow_result)/(x_int_XHII[m_xHII_high] - x_int_XHII[m_xHII_low]);
  final_result *= (xHII_call - x_int_XHII[m_xHII_low]);
  final_result += elow_result;
  
  return final_result;
}

// Function to compute nHeI for an interacting electron, given energy En and IGM ionized fraction
// xHII_call; note we assume xHeI=xHI and all the rest of the helium is in HeII (though that makes
// almost no difference in practice).
// Function returns the number of HeI ionizations produced as a function of the secondary electron
// energy, thus it does NOT include the primary photo-ionization energy loss
// Note that if En>highest element of array, it just uses that high value.  Similarly, if xHII_call
// is less than or smaller than the limits of the ionized fraction array (10^-4 and 0.999),
// it just uses those values.
float interp_nion_HeI(float En, float xHII_call)
{
  int n_low,n_high;
  int m_xHII_low,m_xHII_high;

  float elow_result,ehigh_result,final_result;

  // Check if En is inside interpolation boundaries
  if (En > 0.999*x_int_Energy[x_int_NENERGY-1]) {
    // If it is above the upper limit, we just assume that it is near the upper limit, which
    // has anyway reached the asymptotic limit
    En = x_int_Energy[x_int_NENERGY-1]*0.999;
  } else if (En < x_int_Energy[0]) {
    return 0.0;
  }

  // Check if ionized fraction is within boundaries; if not, adjust to be within
  if (xHII_call > x_int_XHII[x_int_NXHII-1]*0.999) {
    xHII_call = x_int_XHII[x_int_NXHII-1]*0.999;
  } else if (xHII_call < x_int_XHII[0]) {
    xHII_call = 1.001*x_int_XHII[0];
  }

  n_low = locate_energy_index(En);
  n_high = n_low + 1;
  
  m_xHII_low = locate_xHII_index(xHII_call);
  m_xHII_high = m_xHII_low + 1;
  //    printf("%f %f %f\n",xHII_call,XHII[m_xHII_low],XHII[m_xHII_high]);
  
  // First linear interpolation in energy
  elow_result = ((x_int_nion_HeI[m_xHII_low][n_high]-x_int_nion_HeI[m_xHII_low][n_low])/
		 (x_int_Energy[n_high]-x_int_Energy[n_low]));
  elow_result *= (En - x_int_Energy[n_low]);
  elow_result += x_int_nion_HeI[m_xHII_low][n_low];
  
  // Second linear interpolation in energy
  ehigh_result = ((x_int_nion_HeI[m_xHII_high][n_high]-x_int_nion_HeI[m_xHII_high][n_low])/
		  (x_int_Energy[n_high]-x_int_Energy[n_low]));
  ehigh_result *= (En - x_int_Energy[n_low]);
  ehigh_result += x_int_nion_HeI[m_xHII_high][n_low];
  
  // Final interpolation over the ionized fraction
  final_result = (ehigh_result - elow_result)/(x_int_XHII[m_xHII_high] - x_int_XHII[m_xHII_low]);
  final_result *= (xHII_call - x_int_XHII[m_xHII_low]);
  final_result += elow_result;
  
  return final_result;
}

// Function to compute nHeII for an interacting electron, given energy En and IGM ionized fraction
// xHII_call; note we assume xHeI=xHI and all the rest of the helium is in HeII (though that makes
// almost no difference in practice).
// Function returns the number of HeII ionizations produced as a function of the secondary
// electron energy, thus it does NOT include the primary photo-ionization energy loss
// Note that if En>highest element of array, it just uses that high value.  Similarly, if xHII_call
// is less than or smaller than the limits of the ionized fraction array (10^-4 and 0.999),
// it just uses those values.
float interp_nion_HeII(float En, float xHII_call)
{
  int n_low,n_high;
  int m_xHII_low,m_xHII_high;

  float elow_result,ehigh_result,final_result;

  // Check if En is inside interpolation boundaries
  if (En > 0.999*x_int_Energy[x_int_NENERGY-1]) {
    // If it is above the upper limit, we just assume that it is near the upper limit, which
    // has anyway reached the asymptotic limit
    En = x_int_Energy[x_int_NENERGY-1]*0.999;
  } else if (En < x_int_Energy[0]) {
    return 0.0;
  }

  // Check if ionized fraction is within boundaries; if not, adjust to be within
  if (xHII_call > x_int_XHII[x_int_NXHII-1]*0.999) {
    xHII_call = x_int_XHII[x_int_NXHII-1]*0.999;
  } else if (xHII_call < x_int_XHII[0]) {
    xHII_call = 1.001*x_int_XHII[0];
  }

  n_low = locate_energy_index(En);
  n_high = n_low + 1;
  
  m_xHII_low = locate_xHII_index(xHII_call);
  m_xHII_high = m_xHII_low + 1;
  //    printf("%f %f %f\n",xHII_call,XHII[m_xHII_low],XHII[m_xHII_high]);
  
  // First linear interpolation in energy
  elow_result = ((x_int_nion_HeII[m_xHII_low][n_high]-x_int_nion_HeII[m_xHII_low][n_low])/
		 (x_int_Energy[n_high]-x_int_Energy[n_low]));
  elow_result *= (En - x_int_Energy[n_low]);
  elow_result += x_int_nion_HeII[m_xHII_low][n_low];
  
  // Second linear interpolation in energy
  ehigh_result = ((x_int_nion_HeII[m_xHII_high][n_high]-x_int_nion_HeII[m_xHII_high][n_low])/
		  (x_int_Energy[n_high]-x_int_Energy[n_low]));
  ehigh_result *= (En - x_int_Energy[n_low]);
  ehigh_result += x_int_nion_HeII[m_xHII_high][n_low];
  
  // Final interpolation over the ionized fraction
  final_result = (ehigh_result - elow_result)/(x_int_XHII[m_xHII_high] - x_int_XHII[m_xHII_low]);
  final_result *= (xHII_call - x_int_XHII[m_xHII_low]);
  final_result += elow_result;
  
  return final_result;
}

// Function to find bounding indices on the energy array, for an input energy En.
// Note it is done exactly for speed, since all the energy arrays have the same structure.
int locate_energy_index(float En) 
{
  int n_low;

  // Find energy table location analytically!
  if (En < 1008.88) {
    n_low = (int)(log(En/10.0)/1.98026273e-2);
  } else {
    n_low = 232 + (int)(log(En/1008.88)/9.53101798e-2);
  }

  return n_low;
}

// Function to find bounding indices on the ionized fraction array, for an input fraction
// xHII_call.  This is done by comparing to each element, since there are only 14 elements.
int locate_xHII_index(float xHII_call) 
{
  int m_xHII_low;

  // Determine relevant ionized fractions to interpolate; do iteratively because few elements
  m_xHII_low = x_int_NXHII - 1;
  while (xHII_call < x_int_XHII[m_xHII_low]) {
    m_xHII_low--;
  }
  return m_xHII_low;
}

