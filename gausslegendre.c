#define fnID 5
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "parameters.h"
#include "proto.h"
#include "global_vars.h"


/***************************************************************************/

/* Calculate Gauss-Legendre weights */

void gausslegendre()
{
  InitGLMemory();

  absc_and_weight();
}


/***************************************************************************/

void absc_and_weight()
{
  double x, p1,p2,pp,x1;
  double eps = 3.0e-14;
  
  int i, m;   
  
  m = (NWEIGHTS+1)/2;
  
  for(i=1;i<=m;i++)
    {
  
      x = cos(M_PI*(i - 0.25)/(NWEIGHTS + 0.5));
      
      do
	{
	  p1 = plgndr(NWEIGHTS, 0, x);
	  p2 = plgndr(NWEIGHTS-1, 0, x);
	  
	  pp = NWEIGHTS*(x*p1 - p2)/(x*x - 1.0);
	  x1 = x;
	  x  = x1-p1/pp;
	}
      while(fabs(x-x1) >= eps); 
      
      absc[i-1] = -1.0 * x;
      absc[NWEIGHTS-i]   = x;
      weight[i-1]        = 2.0/((1.0 - x*x)*pp*pp);
      weight[NWEIGHTS-i] = weight[i-1];
    }
}

/***************************************************************************/

/*  Computes the associated Legendre polynomial P_l^m(x).  Here m and
    l are integers satisfying 0<=m<=l, while x lies in the range
    -1<=x<=1. */

double plgndr(int l, int m, double x)
{
  double fact,pll,pmm,pmmp1,somx2;
  int i,ll;
  
  if (m < 0 || m > l || fabs(x) > 1.0)
    printf("Bad arguments in routine PLGNDR\n");
  
  pmm = pll = 1.0;
  if (m > 0) 
    {
      /* Compute P_m^m */
      somx2 = sqrt((1.0 - x) * (1.0 + x));
      fact = 1.0;
      for(i=1; i<=m; i++) 
	{
	  pmm *= -fact*somx2;
	  fact += 2.0;
	}
    }
  
  if (l == m)
    return pmm;
  else 
    {
      /* Compute P_m^m+1 */
      pmmp1=x*(2*m+1)*pmm;
      
      if (l == (m+1))
	return pmmp1;
      else 
	{
	  /* Compute P_l^m, l>m+1 */
	  for (ll=(m+2); ll<=l; ll++) 
	    {
	      pll   = (x*(2*ll-1)*pmmp1-(ll+m-1)*pmm)/(ll-m);
	      pmm   = pmmp1;
	      pmmp1 = pll;
	    }
	  return pll;
	}
    }
  }
  

/***************************************************************************/

void InitGLMemory(void)
{
  
  weight = (double *) calloc(NWEIGHTS, sizeof(double));
  if(NULL==weight){free(weight);printf("Memory allocation failed for weights.\n");
    endrun(fnID);}
  
  absc = (double *) calloc(NWEIGHTS, sizeof(double));
  if(NULL==absc){free(absc);printf("Memory allocation failed for absc.\n");
    endrun(fnID);}
}

/***************************************************************************/
