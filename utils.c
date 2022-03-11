#define fnID 3
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


#include "parameters.h"
#include "proto.h"

/*  This function aborts the simulation */
void endrun(int ierr)
{
  if(ierr)
    {
      fprintf(stdout,"Run terminated at level %d\n\n\n", ierr);
      exit(1);
    }
  exit(0);
}

/* returns the maximum of two doubles
 */
double dmax(double x, double y)
{
  if(x > y)
    return x;
  else
    return y;
}

/* returns the minimum of two doubles */
double dmin(double x, double y)
{
  if(x < y)
    return x;
  else
    return y;
}
