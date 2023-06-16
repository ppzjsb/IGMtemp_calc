
/**************************************************************************/

/**! \file utils.c
 *
 * \brief Some general utility functions
 *
 */

/***************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "parameters.h"
#include "proto.h"



/** \brief Returns the maximum of two doubles
 *
 * \param x 
 * \param y
 *
 * \return Maximum of x or y
 */

double dmax(double x, double y)
{
  if(x > y)
    return x;
  else
    return y;
}



/** \brief Returns the minimum of two doubles
 *
 * \param x 
 * \param y
 *
 * \return Minimum of x or y
 */

double dmin(double x, double y)
{
  if(x < y)
    return x;
  else
    return y;
}

