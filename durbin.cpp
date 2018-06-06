#include <math.h>
#include "durbin.h"

/*********************************************************************
*                                                                    *
*   Levinson-Durbin algorithm for LPC coefficients                   *
*                                                                    *
*   input:                                                           *
*   n  = LPC order                                                   *
*   r -> r[0] to r[n]    autocorrelation values                      *
*                                                                    *
*   output:                                                          *
*   a -> a[0] to a[n]    LPC coefficients, a[0] = 1.0                *
*   k -> k[0] to k[n]    reflection coefficients, k[0] = unused      *
*                                                                    *
*                                                                    *
*   Author:                                                          *
*   Robert Bristow-Johnson, comp.dsp, 04.01.2011                     *
*   groups.google.com/forum/#!topic/comp.dsp/s4s1F_w5c30             *
*                                                                    *
*********************************************************************/

#define N 128

float durbin(float *r, float *a, int n, float k_max)
{
    double a_temp[N], alpha, epsilon, ki;       /* n <= N = constant  */
    int i, j;

    if (n > N)
    {
        return 0;
    }

    for (i = 0; i < n; i++)
    {
        a[i] = 0;
    }

    alpha = r[0];

    for (i = 0; i < n; i++)
    {
        epsilon = r[i+1];
        for (j = 0; j<i; j++)
        {
            epsilon += a[j] * r[i - j];
        }

        ki = -epsilon / alpha;

        if (fabs(ki) > k_max)
        {
            return (float) alpha;
        }

        a[i] = (float) ki;

        alpha = alpha*(1.0 - ki * ki);

        for (j = 0; j<i; j++)
        {
            a_temp[j] = a[j] + ki * a[i - j - 1];   /* update a[] array into
                                                       temporary array */
        }

        for (j = 0; j<i; j++)
        {
            a[j] = (float) a_temp[j];               /* update a[] array */
        }
    }

    return (float) alpha;
}

/********************************************************************/
