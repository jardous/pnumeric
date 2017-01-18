/*
  $Id: fft.c,v 1.1 2007/01/21 01:38:35 jp Exp $

  fft.c
     This is implementation of Fast Fourier Transform in C.
     Created as part of pnumeric Python module.

  Copyright (c) 2007 Jiří Popek <jiri.popek@gmail.com>

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License along
  with this program; if not, write to the Free Software Foundation, Inc.,
  51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define    PI 3.1415926535897932
#define TWOPI 6.2831853071795864

#include "m2/m2.h"

void FFT(Float *x, // data vector - will be replaced by FFT
         const int n)    // must be a power of 2)
{
    int i, j, k, n1, n2, m;
    Float c, s, e, a, t1, t2;
    Float *y;

    n2 = n/2;
    y = (Float*)malloc(n*sizeof(Float));

    for(i=0; i<n; i++) y[i]=0;

    // now get the HSB position
    i = 0;
    n1 = n-1;
    while(n1>0) {
        n1 = n1 >> 1;
        i++;
    }
    m = i;
    j = 0;  /* bit-reverse */
    for (i=1; i < n - 1; i++)
    {
        n1 = n2;
        while ( j >= n1 )
        {
            j = j - n1;
            n1 = n1/2;
        }
        j = j + n1;

        if (i < j)
        {
            t1 = x[i];
            x[i] = x[j];
            x[j] = t1;
            t1 = y[i];
            y[i] = y[j];
            y[j] = t1;
        }
    }

    n1 = 0;  /* FFT */
    n2 = 1;
    for (i=0; i<m; i++)
    {
        n1 = n2;
        n2 = n2 + n2;
        e = -TWOPI/n2;
        a = 0.0;

        for (j=0; j < n1; j++)
        {
            c = cos(a);
            s = sin(a);
            a = a + e;
            for (k=j; k < n; k=k+n2)
            {
                t1 = c*x[k+n1] - s*y[k+n1];
                t2 = s*x[k+n1] + c*y[k+n1];
                x[k+n1] = x[k] - t1;
                y[k+n1] = y[k] - t2;
                x[k] = x[k] + t1;
                y[k] = y[k] + t2;
            }
        }
    }

    free(y);
    return;
}
