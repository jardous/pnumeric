/*
  $Id:

  window.c
     This is implementation of windowing in C.
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

#include "m2/m2.h"
#include "window.h"
#include <math.h>

#define    PI 3.1415926535897932
#define TWOPI 6.2831853071795864

// check http://en.wikipedia.org/wiki/Window_function for more window types

/*
 * Rectangular window
 */
void rect(const int N, Float* v)
{
    int i;
    
    for(i=0; i<N; i++) {
        *(v+i) = 1.0;
    }
}



/*
 * Hann window
 */
void hanning(const int N, Float* v)
{
    int i;
    
    if (N == 1) {
        *v = 1.0;
        return;
    }
    
    for(i=0; i<N; i++) {
        v[i] = 0.5 * (1 - cos(TWOPI*i/(N-1)));
    }
}



/*
 * Hamming window
 */
void hamming(const int N, Float* v)
{
    int i;
    
    if (N == 1) {
        *v = 1.0;
        return;
    }
    
    for(i=0; i<N; i++) {
        v[i] = 0.54 - 0.46 * cos(TWOPI*i/(N-1));
    }
}

