/*
  $Id: fft.h,v 1.1 2007/01/21 01:38:35 jp Exp $

  fft.h
     This is declaration of Fast Fourier Transform in C.
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

#ifndef __FFT_H__
#define __FFT_H__

#include "m2/m2.h"

void FFT(Float *x, const int N);

#endif /* fft.h */
