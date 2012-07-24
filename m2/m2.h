/*
  $Id: m2.h,v 1.1 2007/01/21 01:38:35 jp Exp $

  m2.h
      Matrix manipulation routine declarations for C.

  Copyright (c) 2007 Jiří Popek <jiri.popek@gmail.com>
           thanks to David Jedelský <david.jedelsky@gmail.com>
  
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

#ifndef __M2_H__
#define __M2_H__

#if 1

#define Float double
#define Fabs  fabs
#define EPS     1e-16
#define NONZERO 1e-100

#else

#define Float float
#define Fabs  fabsf
#define EPS     1e-7
#deifne NONZERO 1e-16

#endif

extern int m_err;

Float* m_new(int rows, int cols);
Float* m_dup(Float* A, int nrow, int ncol);
Float* m_copy(Float *B, Float *A, int nrow, int ncol);
void   m_free(Float* m);
void   m_print(Float *a, int nrow, int ncol, char *s);
void   m_set0(Float *m, int rows, int cols);
void   m_set1(Float *m, int rows, int cols);
void   m_descramble(Float *A, Float *B, int nrow, int ncol, int *pr, int *pc);
void   m_descramble_rows(Float *A, Float *B, int nrow, int ncol, int *pr);
void   m_descramble_cols(Float *A, Float *B, int nrow, int ncol, int *pc);
int    m_LU(Float *A, int nrow, int *pr, int *per);
int    m_inversion(Float *A, Float *B, int nrow);
void   m_transpose(Float *A, Float *B, int nrow, int ncol);
Float  m_tr(Float *A, int nrow);
Float  m_prod(Float *A, int nrow);
Float  m_det(Float *A, int nrow);
void   m_add(Float *A, Float *B, Float *C, int nrow, int ncol);
void   m_add_scalar(Float a, Float *A, Float *B, int nrow, int ncol);
void   m_sub(Float *A, Float *B, Float *C, int nrow, int ncol);
void   m_scale(Float a, Float *A, int nrow, int ncol);
void   m_mul(Float *A, Float *B, Float *C, int nrowA, int ncolA, int ncolB);

void m_eye(Float *A, int nrow, int ncol);

#endif
