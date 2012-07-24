/*
  $Id: m2_test.c,v 1.1 2007/01/21 01:38:35 jp Exp $

  m2_test.c
      Matrix manipulation routines test program.

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

#include <stdio.h>
#include "m2.h"

void
printper(int *p, int n, char *s)
{
  int i;
  printf("%s:\n",s);

  for(i=0; i<n; i++) printf("  %i",p[i]);
  printf("\n");
}

int main(int argc, char **argv)
{
  Float *A, *B, *C;
  int n;
  int s,per;
  int pr[20];
/*
  n=3;

  A = m_new(n,n);
  B = m_new(n,n);
  C = m_new(n,n);
  A[0]= 0.1; A[1]=0.2;
  A[2]=-0.2; A[3]=0.1;

  A[0]= 0.6; A[1]=0.4; A[2]= 0.2;
  A[3]=-0.2; A[4]=0.5; A[5]= 0.2;
  A[6]= 0.3; A[7]=0.2; A[8]=-0.4;

  m_print(A,n,n,"A");
  B = m_dup(A,n,n);
  s = m_LU(B, n, pr, &per);

  m_print(B,n,n,"LU");
  printper(pr,n,"pr");
  printf("s=%i  per=%i\n",s,per);

  C = m_dup(B,n,n);

  B[0]=1;B[1]=0;B[2]=0;
  C[3]=0;B[4]=1;B[5]=0;
  C[6]=0;C[7]=0;B[8]=1;

  m_mul(B,C,A,n,n,n);
  m_print(A,n,n,"L.U");

  m_inversion(A,B,n);
  m_print(B,n,n,"inv(A)");

  m_mul(A,B,C,n,n,n);
  m_print(C,n,n,"A.inv(A)");

  m_free(A);
  m_free(B);
  m_free(C);
*/
  n=4;
  A=m_new(n,n);
  B=m_new(n,n);
  C=m_new(n,n);


  A[ 0]=  4;  A[ 1]=  3;  A[ 2]= 2;  A[ 3]=-1;
  A[ 4]=  8;  A[ 5]=  7;  A[ 6]=-6;  A[ 7]= 5;
  A[ 8]= 12;  A[ 9]=-11;  A[10]=10;  A[11]= 9;
  A[12]=-16;  A[13]= 15;  A[14]=14;  A[15]=13;

  m_print(A,n,n,"A");

  B = m_dup(A,n,n);
  s = m_LU(B, n, pr, &per);
  m_print(B,n,n,"LU");
  printper(pr,n,"pr");
  printf("s=%i  per=%i\n",s,per);

  printf("Det(A)=%g\n",m_det(A,n));
  m_inversion(A,B,n);
  m_print(B,n,n,"inv(A)");

  printf("Det(B)=%g\n",m_det(B,n));


  A[ 0]= 0.5;  A[ 1]= 0.5;  A[ 2]= 0.0; A[ 3]= 0.5;
  A[ 4]=-0.25; A[ 5]=-0.25; A[ 6]= 0.5; A[ 7]=-0.75;
  A[ 8]=-1.5;  A[ 9]=-0.5;  A[10]= 0.0; A[11]=-2.5;
  A[12]= 0.25; A[13]= 0.25; A[14]= 0.5; A[15]= 0.75;

  m_print(A,n,n,"A");
  printf("Det(A)=%g\n",m_det(A,n));
  m_inversion(A,B,n);
  m_print(B,n,n,"inv(A)");

  printf("Det(B)=%g\n",m_det(B,n));

  m_free(A);
  m_free(B);

  A = m_eye(3,5);
  m_print(A, 3, 5, "eye 3, 5");
  m_free(A);

  return 0;
}

