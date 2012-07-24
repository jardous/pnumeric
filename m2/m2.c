/*
  $Id: m2.c,v 1.1 2007/01/21 01:38:35 jp Exp $

  m2.c
      Matrix manipulation routines for C.

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
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "m2.h"

#define m_Malloc malloc
#define m_Free   free



/*
 Matrix - array of numbers


 element:  a11 | a12 | ... | a1n | a21 | a22 | ...
   index:   0     1          n-1    n    n+1   ...


   for matrix with m rows and n cols:

      Aij = a[(i-1)*n+(j-1)],     i=1..m, j=1..n

*/

/* global error indicator */
int m_err=0;

/*
 allocates a new matrix
 */
Float*
m_new(int rows, int cols)
{
  return (Float*)m_Malloc(sizeof(Float)*rows*cols);
}


Float*
m_dup(Float* A, int nrow, int ncol)
{
  Float *B = m_new(nrow, ncol);

  memcpy(B,A,nrow*ncol*sizeof(Float));

  return B;
}


Float*
m_copy(Float *B, Float *A, int nrow, int ncol)
{
  memcpy(B,A,nrow*ncol*sizeof(Float));

  return B;
}


void
m_free(Float* m)
{
  m_Free(m);
}


void
m_print(Float *a, int nrow, int ncol, char *s)
{
  int i,j;

  printf("%s:\n",s);
  for (i=0; i<nrow; i++){
    for (j=0; j<ncol; j++,a++)
      printf("  %4g", (Fabs(*a)>EPS)?*a:0);
    printf("\n");
  }
}


void
m_set0(Float *m, int rows, int cols)
{
  int i;

  for(i=rows*cols; i>0; i--)
   *m++=0;
}


void
m_set1(Float *m, int rows, int cols)
{
  int i;

  for(i=rows*cols; i>0; i--)
   *m++=1;
}


void
m_eye(Float *A, int nrow, int ncol)
{
  int i;
  
  m_set0(A, nrow, ncol);
  for(i=0; i<nrow; i++)
     *(A + i*ncol + i) = 1.0;
}

/*
  Descramble the matrix A into matrix B acording to
  row (pr) and column (pc) permutation. It means

    B[pr[i]][pc[j]] = A[i][j]

 */
void
m_descramble(Float *A, Float *B, int nrow, int ncol, int *pr, int *pc)
{
  Float *p, *q;
  int i,j;

  p = A;
  for(i=0; i<nrow; i++,p+=nrow){
    q = B + nrow * pr[i];
    for(j=0; j<ncol; j++){
      q[pc[j]] = p[j];
    }
  }
}


/*
  Descramble the matrix A into matrix B acording to
  row (pr) and column (pc) permutation. It means

    B[pr[i]][pc[j]] = A[i][j]

 */
void
m_descramble_cols(Float *A, Float *B, int nrow, int ncol, int *pc)
{
  Float *p, *q;
  int i,j;

  p = A; q = B;
  for(i=0; i<nrow; i++,p+=nrow,q+=nrow)
    for(j=0; j<ncol; j++)
      q[pc[j]] = p[j];
}


/*
  Descramble the matrix A into matrix B acording to
  row (pr) and column (pc) permutation. It means

    B[pr[i]][pc[j]] = A[i][j]

 */
void
m_descramble_rows(Float *A, Float *B, int nrow, int ncol, int *pr)
{
  Float *p;
  int i, rowsz = nrow*sizeof(Float);

  p = A;
  for(i=0; i<nrow; i++,p+=nrow)
    memcpy(B+nrow*pr[i],p,rowsz);
}


/*
  LU decomposition
 */
int
m_LU(Float *A, int nrow, int *pr, int *per)
{
  Float *p,*q,*sf,*ex;
  Float *Ar,*Ai,*Ak=NULL;
  Float a,b;
  int r,s,i,k=0;
  int rowsz = nrow*sizeof(Float);
  int g=0; /* singularity indicator */

  *per=1;

  sf=m_new(2,nrow);
  if(sf==NULL) return 2;
  ex=sf+nrow;

  /* scale factors for implicit pivoting */
  p=A; q=sf;
  for(r=0; r<nrow; r++){
    if(pr) pr[r]=r;
    a=0;
    for(i=0; i<nrow; i++){
      b=Fabs(*p++);
      if(b>a) a=b;
    }
    if(a==0) return -1; /* hard singularity */
    *q++=1.0/a;
  }

  /* LU */
  for(i=0; i<nrow; i++){                       /* columns */

    Ar=A;                                     /* U elements */
    for(r=0; r<i; r++,Ar+=nrow){
      a=0; p=Ar; q=A+i;
      for(s=0; s<r; s++,q+=nrow) a+=*p++**q;
      *q-=a;
    }


    b=0; Ai=Ar;                           /* L elements without 1/Uii */
    for(/*r=i*/; r<nrow; r++,Ar+=nrow){   /*     + partial pivoting   */
      a=0; p=Ar; q=A+i;
      for(s=0; s<i; s++,q+=nrow) a+=*p++**q;
      *p-=a; a=Fabs(*p)*sf[r];
      if(a>b){ b=a; k=r; Ak=Ar; }
    }

    if(i+1 == nrow) break; /* last col, no more row xchgs and 1/Uii */

    if(i!=k){ /* interchange the pivoting row */
      *per=-*per;
      sf[k]=sf[i];
      memcpy(ex,Ai,rowsz); memcpy(Ai,Ak,rowsz); memcpy(Ak,ex,rowsz);
      if(pr){ s=pr[i]; pr[i]=pr[k]; pr[k]=s; }
    }

    p=Ai+i;
    /* Avoid division by zero to give result even in singular cases */
    if(*p==0.0) { *p=NONZERO; g=1; }
    a=1.0 / *p;

    p+=nrow;
    for(r=i+1; r<nrow; r++,p+=nrow) *p*=a;
  }

  m_Free(sf);

  return g;
}



/*
  compute matrix inversion
  Return: 0 - OK
          1 - Singular
          2 - Alloc error
 */
int
m_inversion(Float *A, Float *B, int nrow)
{
  int per,*pr;
  int i,s,r,k,v;
  Float a,b;
  Float *C,*Bk,*Bi,*Ck,*p,*q,*t,*u;
  int g;

  pr=(int*)malloc(sizeof(int)*nrow);
  if( pr==NULL ) return 2;

  C = m_dup(A,nrow,nrow);
  if( C==NULL) { m_Free(pr); return 2;}

  g = m_LU(C, nrow, pr, &per);                       /* Pivoting + LU */
  if(g) goto error;

  Bk=B; Ck=C;
  for(k=0; k<nrow; k++,Bk+=nrow,Ck+=nrow){           /* inv(L) : C->B */
    u=Ck;
    for(i=k-1; i>=0; i--,u-=nrow){
      a=0;
      v=i+1; p=Bk+v; q=u+i; //q=C[v][i];
      for(; v<k; v++,q+=nrow)
	a+=*p++**q;
      Bk[i]=-*q-a;
    }
  }

  Bk = B; Ck = C;
  for(k=0; k<nrow; k++,Bk+=nrow,Ck+=nrow){           /* inv(U) : C->B */
    b=Bk[k]=1/Ck[k];
    p=Bk-nrow; q=Ck-nrow;
    for(i=k-1; i>=0; i--,q-=nrow,p-=nrow){
      a=0; t=p+i; u=q+k;
      for(v=i; v<k; v++,u+=nrow)
	a+=*t++**u;
      *t=-a*b;
    }
  }

  Bi=B; u=C;
  for(i=0; i<nrow; i++,Bi+=nrow){    /* inv(A) = inv(U)*inv(L) : B->C */
    q=B;
    for(s=0; s<nrow; s++,q+=nrow){
      if(i<=s){
	p=Bi+s; a=*p++; t=q+nrow+s;
	for(r=s+1; r<nrow; r++,t+=nrow)
	  a+=*p++**t;
      }
      else{
	a=0; p=Bi+i; t=Bi+s;
	for(r=i; r<nrow; r++,t+=nrow)
	  a+=*p++**t;
      }
      *u++=a;
    }
  }

  m_descramble_cols(C,B,nrow,nrow,pr); /* Descramble C->B */

error:
  m_Free(C);
  m_Free(pr);

  return g;
}



/*
  Transpose matrix A into B
 */
void
m_transpose(Float *A, Float *B, int nrow, int ncol)
{
  int i,j;
  Float *p,*q,*r;

  p=A;
  r=B;
  for(i=0; i<nrow; i++){
   q=r;
   for(j=0; j<ncol; j++,q+=nrow)
    *q=*p++;
   r++;
  }
}


/*
 trace
 */
Float
m_tr(Float *A, int nrow)
{
  Float *p, a;
  int i,j;

  p=A; a=0;
  j=nrow+1;
  for(i=nrow; i>0; i--,p+=j)
   a+=*p;

  return a;
}


/*
 product
 */
Float
m_prod(Float *A, int nrow)
{
  Float *p, a;
  int i,j;

  p=A; a=1;
  j=nrow+1;
  for(i=nrow; i>0; i--,p+=j)
   a*=*p;

  return a;
}


/*
 Determinant
 */
Float
m_det(Float *A, int nrow)
{
  int per, s;
  Float *B, a;

  B = m_dup(A,nrow,nrow);

  s=m_LU(B, nrow, NULL, &per);

  if(s){
    a=0;
  }else{
    a=m_prod(B,nrow);
    a=(per==1)? a:-a;
  }
  m_Free(B);

  return a;
 }


/*
 C = A + B
 */
void
m_add(Float *A, Float *B, Float *C, int nrow, int ncol)
{
  int i;

  for(i=nrow*ncol; i>0; i--)
    *C++=*A+++*B++;

}

/*
 B = a + A
 */
void
m_add_scalar(Float a, Float *A, Float *B, int nrow, int ncol)
{
  int i;

  for(i=nrow*ncol; i>0; i--)
    *B++=a+*A++;

}

/*
 C = A - B
 */
void
m_sub(Float *A, Float *B, Float *C, int nrow, int ncol)
{
  int i;

  for(i=nrow*ncol; i>0; i--)
    *C++=*A++-*B++;

}


/*
 A = a * A
 */
void
m_scale(Float a, Float *A, int nrow, int ncol)
{
  int i;

  for(i=nrow*ncol; i>0; i--)
    *A++*=a;

}


/*
 C = A * B


 nrowC = nrowA
 ncolC = ncolB
 nrowB = ncolA

 */
void
m_mul(Float *A, Float *B, Float *C, int nrowA, int ncolA, int ncolB)
{
  Float *p,*q,*r,*s,*t;
  Float a;
  int i,j,k;

  p = C; s = A;
  for(i=nrowA; i>0; i--){
    t = B;
    for(j=ncolB; j>0; j--){
      q = s;
      r = t;
      a=0;
      for(k=ncolA; k>0; k--,r+=ncolB)
	a+=*q++**r;
      *p++=a;
      t++;
    }
    s+=ncolA;
  }

}
