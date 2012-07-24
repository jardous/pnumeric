/*
  $Id: kf.c,v 1.1 2007/01/21 01:38:35 jp Exp $

  kf.c
     This is implementation of Kalman filter in C. Created as part
     of pnumeric Python module.

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
#include "m2/m2.h"

// for stanalone console app
#include <sys/types.h>
#include <getopt.h>

#define PROGRAM_NAME "tr"

#include <unistd.h>
#ifndef STDIN_FILENO
# define STDIN_FILENO 0
#endif
#ifndef STDERR_FILENO
# define STDERR_FILENO 2
#endif

#define TRUE  1
#define FALSE 0

//#define DEBUG(...) fprintf(stderr, __VA_ARGS__)
#define DEBUG(...) // __VA_ARGS__

/*
Plant
  inputs                 outputs
    p                       q
         +------------+
   u1 ---|            |---- y1
         |    N=3     |---- y2
   u2 ---|            |---- y3
         +------------+

Example:
          A matrix: 3 x 3
             1   2   3
             4   5   6
             7   8  10
          B matrix: 3 x 2
            0  0
            0  1
            1  0
          C matrix: 3 x 3
            1  0  0
            0  1  0
            0  0  1
          D matrix: 3 x 3
            0  0
            0  0
            0  0

*/

/*
 * KF filter algorithm core - one cycle (tick)
 */
void
tick(Float *A /*NxN*/, Float *B  /*Nxp*/, Float *C/*qxN*/, Float *D/*pxq*/,
        int n, int p, int q,
        Float *yv_k/*qx1*/, Float *u_k/*px1*/,
        Float *x   /*Nx1*/, Float *P  /*NxN*/,
        Float *Q  /*NxN*/, Float *R   /*qxq*/,
        // outputs
        Float *x_est/*Nx1*/, Float *y_est/*qxN*/, Float *P_est/*NxN*/
        )
{
    Float *Ax, *Cx, *Bu, *x_hat, *P_hat, *AT, *AP, *APAT, *CT, *PestCT;
    Float *CPest, *K, *K1, *K2, *K2inv, *Cxest, *Du, *X1, *X2, *KX2, *se, *KC, *P1;
    
    // ESTIMATE (TIME UPDATE)
    // a priori estimate of x
    Ax = m_new(n, 1);
    m_mul(A, x, Ax, n, n, 1);
    Bu = m_new(n, 1);
    m_mul(B, u_k, Bu, n, p, 1);
    x_hat = m_new(n, 1);
    m_add(Ax, Bu, x_hat, n, 1);
    m_free(Ax);
    m_free(Bu);
    // a priori estimate of P
    AT = m_new(n, n); // transpose(A)
    m_transpose(A, AT, n, n);
    AP = m_new(n, n); // A*P
    m_mul(A, P, AP, n, n, n);
    APAT = m_new(n, n); // A*P*transpose(A)
    m_mul(AP, AT, APAT, n, n, n);
    P_hat = m_new(n, n);
    m_add(APAT, Q, P_hat, n, n); // P_hat = A*P*transpose(A) + Q
    m_free(AT);
    m_free(AP);
    m_free(APAT);
    
    // CORRECTION (MEASUREMENT UPDATE)
    // K1
    CT = m_new(n, q); // transpose C(q,N)
    m_transpose(C, CT, q, n);
    PestCT = m_new(n, q); // P_hat*transpose(C) (Nxq)
    m_mul(P_hat, CT, PestCT, n, n, q);
    
    K1 = m_new(q, q); //  C*P_hat*CT
    CPest = m_new(q, n); // C*P_hat
    m_mul(C, P_hat, CPest, q, n, n); // C*P_hat
    m_mul(CPest, CT, K1, q, n, q);   // K1 = C*P_hat*CT
    K2 = m_new(q, q);
    m_add(K1, R, K2, q, q);
    K2inv = m_new(q, q);
    m_inversion(K2, K2inv, q);        //          PestCT
    K = m_new(n, q);                  // Nxq! NxN      Nxq         qxq    qxN qxN    Nxq         qxq
    m_mul(PestCT, K2inv, K, n, q, q); // K = P_hat*transpose(C) * inverse( C*P_hat*transpose(C) + R )
    m_free(CT);
    m_free(PestCT);
    m_free(CPest);
    m_free(K1);
    m_free(K2);
    m_free(K2inv);
    
    Cxest = m_new(q, 1); // C*x_hat
    m_mul(C, x_hat, Cxest, q, n, 1);
    Du = m_new(q, 1); // D*u_k
    m_mul(D, u_k, Du, q, p, 1);
    
    X1 = m_new(q, 1); // yv_k - C*x_hat
    m_sub(yv_k, Cxest, X1, q, 1);
    X2 = m_new(q, 1); // yv_k - C*x_hat - D*u_k
    m_sub(X1, Du, X2, q, 1);
    m_free(X1);
    m_free(Cxest);
    
    KX2 = m_new(n, 1);          //              KX2
    m_mul(K, X2, KX2, n, 1, 1); // x = x_hat + K*( transpose(yv_k) - C*x_hat - D*u_k )
    m_free(X2);
    // set new x
    m_add(x_hat, KX2, x_est, n, 1);
    m_free(KX2);
    
    // P = (eye(states) - K*C) * P_hat
    se = m_new(n, n);
    m_eye(se, n, n);
    KC = m_new(n, n);
    P1 = m_new(n, n);
    m_mul(K, C, KC, n, q, n);
    m_sub(se, KC, P1, n, n);
    m_mul(P1, P_hat, P_est, n, n, n);
    m_free(se);
    m_free(P1);
    m_free(x_hat);
    m_free(P_hat);
    
    // y_est = C*x + D*u_k
    Cx = m_new(q, 1);
    m_mul(C, x_est, Cx, q, n, 1);
    m_add(Cx, Du, y_est, q, 1);
    m_free(Cx);
    m_free(Du);
}


void
process(Float *A /*NxN*/, Float *B  /*Nxp*/, Float *C/*qxN*/, Float *D/*pxq*/,
        int n,         // number of states N
        int p,         // number of inputs p
        int q,         // number of outputs q
        Float *yv,     // output values (qxlength)
        Float *u,      // input values  (pxlength)
        Float *x0,     // initial state (Nx1)
        Float *P0,     // initial covariance (NxN)
        Float *Q,      // (NxN)
        Float *R,      // (qxq)
        int length,    // length of input/ouput values vector
        // outputs
        Float *x_est,  // length * (Nx1)
        Float *y_est,  // length * (qx1)
        Float *P_est   // length * (NxN)
        )
{
    int i;
    Float *P, *x;
    
    x = m_new(n, 1);
    P = m_new(n, n);
    m_copy(x, x0, n, 1);
    m_copy(P, P0, n, n);
    
    for (i=0; i<length; i++) {
        tick(A, B, C, D, n, p, q, yv+i*q, u+i*p, x, P, Q, R, x_est+i*n, y_est+i*q, P_est+i*n*n);
        m_copy(x, x_est+i*n, n, 1);
        m_copy(P, P_est+i*n*n, n, n);
    }
}

/* testing */
int 
test(void)
{
    int i, len;
    len = 16;
    int n=1, p=1, q=1;
    Float A=1., B=0., C=1., D=0., *x_est, *y_est, *P_est, *u=NULL, *y, x0=1., P0=1., Q=1., R=1.;
    
    x_est = m_new(n, len);
    y_est = m_new(n, len);
    P_est = m_new(n, n*len);
    
    u = m_new(1, len);
    for (i=0; i<len; i++)
        *(u+i) = 0.0;
    
    y = m_new(1, len);
    for (i=0; i<len; i++) {
        *(y+i) = (Float)rand()/RAND_MAX;
    }
    
    process(&A, &B, &C, &D, n, p, q, y, u, &x0, &P0, &Q, &R, len, x_est, y_est, P_est);
    
    m_print(y, n, len, "y");
    
    m_print(y_est, n, len, "y_est");
    
    m_print(x_est, n, len, "x_est");
    m_print(P_est, n, len, "P_est");
    
    m_free(u);
    m_free(y);
    
    m_free(x_est);
    m_free(y_est);
    m_free(P_est);
    
    return 0;
}




static char eolchar = '\n';


struct option const long_options[] =
{
    {"initial state estimate x0", no_argument, NULL, 'x'},
    {"initial error covariance estimate P0", no_argument, NULL, 'P'},
    {"A", required_argument, 0, 'A'},
    {"B", required_argument, 0, 'B'},
    {"C", required_argument, 0, 'C'},
    {"D", required_argument, 0, 'D'},
    {"Q", required_argument, 0, 'Q'},
    {"R", required_argument, 0, 'R'},
  //{GETOPT_HELP_OPTION_DECL},
  //{GETOPT_VERSION_OPTION_DECL},
    {0, 0, 0, 0}
};



void
usage (int status)
{
    if (status != EXIT_SUCCESS)
        fprintf (stderr, "Try `%s --help' for more information.\n", PROGRAM_NAME);
    else
    {
        printf ("Usage: %s [OPTION]... SET1 [SET2]\n", PROGRAM_NAME);
        fputs ("\
                Process input file by Kalman filter from standard input,\n\
                writing to standard output.\n\
                \n\
                -x     initial state estimate (x0)\n\
                -P     initial error covariance estimate (P0)\n\
                -A     A system parameter\n\
                -B     B system parameter\n\
                -C     C system parameter\n\
                -D     D system parameter\n\
                -Q     Q covariance parameter\n\
                -R     R covariance parameter\n\
                ", stderr);
    }
}


int
main(int argc, char **argv)
{
    char c;//, read;
    FILE * file_in = NULL;
    const short BUFLEN=160;
    Float x=0., P=1., A=1., B=0., C=1., D=0., Q=1., R=100.;
    char eol = eolchar;
    Float decoded;
    char buf[BUFLEN];
    int idx=0;
    Float u_k=0., x_est, y_est, P_est;
    
    /* getopt_long stores the option index here. */
    int option_index = 0;
    while ((c = getopt_long(argc, argv, "x:P:A:B:C:D:Q:R:", long_options, &option_index)) != -1) {
        switch (c)
        {
            case 'x':
                x = atof(optarg);
                break;
            case 'P':
                P = atof(optarg);
                break;
            case 'A':
                A = atof(optarg);
                break;
            case 'B':
                B = atof(optarg);
                break;
            case 'C':
                C = atof(optarg);
                break;
            case 'D':
                D = atof(optarg);
                break;
            case 'Q':
                Q = atof(optarg);
                break;
            case 'R':
                R = atof(optarg);
                break;
            default:
                usage (EXIT_FAILURE);
                break;
        }
    }
    
    DEBUG("A: %g, B: %g, C: %g, D: %g\n", A, B, C, D);
    DEBUG("Q: %g, R: %g\n", Q, R);
    DEBUG("x: %g, P: %g\n", x, P);
    
    file_in = fdopen(STDIN_FILENO, "rt");
    fread (&c, 1, 1, file_in);
    while (1) {
        if (idx == BUFLEN) {
            fprintf(stderr, "error reading file\n");
            exit(1);
        }
        if (c == eol) {
            buf[idx] = "\0";
            decoded = atof(buf);
        
            tick(&A /*NxN*/, &B  /*Nxp*/, &C/*qxN*/, &D/*pxq*/,
                 1 /*n*/, 1 /*p*/, 1 /*q*/,
                 &decoded /*yv_k qx1*/, &u_k /*px1*/,
                 &x   /*Nx1*/, &P  /*NxN*/,
                 &Q  /*NxN*/, &R   /*qxq*/,
                 // outputs
                 &x_est /*Nx1*/, &y_est /*qxN*/, &P_est /*NxN*/
            );
            
            x = x_est;
            P = P_est;
            
            fprintf(stdout, "%f\n", x_est);
            
            idx = 0;
        } else {
            buf[idx] = c;
            idx++;
        }
        
        fread (&c, 1, 1, file_in);
        if (feof(file_in)) {break; }
    }
    
    fclose(file_in);
    
    exit(0);
}

