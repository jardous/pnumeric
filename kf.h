/*

  $Id: kf.h,v 1.1 2007/01/21 01:38:35 jp Exp $

  kf.h
     This declaration of Kalman filter functions in C. Created as part
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

#ifndef __KF_H__
#define __KF_H__


void
tick(Float *A /*NxN*/, Float *B  /*Nxp*/, Float *C/*qxN*/, Float *D/*pxq*/,
        int n, int p, int q,
        Float *yv_k/*qx1*/, Float *u_k/*px1*/,
        Float *x   /*Nx1*/, Float *P  /*NxN*/,
        Float *Q  /*NxN*/, Float *R   /*qxq*/,
        // outputs
        Float *x_est/*Nx1*/, Float *y_est/*qxN*/, Float *P_est/*NxN*/
);



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
);


#endif

