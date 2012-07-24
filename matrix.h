/*
  $Id: matrix.h,v 1.1 2007/01/21 01:38:35 jp Exp $

  matrix.h
      Matrix datatype related function declarations for pnumeric
      Python module.

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

#ifndef __MATRIX_H__
#define __MATRIX_H__

#include "Python.h"
#include "m2/m2.h"

typedef struct {
    PyObject_HEAD
    Float *data;
    int rows;
    int cols;
} MatrixObject;

PyAPI_DATA(PyTypeObject) MatrixType;

#define Matrix_Check(op) PyObject_TypeCheck(op, &MatrixType)

// create a new emtpy matrix object
PyAPI_FUNC(MatrixObject *) matrix_new(int rows, int cols);
// deallocating matrix form memory
PyAPI_FUNC(void) matrix_dealloc(MatrixObject *matrix);
// repr, printing the matrix
PyAPI_FUNC(PyObject *) matrix_repr(MatrixObject *v);
// allocate memory array for matrix data
PyAPI_FUNC(int) matrix_alloc(MatrixObject *mo, int rows, int cols);
// returning matrix row
PyAPI_FUNC(PyObject *) matrix_item(MatrixObject *a, int i);
// comparison of two matrixes
PyAPI_FUNC(int) matrix_cmp(MatrixObject *self, MatrixObject *other);
// matrix inversion
PyAPI_FUNC(PyObject *) matrix_inv(MatrixObject *self);
// return determinant
PyAPI_FUNC(PyObject *) matrix_det(MatrixObject *self);
// matrix add
PyAPI_FUNC(PyObject *) matrix_add(MatrixObject *self, MatrixObject *other);
// matrix negative
PyAPI_FUNC(PyObject *) matrix_negative(MatrixObject *self);
// matrix subtraction
PyAPI_FUNC(PyObject *) matrix_sub(MatrixObject *self, MatrixObject *other);
// matrix multiply
PyAPI_FUNC(PyObject *) matrix_mul(MatrixObject *self, MatrixObject *other);
// create new MatrixObject from list (tuple)
PyAPI_FUNC(PyObject *) MatrixObject_New(PyTypeObject *type, PyObject *args);
// return attribute value
PyAPI_FUNC(PyObject *) matrix_getattr(MatrixObject *self, char *name);
// a little hack - make a 0x0 matrix data set to value
MatrixObject * PyFloat2matrix(const Float value);
// to decide, if we can make arithmetic operation on current datatypes
int matrix_coerce(MatrixObject **v, PyObject **w);
// return eye matrix of given dimensions. Example of 4-dim eye:
PyAPI_FUNC(MatrixObject *) matrix_eye(PyObject *self, PyObject *args);
// retun new matrix object with ones at all positions
PyAPI_FUNC(MatrixObject *) matrix_ones(PyObject *self, PyObject *args);
// retun new matrix object with zeros at all positions
PyAPI_FUNC(MatrixObject *) matrix_zeros(PyObject *self, PyObject *args);


#endif /* matrix.h */

