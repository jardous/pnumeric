/*
  $Id: vector.h,v 1.1 2007/01/21 01:38:35 jp Exp $

  vector.h
      Vector datatype related function declarations for pnumeric
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

#ifndef __VECTOR_H__
#define __VECTOR_H__

#include "Python.h"
#include "m2/m2.h"

typedef struct {
    PyObject_HEAD

    // set this if vector stands for matrix row or a vector slice
    PyObject *object; // pointer to another object (matrix, vector slice)
    Float *p_data; // data pointer into existing data (if vector stands for matrix row of another vector)

    // this is useless if vector is matrix row
    int length;
    Float *data;
} VectorObject;

PyAPI_DATA(PyTypeObject) VectorType;

#define Vector_Check(op) PyObject_TypeCheck(op, &VectorType)

// create new vector object
PyAPI_FUNC(VectorObject *) matrixrow2vector(PyObject *object, Float *p_data);
// create new vector object
PyAPI_FUNC(VectorObject *) vector_new(int length);
// return vector length
PyAPI_FUNC(Py_ssize_t) vector_length(VectorObject *v);
// return vector data pointer
Float * vector_dataptr(VectorObject *v);
// deallocating vector object form memory
PyAPI_FUNC(void) vector_dealloc(VectorObject *v);
// repr, printing of vector or vector row
PyAPI_FUNC(PyObject *) vector_repr(VectorObject *a);
// return vector + number
PyAPI_FUNC(PyObject *) vector_add(VectorObject *v, PyObject *val);
// return vector - number
PyAPI_FUNC(PyObject *) vector_sub(VectorObject *v, PyObject *val);
// return vector item
PyAPI_FUNC(PyObject *) vector_item(VectorObject *a, int i);
// set value to vector row item
PyAPI_FUNC(int) vector_ass_item(VectorObject *self, int idx, PyObject *val);
// comparison of two vectors
PyAPI_FUNC(int) vector_cmp(VectorObject *self, VectorObject *other);
// create new VectorObject from list (tuple)
PyAPI_FUNC(PyObject *) VectorObject_New(PyTypeObject *type, PyObject *args);
// to decide, if we can make arithmetic operation on current datatypes
int vector_coerce(PyObject **v, PyObject **w);
// returns ref (not copy) to current vector object
PyAPI_FUNC(PyObject *) vector_slice(VectorObject *self, int ilow, int ihigh);
// returns a range vector
PyAPI_FUNC(VectorObject *) vector_range(PyObject *self, PyObject *args, PyObject *kws);

#endif /* vector.h */
