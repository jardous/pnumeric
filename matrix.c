/*
  $Id: matrix.c,v 1.1 2007/01/21 01:38:35 jp Exp $

  matrix.c
      Matrix datatype related functions for pnumeric Python module.

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

#include "Python.h"
#include "math.h"  // for EPS
#include "m2/m2.h" // matrixes stuff

#include "pnumeric.h"



/*
 * deallocating matrix form memory
 */
PyAPI_FUNC(void) matrix_dealloc(MatrixObject *matrix)
{
    m_free(matrix->data);
    PyObject_Del(matrix);
}



/*
 * repr, printing the matrix
 */
PyAPI_FUNC(PyObject *)
matrix_repr(MatrixObject *v)
{
    int i, j;
    Float d;
    char fstr[64];
    PyObject *s;
    PyObject *result;
    
    if (v->rows == 0 || v->cols == 0) {
        return PyString_FromString("Matrix([])");
    }
    
    result = PyString_FromString("Matrix([\n");
    for (i = 0; i < v->rows; i++) {
        for (j = 0; j < v->cols; j++) {
            d = v->data[i*(v->cols) + j];
            sprintf(fstr, "   %6g", d);
            PyString_ConcatAndDel(&result, PyString_FromString(fstr));
        }
        PyString_ConcatAndDel(&result, PyString_FromString("\n"));
    }
    
    PyString_ConcatAndDel(&result, PyString_FromString("])"));
    return result;
}



/*
 * allocate memory array for matrix data
 * return 0 - success, 1 - failed
 */
PyAPI_FUNC(int)
matrix_alloc(MatrixObject *mo, int rows, int cols)
{
    mo->data = m_new(rows, cols);
    if (mo->data == NULL) {
        PyErr_NoMemory();
        return 1;
    }
    
    mo->rows = rows;
    mo->cols = cols;
    
    return 0;
}



/*
 * returning matrix row
 */
PyAPI_FUNC(PyObject *)
matrix_item(MatrixObject *a, int i)
{
    PyObject *out = NULL;
    if (i < 0 || i >= a->rows) {
        PyErr_SetObject(PyExc_IndexError, PyString_FromString("row index out of range"));
        return NULL;
    }
    
    out = (PyObject*)matrixrow2vector(a, (a->data) + i*(a->cols));
    //TODO: is this necessary?
    Py_INCREF(out);
    
    return out;
}


/*
 * comparison of two matrixes
 */
PyAPI_FUNC(int)
matrix_cmp(MatrixObject *self, MatrixObject *other)
{
    Float *s, *o, sub, max;
    Float eps = 1e-8; // TODO: hope this will be enough
    int i;
    
    if (self->rows != other->rows || self->cols != other->cols) {
        return 1;
    }
    
    s = self->data;
    o = other->data;
    
    for (i=0; i<(self->rows*self->cols); i++) {
        sub = *s-*o;
        max = MAX(*s, *o);
        
        if (max != 0) {
            sub = sub / max;
        }
        
        s++;
        o++;
        
        if (Fabs(sub) > eps) {
            return 1;
        }
    }
    
    // TODO: which matrix is bigger? (to return -1 instead of 1)
    return 0;
}

/*
 * matrix inversion
 */
PyAPI_FUNC(PyObject *)
matrix_inv(MatrixObject *self)
{
    MatrixObject *out;
    
    if (self->rows != self->cols) {
        PyErr_SetString(PyExc_ValueError, "not a square matrix");
        return NULL;
    }
    
    out = matrix_new(self->rows, self->cols);
    m_inversion(self->data, out->data, self->rows);
    Py_INCREF(out);
    
    return (PyObject *)out;
}

/*
 * return determinant
 */
PyAPI_FUNC(PyObject *)
matrix_det(MatrixObject *self)
{
    Float det;
    
    if (self->rows != self->cols) {
        PyErr_SetString(PyExc_ValueError, "not a square matrix");
        return NULL;
    }
    
    det = m_det(self->data, self->rows);
    return PyFloat_FromDouble(det);
}



/*
 * matrix add
 */
PyAPI_FUNC(PyObject *)
matrix_add(MatrixObject *self, MatrixObject *other)
{
    MatrixObject *out;
    
    if (self->rows==0 && self->cols==0) { // first matrix is scalar
        out = matrix_new(other->rows, other->cols);
        m_add_scalar(*(self->data), other->data, out->data, other->rows, other->cols);
    } else if (other->rows==0 && other->cols==0) { // second matrix is scalar
        out = matrix_new(self->rows, self->cols);
        m_add_scalar(*(other->data), self->data, out->data, self->rows, self->cols);
    } else {
        if (self->rows != other->rows && self->cols != other->cols) {
            PyErr_SetString(PyExc_ValueError, "Matrixes are not aligned");
            return NULL;
        }
        out = matrix_new(self->rows, other->cols);
        m_add(self->data, other->data, out->data, self->rows, self->cols);
    }
    
    Py_INCREF(out);
    return (PyObject*)out;
}



/*
 * matrix negative
 */
PyAPI_FUNC(PyObject *)
matrix_negative(MatrixObject *self)
{
    MatrixObject *out;
    
    out = matrix_new(self->rows, self->cols);
    m_copy(out->data, self->data, self->rows, self->cols);
    m_scale(-1, out->data, out->rows, out->cols);
    Py_INCREF(out);
    return (PyObject*)out;
}



/*
 * matrix subtraction
 */
PyAPI_FUNC(PyObject *)
matrix_sub(MatrixObject *self, MatrixObject *other)
{
    MatrixObject *out;
    Float *tmp;
    
    if (self->rows==0 && self->cols==0) { // first matrix is scalar
        out = matrix_new(other->rows, other->cols);
        tmp = m_new(other->rows, other->cols);
        m_copy(tmp, other->data, other->rows, other->cols);
        m_scale(-1, tmp, other->rows, other->cols);
        m_add_scalar(*(self->data), tmp, out->data, other->rows, other->cols);
        m_free(tmp);
    } else if (other->rows==0 && other->cols==0) { // second matrix is scalar
        out = matrix_new(self->rows, self->cols);
        m_add_scalar(-*(other->data), self->data, out->data, self->rows, self->cols);
    } else {
        if (self->rows != other->rows && self->cols != other->cols) {
            PyErr_SetString(PyExc_ValueError, "Matrixes are not aligned");
            return NULL;
        }
        
        out = matrix_new(self->rows, other->cols);
        m_sub(self->data, other->data, out->data, self->rows, self->cols);
    }
    
    Py_INCREF(out);
    return (PyObject*)out;
}



/*
 * matrix multiply
 */
PyAPI_FUNC(PyObject *)
matrix_mul(MatrixObject *self, MatrixObject *other)
{
    MatrixObject *out;
    DEBUG(" Matrix multiply\n");
    
    if (self->rows==0 && self->cols==0) { // first matrix is scalar
        out = matrix_new(other->rows, other->cols);
        m_copy(out->data, other->data, other->rows, other->cols);
        m_scale(*(self->data), out->data, other->rows, other->cols);
    } else if (other->rows==0 && other->cols==0) { // second matrix is scalar
        out = matrix_new(self->rows, self->cols);
        m_copy(out->data, self->data, self->rows, self->cols);
        m_scale(*(other->data), out->data, self->rows, self->cols);
    } else {
        out = matrix_new(self->rows, other->cols);
        if (self == NULL) {
            PyErr_NoMemory();
            return NULL;
        }
        if (self->rows != other->cols) {
            PyErr_SetString(PyExc_ValueError, "Matrixes not aligned");
            return NULL;
        }
        m_mul(self->data, other->data, out->data, self->rows, self->cols, other->cols);
    }
    
    Py_INCREF(out);
    return (PyObject*)out;
}



PyMethodDef MatrixObject_methods[] = {
    {"inv", (PyCFunction)matrix_inv, METH_NOARGS, "matrix inversion"},
    {"det", (PyCFunction)matrix_det, METH_NOARGS, "determinant of matrix"},
    {NULL}  /* Sentinel */
};



/*
 * return attribute value
 */
PyAPI_FUNC(PyObject *)
matrix_getattr(MatrixObject *self, char *name)
{
    PyObject *result=NULL;
    
    if (strcmp(name, "shape") == 0) {
        result = PyTuple_New(2);
        PyTuple_SetItem(result, 0, PyInt_FromLong(self->rows));
        PyTuple_SetItem(result, 1, PyInt_FromLong(self->cols));
    }
    else if (strcmp(name, "rows") == 0) {
        result = PyInt_FromLong(self->rows);
    }
    else if (!strcmp(name, "cols")) {
        result = PyInt_FromLong(self->cols);
    }
    else  if (strcmp(name, "__members__") == 0) {
        result = PyList_New(3);
        if (result) {
            PyList_SetItem(result, 0, PyString_FromString("shape"));
            PyList_SetItem(result, 1, PyString_FromString("rows"));
            PyList_SetItem(result, 2, PyString_FromString("cols"));
            if (PyErr_Occurred()) {
                Py_DECREF(result);
                result = NULL;
            }
        }
    }
    
    if (result == NULL)
        return Py_FindMethod(MatrixObject_methods, (PyObject *)self, name);
    else
        return result;
}



/*
 * a little hack - make a 0x0 matrix data set to value
 */
MatrixObject *
PyFloat2matrix(const Float value)
{
    MatrixObject *out;
    
    out = matrix_new(1, 1);
    *(out->data) = value;
    out->rows = 0;
    out->cols = 0;
    
    return out;
}



/*
 * to decide, if we can make arithmetic operation on current datatypes
 */
int
matrix_coerce(MatrixObject **v, PyObject **w)
{
    if (PyInt_Check(*w)) {
      *w = (PyObject*)PyFloat2matrix((Float)PyInt_AsLong(*w));
      Py_INCREF(*w);
      return 0;
    }
    
    if (PyFloat_Check(*w)) {
      *w = (PyObject*)PyFloat2matrix((Float)PyFloat_AsDouble(*w));
      Py_INCREF(*w);
      return 0;
    }
    
    return 1; // conversion not possible
}

/*
 * create a new emtpy matrix object
 */
PyAPI_FUNC(MatrixObject *)
matrix_new(int rows, int cols)
{
    MatrixObject *y=NULL;
    
    y = (MatrixObject *) PyMem_Malloc(sizeof(MatrixType));
    if (y == NULL) {
        PyErr_SetString(PyExc_MemoryError, "Can't allocate memory for Matrix");
        return NULL;
    }
    
    PyObject_Init((PyObject *)y, &MatrixType);
    if (matrix_alloc(y, rows, cols)) {
        PyErr_SetString(PyExc_MemoryError, "Can't allocate memory for Matrix data");
        Py_DECREF(y);
        return NULL;
    }
    
    y->rows = rows;
    y->cols = cols;
    
    return y;
}



/*
 * return eye matrix of given dimensions. Example of 4-dim eye:
 * 1 0 0 0
 * 0 1 0 0
 * 0 0 1 0
 * 0 0 0 1
 */
PyAPI_FUNC(MatrixObject *)
matrix_eye(PyObject *self, PyObject *args)
{
    int i;
    long size;
    MatrixObject *result;
    
    if (!PyArg_ParseTuple(args, "l", &size)) {
        return (MatrixObject*)Py_None;
    }
    
    result = matrix_new(size, size);
    if (!result) return (MatrixObject*)Py_None;
    m_set0(result->data, size, size);
    
    for (i=0; i<size; i++) {
        *(result->data + i*size + i) = 1.0;
    }
    
    return result;
}



/*
 * retun new matrix object with ones at all positions
 */
PyAPI_FUNC(MatrixObject *)
matrix_ones(PyObject *self, PyObject *args)
{
    long size;
    MatrixObject *result;
    
    if (!PyArg_ParseTuple(args, "l", &size)) {
        return (MatrixObject*)Py_None;
    }
    
    result = matrix_new(size, size);
    
    if (!result)
        return (MatrixObject*)Py_None;
    
    m_set1(result->data, size, size);
    return result;
}



/*
 * create new MatrixObject from list (tuple)
 */
PyAPI_FUNC(PyObject *)
MatrixObject_New(PyTypeObject *type, PyObject *args)
{
    MatrixObject *self;
    PyObject *item=NULL, *w=NULL;
    int i, cols=0, rows=0;
    PyObject * (*getitem)(PyObject *, int);
    PyObject *listObject;
    
    if (!PyArg_ParseTuple(args, "O", &listObject)) {
        PyErr_BadArgument();
        return NULL;
    }
    
    // [0]
    if (PyList_Check(listObject)) {
        rows = PyList_Size(listObject);
        getitem = PyList_GetItem;
    } else if (PyTuple_Check(listObject)) {
        rows = PyTuple_Size(listObject);
        getitem = PyTuple_GetItem;
    } else if (listObject->ob_type == &MatrixType) {
        Py_INCREF(listObject);
        return listObject;
    } else {
        PyErr_BadArgument();
        return NULL;
    }
    
    if (rows == 0) {
        return Py_None;
    }
    
    // [0][0]
    item = (*getitem)(listObject, 0);
    if (PyList_Check(item)) {
        cols = PyList_Size(item);
    } else if (PyTuple_Check(item)) {
        cols = PyTuple_Size(item);
    } else {
        PyErr_BadArgument();
        return NULL;
    }
    
    self = matrix_new(rows, cols);
    if (self == NULL) {
        PyErr_NoMemory();
        return NULL;
    }
    
    if (cols > 0) {
        for (i=0; i<rows; i++) {
            w = (*getitem)(listObject, i);
            if (!PyArg_GetDoubleArray(w, 1, 0, cols, self->data+cols*i)) {
                PyErr_BadArgument();
                return NULL;
            }
        }
    }
    
    Py_INCREF(self);
    return (PyObject*)self;
}



/*
 * retun new matrix object with zeros at all positions
 */
PyAPI_FUNC(MatrixObject *)
matrix_zeros(PyObject *self, PyObject *args)
{
    long size;
    MatrixObject *result;
    
    if (!PyArg_ParseTuple(args, "l", &size)) {
        return (MatrixObject*)Py_None;
    }
    
    result = matrix_new(size, size);
    
    if (!result)
        return (MatrixObject*)Py_None;
    
    m_set0(result->data, size, size);
    return result;
}



PyAPI_DATA(PySequenceMethods) matrix_as_sequence = {
    0,                     /* sq_length */
    0,//(binaryfunc)row_concat,               /*sq_concat - numerical add */
    0,//(intargfunc)matrix_mul,      /*sq_repeat -  numerical multiply TODO: - how to multiply by float? */
    (ssizeargfunc)matrix_item,      /* sq_item */
    0,     /* sq_slice */
    0,//list_ass_item, /* sq_ass_item */
    0,//list_ass_slice,   /* sq_ass_slice */
    0,          /* sq_contains */
    0,    /* sq_inplace_concat */
    0,  /* sq_inplace_repeat */
};



PyAPI_DATA(PyNumberMethods) matrix_as_number = {
    (binaryfunc)matrix_add,          /*nb_add*/
    (binaryfunc)matrix_sub,          /*nb_subtract*/
    (binaryfunc)matrix_mul,          /*nb_multiply*/
    0,//(binaryfunc)array_divide,    /*nb_divide*/
    0,//(binaryfunc)array_remainder, /*nb_remainder*/
    0,//(binaryfunc)array_divmod,    /*nb_divmod*/
    0,//(ternaryfunc)array_power,    /*nb_power*/
    (unaryfunc)matrix_negative,       /*nb_neg*/
    0,//(unaryfunc)_array_copy_nice, /*nb_pos*/
    0,//(unaryfunc)array_absolute,   /*(unaryfunc)array_abs,*/
    0,//(inquiry)_array_nonzero,     /*nb_nonzero*/
    0,//(unaryfunc)array_invert,     /*nb_invert*/
    0,//proxy_lshift,               /* nb_lshift */
    0,//proxy_rshift,               /* nb_rshift */
    0,//proxy_and,              /* nb_and */
    0,//proxy_xor,              /* nb_xor */
    0,//proxy_or,               /* nb_or */
    // more on coercing: http://www.ragestorm.net/tutorials/25/pyextnum.html
    (coercion)matrix_coerce,               /* nb_coerce */
    0,//matrix_int,              /* nb_int */
    0,//matrix_long,             /* nb_long */
    0,//matrix_float,            /* nb_float */
    0,//matrix_oct,              /* nb_oct */
    0,//matrix_hex,              /* nb_hex */
    /*This code adds augmented assignment functionality*/
    /*that was made available in Python 2.0*/
    0,//(binaryfunc)matrix_add,          /*inplace_add*/
    0,//(binaryfunc)array_inplace_subtract,     /*inplace_subtract*/
    0,//(binaryfunc)matrix_inplace_multiply,     /*inplace_multiply*/
    0,//(binaryfunc)array_inplace_divide,       /*inplace_divide*/
    0,//(binaryfunc)array_inplace_remainder,    /*inplace_remainder*/
    0,//(ternaryfunc)array_inplace_power,       /*inplace_power*/
    0,//(binaryfunc)array_inplace_lshift,       /*inplace_lshift*/
    0,//(binaryfunc)array_inplace_rshift,       /*inplace_rshift*/
    0,//(binaryfunc)array_inplace_bitwise_and,  /*inplace_and*/
    0,//(binaryfunc)array_inplace_bitwise_xor,  /*inplace_xor*/
    0,//(binaryfunc)array_inplace_bitwise_or,   /*inplace_or*/
};



PyAPI_DATA(PyTypeObject) MatrixType = {
    PyObject_HEAD_INIT(NULL)
    0,                          /*ob_size*/
    "pnumeric.Matrix",          /*tp_name*/
    sizeof(MatrixObject),       /*tp_basicsize*/
    0,                          /*tp_itemsize*/
    (destructor)matrix_dealloc, /*tp_dealloc*/
    0,                          /*tp_print*/
    (getattrfunc)matrix_getattr,/*tp_getattr*/
    0,                          /*tp_setattr*/
    (cmpfunc)matrix_cmp,        /*tp_compare*/
    (reprfunc)matrix_repr,      /*tp_repr*/
    &matrix_as_number,          /*tp_as_number*/
    &matrix_as_sequence,        /*tp_as_sequence*/
    0,                          /*tp_as_mapping*/
    0,                          /*tp_hash */
    0,                          /*tp_call*/
    (reprfunc)matrix_repr,       /*tp_str*/
    0,//(getattrofunc)matrix_getattr,/*tp_getattro*/
    0,//(setattrofunc)0,            /*tp_setattro*/
    0,                          /*tp_as_buffer*/
    Py_TPFLAGS_DEFAULT,         /*tp_flags*/
    "Matrix object",            /* tp_doc */
    0,                          /* tp_traverse */
    0,                          /* tp_clear */
    0,                          /* tp_richcompare */
    0,                          /* tp_weaklistoffset */
    0,                          /* tp_iter */
    0,                          /* tp_iternext */
    MatrixObject_methods,       /* tp_methods */
    0, //MatrixObject_members,       /* tp_members */
    0,                         /* tp_getset */
    0,                         /* tp_base */
    0,                         /* tp_dict */
    0,                         /* tp_descr_get */
    0,                         /* tp_descr_set */
    0,                         /* tp_dictoffset */
    0,                         /* tp_init */
    0,                         /* tp_alloc */
    0,//MatrixObject_New,      /* tp_new */
};

