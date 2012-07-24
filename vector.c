/*
  $Id: vector.c,v 1.1 2007/01/21 01:38:35 jp Exp $

  vector.c
      Vector datatype related functions for pnumeric Python module.

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
 * create new empty matrix row object
 */
PyAPI_FUNC(VectorObject *)
matrixrow2vector(PyObject *object, Float *p_data)
{
    VectorObject *y=NULL;
    y = (VectorObject*) PyMem_Malloc(sizeof(VectorType));
    if (y == NULL) {
        PyErr_SetString(PyExc_MemoryError, "Can't allocate memory for vector object");
        return NULL;
    }
    PyObject_Init((PyObject *)y, &VectorType);
    
    DEBUG("matrixrow2vector\n");
    if (object!=NULL) {
        DEBUG("   object NULL\n");
        if (Matrix_Check(object)) {
            DEBUG("   is Matrix\n");
            y->object = object;
        }
    }
    
    y->data = NULL;
    y->p_data = p_data;
    y->length = -1;
    Py_INCREF(object);
    
    return y;
}



/*
 * create new empty vector
 */
PyAPI_FUNC(VectorObject *)
vector_new(int length)
{
    DEBUG("vector_new(%d)", length);
    VectorObject *y = NULL;
    y = (VectorObject*) PyMem_Malloc(sizeof(VectorType));
    if (y == NULL) {
        Py_DECREF(y);
        PyErr_SetString(PyExc_MemoryError, "Can't allocate memory for vector object");
        return NULL;
    }
    
    PyObject_Init((PyObject *)y, &VectorType);
    y->object = NULL;
    y->data = m_new(1, length);
    if (y->data == NULL) {
        Py_DECREF(y);
        PyErr_SetString(PyExc_MemoryError, "Can't allocate memory for vector data");
        return NULL;
    }
    y->p_data = NULL;
    y->length = length;
    
    return y;
}



/*
 * return vector length
 */
PyAPI_FUNC(Py_ssize_t)
vector_length(VectorObject *v)
{
    if (v->length >= 0) {
        DEBUG("vector length - NOT Matrix\n");
        return v->length;
    }
    
    if ((v->object != NULL) && (v->object->ob_type==&MatrixType)) {
        DEBUG("vector length - is Matrix\n");
        return ((MatrixObject*)(v->object))->cols;
    }
    
    return 0;
}



/*
 * return vector data pointer
 */
PyAPI_FUNC(Float *)
vector_dataptr(VectorObject *v)
{
    if (v->length >= 0)
        return v->data;
    
    if (v->object != NULL)
        return v->p_data;
    
    return NULL;
}



/*
 * deallocating vector object form memory
 */
PyAPI_FUNC(void)
vector_dealloc(VectorObject *v)
{
    if (v->object != NULL) {
        Py_DECREF(v->object); // decrement the matrix pointer if exists
    }
    
    PyObject_Del(v);
}



/*
 * repr, printing of vector or matrix row
 */
PyAPI_FUNC(PyObject *)
vector_repr(VectorObject *a)
{
    //char buf[256];
    PyObject *comma;
    int i, len;
    Float *data;
    char fstr[64];
    PyObject *result;
    
    data = vector_dataptr(a);
    len = vector_length(a);
    DEBUG("vector_repr - length:%d\n", len);
    
    if (len == 0) {
        return PyString_FromString("Vector([])");
    }
    
    result = PyString_FromString("Vector([");
    comma = PyString_FromString(", ");
    
    for (i = 0; i < len && !PyErr_Occurred(); i++) {
        if (i > 0)
            PyString_Concat(&result, comma);
        
        sprintf(fstr, "%g", data[i]);
        PyString_ConcatAndDel(&result, PyString_FromString(fstr));
    }
    
    Py_XDECREF(comma);
    PyString_Concat(&result, PyString_FromString("])"));
    
    return result;
}



/*
 * return vector + number
 */
PyAPI_FUNC(PyObject *)
vector_add(VectorObject *v, PyObject *val)
{
    VectorObject *out;
    Float x, *data;
    int len=0;
    
    DEBUG("vector_add\n");
    if PyLong_Check(val)
        x = PyLong_AsDouble(val);
    else if PyInt_Check(val)
        x = PyInt_AsLong(val);
    else if PyFloat_Check(val)
        x = PyFloat_AsDouble(val);
    else {
        PyErr_SetString(PyExc_ValueError, "left must be number");
        return NULL;
    }
    
    data = vector_dataptr(v);
    len = vector_length(v);
    
    out = vector_new(len);
    m_add_scalar(x, data, vector_dataptr(out), 1, len);
    Py_INCREF(out);
    
    return (PyObject*)out;
}



/*
 * return vector - number
 */
PyAPI_FUNC(PyObject *)
vector_sub(VectorObject *v, PyObject *val)
{
    VectorObject *out;
    Float x, *data;
    int len=0;
    
    DEBUG("vector_sub\n");
    if PyLong_Check(val)
        x = PyLong_AsDouble(val);
    else if PyInt_Check(val)
        x = PyInt_AsLong(val);
    else if PyFloat_Check(val)
        x = PyFloat_AsDouble(val);
    else {
        PyErr_SetString(PyExc_ValueError, "left must be number");
        return NULL;
    }
    
    data = vector_dataptr(v);
    len = vector_length(v);
    
    out = vector_new(len);
    m_add_scalar(-x, data, vector_dataptr(out), 1, len);
    Py_INCREF(out);
    
    return (PyObject*)out;
}



/*
 * return multiplied vector (Vector*[Int|Float|Vector]) 
 */
PyAPI_FUNC(VectorObject *)
vector_mul(VectorObject *self, VectorObject *other)
{
    VectorObject *out;
    Float *data_self, *data_other;
    int len_self = 0, len_other = 0, i;
    
    DEBUG("vector_mul\n");
    
    data_self = vector_dataptr(self);
    len_self = vector_length(self);
    
    data_other = vector_dataptr(other);
    len_other = vector_length(other);
    
    if (len_self == 1) { // first vector is scalar
        DEBUG("  first vector is scalar\n");
        out = vector_new(len_other);
        m_copy(out->data, data_other, len_other, 1);
        m_scale(*(data_self), out->data, len_other, 1);
    } else if (len_other == 1) { // second vector is scalar
        DEBUG("  second vector is scalar\n");
        out = vector_new(len_self);
        m_copy(out->data, data_self, len_self, 1);
        m_scale(*(data_other), out->data, len_self, 1);
    } else if (len_self == len_other) {
        out = vector_new(len_self);
        DEBUG("  same length vectors\n");
        for (i=0; i<len_self; i++) {
            *(out->data + i) = *(data_self + i) * *(data_other + i);
        }
    } else {
        PyErr_SetObject(PyExc_ValueError, PyString_FromString("Vectors must be the same length"));
        return NULL;
    }
    
    Py_INCREF(out);
    return out;
}



/*
 * return absolute values of vector items
 */
PyAPI_FUNC(VectorObject *)
vector_absolute(VectorObject *self)
{
    VectorObject *out;
    Float *data, *dst;
    int len=0, i;
    
    DEBUG("vector_abs\n");
    
    data = vector_dataptr(self);
    len = vector_length(self);
    
    out = vector_new(len);
    
    dst = vector_dataptr(out);
    i = 0;
    while(dst < vector_dataptr(out) + len) {
        *(dst++) = Fabs(*(data++));
    }
    
    Py_INCREF(out);
    return out;
}



/*
 * return vector item
 */
PyAPI_FUNC(PyObject *)
vector_item(VectorObject *a, int i)
{
    PyObject *out;
    int len=0;
    Float *data=NULL;
    
    data = vector_dataptr(a);
    len = vector_length(a);
    DEBUG("vector item\n    len=%d\n", len);
    
    if (i < 0 || i >= len) {
        PyErr_SetObject(PyExc_IndexError, PyString_FromString("index out of range"));
        return NULL;
    }
    
    out = (PyObject*)PyFloat_FromDouble(*(data + i));
    return out;
}



/*
 * set value to matrix row item
 */
PyAPI_FUNC(int)
vector_ass_item(VectorObject *self, int idx, PyObject *val)
{
    Float x;
    int len;
    Float *p_data=NULL;
    
    DEBUG("\nass item index %d\n", idx);
    len = vector_length(self);
    DEBUG("   length %d\n", len);
    
    if (idx < 0 || idx >= len) {
        PyErr_SetObject(PyExc_IndexError, PyString_FromString("index out of range"));
        return -1;
    }
    
    if PyLong_Check(val)
        x = PyLong_AsDouble(val);
    else if PyInt_Check(val)
        x = PyInt_AsLong(val);
    else if PyFloat_Check(val)
        x = PyFloat_AsDouble(val);
    else {
        PyErr_SetString(PyExc_ValueError, "item must be number");
        return -1;
    }
    
    p_data = vector_dataptr(self);
    *(p_data + idx) = x;
    
    return 0;
}



/*
 * comparison of two vectors
 * http://docs.python.org/extending/newtypes.html#object-comparison
 */
PyAPI_FUNC(int)
vector_cmp(VectorObject *self, VectorObject *other)
{
    Float *s, *o, sub, max;
    Float eps = 1e-8; // TODO: hope this will be enough
    long i, length, result = 0;
    DEBUG("vector_cmp\n");
    
    if (vector_length(self) > vector_length(other)) {
        result = 1;
    } else if (vector_length(self) < vector_length(other)) {
        result = -1;
    } else { // same length
        length = vector_length(self);
        s = vector_dataptr(self);
        o = vector_dataptr(other);
        DEBUG(" self length: %d\n", length);
        
        for (i=0; i<length; i++) {
            DEBUG(" compare step %d: %g - %g = %g\n", i, *s, *o, *s-*o);
            sub = *s - *o;
            max = MAX(*s, *o);
            
            if (max != 0) {
                sub = sub / max;
            }
            
            s++;
            o++;
            
            if (Fabs(sub) > eps) {
                DEBUG("  FABS");
                if (sub > 0) {
                    DEBUG(" 1\n");
                    result = 1;
                } else {
                    DEBUG(" -1\n");
                    result = -1;
                }
            }
        }
    }

    PyObject *err_type, *err_value, *err_traceback;

    if (PyErr_Occurred()) {
        PyErr_Fetch(&err_type, &err_value, &err_traceback);
        DEBUG(" ERROR!!!! %s\n", PyString_AsString(err_value));
    }
    
    return result;
}



/*
 * make an one-item Vector from the scalar
 */
VectorObject *
PyFloat2VectorObject(const Float value)
{
    VectorObject *out;
    
    out = vector_new(1);
    *(vector_dataptr(out)) = value;
    
    return out;
}



/*
 * to decide, if we can make arithmetic operation on current datatypes
 */
int
vector_coerce(PyObject **v, PyObject **w)
{
    DEBUG("Coercing \n");
    
    if (PyInt_Check(*w)) {
        DEBUG(" - vector and Int number\n");
        *w = (PyObject*)PyFloat2VectorObject((Float)PyInt_AsLong(*w));
        Py_INCREF(*v);
        return 0;
    }
    
    if (PyFloat_Check(*w)) {
        DEBUG(" - vector and Float number\n");
        *w = (PyObject*)PyFloat2VectorObject((Float)PyFloat_AsDouble(*w));
        Py_INCREF(*v);
        return 0;
    }
    
    DEBUG("conversion not possible\n");
    
    return 1; // conversion not possible
}



/*
 * returns ref (not copy) to current vector object
 */
PyAPI_FUNC(PyObject *)
vector_slice(VectorObject *self, int ilow, int ihigh)
{
    VectorObject *np;
    Float *data;
    int len;
    
    if (ilow < 0)
            ilow = 0;
    else if (ilow > vector_length(self))
            ilow = vector_length(self);
    if (ihigh < ilow)
            ihigh = ilow;
    else if (ihigh > vector_length(self))
            ihigh = vector_length(self);
    
    len = ihigh - ilow;
    DEBUG("vector slice\n");
    DEBUG("   len: %d\n", len);
    DEBUG("   ilow: %d\n", ilow);
    
    data = vector_dataptr(self);
    //np = vector_new((PyObject*)self, data + ilow);
    np = vector_new(len);
    m_copy(vector_dataptr(np), data + ilow, 1, len);
    
    if (np == NULL)
        return NULL;
    
    Py_INCREF(np);
    return (PyObject*)np;
}



/*
 * create new VectorObject from list (tuple)
 */
PyAPI_FUNC(PyObject *)
VectorObject_New(PyTypeObject *type, PyObject *args)
{
    VectorObject *self;
    Py_ssize_t length=0;
    PyObject *listObject;
    
    if (!PyArg_ParseTuple(args, "O", &listObject)) {
        PyErr_BadArgument();
        return NULL;
    }
    
    // [0]
    if (PyList_Check(listObject)) {
        length = PyList_Size(listObject);
    } else if (PyTuple_Check(listObject)) {
        length = PyTuple_Size(listObject);
    } else if (listObject->ob_type == &VectorType) {
        Py_INCREF(listObject);
        return listObject;
    } else {
        PyErr_BadArgument();
        return NULL;
    }
    
    if (length == 0) {
        return Py_None;
    }
    
    self = vector_new(length);
    if (self == NULL) {
        PyErr_NoMemory();
        return NULL;
    }
    
    if (!PyArg_GetDoubleArray(listObject, 1, 0, length, vector_dataptr(self))) {
        PyErr_BadArgument();
        return NULL;
    }
    
    Py_INCREF(self);
    return (PyObject*)self;
}



/*
 * returns a range 
 */
PyAPI_FUNC(VectorObject *)
vector_range(PyObject *self, PyObject *args, PyObject *kws)
{
    Float start=0, stop=0, step=1;
    int size, i;
    PyObject *o_start = NULL, *o_stop = NULL, *o_step = NULL;
    static char *kwd[]= {"start", "stop", "step", NULL};
    VectorObject *result;
    DEBUG("   vector_range\n");
    
    if(!PyArg_ParseTupleAndKeywords(args, kws, "O|OO", kwd, &o_start, &o_stop, &o_step)) {
        return NULL;
    }
    
    start = PyFloat_AsDouble(o_start);
    DEBUG("   start: %f, stop:%f, step:%f\n", start, stop, step);
    
    if (o_step && o_step != Py_None) {
        step = PyFloat_AsDouble(o_step);
    }
    
    if (o_stop && o_stop != Py_None) {
        stop = PyFloat_AsDouble(o_stop);
    } else {
        stop = start;
        start = 0;
    }
    
    DEBUG("   start: %f, stop:%f, step:%f\n", start, stop, step);
    
    size = ceil((stop - start) / step);
    DEBUG("   size: %d (%f)\n", size, (stop - start) / step);
    
    if (size < 1)
        return (VectorObject*)Py_None;
    
    result = vector_new(size);
    
    if (!result)
        return (VectorObject*)Py_None;
    
    for (i=0; i<size; i++) {
        result->data[i] = start + i*step;
        DEBUG("    LOOP: #%d %f\n", i, result->data[i]);
    }
    
    DEBUG("   result size: %d\n", vector_length(result));
    
    Py_INCREF(result);
    return result;
}



PyAPI_DATA(PyNumberMethods) vector_as_number = {
    (binaryfunc)vector_add,          /*nb_add*/
    (binaryfunc)vector_sub,          /*nb_subtract*/
    (binaryfunc)vector_mul,          /*nb_multiply*/
    0,//(binaryfunc)array_divide,    /*nb_divide*/
    0,//(binaryfunc)array_remainder, /*nb_remainder*/
    0,//(binaryfunc)array_divmod,    /*nb_divmod*/
    0,//(ternaryfunc)array_power,    /*nb_power*/
    0,//(unaryfunc)_negative,        /*nb_neg*/
    0,//(unaryfunc)_array_copy_nice, /*nb_pos*/
    (unaryfunc)vector_absolute,      /*(unaryfunc)array_abs,*/
    0,//(inquiry)_array_nonzero,     /*nb_nonzero*/
    0,//(unaryfunc)array_invert,     /*nb_invert*/
    0,//proxy_lshift,               /* nb_lshift */
    0,//proxy_rshift,               /* nb_rshift */
    0,//proxy_and,              /* nb_and */
    0,//proxy_xor,              /* nb_xor */
    0,//proxy_or,               /* nb_or */
    // more on coercing: http://www.ragestorm.net/tutorials/25/pyextnum.html
    (coercion)vector_coerce,               /* nb_coerce */
};



PyAPI_DATA(PySequenceMethods) vector_as_sequence = {
    (lenfunc)vector_length,             /*sq_length*/
    0,//(binaryfunc)vector_add,         /*sq_concat - numerical add */
    0,//(intargfunc)vector_mul,         /*sq_repeat -  numerical multiply */
    (ssizeargfunc)vector_item,          /*sq_item*/
    (ssizessizeargfunc)vector_slice,    /*sq_slice*/
    (ssizeobjargproc)vector_ass_item,   /*sq_ass_item*/
    0,//(intintobjargproc)row_ass_slice,/*sq_ass_slice*/
};



PyAPI_DATA(PyTypeObject) VectorType = {
    PyObject_HEAD_INIT(NULL)
    0,                      /*ob_size*/
    "Vector",               /*tp_name*/
    sizeof(VectorObject),   /*tp_basicsize*/
    0,                      /*tp_itemsize*/
    (destructor)vector_dealloc,     /* tp_dealloc */
    0,//(printfunc)row_print,            /* tp_print */
    0,//(getattrfunc)row_getattr,        /* tp_getattr */
    0,                  /* tp_setattr */
    (cmpfunc)vector_cmp,  /* tp_compare */
    //http://docs.python.org/extending/newtypes.html#object-presentation
    (reprfunc)vector_repr,/* tp_repr */
    &vector_as_number,  /* tp_as _number*/
    &vector_as_sequence,/* tp_as _sequence*/
    0,                  /* tp_as _mapping*/
    0,                  /* tp_hash */
    0,                  /* tp_call */
    0,                  /* tp_str */
    (getattrofunc)0,                  /* tp_getattro */
    (setattrofunc)0,                  /* tp_setattro */
    0,//&row_as_buffer,          /* tp_as_buffer*/
    Py_TPFLAGS_DEFAULT, /* tp_flags */
    0,//rowtype_doc,             /* tp_doc */
    0,                  /* tp_traverse */
    0,                  /* tp_clear */
    0,//(richcmpfunc)vector_cmp//row_richcompare,         /* tp_richcompare */
};

