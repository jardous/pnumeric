/*
  $Id: pnumeric.c,v 1.1 2007/01/21 01:38:35 jp Exp $

  pnumeric.c
     This is pnumeric module for Python

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
#include "kf.h"    // Kalman filter
#include "fft.h"   // Fast Fourier Transform
#include "window.h"
//#include "hpspectrum.h"

#include "pnumeric.h"



/*
 * Kalman filter process
 */
static PyObject *
kf_process(PyObject *self, PyObject *args, PyObject *kws)
{
    MatrixObject *A, *B, *C, *D, *x0, *P0, *Q, *R, *x, *P;
    Float *x_est, *y_est, *P_est;
    PyObject *out, *tmp, *tmp1, *x_est_out, *y_est_out, *P_est_out, *result;
    PyObject *mupdate_callback = NULL, *arglist;
    int i, j, k, n, p, q, datalength;
    MatrixObject *y, *u;

    static char *kwlist[] = {"A", "B", "C", "D", "y", "u", "x0", "P0", "Q", "R", "mupdate_callback", NULL};

    if (!PyArg_ParseTupleAndKeywords(args, kws, "O!O!O!O!O!O!O!O!O!O!|O:set_callback", kwlist,
            &MatrixType, &A,
            &MatrixType, &B,
            &MatrixType, &C,
            &MatrixType, &D,
            &MatrixType, &y,
            &MatrixType, &u,
            &MatrixType, &x0,
            &MatrixType, &P0,
            &MatrixType, &Q,
            &MatrixType, &R,
            &mupdate_callback))
        return NULL;

    if (A->rows != A->cols) {
        PyErr_SetString(PyExc_ValueError, "A must be a square matrix");
        return NULL;
    }

    n = A->rows;
    p = B->cols;
    q = C->rows;

    if (B->rows != n) {
        PyErr_SetString(PyExc_ValueError, "B must be Nxp matrix");
        return NULL;
    }
    if (C->cols != n) {
        PyErr_SetString(PyExc_ValueError, "C must be qxN matrix");
        return NULL;
    }
    if (D->rows != p || D->cols != q) {
        PyErr_SetString(PyExc_ValueError, "D must be pxq matrix");
        return NULL;
    }
    datalength = y->cols;
    if (y->cols != u->cols) {
        PyErr_SetString(PyExc_ValueError, "y and u data lengths does not match");
        return NULL;
    }
    if (y->rows != q) {
        PyErr_SetString(PyExc_ValueError, "y must be qxlength matrix");
        return NULL;
    }
    if (u->rows != p) {
        PyErr_SetString(PyExc_ValueError, "u must be pxlength matrix");
        return NULL;
    }
    if (x0->rows != n || x0->cols != 1) {
        PyErr_SetString(PyExc_ValueError, "x0 must be Nx1 matrix");
        return NULL;
    }
    if (P0->rows != n || P0->cols != n) {
        PyErr_SetString(PyExc_ValueError, "P0 must be NxN matrix");
        return NULL;
    }
    if (Q->rows != n || Q->cols != n) {
        PyErr_SetString(PyExc_ValueError, "Q must be NxN matrix");
        return NULL;
    }
    if (R->rows != q || R->cols != q) {
        PyErr_SetString(PyExc_ValueError, "R must be qxq matrix");
        return NULL;
    }
    if (!PyCallable_Check(mupdate_callback)) {
        if (mupdate_callback != NULL) {
            PyErr_SetString(PyExc_TypeError, "parameter must be callable");
            return NULL;
        }
    }

    x_est = m_new(n, datalength);
    y_est = m_new(n, datalength);
    P_est = m_new(n, n*datalength);

    if (mupdate_callback != NULL) {
        Py_XINCREF(mupdate_callback);  // Add a reference to new callback
        x = matrix_new(n, 1);
        P = matrix_new(n, n);
        // initialize x and P
        m_copy(x->data, x0->data, n, 1);
        m_copy(P->data, P0->data, n, n);

        for (i=0; i<datalength; i++) {
            tick(A->data, B->data, C->data, D->data,
                n, p, q,
                y->data+i*q, u->data+i*p,
                x->data, P->data,
                Q->data, R->data,
                x_est+i*n, y_est+i*q, P_est+i*n*n);
            // copy x and P values into output array
            m_copy(x->data, x_est+i*n, n, 1);
            m_copy(P->data, P_est+i*n*n, n, n);
            // update the matrixes
            arglist = Py_BuildValue("(i, O, O, O, O, O)", i, A, B, C, D, x);
            result = PyEval_CallObject(mupdate_callback, arglist);
            Py_DECREF(arglist);
        }
        Py_DECREF(x);
        Py_DECREF(P);
    } else { // model update not needed
        process(A->data, B->data, C->data, D->data,
                n, p, q, y->data, u->data, x0->data, P0->data, Q->data, R->data, datalength,
                x_est, y_est, P_est);
    }

    // create lists from matrixes
    x_est_out = PyList_New(datalength);
    for (i=0; i<datalength; i++) {
        tmp = PyList_New(n);
        for (j=0; j<n; j++) {
            PyList_SetItem(tmp, j, PyFloat_FromDouble( *(x_est+i*n+j) ));
        }
        PyList_SetItem(x_est_out, i, tmp);
    }
    y_est_out = PyList_New(datalength);
    for (i=0; i<datalength; i++) {
        tmp = PyFloat_FromDouble( *(y_est+i) );
        PyList_SetItem(y_est_out, i, tmp);
    }
    P_est_out = PyList_New(datalength);
    for (i=0; i<datalength; i++) {
        tmp = PyList_New(n);
        for (j=0; j<n; j++) {
            tmp1 = PyList_New(n);
            for (k=0; k<n; k++) {
                PyList_SetItem(tmp1, k, PyFloat_FromDouble( *(x_est+i*n+j*n+k) ));
            }
            PyList_SetItem(tmp, j, tmp1);
        }
        tmp = PyFloat_FromDouble( *(P_est+i) );
        PyList_SetItem(P_est_out, i, tmp);
    }

    m_free(x_est);
    m_free(y_est);
    m_free(P_est);

    out = PyTuple_New(3);
    PyTuple_SetItem(out, 0, x_est_out);
    PyTuple_SetItem(out, 1, y_est_out);
    PyTuple_SetItem(out, 2, P_est_out);
    Py_INCREF(out); // TODO add Py_INCREF every time when returning object from fnc

    return out;
}



/*
 * FFT
 */
static PyObject *
fft_process(PyObject *self, VectorObject *v)
{
    VectorObject *out;

    if (v->ob_type == &VectorType) {
        DEBUG("length %d\n", v->length);
        out = vector_new(v->length);
        m_copy(out->data, v->data, 1, v->length);
        FFT(out->data, v->length);
        return (PyObject*)out;
    } else {
        PyErr_SetString(PyExc_ValueError, "argument must be pnumeric.Vector type");
        return NULL;
    }
}



/*
 * hpspectrum
 */
/*
static PyObject *
py_hpspectrum(PyObject *self, PyObject *args)
{
    VectorObject *out, *v;
    long N, R;

    if (!PyArg_ParseTuple(args, "Oii", &v, &N, &R)) {
        PyErr_SetString(PyExc_ValueError, "argument must be a positive Integer");
        return NULL;
    }

    if (v->ob_type != &VectorType) {
        PyErr_SetString(PyExc_ValueError, "argument must be pnumeric.Vector type");
        return NULL;
    }

    DEBUG("py_hpspectrum:\n");

    DEBUG("    length %d\n", v->length);

    out = vector_new(0);

    //m_copy(out->data, v->data, 1, v->length);
    int K = 0;
    hpspectrum(out->data, &K, v->data, v->length, N, R);
    //printf("   returned length K=%d\n", K);
    out->length = K;

    Py_INCREF(out);
    return (PyObject*)out;
}
*/


/*
 * compute Root Mean Square
*/
Float
rms(int length, Float *array)
{
    int i;
    Float res, X=0, mean=0;

    for(i=0; i<length; i++) {
        X = X + pow(*array, 2);
        mean = mean + *array;
        array++;
    }

    mean = pow(mean, 2) / length;

    res = (X-mean) / length;
    res = sqrt(res);
    return res;
}



/*
 * RMS
 */
static PyObject *
py_rms(PyObject *self, VectorObject *v)
{
    Float res, *data;
    PyObject *out;
    int length;

    DEBUG("\npy_rms:\n");

    if (v->ob_type != &VectorType) {
        PyErr_SetString(PyExc_ValueError, "argument must be pnumeric.Vector type");
        return NULL;
    }

    data = vector_dataptr(v);
    length = vector_length(v);
    DEBUG("   length %d\n", length);

    res = rms(length, data);
    DEBUG("   res %g\n", res);

    out = Py_BuildValue("f", res);
    return out;
}



/*
 * mean value
 */
static PyObject *
py_mean(PyObject *self, VectorObject *v)
{
    Float res=0, *data;
    PyObject *out;
    int length, i;

    if (v->ob_type != &VectorType) {
        PyErr_SetString(PyExc_ValueError, "argument must be pnumeric.Vector type");
        return NULL;
    }

    data = vector_dataptr(v);
    length = vector_length(v);

    for(i=0; i<length; i++) {
        res = res + *data++;
    }

    res = res / length;
    out = Py_BuildValue("f", res);
    return out;
}



/*
 * Rectangular window
 */
static VectorObject *
py_rect(PyObject *self, PyObject *args)
{
    int N;
    VectorObject *out;

    if (!PyArg_ParseTuple(args, "i", &N) || N < 1) {
        PyErr_SetString(PyExc_ValueError, "argument must be a positive Integer");
        return NULL;
    }

    out = vector_new(N);
    rect(N, out->data);

    return out;
}



/*
 * Hanning window
 */
static VectorObject *
py_hanning(PyObject *self, PyObject *args)
{
    int N;
    VectorObject *out;

    if (!PyArg_ParseTuple(args, "i", &N) || N < 1) {
        PyErr_SetString(PyExc_ValueError, "argument must be a positive Integer");
        return NULL;
    }

    out = vector_new(N);
    hanning(N, out->data);

    return out;
}



/*
 * Hamming window
 */
static VectorObject *
py_hamming(PyObject *self, PyObject *args)
{
    int N;
    VectorObject *out;

    if (!PyArg_ParseTuple(args, "i", &N) || N < 1) {
        PyErr_SetString(PyExc_ValueError, "argument must be a positive Integer");
        return NULL;
    }

    out = vector_new(N);
    hamming(N, out->data);

    return out;
}



static PyMethodDef pnumeric_methods[] = {
    /*
     * The cast of the function is necessary since PyCFunction values
     * only take two PyObject* parameters, and keywdarg_parrot() takes
     * three.
     */
    {"zeros",   (PyCFunction)matrix_zeros, METH_VARARGS, "returns zeros matrix"},
    {"ones",   (PyCFunction)matrix_ones, METH_VARARGS, "returns ones matrix"},
    {"eye",   (PyCFunction)matrix_eye, METH_VARARGS, "returns eye matrix"},
    {"vrange", (PyCFunction)vector_range, METH_VARARGS | METH_KEYWORDS, "generate a range vector"},
    {"kf_process", (PyCFunction)kf_process, METH_VARARGS | METH_KEYWORDS, "Kalman filter"},
    //{"resample", (PyCFunction)py_resample, METH_VARARGS, "Resample a vector decimation or interpolation"},
    {"fft", (PyCFunction)fft_process, METH_O, "Fast Fourier Transform"},
    //{"hpspectrum", (PyCFunction)py_hpspectrum, METH_VARARGS, "Harmonic product spectrum of a vector"},
    {"rms", (PyCFunction)py_rms, METH_O, "Root Mean Square"},
    {"mean", (PyCFunction)py_mean, METH_O, "Mean value"},
    {"rect", (PyCFunction)py_rect, METH_VARARGS, "Rectangular window"},
    {"hanning", (PyCFunction)py_hanning, METH_VARARGS, "Hanning window"},
    {"hann", (PyCFunction)py_hanning, METH_VARARGS, "Hanning window"},
    {"hamming", (PyCFunction)py_hamming, METH_VARARGS, "Hamming window"},
    {NULL, NULL, 0, NULL}   /* sentinel */
};



DL_EXPORT(void)
initpnumeric(void)
{
    PyObject *m;

    //VectorType.ob_type = &PyType_Type;
    VectorType.tp_new = VectorObject_New;
    if (PyType_Ready(&VectorType) < 0)
        return;

    MatrixType.tp_new = MatrixObject_New;
    if (PyType_Ready(&MatrixType) < 0)
        return;

    // Create the module and add the functions
    m = Py_InitModule("pnumeric", pnumeric_methods);

    Py_INCREF(&VectorType);
    Py_INCREF(&MatrixType);
    PyModule_AddObject(m, "Vector", (PyObject *)&VectorType);
    PyModule_AddObject(m, "Matrix", (PyObject *)&MatrixType);
}

