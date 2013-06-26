/* -*- mode: C; tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 4 -*-
 *      vim: sw=4 ts=8 et tw=80
 */

/* lookup.c
 *
 * This file contains rountines to help lookup time-ordered data such as
 *  APEX weather readings.
 */

#include "pyactpol.h"
#include "myassert.h"

/* #include <stdbool.h> */


PyDoc_STRVAR( get_nearest__doc__,
"get_nearest(data_array, time)\n"
"\n"
"The data_array must be a numpy 2d array of float64.  The first row\n"
"of the matrix must be strictly increasing.  This function performs\n"
"the equivalent of\n"
"\n"
"       nearest_index = np.argmin(abs(data[:,0]-time))\n"
"       dt = data[nearest_index,0] - time\n"
"       vals = data[nearest_index,1:]\n"
"       return np.hstack((vals, dt))\n"
    );


static PyObject *get_nearest(PyObject *self, PyObject *args)
{
    PyArrayObject *data_ar;
    double at_ctime;

    if (!PyArg_ParseTuple(args, "O!d",
                          &PyArray_Type, &data_ar,
                          &at_ctime
            ))
        po_raise("invalid arguments");

    ASSERT_CARRAY_TYPE_NDIM(data_ar, NPY_FLOAT64, 2);
    int n_rows = PyArray_DIMS(data_ar)[1];
    double *times = PyArray_DATA(data_ar);
    double *data = times + n_rows;

    int n_vals = PyArray_DIMS(data_ar)[0] - 1;

    // Trivial cases
    int out_index;
    if (at_ctime < times[0]) {
        out_index = 0;
    } else if (at_ctime > times[n_rows-1]) {
        out_index = n_rows -1 ;
    } else {
        // Binary search.
        int left = 0, right = n_rows - 1;
        while (right - left > 1) {
            int c = (left + right) / 2;
            if (times[c] > at_ctime)
                right = c;
            else
                left = c;
        }
        if (at_ctime - times[left] < times[right] - at_ctime)
            out_index = left;
        else
            out_index = right;
    }

    PyTupleObject *out = (PyTupleObject*)PyTuple_New(n_vals + 1);
    for (int i=0; i<n_vals; i++)
        PyTuple_SET_ITEM(out, i, PyFloat_FromDouble(data[i*n_rows + out_index]));
    PyTuple_SET_ITEM(out, n_vals, PyFloat_FromDouble(times[out_index] - at_ctime));
    
    return (PyObject*)out;
}


PyMethodDef pyactpol_lookup_methods[] = {
    {"get_nearest", get_nearest, METH_VARARGS, get_nearest__doc__},
    {NULL, NULL, 0, NULL}        /* Sentinel */
};


