/* -*- mode: C; tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 4 -*-
 *      vim: sw=4 ts=8 et tw=80
 */

#include "pyactpol.h"

/*
 * Some routines for dealing with MCE data quickly.  E.g. extracting
 * data fields from mixed mode detector data, reconstructing windowed
 * data, etc.
 */

/*
 * extract_dm<mode>_<field>
 *
 * Input:
 *   src  (array of shape (n) and dtype int32
 *   dest (array of shape (n) and dtype float32)
 *
 * Returns:
 *   True
 */

static PyObject *extract_dm10_filt(PyObject *self, PyObject *args)
{
    PyArrayObject *src_array;
    PyArrayObject *dest_array;
    if (!PyArg_ParseTuple(args, "O!O!",
                          &PyArray_Type, &src_array,
                          &PyArray_Type, &dest_array))
        Py_RETURN_NONE;

    long n = PyArray_DIMS(src_array)[0];
    float *dest = PyArray_DATA(dest_array);
    int *src = PyArray_DATA(src_array);
    int i;

    for (i=0; i<n; i++)
        dest[i] = (float)(src[i] >> 7) * 8;

    Py_RETURN_TRUE;
}


PyMethodDef pyactpol_mce_methods[] = {
    {"extract_dm10_filt", extract_dm10_filt, METH_VARARGS,
     ""},
    {NULL, NULL, 0, NULL}        /* Sentinel */
};
