#pragma once

#include <Python.h>

/* Numpy stuff, to traverse multiple source files. */

#define PY_ARRAY_UNIQUE_SYMBOL PYACTPOL_Array_API
#ifndef PYACTPOL_MAIN /* Only defined by the module that calls import_array */
#define NO_IMPORT_ARRAY /* Disable import_array warnings from numpy */
#endif
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <numpy/arrayobject.h>

#include <stddef.h>

#include "myassert.h"

extern PyMethodDef pyactpol_cuts_methods[];
extern PyMethodDef pyactpol_dirfile_methods[];
extern PyMethodDef pyactpol_fft_methods[];
extern PyMethodDef pyactpol_filter_methods[];
extern PyMethodDef pyactpol_linalg_methods[];
extern PyMethodDef pyactpol_lookup_methods[];
extern PyMethodDef pyactpol_mce_methods[];
extern PyMethodDef pyactpol_pointing_methods[];
extern PyMethodDef pyactpol_sync_methods[];
extern PyMethodDef pyactpol_waterfall_methods[];

void pyactpol_dirfile_init(void);


typedef struct {
    PyObject_HEAD
    void *p;
} pyactpol_ptrobj;

int pyactpol_ptrobj_decode(PyObject *o, void **dest);
pyactpol_ptrobj *pyactpol_ptrobj_new(void *item);

int pyactpol_print_error(const char *format, ... );

int *pyactpol_cuts_vector_from_mask(char *mask_array, int n_mask, int *n_cuts);


struct cuts_set_data {
    int n;
    int *ncut;
    int **cuts;
};

int pyactpol_unpack_cuts_data(struct cuts_set_data *cuts_set,
                              PyObject *pycutsset);
void pyactpol_destroy_cuts_data(struct cuts_set_data *cuts_set);

/* This stats_float_select is just gsl_stats_float_select, first
 * included in GSL 2.5.  Since this is a fairly recent addition to the
 * GSL, we'll include it explicitly for now. */

float moby2_stats_float_select(float *arr, 
                               const size_t stride,
                               const size_t n,
                               const size_t k);
    
static float pyactpol_sselect(unsigned long k, unsigned long n, float *arr)
{
    return moby2_stats_float_select(arr, 1, n, k);
}

static float compute_median(unsigned long n, float *arr)
{
    return moby2_stats_float_select(arr, 1, n, n/2);
}


/* Python 2/3 compatibility */
#if PY_MAJOR_VERSION >= 3
#  define PyString_AsString PyUnicode_AsUTF8
#  define PyString_Check    PyUnicode_Check
#  define PyInt_AsLong      PyLong_AsLong
#  define PyInt_Check       PyLong_Check
#  define PyInt_AS_LONG     PyLong_AS_LONG
#endif
