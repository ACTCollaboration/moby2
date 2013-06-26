/* -*- mode: C; tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 4 -*-
 *      vim: sw=4 ts=8 et tw=80
 */

/* linalg.c
 *
 * This file contain some simple accelerated array processing
 * routines.  Like means and medians and dot products and stuff.
 *
 */


#include "pyactpol.h"
#include "myassert.h"

#include <stdbool.h>

static PyObject *data_median_axis0(PyObject *self, PyObject *args)
{
//    mbTOD *tod, const int *dets, int ndet, double *output, int ndata

    PyArrayObject *tod_array;
    PyArrayObject *det_idx_array;

    if (!PyArg_ParseTuple(args, "O!O!",
                          &PyArray_Type, &tod_array,
                          &PyArray_Type, &det_idx_array
            ))
        po_raise("invalid arguments");

    ASSERT_CARRAY_TYPE_NDIM(tod_array, NPY_FLOAT32, 2);
    int ndata = PyArray_DIMS(tod_array)[1];
    float *data = PyArray_DATA(tod_array);
    
    ASSERT_CARRAY_TYPE_NDIM(det_idx_array, NPY_INT32, 1);
    int ndet = PyArray_DIMS(det_idx_array)[0];
    int *dets = PyArray_DATA(det_idx_array);

    PyObject *output_array =
        PyArray_SimpleNew(1, PyArray_DIMS(tod_array)+1, NPY_FLOAT64);
    double *output = PyArray_DATA((PyArrayObject*)output_array);

#pragma omp parallel shared(ndet, ndata, data, dets, output) default(none)
    {
        float *mymedian=(float *)malloc(ndet*sizeof(*mymedian));
        for (int i = 0; i < ndet; i++) mymedian[i]=0;

#pragma omp for
        for (int i = 0; i < ndata; i++)  {
            for (int j = 0; j < ndet; j++)
                mymedian[j] = data[dets[j]*ndata + i];
            output[i] = (double)compute_median(ndet, mymedian);
        }
        free(mymedian);
    }

    return (PyObject*)output_array;
}


static PyObject *data_norm_axis1(PyObject *self, PyObject *args)
{
    PyArrayObject *tod_array;
    if (!PyArg_ParseTuple(args, "O!",
                          &PyArray_Type, &tod_array
            ))
        po_raise("invalid arguments");

    ASSERT_CARRAY_TYPE_NDIM(tod_array, NPY_FLOAT32, 2);
    const int ndet = PyArray_DIMS(tod_array)[0];
    const int ndata = PyArray_DIMS(tod_array)[1];
    float *data = PyArray_DATA(tod_array);
    
    PyObject *output_array =
        PyArray_SimpleNew(1, PyArray_DIMS(tod_array), NPY_FLOAT64);
    double *output = PyArray_DATA((PyArrayObject*)output_array);

#pragma omp parallel for
    for (int j = 0; j < ndet; j++) {
        output[j] = 0.;
        float *det_data = data + j*ndata;
        for (int i = 0; i < ndata; i++)
            output[j] += det_data[i]*det_data[i];
        output[j] = sqrt(output[j] / ndata);
    }

    return (PyObject*)output_array;
}

/* data_moments_chunks
 *
 * Compute arbitrary raw moments on chunks of an array.  Returns
 * 3d array; (n_moment, n_det, n_chunk).
 */

static PyObject *data_moments_chunks(PyObject *self, PyObject *args)
/* void mbPseudoMomentsTODChunks(mbTOD *tod, */
/*                               double *var, int nvar, int nvch, */
/*                               double *skew, int nskew, int nsch, */
/*                               double *kurt, int nkurt, int nkch, */
/*                               int pivot, */
/*                               int T) */
{
    PyArrayObject *tod_array;
    PyObject *moments_list;
    int offset;
    int chunk_size;

    if (!PyArg_ParseTuple(args, "O!Oii",
                          &PyArray_Type, &tod_array,
                          &moments_list,
                          &offset,
                          &chunk_size
            ))
        po_raise("invalid arguments");

    /* Check TOD data */
    ASSERT_CARRAY_TYPE_NDIM(tod_array, NPY_FLOAT32, 2);
    int ndet = PyArray_DIMS(tod_array)[0];
    int ndata = PyArray_DIMS(tod_array)[1];
    float *data = PyArray_DATA(tod_array);
    int nchunk = (ndata-offset) / chunk_size;

    /* Check the moments list */
    po_assert(PySequence_Check(moments_list));
    int n_moments = PySequence_Length(moments_list);
    int *moments = malloc(n_moments * sizeof(int));
    int last_moment = 0;
    for (int i=0; i<PySequence_Length(moments_list); i++) {
        PyObject* o = PySequence_GetItem(moments_list, i);
        if (!PyInt_Check(o))
            last_moment = -1;
        moments[i] = PyInt_AS_LONG(o);
        if (moments[i] <= last_moment)
            last_moment = -1;
        last_moment = moments[i];
        Py_DECREF(o);
    }
    if (last_moment == -1)
        po_raise("moments list was not increasing sequence of positive ints");

    /* Create output array */
    npy_intp dims[3] = {n_moments, ndet, nchunk};
    PyObject *output_array = PyArray_SimpleNew(
        3, dims, NPY_FLOAT64);
    double *output = PyArray_DATA((PyArrayObject*)output_array);
    
#pragma omp parallel for shared(output, offset, chunk_size, nchunk)
    for (int d = 0; d < ndet; d++) {
        for (int i = 0; i < nchunk; i++) {
            for (int im=0; im<n_moments; im++)
                output[(im*ndet+d)*nchunk + i] = 0;
            for (int j = offset+chunk_size*i; j < offset+chunk_size*(i+1); j++) {
                int last_moment = 1;
                double x = data[d*ndata+j];
                double m = x;
                for (int im=0; im<n_moments; im++) {
                    for (; last_moment < moments[im]; last_moment++)
                        m *= x;
                    output[(im*ndet+d)*nchunk + i] += m;
                }
            }
            for (int im=0; im<n_moments; im++)
                output[(im*ndet+d)*nchunk + i] /= chunk_size;
        }
    }

    /* And out */
    free(moments);
    return output_array;


}



/* data_mean_axis0
 *
 * Input:
 *   data: (ndet,nsamps) time-ordered data to filter (modified in place).
 *   dets: (ndet_used) index of data vectors to actually filter.
 *   
 */

static PyObject *data_mean_axis0(PyObject *self, PyObject *args)
{
    PyArrayObject *tod_array;
    PyArrayObject *det_idx_array;

    if (!PyArg_ParseTuple(args, "O!O!",
                          &PyArray_Type, &tod_array,
                          &PyArray_Type, &det_idx_array
            ))
        po_raise("invalid arguments.");


    // Those other vectors
    po_assert((PyArray_TYPE(tod_array) == NPY_FLOAT32) &&
              (PyArray_NDIM(tod_array) == 2));
    po_assert(PyArray_FLAGS(tod_array) & NPY_ARRAY_CARRAY);

    float *tod = PyArray_DATA(tod_array);
    ASSERT_CARRAY_TYPE_NDIM(det_idx_array, NPY_INT, 1);
    int *det_idx = PyArray_DATA(det_idx_array);

    int ndata = PyArray_DIMS(tod_array)[1];
    int ndets = PyArray_SIZE(det_idx_array);

    npy_intp d = ndata;
    PyArrayObject* output_array = (PyArrayObject*)
        PyArray_SimpleNew(1, &d, NPY_FLOAT32);
    float *output = PyArray_DATA(output_array);

#pragma omp parallel for
    for (int i = 0; i < ndata; i++)  {
        double m = 0.0;
        for (int j = 0; j < ndets; j++)
            m += tod[det_idx[j]*ndata+i];
        output[i] = m/ndets;
    }

    return (PyObject*)output_array;
}


PyDoc_STRVAR(remove_modes__doc__,
	     "remove_modes[64](data, dets, modes, amps)\n"
	     "\n"
	     "The arrays must be C-ordered with dimensions and types:\n"
	     "     data [     *,n_data]   (float32)\n"
	     "     dets [n_det]           (int32)\n"
             "     modes[n_mode,n_data]   (float32[float64])\n"
	     "     amps [n_mode,n_det]    (float64)\n"
	     "\n"
	     "The dets array gives indices into the first dimension of data, "
	     "so the largest value in dets must be smaller than data.shape[0].\n"
	     "\n"
	     "Anyway, this computes \n"
	     "        data[dets,:] -= (amp[:,:,None] * mode[:,None,:]).sum(axis=0)\n"
	     "but probably faster and probably more efficiently.\n"
    );

static PyObject *remove_modes(PyObject *self, PyObject *args)
{
    PyArrayObject *data_array;
    PyArrayObject *dets_array;
    PyArrayObject *modes_array;
    PyArrayObject *amps_array;

    if (!PyArg_ParseTuple(args, "O!O!O!O!",
                          &PyArray_Type, &data_array,
                          &PyArray_Type, &dets_array,
                          &PyArray_Type, &modes_array,
                          &PyArray_Type, &amps_array
            ))
        po_raise("invalid arguments.");

    // The easy one...
    ASSERT_CARRAY_TYPE_NDIM(dets_array, NPY_INT32, 1);

    // Types and ordering... the ordering is stricter than we really need.
    ASSERT_CARRAY_TYPE_NDIM(data_array, NPY_FLOAT32, 2);
    ASSERT_CARRAY_TYPE_NDIM(modes_array, NPY_FLOAT32, 2);
    ASSERT_CARRAY_TYPE_NDIM(amps_array, NPY_FLOAT64, 2);

    int ndata = PyArray_DIMS(data_array)[1];
    int ndet = PyArray_DIMS(dets_array)[0];
    int nmode = PyArray_DIMS(modes_array)[0];

    po_assert(PyArray_DIMS(amps_array)[0] == nmode);
    po_assert(PyArray_DIMS(amps_array)[1] == ndet);
    po_assert(PyArray_DIMS(modes_array)[1] == ndata);

    float *data = PyArray_DATA(data_array);
    float *modes = PyArray_DATA(modes_array);
    double *amps = PyArray_DATA(amps_array);
    int *dets = PyArray_DATA(dets_array);

    // We've been so thorough, it would be a shame to segfault now.
    for (int i=0; i<ndet; i++)
        po_assert(dets[i] >= 0 &&
                  dets[i] < PyArray_DIMS(data_array)[0]);
	 
#pragma omp parallel for
    for (int deti=0; deti<ndet; deti++) {
        float *this_data = data + ndata*dets[deti];
        for (int modei=0; modei<nmode; modei++) {
            float *this_mode = modes + modei*ndata;
            float this_amp = amps[modei*ndet+deti];
            for (int i=0; i<ndata; i++)
                this_data[i] -= this_amp * this_mode[i];
        }
    }

    Py_RETURN_NONE;
}

static PyObject *remove_modes64(PyObject *self, PyObject *args)
{
    PyArrayObject *data_array;
    PyArrayObject *dets_array;
    PyArrayObject *modes_array;
    PyArrayObject *amps_array;

    if (!PyArg_ParseTuple(args, "O!O!O!O!",
                          &PyArray_Type, &data_array,
                          &PyArray_Type, &dets_array,
                          &PyArray_Type, &modes_array,
                          &PyArray_Type, &amps_array
            ))
        po_raise("invalid arguments.");

    // The easy one...
    ASSERT_CARRAY_TYPE_NDIM(dets_array, NPY_INT32, 1);

    // Types and ordering... the ordering is stricter than we really need.
    ASSERT_CARRAY_TYPE_NDIM(data_array, NPY_FLOAT32, 2);
    ASSERT_CARRAY_TYPE_NDIM(modes_array, NPY_FLOAT64, 2);
    ASSERT_CARRAY_TYPE_NDIM(amps_array, NPY_FLOAT64, 2);

    int ndata = PyArray_DIMS(data_array)[1];
    int ndet = PyArray_DIMS(dets_array)[0];
    int nmode = PyArray_DIMS(modes_array)[0];

    po_assert(PyArray_DIMS(amps_array)[0] == nmode);
    po_assert(PyArray_DIMS(amps_array)[1] == ndet);
    po_assert(PyArray_DIMS(modes_array)[1] == ndata);

    float *data = PyArray_DATA(data_array);
    double *modes = PyArray_DATA(modes_array);
    double *amps = PyArray_DATA(amps_array);
    int *dets = PyArray_DATA(dets_array);

    // We've been so thorough, it would be a shame to segfault now.
    for (int i=0; i<ndet; i++)
        po_assert(dets[i] >= 0 &&
                  dets[i] < PyArray_DIMS(data_array)[0]);
	 
#pragma omp parallel
{
    double *tmp = malloc(ndata * sizeof(*tmp));

#pragma omp for
    for (int deti=0; deti<ndet; deti++) {
        memset(tmp, 0, ndata * sizeof(*tmp));
        for (int modei=0; modei<nmode; modei++) {
            // Loop over samples...
            for (int i=0; i<ndata; i++)
                tmp[i] += amps[modei*ndet+deti] * modes[modei*ndata+i];
        }
        float *this_data = data + ndata*dets[deti];
        for (int i=0; i<ndata; i++)
            this_data[i] -= tmp[i];
    }

    free(tmp);
}
/* #pragma omp parallel for */
/*     for (int deti=0; deti<ndet; deti++) { */
/*         float *this_data = data + ndata*dets[deti]; */
/*         // Loop over samples... */
/*         for (int i=0; i<ndata; i++) { */
/*             // And over modes. */
/*             double to_remove = 0; */
/*             /\*     double *this_mode = modes + modei*ndata; *\/ */
/*             /\* double this_amp = amps[modei*ndet+deti]; *\/ */
/*             for (int modei=0; modei<nmode; modei++) { */
/*                 to_remove += amps[modei*ndet+deti] * modes[modei*ndata+i]; */
/*             } */
/*             this_data[i] -= to_remove; */
/*         } */
/*     } */

    Py_RETURN_NONE;
}

/*
 * data_dot_modes replaces mbGetModeRenorm and mbDotProductTOD, which
 * also replaces mbGetModeRenorm.  This is why we're doing this.
 */

PyDoc_STRVAR(data_dot_modes__doc__,
	     "data_dot_modes(data, dets, modes, cuts)\n"
	     "\n"
	     "All arguments are two-dimensional numpy arrays of floats, except "
	     "for dets, which is 1d array of integers.\n"
	     "\n"
	     "The arrays must be C-ordered with dimensions like:\n"
	     "     data [     *,n_data]\n"
	     "     dets [n_det]\n"
             "     modes[n_mode,n_data]\n"
	     "\n"
	     "The dets array gives indices into the first dimension of data, "
	     "so the largest value in dets must be smaller than data.shape[0].\n"
	     "\n"
	     "Anyway, this computes \n"
	     "        data[dets,:] -= (amp[:,:,None] * mode[:,None,:]).sum(axis=0)\n"
	     "but probably faster and probably more efficiently.\n"
    );

static PyObject *data_dot_modes(PyObject *self, PyObject *args)
{
    PyArrayObject *data_array;
    PyArrayObject *dets_array;
    PyArrayObject *modes_array;
    PyObject *cuts_data;

    if (!PyArg_ParseTuple(args, "O!O!O!O",
                          &PyArray_Type, &data_array,
                          &PyArray_Type, &dets_array,
                          &PyArray_Type, &modes_array,
                          &cuts_data
            ))
        po_raise("invalid arguments");

    // The easy one...
    ASSERT_CARRAY_TYPE_NDIM(dets_array, NPY_INT32, 1);

    struct cuts_set_data cuts_set = {0,NULL,NULL};
    if (cuts_data != Py_None) {
        if (!pyactpol_unpack_cuts_data(&cuts_set, cuts_data))
            po_raise("could not unpack cuts data.");
    }

    // Types and ordering... the ordering is stricter than we really need.
    ASSERT_CARRAY_TYPE_NDIM(data_array, NPY_FLOAT32, 2);
    ASSERT_CARRAY_TYPE_NDIM(modes_array, NPY_FLOAT32, 2);

    int nsamps = PyArray_DIMS(data_array)[1];
    int ndet = PyArray_DIMS(dets_array)[0];
    int nmode = PyArray_DIMS(modes_array)[0];

    po_assert(PyArray_DIMS(modes_array)[1] == nsamps);

    float *data = PyArray_DATA(data_array);
    float *modes = PyArray_DATA(modes_array);
    int *dets = PyArray_DATA(dets_array);

    // We've been so thorough, it would be a shame to segfault now.
    for (int i=0; i<ndet; i++)
        po_assert(dets[i] >= 0 &&
                  dets[i] < PyArray_DIMS(data_array)[0]);

    npy_intp dims[2];
    dims[0] = nmode;
    dims[1] = ndet;
    PyArrayObject *coeff_array =
        (PyArrayObject*)PyArray_SimpleNew(2, dims, NPY_FLOAT64);
    double *coeff = PyArray_DATA(coeff_array);
	 
#pragma omp parallel for
    for (int deti=0; deti<ndet; deti++) {
        float *this_data = data + nsamps*dets[deti];
        for (int modei=0; modei<nmode; modei++) {
            float *this_mode = modes + modei*nsamps;
            int ind = modei*ndet + deti;
            coeff[ind] = 0;

            // Loop through uncut regions.
            int j = 0;
            int seg_idx = 0;
            int seg_end;
            while (j < nsamps) {
                if (cuts_data!=Py_None && seg_idx < cuts_set.ncut[deti])
                    seg_end = cuts_set.cuts[deti][seg_idx*2];
                else
                    seg_end = nsamps;

                //Operation
                for (; j<seg_end; j++)
                    coeff[ind] += this_data[j]*this_mode[j];
                
                // Cue up next region
                if (cuts_data!=Py_None && seg_idx < cuts_set.ncut[deti]) {
                    j = cuts_set.cuts[deti][seg_idx*2+1];
                    seg_idx++;
                }
            }                
        }
    }

    if (cuts_data!=Py_None)
        pyactpol_destroy_cuts_data(&cuts_set);

    return (PyObject*)coeff_array;
}


PyDoc_STRVAR(apply_calibration__doc__,
	     "apply_calibration(data, dets, cal)\n"
	     "\n"
	     "The dets array gives indices into the first dimension of data, "
	     "so the largest value in dets must be smaller than data.shape[0].\n"
	     "\n"
	     "Applies\n"
	     "        data[dets,:] *= cal[:,None]\n"
    );

static PyObject *apply_calibration(PyObject *self, PyObject *args)
{
    PyArrayObject *data_array;
    PyArrayObject *dets_array;
    PyArrayObject *amps_array;

    if (!PyArg_ParseTuple(args, "O!O!O!",
                          &PyArray_Type, &data_array,
                          &PyArray_Type, &dets_array,
                          &PyArray_Type, &amps_array
            ))
        po_raise("invalid arguments.");

    // The easy ones...
    ASSERT_CARRAY_TYPE_NDIM(dets_array, NPY_INT32, 1);
    ASSERT_CARRAY_TYPE_NDIM(amps_array, NPY_FLOAT32, 1);

    // Types and ordering... the ordering is stricter than we really need.
    ASSERT_CARRAY_TYPE_NDIM(data_array, NPY_FLOAT32, 2);

    int nsamps = PyArray_DIMS(data_array)[1];
    int ndet = PyArray_DIMS(dets_array)[0];

    float *data = PyArray_DATA(data_array);
    float *amps = PyArray_DATA(amps_array);
    int *dets = PyArray_DATA(dets_array);

    // We've been so thorough, it would be a shame to segfault now.
    for (int i=0; i<ndet; i++)
        po_assert(dets[i] >= 0 &&
                  dets[i] < PyArray_DIMS(data_array)[0]);
	 
#pragma omp parallel for
    for (int deti=0; deti<ndet; deti++) {
        float *this_data = data + nsamps*dets[deti];
        float this_amp = amps[deti];
        for (int i=0; i<nsamps; i++)
            this_data[i] *= this_amp;
    }

    Py_RETURN_NONE;
}


/* bin_data(data, index, binned)

   This function is to assist with histogramming, especially for
   things like power spectra, where the mapping from each detector's
   frequency vector into the binned space is the same for all
   detectors.

   Given the data[n_det, n_freqs] and the bin indices index[n_freqs],
   the function accumulates binned[n_det, n_bins].  Caller is
   responsible to enforce 0 <= index < n_bins.
  */

static PyObject *bin_data(PyObject *self, PyObject *args)
{
    PyArrayObject *index_array;
    PyArrayObject *data_array;
    PyArrayObject *binned_array;

    if (!PyArg_ParseTuple(args, "O!O!O!",
                          &PyArray_Type, &index_array,
                          &PyArray_Type, &data_array,
                          &PyArray_Type, &binned_array
            ))
        po_raise("invalid arguments.");

    ASSERT_CARRAY_TYPE_NDIM(index_array, NPY_INT32, 1);
    ASSERT_CARRAY_TYPE_NDIM(data_array, NPY_FLOAT32, 2);
    ASSERT_CARRAY_TYPE_NDIM(binned_array, NPY_FLOAT64, 2);

    int nsamps = PyArray_DIMS(data_array)[1];
    int ndets = PyArray_DIMS(data_array)[0];
    int nbins = PyArray_DIMS(binned_array)[1];

    float *data = PyArray_DATA(data_array);
    int *index = PyArray_DATA(index_array);
    double *out = PyArray_DATA(binned_array);

#pragma omp parallel for
    for (int deti=0; deti<ndets; deti++) {
        float *d = data + deti*nsamps;
        double *outp = out + deti*nbins;
        for (int i=0; i<nsamps; i++)
            if (index[i]>=0)
                outp[index[i]] += *(d+i);
    }

    Py_RETURN_NONE;
}


static PyObject *bin_data_complex(PyObject *self, PyObject *args)
{
    PyArrayObject *index_array;
    PyArrayObject *data_array;
    PyArrayObject *binned_array;

    if (!PyArg_ParseTuple(args, "O!O!O!",
                          &PyArray_Type, &index_array,
                          &PyArray_Type, &data_array,
                          &PyArray_Type, &binned_array
            ))
        po_raise("invalid arguments.");

    ASSERT_CARRAY_TYPE_NDIM(index_array, NPY_INT32, 1);
    ASSERT_CARRAY_TYPE_NDIM(data_array, NPY_COMPLEX64, 2);
    ASSERT_CARRAY_TYPE_NDIM(binned_array, NPY_COMPLEX128, 2);

    int nsamps = PyArray_DIMS(data_array)[1];
    int ndets = PyArray_DIMS(data_array)[0];
    int nbins = PyArray_DIMS(binned_array)[1];

    float *data = PyArray_DATA(data_array);
    int *index = PyArray_DATA(index_array);
    double *out = PyArray_DATA(binned_array);

#pragma omp parallel for
    for (int deti=0; deti<ndets; deti++) {
        float *d = data + deti*nsamps*2;
        double *outp = out + deti*nbins*2;
        for (int i=0; i<nsamps; i++)
            if (index[i]>=0) {
                outp[index[i]*2  ] += *(d+i*2  );
                outp[index[i]*2+1] += *(d+i*2+1);
            }
    }

    Py_RETURN_NONE;
}


PyMethodDef pyactpol_linalg_methods[] = {
    {"data_mean_axis0", data_mean_axis0, METH_VARARGS, ""},
    {"data_median_axis0", data_median_axis0, METH_VARARGS, ""},
    {"data_norm_axis1", data_norm_axis1, METH_VARARGS, ""},
    {"data_dot_modes", data_dot_modes, METH_VARARGS,
     data_dot_modes__doc__},
    {"remove_modes", remove_modes, METH_VARARGS,
     remove_modes__doc__},
    {"remove_modes64", remove_modes64, METH_VARARGS,
     remove_modes__doc__},
    {"apply_calibration", apply_calibration, METH_VARARGS,
     apply_calibration__doc__},
    {"data_moments_chunks", data_moments_chunks, METH_VARARGS,
     ""}, //data_moments_chunks__doc__},
    {"bin_data", bin_data, METH_VARARGS,
     ""},
    {"bin_data_complex", bin_data_complex, METH_VARARGS,
     ""},
    {NULL, NULL, 0, NULL}        /* Sentinel */
};


