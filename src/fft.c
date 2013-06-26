/* -*- mode: C; tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 4 -*-
 *      vim: sw=4 ts=8 et tw=80
 */

/* fft.c
 *
 * Provide OMP'd FFTW3 wrappers, using the real<->complex transforms
 * (that seem to be twice as fast as numpy, in addition to
 * parallelization).
 */


#include "pyactpol.h"
#include "myassert.h"

#include <stdbool.h>

#include <fftw3.h>

/*
 *  The plan is to name routines like:
 *
 *    [i]fft_T_C
 *
 *  where T is the type (f for float) and C is the complexity (r for real).
 */

PyDoc_STRVAR(fft_f_r__doc__,
	     "fft_f_r(data, dets, full_freq)\n"
	     "\n"
	     "Take 32-bit FFT of data, parallelized.  Arrays must be C-ordered\n"
	     "with dimensions\n"
	     "and types:\n"
	     "     data [     *,n_data]   (float32)\n"
	     "     dets [n_det]           (int64)\n"
	     "\n"
	     "The dets array gives indices into the first dimension of data, "
	     "so the largest value in dets must be smaller than data.shape[0].\n"
	     "\n"
	     "If full_freq!=0, the full FFT is returned (like numpy.fft.fft),\n"
	     "with shape [n_det,n_data] and dtype complex64.\n"
	     "\n"
	     "If full_freq==0, only the first n/2+1 frequency components are\n"
	     "returned.  This saves a bit of space but is not much faster.\n"
    );

static PyObject *fft_f_r(PyObject *self, PyObject *args)
{
    PyArrayObject *tod_array;
    PyArrayObject *dets_array;
    int full_freq;

    PyArrayObject *tform_array;

    if (!PyArg_ParseTuple(args, "O!O!i",
                          &PyArray_Type, &tod_array,
                          &PyArray_Type, &dets_array,
			  &full_freq
            ))
        po_raise("invalid arguments");

    ASSERT_CARRAY_TYPE_NDIM(tod_array, NPY_FLOAT32, 2);
    ASSERT_CARRAY_TYPE_NDIM(dets_array, NPY_INT64, 1);

    long int ndata = PyArray_DIMS(tod_array)[1];
    int ndet = PyArray_DIMS(dets_array)[0];
    long int *dets = PyArray_DATA(dets_array);
    float *data = PyArray_DATA(tod_array);

    npy_intp dims[2];
    dims[0] = ndet;
    if (full_freq)
	dims[1] = ndata;
    else
	dims[1] = ndata/2+1;

    tform_array = (PyArrayObject *)PyArray_SimpleNew(2, dims, NPY_CFLOAT);
    fftwf_complex *tform = (fftwf_complex *)PyArray_DATA(tform_array);

    // Make a plan to share.
    float *test_data = fftwf_malloc(ndata * sizeof(float));
    fftwf_plan f_plan = fftwf_plan_dft_r2c_1d(ndata, test_data, tform,
						 FFTW_ESTIMATE);
    fftwf_free(test_data);
    
#pragma omp parallel for shared(ndata, ndet, dets, stdout, data, tform, f_plan, \
				dims, full_freq) default (none)
    for (int d=0; d<ndet; d++)
    {
	//Source
	float *src = data + dets[d]*ndata;
	fftwf_complex *dest = tform + d*dims[1];
	fftwf_execute_dft_r2c(f_plan, src, dest);
	// Unload conjugate.  This is not very expensive.
	if (full_freq) {
	    for (int i=1; i<ndata/2; i++) {
		dest[ndata-i][0] =  dest[i][0];
		dest[ndata-i][1] = -dest[i][1];
	    }
	}
    }
    // End parallel block
    fftwf_destroy_plan(f_plan);

    return PyArray_Return(tform_array);
}

PyDoc_STRVAR(ifft_f_r__doc__,
	     "ifft_f_r(fdata, dets, ndata)\n"
	     "\n"
	     "Take 32-bit IFFT of fdata, parallelized.  Arrays must be C-ordered\n"
	     "with dimensions\n"
	     "and types:\n"
	     "     fdata[     *,n_freq]   (complex64)\n"
	     "     dets [n_det]           (int64)\n"
	     "\n"
	     "The dets array gives indices into the first dimension of fdata, "
	     "so the largest value in dets must be smaller than fdata.shape[0].\n"
	     "\n"
	     "If ndata<=0, then fdata is assumed to contain a full set of\n"
	     "positive and negative frequencies, so that the output vectors\n"
	     "have shape (n_det,ndata).\n"
	     "\n"
	     "Otherwise, fdata is assumed to contain the minimal set of\n"
	     "positive frequencies, and ndata must be set such that\n"
	     "ndata/2 + 1 = fdata.shape[0].  This is the case to use when\n"
	     "you passed full_freq==0 to the fft_f_r call.\n"
    );

static PyObject *ifft_f_r(PyObject *self, PyObject *args)
{
    PyArrayObject *fft_array;
    PyArrayObject *dets_array;
    int ndata;
    PyArrayObject *itform_array;

    if (!PyArg_ParseTuple(args, "O!O!i",
                          &PyArray_Type, &fft_array,
                          &PyArray_Type, &dets_array,
			  &ndata
            ))
        po_raise("invalid arguments");

    ASSERT_CARRAY_TYPE_NDIM(fft_array, NPY_COMPLEX64, 2);
    ASSERT_CARRAY_TYPE_NDIM(dets_array, NPY_INT64, 1);

    long int nfreq = PyArray_DIMS(fft_array)[1];
    if (ndata <= 0)
	ndata = nfreq;
    int ndet = PyArray_DIMS(dets_array)[0];
    long int *dets = PyArray_DATA(dets_array);
    fftwf_complex *fdata = PyArray_DATA(fft_array);

    npy_intp dims[2];
    dims[0] = ndet;
    dims[1] = ndata;
    itform_array = (PyArrayObject *)PyArray_SimpleNew(2, dims, NPY_FLOAT32);
    float *itform = (float*)PyArray_DATA(itform_array);

    // Make a plan to share.
    fftwf_complex *test_data = fftwf_malloc(nfreq * sizeof(fftwf_complex));
    fftwf_plan f_plan = fftwf_plan_dft_c2r_1d(ndata, test_data, itform,
						 FFTW_ESTIMATE);
    fftwf_free(test_data);

    /* Note that IFFT requires making a copy of the input, otherwise
     * FFTW may trash it.  Also, only _execute is threadsafe, so
     * protect allocations of temporary storage with omp critical. */
    
#pragma omp parallel shared(ndata, nfreq, ndet, dets, stdout, fdata, itform, f_plan) \
    default (none)
    {
	// Don't allocate until we know the thread has work to do.
	fftwf_complex *src_copy = NULL;
#pragma omp for
	for (int d=0; d<ndet; d++)
	{
	    if (src_copy == NULL) {
#pragma omp critical
		{
		    src_copy = fftwf_malloc((ndata/2+1) * sizeof(*src_copy));
		}
	    }
	    // Copy input, and scale by 1/n since FFTW does not bother
	    fftwf_complex *src = fdata + dets[d]*nfreq;
	    for (int i=0; i<ndata/2+1; i++) {
		src_copy[i][0] = src[i][0] / ndata;
		src_copy[i][1] = src[i][1] / ndata;
	    }
	    float *dest = itform + d*ndata;
	    fftwf_execute_dft_c2r(f_plan, src_copy, dest);
	}
	if (src_copy != NULL) {
#pragma omp critical
	    {
		fftwf_free(src_copy);
	    }
	}
    }
    // End parallel block
    fftwf_destroy_plan(f_plan);

    return PyArray_Return(itform_array);
}

PyMethodDef pyactpol_fft_methods[] = {
    {"fft_f_r", fft_f_r, METH_VARARGS, fft_f_r__doc__},
    {"ifft_f_r", ifft_f_r, METH_VARARGS, ifft_f_r__doc__},
    {NULL, NULL, 0, NULL}        /* Sentinel */
};


