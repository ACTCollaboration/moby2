/*!
 * \file todUtils.c
 * Utilities for actSample and actTOD objects.
 */

#include "pyactpol.h"
#include <assert.h>
#include "fftw3.h"


static PyObject *freq_space_waterfall(PyObject *self, PyObject *args);
static float extractPSValue(fftwf_complex *fftdata, int nn, float f, float fmin, 
                            float fmax, float df);
static PyObject *time_space_waterfall(PyObject *self, PyObject *args);


/*--------------------------------------------------------------------------------*/
/// \brief Obtain frequency-space waterfall matrix for TOD
static PyObject *freq_space_waterfall(PyObject *self, PyObject *args)
{
    PyArrayObject *tod_array;
    PyArrayObject *freq_array;
    PyArrayObject *mat_array;
    float tsamp;

    if (!PyArg_ParseTuple(args, "O!O!f",
                          &PyArray_Type, &tod_array,
                          &PyArray_Type, &freq_array,
                          &tsamp
            ))
        po_raise("invalid arguments");

    int ndet = PyArray_DIMS(tod_array)[0];
    int n = PyArray_DIMS(tod_array)[1];
    int nfreq = PyArray_DIMS(freq_array)[0];
    double *freq = PyArray_DATA(freq_array);
    po_assert(PyArray_TYPE(tod_array) == NPY_FLOAT32);
    float *tod = PyArray_DATA(tod_array);

    npy_intp dims[2];
    dims[0] = ndet; dims[1] = nfreq;
    mat_array = (PyArrayObject *)PyArray_SimpleNew(2, dims, NPY_DOUBLE);
    double *mat = (double *)PyArray_DATA(mat_array);

    int nn = n/2 + 1;
    float df = 1.0/(float)n/tsamp;


#pragma omp parallel shared(tod,mat,ndet,nfreq,n,nn,df,freq) default (none)
    {
        fftwf_complex *fftdata = fftwf_malloc(nn * sizeof(fftwf_complex));
        float *data = fftwf_malloc(n * sizeof(float));
        fftwf_plan p_forward;
#pragma omp critical
        {
            p_forward = fftwf_plan_dft_r2c_1d(n, data, fftdata, FFTW_ESTIMATE);
        }

#pragma omp for
        for (int i = 0; i < ndet; i++) {
            for (int j = 0; j < n; j++) {
                data[j] = tod[i*n + j];
            }
            fftwf_execute(p_forward);
            float f0 = freq[0] - (freq[1]-freq[0])/2;
            float f1 = (freq[1]+freq[0])/2;
            mat[i*nfreq] = (double)extractPSValue(fftdata, nn, freq[0], f0, f1, df)*2.0/df/(double)n/(double)n;
            for (int j = 1; j < nfreq-1; j++) {
                f0 = (freq[j]+freq[j-1])/2;
                f1 = (freq[j+1]+freq[j])/2;
                mat[i*nfreq + j] = (double)extractPSValue(fftdata, nn, freq[j], f0, f1, df)*2.0/df/(double)n/(double)n;
            }
            f0 = (freq[nfreq-1] + freq[nfreq-2])/2;
            f1 = freq[nfreq-1] + (freq[nfreq-1]-freq[nfreq-2])/2;
            mat[(i+1)*nfreq - 1] = (double)extractPSValue(fftdata, nn, freq[nfreq-1], f0, f1, df)*2.0/df/(double)n/(double)n;
        }

#pragma omp critical
        {
            fftwf_destroy_plan(p_forward);
            fftwf_free(fftdata);
            fftwf_free(data);
        }
    }
    return PyArray_Return(mat_array);
}

static float extractPSValue(fftwf_complex *fftdata, int nn, float f, float fmin, float fmax, float df)
{
  float value, p1, p2;
  int ini = (int)ceil(fmin/df);
  int end = (int)(fmax/df);
  if (end > nn) end = nn;
  if (end - ini < 2) {
    ini = (int)(f/df);
    end = (int)ceil(f/df);
    p1 = fftdata[ini][0]*fftdata[ini][0] + fftdata[ini][1]*fftdata[ini][1];
    if (end > nn) {
      value = p1;
      return value;
    }
    else {
      p2 = fftdata[end][0]*fftdata[end][0] + fftdata[end][1]*fftdata[end][1];
      value = (p1*(end*df-f) + p2*(f-ini*df))/df;
      return value;
    }
  } else {
    value = 0.0;
    for (int i = ini; i <= end; i++) {
      value += (fftdata[i][0]*fftdata[i][0] + fftdata[i][1]*fftdata[i][1])/(end-ini+1);
    }
  }
  return value;
}



///--------------------------------------------------------------------------------
/// \brief Casting function for mbTimeSpaceWaterfall for using from Python. 
static PyObject *time_space_waterfall(PyObject *self, PyObject *args)
{
    PyArrayObject *tod_array;
    PyArrayObject *time_array;
    PyArrayObject *mat_array;
    float dt;

    if (!PyArg_ParseTuple(args, "O!O!f",
                          &PyArray_Type, &tod_array,
                          &PyArray_Type, &time_array,
                          &dt
            ))
        po_raise("invalid arguments");

    int ndet = PyArray_DIMS(tod_array)[0];
    int n = PyArray_DIMS(tod_array)[1];
    int ntime = PyArray_DIMS(time_array)[0];
    double *time = PyArray_DATA(time_array);
    float step = (time[ntime-1]-time[0])/(ntime-1);

    po_assert(PyArray_TYPE(tod_array) == NPY_FLOAT32);
    float *tod = PyArray_DATA(tod_array);

    npy_intp dims[2];
    dims[0] = ndet; dims[1] = ntime;
    mat_array = (PyArrayObject *)PyArray_SimpleNew(2, dims, NPY_DOUBLE);
    double *mat = (double *)PyArray_DATA(mat_array);

    assert(time[0] >= 0.0);
    assert(time[ntime-1] <= dt*n);

#pragma omp parallel shared(tod,mat,ndet,ntime,n,dt,time,step) default (none)
    {
#pragma omp for
        for (int i = 0; i < ndet; i++) {
            for (int j = 0; j < ntime; j++) {
                float t0 = time[j] - step/2;
                float t1 = time[j] + step/2;
                int ki = (int)ceil(t0/dt);
                int ke = (int)(t1/dt)+1;
                if (ki < 0) ki = 0;
                if (ke <= ki) ke = ki + 1;
                if (ki == n) ki -= 1;
                if ((ke > n) && (i<10)) ke = n;
                mat[i*ntime + j] = 0.0;
                for (int k = ki; k < ke; k++) mat[i*ntime + j] += tod[i*n + k];
                mat[i*ntime + j] /= ke-ki;
            }
        }
    }
    return PyArray_Return(mat_array);
}


PyMethodDef pyactpol_waterfall_methods[] = {
    {"freq_space_waterfall", freq_space_waterfall, METH_VARARGS, ""},
    {"time_space_waterfall", time_space_waterfall, METH_VARARGS, ""},
    {NULL, NULL, 0, NULL}        /* Sentinel */
};
