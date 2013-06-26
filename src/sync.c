#include "pyactpol.h"
#include <gsl/gsl_linalg.h>

/* sync.c
 *
 * scan-synchronous mode fitting and removal.
 */


static int fit_poly(float *x, double *y, int n, int order, int samp,
                    double *fit);
static int fit_sine(double *t, double *y, int n, double f,
                  double *fit, double *coeff);


/* Find the amplitude and phase of the synchronous signal in a set of
 * detectors together with the phase of the azimuth scan.
 *
 * az should already have the mean removed.  ctime should start at 0
 * to minimize rounding errors. nwin should be the number of samples
 * to use for the amplitude fit; i.e. an round number of azimuth scans.
 */

static int mbSyncFit(double* az, double* ctime,
		     float* data, int ndata,
		     int *dets, int ndets,
		     float *sync_amp, float* sync_phase,
		     float *sync_rms,
		     float *phaseAz,
		     float f,
		     int order, int samp, int nT)
{
  int n, nwin;
  double sampRate, thetaAz;

  // Calculate size of data window to do analysis
  sampRate = (double)ndata / (ctime[ndata-1] - ctime[0]);
  nwin = floor( (double)nT / f * sampRate );
  if (nwin > ndata) {
    nwin = ndata;
  }

  double *fit_az = malloc(nwin * sizeof(double));
  double *coeff_az = malloc(2 * sizeof(double));

  // Obtain phase of azimuth scan
  if (fit_sine(ctime, az, nwin, f, fit_az, coeff_az) != 0)
       return 1;
  if (coeff_az[0] > 0) thetaAz = atan(coeff_az[1]/coeff_az[0]);
  else thetaAz = atan(coeff_az[1]/coeff_az[0]) + M_PI;
  *phaseAz = thetaAz;

  free(fit_az);
  free(coeff_az);

  float *tcopy = malloc(nwin * sizeof(float));
  for (int j=0; j<nwin; j++)
       tcopy[j] = ctime[j] - ctime[0];

  // Process all selected detectors
#pragma omp parallel private(n)
  {
    double *y = malloc(nwin * sizeof(double));
    double *fit_1f = malloc(nwin * sizeof(double));
    double *fit_syn = malloc(nwin * sizeof(double));
    double *coeff_syn = malloc(2 * sizeof(double));
    
#pragma omp for
    for (n = 0; n < ndets; n++) {
      int i;
      double A, theta, rms;

      float *det_data = data + ndata * dets[n];

      // Obtain copy of data
      for (i = 0; i < nwin; i++)
	   y[i] = det_data[i];

      // Fit and subtract polynomial
      fit_poly(tcopy, y, nwin, order, samp, fit_1f);
      for (i = 0; i < nwin; i++)
	y[i] -= fit_1f[i];
 
      // Fit sinusoidal to obtain amplitude and phase
      fit_sine(ctime, y, nwin, f, fit_syn, coeff_syn);
      A = sqrt(coeff_syn[0]*coeff_syn[0] + coeff_syn[1]*coeff_syn[1]);
      if (A == 0) theta = 0.0;
      else if (coeff_syn[0] > 0) theta = atan(coeff_syn[1]/coeff_syn[0]);
      else theta = atan(coeff_syn[1]/coeff_syn[0]) + M_PI;
      sync_amp[n] = A;
      sync_phase[n] = theta;
      
      // Find fit RMS error
      rms = 0.0;
      for (i = 0; i < nwin; i++)
	   rms += (y[i] - fit_syn[i])*(y[i] - fit_syn[i]);
      sync_rms[n] = sqrt(rms/nwin);
    } /* pragma omp for */

    free(y);
    free(fit_1f);
    free(fit_syn);
    free(coeff_syn);
  } /* pragma omp parallel */


  return 0;
}

/// Fit sinusoidal of a given frequency to data.
/// \param t      Independent variable (time).
/// \param y      Dependent variable (data).
/// \param n      Number of elements in data.
/// \param f      Frequency of the sinusoid.
/// \param fit    Vector of fitted data (sinusoid).
/// \param coeff  Coefficients of sinusoid (A*cos(th) + B*sin(th)).
static int fit_sine(double *t, double *y, int n, double f, 
	          double *fit, double *coeff) {
  int k;

  double temp_c, temp_s;
  double cos2, sin2, sincos, ycos, ysin;

  cos2 = sin2 = sincos = ycos = ysin = 0.0; 
  for (k = 0; k < n; k++) {
    temp_c = cos(2*M_PI*f*t[k]);
    temp_s = sin(2*M_PI*f*t[k]);
    cos2 += temp_c*temp_c;
    sin2 += temp_s*temp_s;
    sincos += temp_s*temp_c;
    ycos += y[k]*temp_c;
    ysin += y[k]*temp_s;
  }

  coeff[0] = (sin2*ycos - sincos*ysin) / (cos2*sin2 - sincos*sincos);
  coeff[1] = (cos2*ysin - sincos*ycos) / (cos2*sin2 - sincos*sincos);

  // Evaluate result
  for (k = 0; k < n; k++) {
    fit[k] = cos(2*M_PI*f*t[k])*coeff[0] + sin(2*M_PI*f*t[k])*coeff[1];
  }

  return 0;
}


/// Fit polynomial to data to remove common mode.
/// \param x       Independent variable.
/// \param y       Dependent variable (data).
/// \param n       Number of elements in data.
/// \param order   Order of polinomial to fit.
/// \param samp    Number of samples to take from data to do the fit.
/// \param fit     vector with evaluated polynomial

static int fit_poly(float *x, double *y, int n, int order, int samp, 
                    double *fit) {
  int i, j, k;
  int N;
  double *coeff;
  double **A, *pows;
  double temp;

  gsl_vector *g_Ab = gsl_vector_calloc(order+1);
  gsl_matrix *g_AA = gsl_matrix_alloc(order+1,order+1);
  gsl_vector *g_x = gsl_vector_alloc(order+1);
  gsl_permutation *g_p = gsl_permutation_alloc(order+1);
  int signum = 0;

  int status = 1; // failure.

  if (order <= 0) {
       print_error("Invalid 'order=%i' in fit_poly\n", order);
       return 1;
  }

  // Initializa Matrices and Vectors
  coeff = malloc((order+1) * sizeof(double));
  N = n/samp;
  if (N < order) {
       print_error("Fitting 'order=%i' with 'N=%i' data points in fit_poly\n", order, N);
       return 1;
  }
  A = (double **)malloc((order+1) * sizeof(double *));
  A[0] = (double *)malloc(N * (order+1) * sizeof(double));
  for (i = 0; i < order+1; i++){
    A[i] = A[0] + N * i;
    for( j = 0; j < N; j++ ) A[i][j] = 0.0;
  }
  
  pows = (double *)malloc((2*order+1) * sizeof(double));
  for (i = 0; i < 2*order+1; i++)
    pows[i] = 0.0;

  // Generate A matrix
  for (i = 0; i < order+1; i++)
    for (k = 0; k < N; k++) {
      A[i][k] = 1.0;
      for (j = 0; j < i; j++)
        A[i][k] *= x[k*samp + samp/2];
    }

  // Generate AA as transpose(A)*A
  for (i = 2*order; i > 0; i -= 2) {
    for (k = 0; k < N; k++) {
      pows[i] += A[i/2][k]*A[i/2][k];
      pows[i-1] += A[i/2][k]*A[i/2-1][k];
    }
  }
  for (k = 0; k < N; k++)
    pows[0]++;
  
  for (i = 0; i < order+1; i++)
    for (j = 0; j < order+1; j++)
      g_AA->data[i*g_AA->tda+j] = pows[i+j];
	
/*  for (i = 0; i < order+1; i++) {
	for (j = 0; j <= i; j++) {
	  for (k = 0; k < N; k++)
		AA[i][j] += A[j][k]*A[i][k];
	}
  }*/

  // Generate Ab as transpose(A)*y
  for (i = 0; i < order+1; i++) {
    for (k = 0; k < N; k++) {
      temp = 0.0;
      for (j = 0; j < samp; j++)
        temp += y[k*samp + j];
      g_Ab->data[i] += A[i][k] * temp / (float)samp;
    }
  }

  // Solve system
  if ((status = gsl_linalg_LU_decomp(g_AA, g_p, &signum))!=0)
      goto exit_now;
  if ((status = gsl_linalg_LU_solve(g_AA, g_p, g_Ab, g_x))!=0)
      goto exit_now;

  // Copy out...
  for (i = 0; i < order+1; i++)
      coeff[i] = g_x->data[i];
  
  // Evaluate result
  for (k = 0; k < n; k++) {
    fit[k] = 0.0;
    for (i = 0; i < order+1; i++) {
      temp = 1.0;
      for (j = 0; j < i; j++)
        temp *= x[k];
      fit[k] += temp * coeff[i];
    }
  }

/*
  for (k = 0; k < n; k++)
	fit[k] = 0.0;
  for (i = 0; i < order+1; i++)
	for (k = 0; k < n; k++)
	  fit[k] += A[i][k] * Ab[i];
*/

  status = 0;

exit_now:
  free(A[0]);
  
  free(A);
  free(pows);
  free(coeff);

  gsl_vector_free(g_Ab);
  gsl_vector_free(g_x);
  gsl_matrix_free(g_AA);
  gsl_permutation_free(g_p);

  return status;
}
	

PyDoc_STRVAR(get_sync_amps__doc__,
	     "get_sync_amps(data, dets, az, ctime)\n"
	     "\n"
	     "The arrays must be C-ordered with dimensions like:\n"
	     "     data [     *,n_data]    (float)\n"
	     "     dets [n_det]    	   (int)\n"
             "     az   [n_data]   	   (double)\n"
	     "     ctime[n_data]   	   (double)\n"
	     "\n"
	     "Returns (az_phase, amps_cos, amps_sin)."
     );

static PyObject *get_sync_amps(PyObject *self, PyObject *args)
{
    PyArrayObject *data_array;
    PyArrayObject *dets_array;
    PyArrayObject *az_array;
    PyArrayObject *ctime_array;
    double scan_freq;
    int order, samp, nT;

    if (!PyArg_ParseTuple(args, "O!O!O!O!diii",
                          &PyArray_Type, &data_array,
                          &PyArray_Type, &dets_array,
                          &PyArray_Type, &az_array,
                          &PyArray_Type, &ctime_array,
			  &scan_freq,
			  &order,
			  &samp,
			  &nT
            ))
        po_raise("invalid arguments.");

    // Types and ordering
    ASSERT_CARRAY_TYPE_NDIM(data_array, NPY_FLOAT32, 2);
    ASSERT_CARRAY_TYPE_NDIM(dets_array, NPY_INT32, 1);
    ASSERT_CARRAY_TYPE_NDIM(az_array, NPY_FLOAT64, 1);
    ASSERT_CARRAY_TYPE_NDIM(ctime_array, NPY_FLOAT64, 1);

    int ndata = PyArray_DIMS(data_array)[1];
    int ndet = PyArray_DIMS(dets_array)[0];

    po_assert(PyArray_DIMS(az_array)[0] == ndata);
    po_assert(PyArray_DIMS(ctime_array)[0] == ndata);

    float *data = PyArray_DATA(data_array);
    double *az = PyArray_DATA(az_array);
    double *ctime = PyArray_DATA(ctime_array);
    int *dets = PyArray_DATA(dets_array);

    // We've been so thorough, it would be a shame to segfault now.
    for (int i=0; i<ndet; i++)
	 po_assert(dets[i] >= 0 &&
		   dets[i] < PyArray_DIMS(data_array)[0]);

    // And, places for the results.
    npy_intp ndet_ = ndet;
    PyArrayObject *amp_array = (PyArrayObject*)
	 PyArray_SimpleNew(1, &ndet_, NPY_FLOAT32);
    PyArrayObject *phase_array = (PyArrayObject*)
	 PyArray_SimpleNew(1, &ndet_, NPY_FLOAT32);
    PyArrayObject *rms_array = (PyArrayObject*)
	 PyArray_SimpleNew(1, &ndet_, NPY_FLOAT32);
    
    float phaseAz = -1.;
    mbSyncFit(az, ctime,
	      data, ndata,
	      dets, ndet,
	      PyArray_DATA(amp_array),
	      PyArray_DATA(phase_array),
	      PyArray_DATA(rms_array),
	      &phaseAz,
	      scan_freq, order, samp, nT);

    return Py_BuildValue("NNNf",
			 amp_array,
			 phase_array,
			 rms_array,
			 phaseAz);
}

PyMethodDef pyactpol_sync_methods[] = {
    {"get_sync_amps", get_sync_amps, METH_VARARGS, 
    get_sync_amps__doc__},
    {NULL, NULL, 0, NULL}        /* Sentinel */
};


