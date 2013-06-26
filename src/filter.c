/* -*- mode: C; tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 4 -*-
 *      vim: sw=4 ts=8 et tw=80
 */

/* filter.c
 *
 * This file contains a few time-domain processing routines, copied
 * without much modification from moby.
 *
 * These are mostly needed for filtering during mapping, and for
 * various pathologies processing.
 */


#include "pyactpol.h"
#include "myassert.h"

#include <stdbool.h>

#include <fftw3.h>
#include <gsl/gsl_cblas.h>

static PyObject *remove_mean(PyObject *self, PyObject *args)
{
    PyArrayObject *tod_array;
    PyArrayObject *dets_array;
    PyArrayObject *means_array;

    if (!PyArg_ParseTuple(args, "O!O!",
                          &PyArray_Type, &tod_array,
                          &PyArray_Type, &dets_array
            ))
        po_raise("invalid arguments");

    ASSERT_CARRAY_TYPE_NDIM(tod_array, NPY_FLOAT32, 2);
    ASSERT_CARRAY_TYPE_NDIM(dets_array, NPY_INT64, 1);

    long int ndata = PyArray_DIMS(tod_array)[1];
    int ndet = PyArray_DIMS(dets_array)[0];
    long int *dets = PyArray_DATA(dets_array);
    float *data = PyArray_DATA(tod_array);

    npy_intp dims[1];
    dims[0] = ndet;
    means_array = (PyArrayObject *)PyArray_SimpleNew(1, dims, NPY_DOUBLE);
    double *means = (double *)PyArray_DATA(means_array);

#pragma omp parallel for shared(data,ndata,dets,ndet)
  for (int d=0; d<ndet; d++)
  {
    double m = 0.0;
    for (long int i = dets[d]*ndata; i < (dets[d]+1)*ndata; i++)
      m += (double)data[i];
    m /= ndata;
    means[d] = m;
    for (long int i = dets[d]*ndata; i < (dets[d]+1)*ndata; i++)
      data[i] -= m;
  }
  return PyArray_Return(means_array);
}


static void filter_one_detector(
    fftwf_complex *data,
    const float *real, const float *imag,
    int ndata,
    fftwf_plan p_forward, fftwf_plan p_back, 
    int detrend, int retrend);

static void mbTODFourierFilter(
    float *tod,        ///  2d data (nd,nsamps) to filter
    int nsamps,        ///< data size, fast dimension
    const int *dets,  ///< tod detector numbers -- designates
                       ///< which detectors to filter
    int ndets,         ///< length of the dets vector
    const float *freq, ///> frequency vector (nsamps)
    const float *real, ///< real component of filter, (nsamps)
    const float *imag, ///< imaginary component of filter, (nsamps)
    const float *tau,  ///< time constants to deconvolve (ndets)
    int tau_do_inverse,///< true if deconvolving time constants
    int detrend,       ///< whether to detrend the data before filtering
    int retrend)       ///< whether to retrend the detrended data after filtering
{
    //float *freqs;
    //freqs = genFreqs(tod->ndata, tod->sampleTime);

    int n = nsamps;
    fftwf_complex *data = fftwf_malloc(n * sizeof(fftwf_complex));
    fftwf_plan p_forward;
    fftwf_plan p_back;
    p_forward = fftwf_plan_dft_1d(n, data, data, FFTW_FORWARD, FFTW_ESTIMATE);
    p_back    = fftwf_plan_dft_1d(n, data, data, FFTW_BACKWARD, FFTW_ESTIMATE);
    fftwf_free(data);

#pragma omp parallel shared(tod,real,imag,dets,ndets,freq,tau,tau_do_inverse, \
    detrend,retrend,n,p_forward,p_back) default (none)
    {
        fftwf_complex *fftdata = fftwf_malloc(n * sizeof(fftwf_complex));
            // Customize the filter with a time constant per-detector
        float *tmp_real = malloc(n * sizeof(float));
        float *tmp_imag = malloc(n * sizeof(float));
        if (tau == NULL) {
            for (int j=0; j<n; j++) {
                tmp_real[j] = real[j];
                tmp_imag[j] = imag[j];
            }
        }

#pragma omp for
        for (int i = 0; i < ndets; i++) {
            for (int j = 0; j < n; j++) {
                if (tau != NULL) {
                    float gain, dtc;
                    dtc = 2.0*tau[i]*freq[j]*M_PI;
                    if (tau_do_inverse) {
                        gain = 1.0;
                    } else {
                        gain = 1.0 / (1.0 + dtc*dtc);
                        dtc = -dtc;
                    }
                    tmp_real[j] = gain * (real[j] - imag[j]*dtc);
                    tmp_imag[j] = gain * (real[j]*dtc + imag[j]);
                }
                fftdata[j][0] = tod[dets[i]*n+j];
                fftdata[j][1] = 0.0;
            }
            filter_one_detector(fftdata, tmp_real, tmp_imag, n, p_forward, p_back, detrend, retrend);
            for (int j = 0; j < n; j++) tod[dets[i]*n+j] = fftdata[j][0];
        }

        fftwf_free(fftdata);
        free(tmp_real);
        free(tmp_imag);
    }
#pragma omp critical
    {
        fftwf_destroy_plan(p_forward);
        fftwf_destroy_plan(p_back);
    }
}


/*--------------------------------------------------------------------------------*/
/*!
 * Filter a detector from a TOD with a give Fourier-domain filter.
 * \param data      The data vector to be filtered.  Modified and returns the filtered result.
 * \param real      The real part of the filter.
 * \param imag      The imaginary part of the filter.
 * \param p_forward The FFTW plan for going to the frequency domain.
 * \param p_back    The FFTW plan for going back to the time domain.
 * \param ndata     The length of the data set.
 */
static void filter_one_detector(fftwf_complex *data, const float *real, const float *imag,
                                int ndata, fftwf_plan p_forward, fftwf_plan p_back, 
                                int detrend, int retrend)
{
  float tmp;
  double x0 = 0.0, x1 = 0.0, m = 0, trendSlope = 0;
  int win = 1000;
  if (win > ndata/2) win = ndata/2;
  if (detrend) {
    for (int i = 0; i < win; i++) {
      x0 += data[i][0];
      x1 += data[ndata-1 - i][0];
    }
    x0 /= (double)win;
    x1 /= (double)win;
 
    m = (x0+x1)/2.0;
    trendSlope = (x1-x0)/(double)(ndata-1.0-win);
    
    for (int j = 0; j < ndata; j++)
      data[j][0] -= (x0 - m + trendSlope*((double) (j - win/2)));
  }

  fftwf_execute_dft(p_forward, data, data);
  for (int j = 0; j < ndata; j++) {
    // Multiply by complex filter
    tmp = data[j][0]*real[j] - data[j][1]*imag[j];
    data[j][1] = (data[j][0]*imag[j] + data[j][1]*real[j]) / (float)ndata;
    data[j][0] = tmp / (float)ndata;
  }
  fftwf_execute_dft(p_back, data, data);

  if (retrend && detrend) 
    for (int j = 0; j < ndata; j++)
      data[j][0] += (x0 - m + trendSlope*((double) (j - win/2)));

}


/*--------------------------------------------------------------------------------*/
/*!
 * Apply a pre-calculated smoothing filter to a data vector.  Then do one of
 * (1) Just smooth the data if (do_smooth).
 * (2) Preserve the data, replacing with smoothed data only during large glitches if (apply_glitch).
 * (3) Fill a vector cutvec saying whether to cut each point due to glitches if (do_cuts).
 * You can call any combination that doesn't involve both #1 and #2.
 * If apply_glitch is true, then replace anything flagged as over the cut threshold (defined as a
 * smoothed value differing from raw by at least nsig*median absolute deviation of smoothed values)
 * with the smoothed value.  Since the glitch filter should have zero in the middle, height of
 * cosmic ray doesn't matter.
 */

typedef float actData;

static void glitch_one_detector(
  fftwf_complex *fftdata, ///< Real data vector to filter.  Modified and returned.
  const float *real,      ///< The pre-calculated smoothing filter to use (real part).
  char *cutvec,           ///< Pre-allocated vector of length n.  Returns cut list if (do_cuts).
  int n,                  ///< Length of the input data vector mydat.
  fftwf_plan p_forward,   ///< The FFTW plan for going to the frequency domain.
  fftwf_plan p_back,      ///< The FFTW plan for going back to the time domain.
  actData nsig,           ///< The cut threshold is this factor times the median abs deviation of smooths.
  int minSep)             ///< Min separation between cut detections to consider them different cuts.
{
  // Save a copy of the input data, then apply the smoothing filter only.
  actData *mydat = malloc(sizeof(actData)*n);
  float *imag = malloc(n * sizeof(float));
  for (int i = 0; i < n; i++) {
    mydat[i] = fftdata[i][0];
    imag[i] = 0.0;
  }
  filter_one_detector(fftdata, real, imag, n, p_forward, p_back, true, false);
  free(imag);

  // tmpclean is the absolute difference between smooth and raw vectors.  
  // Find its median (for use in cuts).
  float *tmpclean = malloc(n*sizeof(actData));
  for (int j = 0; j < n; j++)
      tmpclean[j] = fabs(fftdata[j][0]);
  actData thresh = pyactpol_sselect((n*3)/4, n, tmpclean);
  thresh -= pyactpol_sselect(n/4, n, tmpclean);
  thresh *= 0.741*nsig;

  // Find cuts and reset raw data
  int lind = 0;
  memset(cutvec,0,sizeof(*cutvec)*n);
  for (int j = 0; j < n; j++) {
      if (fabs(fftdata[j][0]) > thresh) {
          cutvec[j] = 1;
          if (j-lind < minSep)
              for (int k = lind+1; k < j; k++) cutvec[k] = 1;
          lind = j;
      }
      fftdata[j][0] = mydat[j];
  }

  free(mydat);
  free(tmpclean);
}


/*--------------------------------------------------------------------------------*/
/*!
 * Remove glitches from multiple detector data vectors in a TOD using a pre-computed smoothing (and
 * glitch removal) filter. Also, if (do_cuts), then also apply all glitches as new cuts in the cuts
 * object.
 */

static void mbGlitchC(
    float *tod_data,       ///< The TOD to modify.
    int nsamps,            ///< tod.shape[1]
    const float *filt,     ///< The pre-calculated smoothing filter to use.
    const int *dets,       ///< List of detectors (by index) to clean.
    int ndet,              ///< Number of detectors to clean (length of dets)
    PyObject **cuts,       ///< An array of pointers to PyObjects, where we should put the cuts.
    actData nsig,          ///< The cut threshold is this factor times the median abs deviation of smooths.
    int maxGlitch,         ///< Maximum allowed number of glitches.  If exceeded, cut whole detector.
    int minSep)            ///< Minimum separation between cuts to consider them different cuts.
{

    int **cuts_data = malloc(ndet * sizeof(*cuts_data));
    int *n_cuts = malloc(ndet * sizeof(*n_cuts));
    assert(cuts_data != NULL && n_cuts != NULL);

#pragma omp parallel shared(tod_data,filt,nsamps,dets,ndet,cuts_data,n_cuts,nsig,maxGlitch,minSep) default(none)
    {
        fftwf_plan p_forward;
        fftwf_plan p_back;
    
        const int n = nsamps;
        fftwf_complex *fftdata = fftwf_malloc(sizeof(fftwf_complex)*n);
        char *cutvec = malloc(n*sizeof(int));
    
#pragma omp critical
        {
            p_forward = fftwf_plan_dft_1d(n, fftdata, fftdata, FFTW_FORWARD, FFTW_ESTIMATE);
            p_back    = fftwf_plan_dft_1d(n, fftdata, fftdata, FFTW_BACKWARD, FFTW_ESTIMATE);
        }

#pragma omp for
        for (int i = 0; i < ndet; i++) {
            for (int j = 0; j < n; j++) {
                fftdata[j][0] = tod_data[dets[i]*nsamps+j];
                fftdata[j][1] = 0.0;
            }

            glitch_one_detector(fftdata, filt, cutvec ,n, p_forward, p_back, nsig, minSep);

            for (int j = 0; j < n; j++) tod_data[dets[i]*nsamps+j] = fftdata[j][0];

            /* Convert cuts vector start / stop positions */
            cuts_data[i] = pyactpol_cuts_vector_from_mask(cutvec, nsamps, &n_cuts[i]);
        }
    
#pragma omp critical
        {
            fftwf_destroy_plan(p_forward);
            fftwf_destroy_plan(p_back);
            fftwf_free(fftdata);
            free(cutvec);
        }
    }
    /* Non-parallel: create arrays for the resulting cuts vectors. */
    npy_intp dims[2] = {0, 2};
    for (int i=0; i<ndet; i++) {
        dims[0] = n_cuts[i];
        cuts[i] = PyArray_SimpleNew(2, dims, NPY_INT);
        if (n_cuts[i] > 0) {
            int *dest = PyArray_DATA((PyArrayObject *)cuts[i]);
            memcpy(dest, cuts_data[i], 2*n_cuts[i]*sizeof(*dest));
            free(cuts_data[i]);
        }
    }
}



/* filter_tod_data
 *
 * Apply an arbitrary Fourier filter to time ordered data.  Possibly
 * also convolve or deconvolve the effects of a time constant.
 *
 * Input:
 *   data: (ndet,nsamps) time-ordered data to filter (modified in place).
 *   dets: (ndet_used) index of data vectors to actually filter.
 *   filter: (nsamps) complex filter to apply to the data
 *   taus: (ndet_used) time constants, in seconds, to convolve/deconvolve.
 *   dt: the sample separation, in seconds.
 *   tau_gain: boolean, set to true for time constant deconvolution.
 *   detrend: boolean, set to true to remove a trend prior to filtering.
 *   retrend: boolean, set to true to restore the removed trend after filtering.
 */

PyDoc_STRVAR(filter_tod_data__doc__,
	     "filter_tod_data(data, dets, filter, time_const)\n"
	     "\n"
	     "The dets array gives indices into the first dimension of data, "
	     "so the largest value in dets must be smaller than data.shape[0].\n"
	     "\n"
	     "Applies\n"
	     "        data[dets,:] *= cal[:,None]\n"
    );


static PyObject *filter_tod_data(PyObject *self, PyObject *args)
{
    PyArrayObject *tod_array;
    PyArrayObject *det_idx_array;
    PyArrayObject *filter_array;
    PyObject *taus_object;
    double dt;
    int tau_do_inverse;
    int detrend;
    int retrend;

    if (!PyArg_ParseTuple(args, "O!O!O!Odiii",
                          &PyArray_Type, &tod_array,
                          &PyArray_Type, &det_idx_array,
                          &PyArray_Type, &filter_array,
                          &taus_object,
                          &dt, &tau_do_inverse, &detrend, &retrend
            ))
        po_raise("invalid arguments");

    // Make a frequency vector
    int nsamps = PyArray_DIMS(tod_array)[1];
    float *freq = malloc(nsamps*sizeof(*freq));
    float df = 1.0/nsamps/dt;
    int i=0;
    for (; i<nsamps/2; i++)
        freq[i] = i*df;
    for (; i<nsamps; i++)
        freq[i] = (i-nsamps)*df;

    // Copy filter to real and imaginary pieces
    ASSERT_CARRAY_TYPE_NDIM(filter_array, NPY_COMPLEX64, 1);
    po_assert(PyArray_DIMS(filter_array)[0] == nsamps);
    float *filter = PyArray_DATA(filter_array);
    float *real = malloc(nsamps*sizeof(*real));
    float *imag = malloc(nsamps*sizeof(*imag));
    for (i=0; i<nsamps; i++) {
        real[i] = filter[i*2];
        imag[i] = filter[i*2+1];
    }

    // Those other vectors
    po_assert(PyArray_TYPE(tod_array) == NPY_FLOAT32);
    float *tod = PyArray_DATA(tod_array);
    ASSERT_CARRAY_TYPE_NDIM(det_idx_array, NPY_INT32, 1);
    int *det_idx = PyArray_DATA(det_idx_array);
    float *tau = NULL;
    if (taus_object != Py_None) {
        PyArrayObject *taus_array = (PyArrayObject*)taus_object;
        ASSERT_CARRAY_TYPE_NDIM(taus_array, NPY_FLOAT32, 1);
        po_assert(
                  PyArray_DIMS(det_idx_array)[0]
                  == PyArray_DIMS(taus_array)[0]);
        tau = PyArray_DATA(taus_array);
    }        

    // Go
    mbTODFourierFilter(tod, nsamps, det_idx, PyArray_SIZE(det_idx_array),
                       freq, real, imag, tau,
                       tau_do_inverse, detrend, retrend);

    free(freq);
    free(real);
    free(imag);
    Py_RETURN_NONE;
}


/* deglitch_tod_data
 *
 * Input:
 *   data: (ndet,nsamps) time-ordered data to filter (modified in place).
 *   dets: (ndet_used) index of data vectors to actually filter.
 *   filter: (real!) low-pass filter that will be used to remove low
 *     frequency signal
 *   nsig: float that controls strength of glitch rejection.
 *   maxGlitch: maximum number of allowed glitches; when limit is reached
 *     the whole detector will be cut.
 *   minSep: integer, minimum separation cuts must have to not get
 *     auto-joined.
 *
 * Returns:
 *   cuts: list of cuts vectors, indicating which samples were deglitched.
 *   
 */

static PyObject *get_glitch_cuts(PyObject *self, PyObject *args)
{
    PyArrayObject *tod_array;
    PyArrayObject *det_idx_array;
    PyArrayObject *filter_array;
    double nsig;
    int maxGlitch;
    int minSep;

    if (!PyArg_ParseTuple(args, "O!O!O!dii",
                          &PyArray_Type, &tod_array,
                          &PyArray_Type, &det_idx_array,
                          &PyArray_Type, &filter_array,
                          &nsig, &maxGlitch, &minSep
            ))
        po_raise("invalid arguments.");

    int nsamps = PyArray_DIMS(tod_array)[1];
    ASSERT_CARRAY_TYPE_NDIM(filter_array, NPY_FLOAT32, 1);
    float *filter = PyArray_DATA(filter_array);

    // Those other vectors
    po_assert(PyArray_TYPE(tod_array) == NPY_FLOAT32);
    float *tod = PyArray_DATA(tod_array);
    ASSERT_CARRAY_TYPE_NDIM(det_idx_array, NPY_INT32, 1);
    int *det_idx = PyArray_DATA(det_idx_array);

    // Deglitcher will return some cuts information
    int ndets = PyArray_DIMS(det_idx_array)[0];
    PyObject **cuts = malloc(ndets*sizeof(PyObject*));
    assert(cuts != NULL);

    // Go
    mbGlitchC(tod, nsamps, filter, det_idx, PyArray_SIZE(det_idx_array),
              cuts,
              nsig, maxGlitch, minSep);

    // Turn the cuts into a list.
    PyObject *cuts_list = PyList_New(ndets);
    for (int i=0; i<ndets; i++)
        PyList_SET_ITEM(cuts_list, i, cuts[i]);

    free(cuts);

    return cuts_list;
}



// --------------------------------------------------------------------
// Analyzes the scan to determine scan frequency, amplitude, speed, and 
// identify sections of good and bad scan
/* Error codes:
   0: Good scan
   1: No scan
   2: Bad scan ending
   3: Bad scan beginning
   4: Multiple scan sections
   5: Too many scan sections
*/

typedef struct {
    int min_chunk;
    float min_scan;
    float max_scan;
    float freq;
    float speed;
    float az_max;
    float az_min;
    int *section_limits;
    int nsection;
    int max_sections;
} mbScan;

static int mbAnalyzeScan(double *az, double *ctime, int ndata,
                         mbScan *scan)
{
    // Hack the downsampleLevel to none...
    int downsampleLevel = 1;

    double sampleTime = ctime[1] - ctime[0];
    
    // Consider N_chunk segments of length chunck
    int N_chunk = (ndata*downsampleLevel)/scan->min_chunk;
    if (N_chunk == 0) N_chunk = 1;
    int chunk = ndata/N_chunk;
    int i, k, nSamp;
    double y, Sy, Sxy;
    double az_Hi, az_Low;

    float *az_max = (float *)malloc(N_chunk * sizeof(float));
    float *az_min = (float *)malloc(N_chunk * sizeof(float));
    float *freq = (float *)malloc(N_chunk * sizeof(float));
    float *speeds = (float *)malloc(N_chunk * sizeof(float));

    // For each chunk
    for (int c = 0; c < N_chunk; c++) {
        i = c*chunk + 1;
        /* Determine az_min and az_max by finding the first two
         * turn-arounds in this chunk, if possible. */
        if (az[i] > az[i-1]) {
            while ((az[i] >= az[i-1]) && (i < (c+1)*chunk)) i++;
            az_max[c] = az[i-1];
            if (i > c*chunk + 1 + chunk/2) {
                freq[c] = -1;
                continue;
            }
            i += 100;
            while ((az[i] <= az[i-1]) && (i < (c+1)*chunk)) i++;
            az_min[c] = az[i-1];
        }
        else {
            while ((az[i] <= az[i-1]) && (i < (c+1)*chunk)) i++;
            az_min[c] = az[i-1];
            if (i > c*chunk + 1 + chunk/2) {
                freq[c] = -1;
                continue;
            }
            i += 100;
            while ((az[i] >= az[i-1]) && (i < (c+1)*chunk)) i++;
            az_max[c] = az[i-1];
        }
        /* Reject weird scan amplitudes */
        if ((az_max[c] - az_min[c] < scan->min_scan) || 
            (az_max[c] - az_min[c] > scan->max_scan) ||
            (i == (c+1)*chunk)) {
            freq[c] = -1;
            continue;
        }
        /* Now go through chunk again, getting speed measurements at
         * every sample that lies in the central 20% of the az range.
         * Also get the mean time of each such segment, and basically
         * fit a line to those times to get the scan frequency.  */
        az_Hi = az_max[c]*0.6 + az_min[c]*0.4;
        az_Low = az_max[c]*0.4 + az_min[c]*0.6;
        nSamp = 0;
        Sy = 0;
        Sxy = 0;
        speeds[c] = 0.0;
        i = c*chunk;
        /* Skip any initial segment that is in the central range */
        while ((az[i] > az_Low) && (az[i] < az_Hi) && (i < (c+1)*chunk)) i++;
        while (i < (c+1)*chunk) {
            if ((az[i] > az_Low) && (az[i] < az_Hi)) {
                nSamp++;
                k = 0;
                y = 0.0;
                while ((az[i] > az_Low) && (az[i] < az_Hi) && (i < (c+1)*chunk)) {
                    y += ctime[i] - ctime[0];
                    k++;
                    i++;
                }
                speeds[c] += fabs(az[i-1]-az[i-k-1])/k;
                y /= k;
                Sy += y;
                Sxy += nSamp*y;
            }
            else i++;
        }
        speeds[c] /= nSamp*sampleTime;
        if (nSamp > 1) freq[c] = nSamp*(nSamp-1) / 12.0 / (2*Sxy/(nSamp+1) - Sy);
        else {
            freq[c] = -1;
        }
    }

    /* Now just throw out any bad chunks, collapsing the various
     * arrays. */
    k = 0;
    for (int c = 0; c < N_chunk; c++) {
        /* printf("freq[%d] = %f; speed[%d] = %f; az_max[%d] = %f; az_min[%d] = %f\n", */
        /*        c,freq[c], c, speeds[c], c, az_max[c], c, az_min[c]); */
        if (freq[c] != -1) {
            az_max[k] = az_max[c];
            az_min[k] = az_min[c];
            freq[k] = freq[c];
            speeds[k] = speeds[c];
            k++;
        }
    }
    /* You either of zero, many, or a small number of good chunks.  So
     * compute and store nothing, the median, or the mean,
     * respectively. */ 
    if (k == 0) { // This TOD seems to have no scanning
        return 1;
    }
    if (k > 2) {
        scan->az_max = compute_median(k, az_max);
        scan->az_min = compute_median(k, az_min);
        scan->freq = compute_median(k, freq);
        scan->speed = compute_median(k,speeds);
    }
    else {
        scan->az_max = scan->az_min = scan->freq = scan->speed = 0.0;
        for (i = 0; i < k; i++) {
            scan->az_max += az_max[i]/k;
            scan->az_min += az_min[i]/k;
            scan->freq += freq[i]/k;
            scan->speed += speeds[i]/k;
        }
    }

    free(az_max);
    free(az_min);
    free(freq);
    free(speeds);
    //printf("freq_m = %f; speed_m = %f\n", scan->freq, scan->speed);

/*
 * ACT 2
 *
 * Look for bad turn-arounds.  Construct list of segments over which
 * scanning looks ok.  Return 0 only if the whole scan looks like one
 * good scan session.
 */

    int turn = 0, pivot = 0;
    int hperiod = (int)(1./2/scan->freq/sampleTime);
    float speed, speed300;
    float delta = scan->az_max - scan->az_min;
    float margin1 = delta*0.08;
    float margin2 = delta*0.01;

    scan->nsection = 0;
    i = 3;
    /* Loop over half-scans... hopefully */
    while (i+hperiod < ndata) {
        while (turn == 0) {
            if (scan->nsection == scan->max_sections) {
                return 5; // Too many scan sections
            }
            speed = (az[i] - az[i-3])/sampleTime/3;
            // Case in middle of the scan
            if ((fabs(fabs(speed) - scan->speed)/scan->speed <  0.1) &&
                (az[i] <= scan->az_max - margin1) &&
                (az[i] >= scan->az_min + margin1)) {  // good speed
                if (speed > 0) turn = 1;
                else turn = -1;
                scan->nsection++;
                scan->section_limits[(scan->nsection-1)*2] = i-3;
                i++;
                // Search for pivot
                while ((speed*turn > 0) && (i < ndata)) {
                    speed = (az[i] - az[i-1])/sampleTime;
                    i++;
                }
                pivot = i-1;
                turn *= -1;
                //printf("Add1: az_pivot = %f, turn = %d\n", az[pivot], turn);
                break;
            }
            // Case in the lower or upper turnarounds
            else if (((az[i] > scan->az_max - margin1) &&
                      (az[i] < scan->az_max + margin2)) ||
                     ((az[i] < scan->az_min + margin1) &&
                      (az[i] > scan->az_min - margin2))) {
                if (i >= ndata - 300) {
                    return 2; // Bad scan end
                }
                speed300 = (az[i+300] - az[i+297])/sampleTime/3;
                if ((((speed300 < 0) && (az[i] > (scan->az_max + scan->az_min)/2)) ||
                     ((speed300 > 0) && (az[i] < (scan->az_max + scan->az_min)/2))) &&
                    (fabs(fabs(speed300) - scan->speed)/scan->speed < 0.1)) { //good turnaround
                    while (speed == 0) {
                        i++;
                        speed = (az[i] - az[i-1])/sampleTime;
                    }
                    scan->nsection++;
                    scan->section_limits[(scan->nsection-1)*2] = i-3;
                    if (az[i] > (scan->az_max + scan->az_min)/2) turn = -1;
                    else turn = 1;

                    // Search for pivot
                    while ((speed*turn < 0) && (i < ndata-300)) {
                        i++;
                        speed = (az[i] - az[i-1])/sampleTime;
                    }
                    pivot = i-1;
                    if (((turn == -1) && 
                         (fabs(az[pivot] - scan->az_max) > margin2)) ||
                        ((turn == 1) &&
                         (fabs(az[pivot] - scan->az_min) > margin2))) {
                        while ((speed*turn > 0) && (i < ndata-300)) {
                            i++;
                            speed = (az[i] - az[i-1])/sampleTime;
                        }
                        turn *= -1;
                    }
                    pivot = i-1;
                    //printf("Add2: speed300 = %f; az = %f; az_m = %f\n", 
                    //     speed300, az[i], (scan->az_max + scan->az_min)/2);
                    break;
                }
                else i++;
            }
            else {
                i++;
            }
            if (i > ndata - hperiod) {
                if (scan->nsection > 1)
                    return 4; //multiple scan sections
                else
                    return 2; // Bad scan end
            }
        }

        //printf("pivot = %d; turn = %d; az = %f\n",pivot, turn, az[pivot+hperiod]);
        i = pivot;
        // Check if we are on a turnarond
        if (((turn == 1) && 
             (fabs(az[pivot+hperiod] - scan->az_max) > margin2)) ||
            ((turn == -1) &&
             (fabs(az[pivot+hperiod] - scan->az_min) > margin2))) { // not good
            //printf("az = %f\n",fabs(az[pivot+hperiod]));
            //printf("turn = %d\n",turn);
            turn = 0;
            scan->section_limits[(scan->nsection-1)*2+1] = i;
            i += hperiod;
        }
        else {
            // Reset the pivot every time to avoid numerical errors
            speed = (az[i] - az[i-1])/sampleTime;
            if (speed*turn < 0) {
                while (speed*turn < 0) {
                    i++;
                    speed = (az[i] - az[i-1])/sampleTime;
                }
                pivot = i-1;
            }
            else {
                while (speed*turn > 0) {
                    i--;
                    speed = (az[i] - az[i-1])/sampleTime;
                }
                pivot = i+1;
            }
            pivot += hperiod;
            i = pivot;
            turn *= -1;
        }
    }
    if (turn != 0)
        scan->section_limits[(scan->nsection-1)*2+1] = ndata-1;

    if (scan->nsection > 1)
        return 4; //multiple scan sections
    else if ((scan->section_limits[0] == 0) && (scan->section_limits[1] == ndata-1))
        return 0;  // Good scan
    else
        return 3;  // Bad beginning of scan
}


PyDoc_STRVAR(analyze_scan__doc__,
             "analyze_scan(az, times)\n"
             "\n"
             "Determine basic properties of a triangle-wave azimuth scan. "
             "Each vector should be a simple ndarray of doubles, with az "
             "in radians and times in seconds.\n"
             "\n"
             "Returns a tuple of items:\n"
             "  (status, speed, scan_freq, az_max, az_min, sections_array)\n"
             "\n"
             "Units are combinations of radians and seconds.  status=0 is "
             "good; otherwise there is something weird with the scan, which "
             "can you can figure out from sections_array.\n");

static PyObject *analyze_scan(PyObject *self, PyObject *args)
{
    PyArrayObject *az_array;
    PyArrayObject *ctime_array;

    if (!PyArg_ParseTuple(args, "O!O!",
                          &PyArray_Type, &az_array,
                          &PyArray_Type, &ctime_array
            ))
        po_raise("invalid arguments.");

    double *az = PyArray_DATA(az_array);
    ASSERT_CARRAY_TYPE_NDIM(az_array, NPY_FLOAT64, 1);
    double *ctime = PyArray_DATA(ctime_array);
    ASSERT_CARRAY_TYPE_NDIM(ctime_array, NPY_FLOAT64, 1);
    int nsamps = PyArray_DIMS(ctime_array)[0];

    mbScan scan;
    scan.min_chunk = 24000;
    scan.min_scan = 0.00873; // Half a degree in azimuth
    scan.max_scan = 0.34907; // 20 degrees in azimuth
    scan.freq = 0.0;
    scan.speed = 0.0;
    scan.az_max = 0.0;
    scan.az_min = 0.0;
    scan.nsection = 0;
    scan.max_sections = 5;
    scan.section_limits = (int*)malloc(2*scan.max_sections *
                                       sizeof(int));
    for (int i = 0; i < 2*scan.max_sections; i++)
        scan.section_limits[i] = -1;

    
    int scan_err = mbAnalyzeScan(az, ctime, nsamps, &scan);

    // Pack up the scan sections in a little array...
    npy_intp dims[2] = {scan.nsection, 2};
    PyArrayObject* sections_array = (PyArrayObject*)
        PyArray_SimpleNew(2, dims, NPY_INT32);
    for (int i=0; i<scan.nsection*2; i++)
        ((int*)PyArray_DATA(sections_array))[i] = scan.section_limits[i];
    free(scan.section_limits);

    /* Return everything in a tuple.  Update the docstring ^^^ if you
     * alter this. */
    return Py_BuildValue("iffffN",
                         scan_err,
                         scan.speed,
                         scan.freq,
                         scan.az_max,
                         scan.az_min,
                         sections_array);
}


/// ------------------------------------------------------------------------
/// Simple search for TOD jumps
static PyObject *find_jumps(PyObject *self, PyObject *args)
{
    PyArrayObject *data_array;
    int dsStep, window;

    if (!PyArg_ParseTuple(args, "O!ii",
                          &PyArray_Type, &data_array,
                          &dsStep, &window
            ))
        po_raise("invalid arguments.");


    // Types and ordering... the ordering is stricter than we really need.
    ASSERT_CARRAY_TYPE_NDIM(data_array, NPY_FLOAT32, 2);
    int ndet = PyArray_DIMS(data_array)[0];
    int nsamps = PyArray_DIMS(data_array)[1];
    int ndata = nsamps/dsStep;

    float *tod_data = PyArray_DATA(data_array);

    npy_intp dims[1];
    dims[0] = ndet;
    PyObject *jump_array = PyArray_SimpleNew(1, dims, NPY_FLOAT32);
    float *jump = PyArray_DATA((PyArrayObject*) jump_array);


    // Downsample data
    float *data = (float *)malloc(ndet * ndata * sizeof(float));
    for (int i = 0; i < ndet; i++) {
        float *this_data = tod_data + nsamps * i;
        for (int j = 0; j < ndata; j++) {
            data[i*ndata + j] = this_data[dsStep * j];
        }
    }

    // Find jumps
    int w = window;
    int nj = ndata-2*w+1;
    float jump_j;
    for (int i = 0; i < ndet; i++) {
        int d = i*ndata;
        jump[i] = 0;
        for (int j = 0; j < nj; j++) {
            jump_j = 0;
            for (int k=0; k<w; k++) jump_j += data[d+j+k] - data[d+j+k+w];
            jump_j /= w;
            if (fabs(jump_j) > jump[i]) jump[i] = fabs(jump_j);
        }
    }
    return jump_array;
}


/// ------------------------------------------------------------------------
/// Perform preselection of detectors depending on correlation with common
/// mode and gain with respect to the common mode. The common mode is
/// obtained using the median of the initial set of candidates.
static PyObject *preselect_dets(PyObject *self, PyObject *args)
{
    PyArrayObject *data_array;
    PyArrayObject *sel_array;
    double min_corr, norm_factor;
    int dsStep;

    if (!PyArg_ParseTuple(args, "O!O!ddi",
                          &PyArray_Type, &data_array,
                          &PyArray_Type, &sel_array,
                          &min_corr, &norm_factor, &dsStep
            ))
        po_raise("invalid arguments.");


    // Types and ordering... the ordering is stricter than we really need.
    ASSERT_CARRAY_TYPE_NDIM(data_array, NPY_FLOAT32, 2);
    int ndet = PyArray_DIMS(data_array)[0];
    int nsamps = PyArray_DIMS(data_array)[1];
    int ndata = nsamps/dsStep;

    // Selections array
    ASSERT_CARRAY_TYPE_NDIM(sel_array, NPY_INT32, 1);
    po_assert(PyArray_DIMS(sel_array)[0] == ndet);

    float *tod_data = PyArray_DATA(data_array);
    int *sel = PyArray_DATA(sel_array);


    // Produce list of detector candidates
    int N = 0;
    for (int i = 0; i < ndet; i++) if (sel[i]) N++;
    int *dets = (int *)malloc(N * sizeof(int));
    int j = 0;
    for (int i = 0; i < ndet; i++) 
        if (sel[i]) {dets[j] = i; j++;}

    // Downsample data and obtain norm of every detector TOD
    int NN = N;
    float *norm = (float *)malloc(N * sizeof(float));
    float *norm2 = (float *)malloc(N * sizeof(float));
    float *data = (float *)malloc(N * ndata * sizeof(float));
    for (int i = 0; i < N; i++) {
        float *this_data = tod_data + nsamps * dets[i];
        norm[i] = 0.0;
        for (int j = 0; j < ndata; j++) {
            data[i*ndata + j] = this_data[dsStep * j];
            norm[i] += this_data[dsStep * j]*this_data[dsStep * j];
        }
        norm[i] = sqrt(norm[i]);
        norm2[i] = norm[i];
    }
    float norm_m = compute_median(N,norm2);
    for (int i = 0; i < N; i++) {
        if ((norm[i]/norm_m > norm_factor) || (norm[i]/norm_m < 1./norm_factor)) {
            sel[dets[i]] = 0;
            NN--;
        }
    }

    // Obtain median common mode
    float cm2 = 0.0;
    float *cm = (float *)malloc(ndata * sizeof(float));
    float *temp = (float *)malloc(N * sizeof(int));
    int k;
    for (int s = 0; s < ndata; s++) {
        k = 0;
        for (int j = 0; j < N; j++) {
            if (sel[dets[j]]) {
                temp[k] = data[j*ndata + s]/norm[j];
                k++;
            }
        }
        cm[s] = compute_median(NN, temp);
        cm2 += cm[s] * cm[s];
    }

    // Select according to correlation and gain
    float *dot = (float *)malloc(N * sizeof(float));
    cblas_sgemv(CblasRowMajor, CblasNoTrans,
                N, ndata, 1.0, data, ndata, cm, 1, 0.0, dot, 1);

    // Select according to correlation first
    double c;
    for (int i = 0; i < N; i++) {
        if (sel[dets[i]]) {
            if (norm[i] != 0.0) c = dot[i]/norm[i]/sqrt(cm2);
            else c = 0.0;
            if (c < min_corr) {
                sel[dets[i]] = 0;
                NN--;
            }
        }
    }

    free(norm);
    free(data);
    free(dets);
    free(cm);
    free(temp);
    free(dot);
    Py_RETURN_NONE;
}


PyDoc_STRVAR(preselect_dets__doc__,
	     "preselect_dets(data, sel, min_corr, norm_factor, dsStep)\n"
    );



PyDoc_STRVAR( get_drift_tests__doc__,
              "get_drift_tests(...)\n"
    );

////////////////////////////////////////////////////////////////////////////
//                            CLAPACK WRAPPER                             //
////////////////////////////////////////////////////////////////////////////
static int sgesdd(char JOBZ, int M, int N, float *A, int LDA, float *S, float *U,
           int LDU, float *VT, int LDVT, float *WORK, int LWORK, int *IWORK)
{
    extern void sgesdd_(char *JOBZp, int *Mp, int *Np, float *A, int *LDAp,
                        float *S, float *U, int *LDUp, float *VT, int *LDVTp,
                        float *WORK, int *LWORKp, int *IWORK, int *INFOp);
    int INFO;
    sgesdd_(&JOBZ, &M, &N, A, &LDA, S, U, &LDU, VT, &LDVT, WORK, &LWORK, IWORK,
            &INFO);
    return INFO;
}

static PyObject *get_drift_tests(PyObject *self, PyObject *args)
{
    PyArrayObject *tod_array;
    PyArrayObject *row_array;
    PyArrayObject *col_array;
    PyArrayObject *ldet_array;
    PyArrayObject *ddet_array;
    PyArrayObject *temp_array;
    
    int dsStep, nBlock, nModes, removeDarkCM;

    if (!PyArg_ParseTuple(args, "O!O!O!O!O!O!iiii",
                          &PyArray_Type, &tod_array,
                          &PyArray_Type, &row_array,
                          &PyArray_Type, &col_array,
                          &PyArray_Type, &ldet_array,
                          &PyArray_Type, &ddet_array,
                          &PyArray_Type, &temp_array,
                          &dsStep,
                          &nBlock,
                          &nModes,
                          &removeDarkCM
            ))
        po_raise("invalid arguments.");

    // Check TOD data...
    ASSERT_CARRAY_TYPE_NDIM(tod_array, NPY_FLOAT32, 2);
    int ndet = PyArray_DIMS(tod_array)[0];
    int tod_ndata = PyArray_DIMS(tod_array)[1];
    float *tod_data = PyArray_DATA(tod_array);

    // Check row and col
    ASSERT_CARRAY_TYPE_NDIM(row_array, NPY_INT32, 1);
    ASSERT_CARRAY_TYPE_NDIM(col_array, NPY_INT32, 1);
    int *tod_rows = PyArray_DATA(row_array);
    int *tod_cols = PyArray_DATA(col_array);
    po_assert(PyArray_DIMS(row_array)[0] == ndet);
    po_assert(PyArray_DIMS(col_array)[0] == ndet);
    
    // Check det lists...
    ASSERT_CARRAY_TYPE_NDIM(ldet_array, NPY_INT32, 1);
    ASSERT_CARRAY_TYPE_NDIM(ddet_array, NPY_INT32, 1);
    int *ldet = PyArray_DATA(ldet_array);
    int nldet = PyArray_DIMS(ldet_array)[0];
    int *ddet = PyArray_DATA(ddet_array);
    int nddet = PyArray_DIMS(ddet_array)[0];

    // Check temperature vector
    ASSERT_CARRAY_TYPE_NDIM(temp_array, NPY_FLOAT32, 1);
    po_assert(PyArray_DIMS(temp_array)[0] == tod_ndata);
    int ndata = tod_ndata / dsStep;
    float *temp = PyArray_DATA(temp_array);
    bool valid_temp = false;
    if ((temp[0] != temp[ndata/2]) || (temp[0] != temp[ndata/3])) valid_temp = true;

    // Downsampling -- ndata is the downsampled size.

    // Sanity checks
    assert((nBlock == 1) || (nBlock == 4) || (nBlock = 9) || (nBlock == 16));
    assert(ndet > 0);
    assert(ndata > 0);
    if (nModes > nBlock)
        nModes = nBlock;

    // Ready output arrays.
    npy_intp dims[2];
    /* Time-like... */
    dims[0] = ndata;
    PyObject *cm_array = PyArray_SimpleNew(1, dims, NPY_FLOAT32);
    PyObject *dcm_array = PyArray_SimpleNew(1, dims, NPY_FLOAT32);
    float *cm = PyArray_DATA((PyArrayObject*) cm_array);
    float *dcm = PyArray_DATA((PyArrayObject*) dcm_array);
    /* Det-like... */
    dims[0] = ndet;
    /* Template:
       PyObject *X_array = PyArray_SimpleNew(1, dims, NPY_FLOAT64);
       double *X         = PyArray_DATA((PyArrayObject*) X_array);
    */
    PyObject *norm_array     = PyArray_SimpleNew(1, dims, NPY_FLOAT64);
    double   *norm           = PyArray_DATA((PyArrayObject*) norm_array);
    PyObject *corrLive_array = PyArray_SimpleNew(1, dims, NPY_FLOAT64);
    double   *corrLive       = PyArray_DATA((PyArrayObject*) corrLive_array);
    PyObject *DELive_array   = PyArray_SimpleNew(1, dims, NPY_FLOAT64);
    double   *DELive         = PyArray_DATA((PyArrayObject*) DELive_array);
    PyObject *gainLive_array = PyArray_SimpleNew(1, dims, NPY_FLOAT64);
    double   *gainLive       = PyArray_DATA((PyArrayObject*) gainLive_array);
    PyObject *corrDark_array = PyArray_SimpleNew(1, dims, NPY_FLOAT64);
    double   *corrDark       = PyArray_DATA((PyArrayObject*) corrDark_array);
    PyObject *DEDark_array   = PyArray_SimpleNew(1, dims, NPY_FLOAT64);
    double   *DEDark         = PyArray_DATA((PyArrayObject*) DEDark_array);
    PyObject *gainDark_array = PyArray_SimpleNew(1, dims, NPY_FLOAT64);
    double   *gainDark       = PyArray_DATA((PyArrayObject*) gainDark_array);
    
    /* Resume original moby code... */
  
  int side = (int)sqrt((float)nBlock);
  int nSide = (int)ceil(32.0/side);
  int rmin, cmin, N, zeros = 0;

  int *bdet = (int *)malloc(ndet*sizeof(int));
  int *rows = (int *)malloc(nldet*sizeof(int));
  int *cols = (int *)malloc(nldet*sizeof(int));
  for (int i = 0; i < nldet; i++) {
    rows[i] = tod_rows[ldet[i]];
    cols[i] = tod_cols[ldet[i]];
  }
  
  // Downsample data for every detector TOD
  float *data = (float *)malloc(ndet * ndata * sizeof(float));
  for (int i = 0; i < ndet; i++) {
    double mean = 0.;
    for (int j = 0; j < ndata; j++) {
      data[i*ndata + j] = tod_data[i*tod_ndata + dsStep * j];
      mean += data[i*ndata+j];
    }
    mean /= ndata;
    for (int j = 0; j < ndata; j++) {
      data[i*ndata + j] -= mean;
    }
  }

  float *dot = (float *)malloc(ndet * sizeof(float));
  if (valid_temp) { 
    // Downsample temperature
    float *dstemp = (float *)malloc(ndata * sizeof(float));
    double mean = 0;
    double temp2 = 0;
    for (int j = 0; j < ndata; j++) {
      dstemp[j] = temp[dsStep * j];
      mean += dstemp[j];
    }
    // Remove mean and obtain norm^2 of temp
    mean /= ndata;
    for (int j = 0; j < ndata; j++) {
      dstemp[j] -= mean;
      temp2 += dstemp[j] * dstemp[j];
    }
    if (temp2==0.)
        temp2 = 1.;
 
    // Obtain fit temperature to data
    cblas_sgemv(CblasRowMajor, CblasNoTrans, ndet, ndata, 1.0, data, ndata,
                dstemp, 1, 0.0, dot, 1);
 
    // Remove temperature signal
//#pragma omp parallel shared(ndet, ndata, data, dot, dstemp)
    {
//#pragma omp for
      for (int d = 0; d < ndet; d++) {
        for (int i = 0; i < ndata; i++) data[d*ndata + i] -= dot[d]*dstemp[i]/temp2;
      }
    }
    free(dstemp);
  }


  // Obtain dark common mode
  float cm2 = 0;
  for (int s = 0; s < ndata; s++) {
    dcm[s] = 0;
    for (int j = 0; j < nddet; j++) dcm[s] += data[ddet[j]*ndata + s]/nddet;
    cm2 += dcm[s] * dcm[s];
  }

  // Obtain dark correlation and gain
  cblas_sgemv(CblasRowMajor, CblasNoTrans, ndet, ndata, 1.0, data, ndata,
              dcm, 1, 0.0, dot, 1);
  for (int i = 0; i < ndet; i++) {
    norm[i] = 0.0;
    for (int j = 0; j < ndata; j++)  norm[i] += data[i*ndata+j]*data[i*ndata+j];
    norm[i] = sqrt(norm[i]/ndata);
    if (norm[i] != 0.0) corrDark[i] = dot[i]/norm[i]/sqrt(cm2*ndata);
    else corrDark[i] = 0.0;
    gainDark[i] = dot[i]/cm2;
  }

  // Obtain dark drift error
//#pragma omp parallel shared(DEDark, ndet, ndata, data, gainDark)
  {
//#pragma omp for
    for (int d = 0; d < ndet; d++) {
      float dat;
      // Remove dark common mode
      DEDark[d] = 0.0;
      for (int i = 0; i < ndata; i++) {
        dat = data[d*ndata + i] - gainDark[d]*dcm[i];
        DEDark[d] += dat * dat / ndata;
        if (removeDarkCM==1) data[d*ndata + i] = dat;
      }
      DEDark[d] = sqrt(DEDark[d]);
      if (norm[d] != 0.0) DEDark[d] /= norm[d];
    }
  }
  //Note that DEDark here constitutes, for live detectors, the norm of the 
  //data minus the dark common mode

  // Remove slope
  int Ns = ((ndata+1)/2)*2 - 1;
  float Ns2 = (float)(Ns/2);
  long int sx2 = (long int)(Ns2*Ns2*Ns2/3 + Ns2*Ns2/2 + Ns2/6)*2;
  double a, A;
  A = 0;
  for (int j = 0; j < ndet; j++) {
    a = 0;
    for (int s = 0; s < Ns; s++) a += (s - Ns2)*data[j*ndata + s];
    for (int d = 0; d < nldet; d++) if (j == ldet[d]) A += a;
    for (int s = 0; s < ndata; s++) data[j*ndata + s] -= (s - Ns2)*a/sx2;
  }
  A /= nldet;

  // Obtain live common mode
  cm2 = 0;
  for (int s = 0; s < ndata; s++) {
    cm[s] = 0;
    for (int j = 0; j < nldet; j++) cm[s] += data[ldet[j]*ndata + s]/nldet;
    cm2 += cm[s] * cm[s];
  }
  // Obtain live correlation and gain
  cblas_sgemv(CblasRowMajor, CblasNoTrans, ndet, ndata, 1.0, data, ndata,
              cm, 1, 0.0, dot, 1);
  for (int d = 0; d < ndet; d++) {
    norm[d] = 0.0;
    for (int j = 0; j < ndata; j++)  norm[d] += data[d*ndata+j]*data[d*ndata+j];
    norm[d] = sqrt(norm[d]/ndata);
    if (norm[d] != 0.0) corrLive[d] = dot[d]/norm[d]/sqrt(cm2*ndata);
    else corrLive[d] = 0.0;
    gainLive[d] = dot[d]/cm2;
  }
  
  // Obtain block common modes, weighted by the number of detectors in
  // each block.
  float *cmodes = (float *)malloc(nBlock * ndata * sizeof(float));
  for (int i = 0; i < nBlock; i++) {
    rmin = nSide * (i/side);
    cmin = nSide*i - rmin*side; 
    N = 0;
    for (int j = 0; j < nldet; j++) {
      if (rows[j] >= rmin && rows[j] < rmin + nSide &&
          cols[j] >= cmin && cols[j] < cmin + nSide) {
        bdet[N] = ldet[j];
        N++;
      }
    }
    if (N > 0) {
      for (int s = 0; s < ndata; s++) {
        cmodes[(i-zeros)*ndata + s] = 0;
        for (int j = 0; j < N; j++) 
          cmodes[(i-zeros)*ndata + s] += data[bdet[j]*ndata + s];
      }
    } 
    else {
      zeros++;
    }
  }
  nBlock -= zeros;
  if (nModes > nBlock) nModes = nBlock;

  // Obtain the SVD
    // A (M-by-N) with M rows and N columns
    // JOBZ = 'S'
    // M = ndata
    // N = nBlock
    // LDA = ndata
  float *S = (float *)malloc(nBlock * sizeof(float));
  float *U = (float *)malloc(nBlock * ndata * sizeof(float));
  float *VT = (float *)malloc(nBlock * nBlock * sizeof(float));
  int LWORK = nBlock * ndata;
  float *WORK = (float *)malloc(LWORK * sizeof(float));
  int *IWORK = (int *)malloc(8 * nBlock * sizeof(int));
  // Apply Lapack SVD
  sgesdd('S', ndata, nBlock, cmodes, ndata, S, U, ndata, VT, nBlock, 
         WORK, LWORK, IWORK);

  // Find fit coefficients (mT A). 
  // Note that matrices are transposed in this definition
  float *coeff = (float *)malloc(nBlock * ndet * sizeof(float));
  cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasTrans, ndet, nBlock, ndata, 1.0, 
              data, ndata, U, ndata, 0.0, coeff, nBlock);


  // Rearrange coeff matrix to contain only nModes modes
  for (int i = 0; i < ndet; i++) {
    for (int j = 0; j < nModes; j ++) {
      coeff[i*nModes + j] = coeff[i*nBlock + j];
    }
  }

  // Remove modes (A - m mT A)
  cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, ndet, ndata, nModes, -1.0, 
              coeff, nModes, U, ndata, 1.0, data, ndata);

  // Obtain drift error
//#pragma omp parallel shared(DELive, ndet, ndata, data)
  {
//#pragma omp for
    for (int d = 0; d < ndet; d++) {
      DELive[d] = 0.0;
      for (int i = 0; i < ndata; i++) {
        DELive[d] += data[d*ndata + i] * data[d*ndata + i] / ndata;
      }
      if (norm[d] != 0.0) DELive[d] = sqrt(DELive[d])/norm[d];
    }
  }
  //Note that we have normalized the DriftError to cancel the effect of the gain in 
  //the value, because we want to represent differences in shape with respect to the
  //common mode

  //Re-add the slope to the common-mode
  for (int s = 0; s < ndata; s++) cm[s] += (s - Ns2)*A/sx2;

  free(bdet);
  free(rows);
  free(cols);
  free(data);
  free(dot);
  free(cmodes);
  free(S);
  free(U);
  free(VT);
  free(WORK);
  free(IWORK);
  free(coeff);

    return Py_BuildValue("NNNNNNNNN",
                         cm_array,
                         dcm_array,
                         norm_array,
                         corrLive_array,
                         DELive_array,
                         gainLive_array,
                         corrDark_array,
                         DEDark_array,
                         gainDark_array);

}

PyMethodDef pyactpol_filter_methods[] = {
    {"remove_mean", remove_mean, METH_VARARGS, ""},
    {"filter_tod_data", filter_tod_data, METH_VARARGS,
    filter_tod_data__doc__},
    {"get_glitch_cuts", get_glitch_cuts, METH_VARARGS, ""},
    {"analyze_scan", analyze_scan, METH_VARARGS,
     analyze_scan__doc__},
    {"preselect_dets", preselect_dets, METH_VARARGS,
     preselect_dets__doc__},
    {"find_jumps", find_jumps, METH_VARARGS, ""},
    {"get_drift_tests", get_drift_tests, METH_VARARGS,
     get_drift_tests__doc__},
    {NULL, NULL, 0, NULL}        /* Sentinel */
};


