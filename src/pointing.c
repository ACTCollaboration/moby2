/* -*- mode: C; tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 4 -*-
 *      vim: sw=4 ts=8 et tw=80
 */

#include "pyactpol.h"
#include "myassert.h"

#include <actpol/actpol.h>

#ifndef LIBACTPOL_API
#define LIBACTPOL_API 20130101
#warning "LIBACTPOL_API was not defined."
#endif

#if LIBACTPOL_API <= 20180101
#define FC_RA(fc) fc->ra
#define FC_DEC(fc) fc->dec
#else
#define FC_RA(fc) fc->a
#define FC_DEC(fc) fc->b
#endif

#define REFRAC_UPDATE (0.2 * M_PI / 180)

struct string_table {
    int code;
    char *text;
};

/*
 * libactpol computes a bunch of different numbers for each feedhorn
 * pointing.  The enum and strings below are used to request
 * particular coordinates to be returned.
 */

enum pointing_field_id {
    POINTING_RA,
    POINTING_DEC,
    POINTING_GAMMA,
    POINTING_2GAMMA,
    POINTING_COS_DEC,
    POINTING_SIN_DEC,
    POINTING_COS_GAMMA,
    POINTING_SIN_GAMMA,
    POINTING_COS_2GAMMA,
    POINTING_SIN_2GAMMA,
    POINTING_R_X,
    POINTING_R_Y,
    POINTING_R_Z
};

static struct string_table pointing_field_names[] = {
    {POINTING_RA, "ra"},
    {POINTING_DEC, "dec"},
    {POINTING_GAMMA, "gamma"},
    {POINTING_2GAMMA, "2gamma"},
    {POINTING_COS_DEC, "cos_dec"},
    {POINTING_SIN_DEC, "sin_dec"},
    {POINTING_COS_GAMMA, "cos_gamma"},
    {POINTING_SIN_GAMMA, "sin_gamma"},
    {POINTING_COS_2GAMMA, "cos_2gamma"},
    {POINTING_SIN_2GAMMA, "sin_2gamma"},
    {POINTING_R_X, "x"},
    {POINTING_R_Y, "y"},
    {POINTING_R_Z, "z"},
    {-1, NULL}
};


/*
 * projection_class enumerates the projections for which there exist
 * accelerated wand routines.
 */

enum projection_class_id {
    PROJECTION_TANGENT,
    PROJECTION_RADEC,
    PROJECTION_RASINDEC
};

static struct string_table projection_class_names[] = {
    {PROJECTION_TANGENT, "tangent"},
    {PROJECTION_RADEC, "ra_dec"},
    {PROJECTION_RASINDEC, "ra_sindec"},
    {-1, NULL},
};

static int get_string_code(const char* text, struct string_table *table)
{
    while (table->text != NULL) {
        if (strcmp(table->text, text)==0)
            return table->code;
        table++;
    }
    return -1;
}


/*
 * Assistance functions for unpacking common Python objects
 *
 * These can be used in PyArgs_ParseTuple with the "O&" field code.
 */

static int extract_double_data(PyObject *o, void **dest)
{
    if (!PyArray_Check(o)) {
        print_error("object is not an array\n");
        return 0;
    }
    if (!IS_ARRAY_SIMPLE_1D_TYPE((PyArrayObject*)o, NPY_FLOAT64)) {
        print_error("array has unexpected shape, packing, or dtype\n");
        return 0;
    }
    *dest = (double*)PyArray_DATA((PyArrayObject*)o);
    return 1;
}


static int extract_double_data_or_null(PyObject *o, void **dest)
{
    if (o == Py_None) {
        *dest = NULL;
        return 1;
    }
    return extract_double_data(o, dest);
}


/*
 * wand_data is a local container for WandData
 */

struct wand_data {
    enum projection_class_id proj_id;
    int polarized;
    int has_hwp;
    int n;
    const double *ra;
    const double *dec;
    const double *cosdec;
    const double *sindec;
    const double *cosg;
    const double *sing;
    const double *x;
    const double *y;
    const double *cos2g;
    const double *sin2g;
    const double *cos_2hwp;
    const double *sin_2hwp;
};




/*
 * pixelization here refers to a mapping between an infinite 2-d
 * cartesian coordinate system and a finite 2-d grid of pixels.  These
 * are analagous to the FITS / wcslib "projection plane coordaintes"
 * and "pixel coordinates", respectively, except we don't support
 * arbitrary linear transformations to get from the former to the
 * latter.  Rotations are supported by the libactpol pointing system,
 * where you can drop in an arbitary rotation as part of the
 * transformation to the projection plane coordinates.  You can't do
 * skewness transformations in this system.  But why do you want to?
 */

/* pixelization_data is a local container for SimpleProjection. */

struct pixelization_data {
    int n0;
    int n1;
    int pix0; 
    int pix1;
    double x0;
    double x1;
    double dx0;
    double dx1;
};

static int unpack_pixelization(struct pixelization_data *pix,
                               PyObject *pypix)
{
    if (!PyArg_ParseTuple(pypix, "iiiidddd",
                          &pix->n0, &pix->n1,
                          &pix->pix0, &pix->pix1,
                          &pix->x0, &pix->x1,
                          &pix->dx0, &pix->dx1))
        return 0;
    return 1;
}


static int unpack_wand(struct wand_data *wand, PyObject *pywand)
{
    /* This must be a sequence, with at least one item: a string that
     * describes the wand type. */
    if (!PySequence_Check(pywand) || 
        PySequence_Length(pywand) < 1) {
        print_error("not a valid wand tuple\n");
        return 0;
    }

    // Steal a reference to first item (projection type)
    PyObject *proj_obj = PySequence_GetItem(pywand, 0);
    Py_DECREF(proj_obj);
    if (!PyString_Check(proj_obj)) {
        print_error("first argument is not a string\n");
        return 0;
    }

    const char *proj_type = PyString_AsString(proj_obj);
    wand->proj_id = get_string_code(proj_type, projection_class_names);

    if (wand->proj_id < 0) {
        print_error("invalid format string, %s\n", proj_type);
        return 0;
    }

    int fields_ok = 0;
    switch(wand->proj_id) {
    case PROJECTION_TANGENT:
        fields_ok = PyArg_ParseTuple(pywand, "siiO&O&O&O&O&O&O&O&",
                                     &proj_type,
                                     &wand->polarized,
                                     &wand->has_hwp,
                                     extract_double_data, &wand->cosg,
                                     extract_double_data, &wand->sing,
                                     extract_double_data_or_null, &wand->cos2g,
                                     extract_double_data_or_null, &wand->sin2g,
                                     extract_double_data, &wand->x,
                                     extract_double_data, &wand->y,
                                     extract_double_data_or_null, &wand->cos_2hwp,
                                     extract_double_data_or_null, &wand->sin_2hwp);
        break;
    case PROJECTION_RADEC:
    case PROJECTION_RASINDEC:
        fields_ok = PyArg_ParseTuple(pywand, "siiO&O&O&O&O&O&O&O&O&O&",
                                     &proj_type,
                                     &wand->polarized,
                                     &wand->has_hwp,
                                     extract_double_data, &wand->cosg,
                                     extract_double_data, &wand->sing,
                                     extract_double_data_or_null, &wand->cos2g,
                                     extract_double_data_or_null, &wand->sin2g,
                                     extract_double_data, &wand->ra,
                                     extract_double_data, &wand->dec,
                                     extract_double_data, &wand->cosdec,
                                     extract_double_data, &wand->sindec,
                                     extract_double_data_or_null, &wand->cos_2hwp,
                                     extract_double_data_or_null, &wand->sin_2hwp);
        break;
    }        
    if (!fields_ok) {
        print_error("parsing of %i items failed\n", PySequence_Size(pywand));
        return 0;
    }

    // Oh, and the length.
    PyObject *item = PySequence_GetItem(pywand, 3);
    Py_DECREF(item);
    if (!PyArray_Check(item)) {
        print_error("first vector argument in wand tuple is not an ndarray\n");
        return 0;
    }
    wand->n = PyArray_SIZE((PyArrayObject*)item);

    if ((wand->cos_2hwp == NULL) != (wand->sin_2hwp == NULL)) {
        print_error("HWP data is invalid.\n");
        return 0;
    }

    return 1;
}


/*
 * focal_plane_data is a local container for the FocalPlane object.
 */

struct focal_plane_data {
    int n;
    double *x;
    double *y;
    double *phi;
};

static int unpack_focal_plane(struct focal_plane_data *fplane,
                              PyObject *pyfplane)
{
    if (!PyArg_ParseTuple(pyfplane, "O&O&O&",
                          extract_double_data, &fplane->x,
                          extract_double_data, &fplane->y,
                          extract_double_data, &fplane->phi))
        return 0;

    // Oh, and the length.
    PyObject *item = PySequence_GetItem(pyfplane, 0);
    fplane->n = PyArray_SIZE((PyArrayObject*)item);
    Py_DECREF(item);

    return 1;
}


/* Unpacking of ACTpolWeather objects. */

static int unpack_weather(ACTpolWeather *weather,
                          PyObject *pyweather)
{
    if (!PyArg_ParseTuple(pyweather, "dddd",
                          &weather->temperature_C,
                          &weather->pressure_mbar,
                          &weather->relative_humidity,
                          &weather->tropospheric_lapse_rate_K_per_m
            )) {
        print_error("Failed to parse 4 doubles from tuple.\n");
        return 0;
    }
    return 1;
}

/*
static PyObject *pack_weather(ACTpolWeather *weather)
{
    return Py_BuildValue("dddd",
                         weather->temperature_C,
                         weather->pressure_mbar,
                         weather->relative_humidity,
                         weather->tropospheric_lapse_rate_K_per_m);
}
*/


/* Unpacking of encoded TODCuts data
 *
 * Note that this requires local storage, so the destructor should be
 * called before returning to python.
 */

int pyactpol_unpack_cuts_data(struct cuts_set_data *cuts_set,
                                     PyObject *pycutsset)
{
    // pycutsset is a sequence of numpy arrays.
    if (!PySequence_Check(pycutsset))
        return 0;
    /* For internal use only... steal a reference to each numpy array
     * and store the pointer to its data.  This might not be legal. */
    cuts_set->n = PySequence_Length(pycutsset);
    cuts_set->ncut = malloc(cuts_set->n * sizeof(*cuts_set->ncut));
    cuts_set->cuts = malloc(cuts_set->n * sizeof(*cuts_set->cuts));
    int i;
    for (i=0; i<cuts_set->n; i++) {
        PyObject *o = PySequence_GetItem(pycutsset, i);
        Py_DECREF(o); // Whatever.  We own you, sequence of arrays.
        PyArrayObject *a = (PyArrayObject*)o;
        cuts_set->ncut[i] = PyArray_SIZE(a) / 2;
        cuts_set->cuts[i] = PyArray_DATA(a);
    }
    return 1;
}

void pyactpol_destroy_cuts_data(struct cuts_set_data *cuts_set)
{
    free(cuts_set->ncut);
    free(cuts_set->cuts);
}


/******
 * Main pointing code.  Wrapped calls to libactpol to get exact
 * mapping between celestial sky and focal plane.
 *
 * From the python level call:
 * - get_coords : take focal plane coords to equatorial coords
 * - get_coords_inverse : the inverse
 *
 */


/* get_misc_pointing
 *
 * This is broken off from get_coords to allow for more specialized
 * alternatives, and to have more, shorter functions...
 */

static PyObject *get_misc_pointing(ACTpolState *state,
                                   ACTpolArrayCoords *coords,
                                   ACTpolWeather *weather,
                                   long nsamps,
                                   int nfields,
                                   enum pointing_field_id *fields,
                                   double *ctime,
                                   double *az,
                                   double *alt,
                                   int one_d
    )
{
    int nhorns = coords->array->nhorns;
    npy_intp dims[2];
    int nd;
    if (one_d) {
        po_assert(nhorns == 1);
        dims[0] = nsamps;
        nd = 1;
    } else {
        dims[0] = nhorns;
        dims[1] = nsamps;
        nd = 2;
    }

    PyArrayObject **output_arrays = malloc(nfields * sizeof(void*));
    double **output = malloc(nfields * sizeof(double*));
    
    for (int i=0; i<nfields; i++) {
        output_arrays[i] = (PyArrayObject*)
            PyArray_SimpleNew(nd, dims, NPY_FLOAT64);
        output[i] = PyArray_DATA(output_arrays[i]);
    }

    double ref_alt = alt[0];

    for (int j=0; j<nsamps; j++) {
        // Use ctime=0 as a "no data" indicator.
        if (ctime[j] == 0.) {
            for (int k=0; k<nhorns; k++)
                for (int i=0; i<nfields; i++)
                    *(output[i] + k*nsamps + j) = 0;
            continue;
        }

        if (fabs(ref_alt - alt[j]) > REFRAC_UPDATE) {
            double throw = .1;
            ACTpolScan scan;
            ACTpolScan_init(&scan, alt[j], az[j], throw);
            ACTpolArrayCoords_update_refraction(coords, &scan, weather);
            ref_alt = alt[j];
        }

        ACTpolState_update(state, ctime[j], alt[j], az[j]);
        ACTpolArrayCoords_update(coords, state);
        

        for (int k=0; k<nhorns; k++) {
            ACTpolFeedhornCoords *fc = coords->horn + k;

            for (int i=0; i<nfields; i++) {
                double *dest = output[i] + k*nsamps + j;
                switch(fields[i]) {
                case POINTING_RA:
                    *dest = FC_RA(fc);
                    break;
                case POINTING_DEC:
                    *dest = FC_DEC(fc);
                    break;
                case POINTING_GAMMA:
                    *dest = atan2(fc->singamma, fc->cosgamma);
                    break;
                case POINTING_2GAMMA:
                    *dest = atan2(fc->sin2gamma, fc->cos2gamma);
                    break;
                case POINTING_COS_DEC:
                    *dest = cos(FC_DEC(fc));
                    break;
                case POINTING_SIN_DEC:
                    *dest = sin(FC_DEC(fc));
                    break;
                case POINTING_COS_GAMMA:
                    *dest = fc->cosgamma;
                    break;
                case POINTING_SIN_GAMMA:
                    *dest = fc->singamma;
                    break;
                case POINTING_COS_2GAMMA:
                    *dest = fc->cos2gamma;
                    break;
                case POINTING_SIN_2GAMMA:
                    *dest = fc->sin2gamma;
                    break;
                case POINTING_R_X:
                    *dest = fc->m[2][0];
                    break;
                case POINTING_R_Y:
                    *dest = fc->m[2][1];
                    break;
                case POINTING_R_Z:
                    *dest = fc->m[2][2];
                    break;
                }
            }
        }
    }

    PyObject* tuple_out = PyTuple_New(nfields);
    for (int i=0; i<nfields; i++)
        PyTuple_SET_ITEM(tuple_out, i, (PyObject*)output_arrays[i]);

    free(output_arrays);
    free(output);

    return tuple_out;
}

/* There's probably not much point in this.  If you're just getting
 * ra,dec that will be fast anyway.  And a switch won't be the bottle
 * neck. */

static PyObject *get_ra_dec(ACTpolState *state,
                            ACTpolArrayCoords *coords,
                            ACTpolWeather *weather,
                            long nsamps,
                            double *ctime,
                            double *alt,
                            double *az,
                            int one_d
    )
{
    int nhorns = coords->array->nhorns;
    npy_intp dims[2];
    int nd;
    if (one_d) {
        po_assert(nhorns == 1);
        dims[0] = nsamps;
        nd = 1;
    } else {
        dims[0] = nhorns;
        dims[1] = nsamps;
        nd = 2;
    }
    PyArrayObject *ra_array = (PyArrayObject*)
        PyArray_SimpleNew(nd, dims, NPY_FLOAT64);
    PyArrayObject *dec_array = (PyArrayObject*)
        PyArray_SimpleNew(nd, dims, NPY_FLOAT64);

    double *ra = PyArray_DATA(ra_array);
    double *dec = PyArray_DATA(dec_array);

    int j,k;
    double ref_alt = alt[0];

    for (j=0; j<nsamps; j++) {
        if (ctime[j] == 0) {
            for (k=0; k<nhorns; k++) {
                ra[k*nsamps+j] = 0.;
                dec[k*nsamps+j] = 0.;
            }
            continue;
        }

        if (fabs(ref_alt - alt[j]) > REFRAC_UPDATE) {
            double throw = .1;
            ACTpolScan scan;
            ACTpolScan_init(&scan, alt[j], az[j], throw);
            ACTpolArrayCoords_update_refraction(coords, &scan, weather);
            ref_alt = alt[j];
        }

        ACTpolState_update(state, ctime[j], alt[j], az[j]);
        ACTpolArrayCoords_update(coords, state);
        
        for (k=0; k<nhorns; k++) {
            ACTpolFeedhornCoords *fc = coords->horn + k;
            ra[k*nsamps+j] = FC_RA(fc);
            dec[k*nsamps+j] = FC_DEC(fc);
        }
    }

    return Py_BuildValue("NN",
                         ra_array, dec_array);
}


static PyObject *get_coords(PyObject *self, PyObject *args, PyObject *kwargs)
{
    PyArrayObject *ctime_array;
    PyArrayObject *az_array;
    PyArrayObject *alt_array;
    PyObject *field_strings = Py_None;
    PyObject *fplane_data = Py_None;
    PyObject *final_rotation = Py_None;
    PyObject *weather_data = Py_None;

    char *kwnames[] = {"ctime",
                       "az",
                       "alt",
                       "fields",
                       "focal_plane",
                       "final_q",
                       "weather",
                       NULL};


    if (!PyArg_ParseTupleAndKeywords(
            args, kwargs,
            "O!O!O!|OOOO", kwnames,
            &PyArray_Type, &ctime_array,
            &PyArray_Type, &az_array,
            &PyArray_Type, &alt_array,
            &field_strings,
            &fplane_data,
            &final_rotation,
            &weather_data
            ))
        return NULL;

    /*
     * Check types and stuff of input ctime, alt, az arrays.
     */
    ASSERT_CARRAY_TYPE_NDIM(ctime_array, NPY_FLOAT64, 1);
    ASSERT_CARRAY_TYPE_NDIM(az_array, NPY_FLOAT64, 1);
    ASSERT_CARRAY_TYPE_NDIM(alt_array, NPY_FLOAT64, 1);

    long n = PyArray_DIMS(ctime_array)[0];
    ASSERT_AXIS_N(az_array, 0, n);
    ASSERT_AXIS_N(alt_array, 0, n);
    
    double *ctime = PyArray_DATA(ctime_array);
    double *az = PyArray_DATA(az_array);
    double *alt = PyArray_DATA(alt_array);

    /*
     * Decode pointing field instructions.
     */
    int nfields = 0;
    enum pointing_field_id *fields = NULL;
    const char *fast_field = NULL;
    if (field_strings == Py_None) {
        fast_field = "ra_dec";
    } else if (PyString_Check(field_strings)) {
        fast_field = PyString_AsString(field_strings);
    } else if (PySequence_Check(field_strings)) {
        nfields = PySequence_Length(field_strings);
        fields = malloc(nfields*sizeof(*fields));
        int i;
        for (i=0; i<nfields; i++) {
            PyObject *o = PySequence_GetItem(field_strings, i);
            Py_DECREF(o);
            po_assert(PyString_Check(o));
            fields[i] = get_string_code(PyString_AsString(o), pointing_field_names);
        }
    } else {
        po_raise("fields is neither string nor list");
    }

    int one_d = (fplane_data == Py_None);
    double dummy_horn = 0.;
    double dummy_pol = M_PI/2;
    struct focal_plane_data fplane = {.n = 1,
                                      .x = &dummy_horn,
                                      .y = &dummy_horn,
                                      .phi = &dummy_pol };
    if (!one_d && !unpack_focal_plane(&fplane, fplane_data)) {
        po_raise("invalid focal plane object.");
    }

    /* Unpack the final rotation quaternion */
    Quaternion *BCRS_to_NSC = NULL;
    if (final_rotation != Py_None) {
        BCRS_to_NSC = PyArray_DATA((PyArrayObject*)final_rotation);
    }

    /* Unpack any weather */
    ACTpolWeather weather;
    ACTpolWeather_default(&weather);
    if (weather_data != Py_None && 
        !unpack_weather(&weather, weather_data))
        po_raise("invalid weather specifier.");

    ACTpolArray *array = ACTpolArray_alloc(fplane.n);
    ACTpolArray_init(array, 150., 0., 0.);
    for (int i=0; i<fplane.n; i++) {
//        printf("%i %lf %lf\n", i, fplane.y[i], fplane.x[i]);
        ACTpolFeedhorn_init(array->horn+i, fplane.y[i], fplane.x[i],
                            -M_PI/2 + fplane.phi[i]);
    }
        
    ACTpolArrayCoords *coords = ACTpolArrayCoords_alloc(array);
#if LIBACTPOL_API < 20180101
    ACTpolArrayCoords_init(coords);
#else
    ACTpolArrayCoords_init(coords, ACTPOL_COORDSYS_GENERAL);
#endif    

    ACTpolState *state = ACTpolState_alloc();
    ACTpolState_init(state);
    if (BCRS_to_NSC != NULL)
        memcpy(state->BCRS_to_NSC_q, BCRS_to_NSC, sizeof(*BCRS_to_NSC));

    double throw = .1;
    ACTpolScan scan;
    ACTpolScan_init(&scan, alt[0], az[0], throw);

    ACTpolArrayCoords_update_refraction(coords, &scan, &weather);

    PyObject *output = NULL;

    if (fields != NULL) {
        output = get_misc_pointing(state, coords, &weather,
                                   n, nfields, fields,
                                   ctime, az, alt, one_d);
        free(fields);
    } else if (strcmp(fast_field, "ra_dec")==0) {
        output = get_ra_dec(state, coords, &weather, n, ctime, alt, az, one_d);
    }

    ACTpolState_free(state);
    ACTpolArrayCoords_free(coords);
    ACTpolArray_free(array);

    if (output != NULL)
        return output;

    po_raise("invalid field specification");
}


static PyObject *get_coords_inverse(PyObject *self, PyObject *args, PyObject *kwargs)
{
    PyArrayObject *ctime_array;
    PyArrayObject *az_array;
    PyArrayObject *alt_array;
    PyArrayObject *ra_array;
    PyArrayObject *dec_array;

    PyObject *final_rotation = Py_None;
    PyObject *weather_data = Py_None;

    char *kwnames[] = {"ctime",
                       "az",
                       "alt",
                       "ra",
                       "dec",
                       "final_q",
                       "weather",
                       NULL};


    if (!PyArg_ParseTupleAndKeywords(
            args, kwargs,
            "O!O!O!O!O!|OO", kwnames,
            &PyArray_Type, &ctime_array,
            &PyArray_Type, &az_array,
            &PyArray_Type, &alt_array,
            &PyArray_Type, &ra_array,
            &PyArray_Type, &dec_array,
            &final_rotation,
            &weather_data
            ))
        return NULL;

    /*
     * Check types and stuff of input ctime, alt, az arrays.
     */
    ASSERT_CARRAY_TYPE_NDIM(ctime_array, NPY_FLOAT64, 1);
    ASSERT_CARRAY_TYPE_NDIM(az_array, NPY_FLOAT64, 1);
    ASSERT_CARRAY_TYPE_NDIM(alt_array, NPY_FLOAT64, 1);
    ASSERT_CARRAY_TYPE_NDIM(ra_array, NPY_FLOAT64, 1);
    ASSERT_CARRAY_TYPE_NDIM(dec_array, NPY_FLOAT64, 1);

    long n = PyArray_DIMS(ctime_array)[0];
    ASSERT_AXIS_N(az_array, 0, n);
    ASSERT_AXIS_N(alt_array, 0, n);
    ASSERT_AXIS_N(ra_array, 0, n);
    ASSERT_AXIS_N(dec_array, 0, n);
    
    double *ctime = PyArray_DATA(ctime_array);
    double *az = PyArray_DATA(az_array);
    double *alt = PyArray_DATA(alt_array);
    double *ra = PyArray_DATA(ra_array);
    double *dec = PyArray_DATA(dec_array);

    double dummy_horn = 0.;
    struct focal_plane_data fplane = {.n = 1,
                                      .x = &dummy_horn,
                                      .y = &dummy_horn,
                                      .phi = &dummy_horn };

    /* Unpack any weather */
    ACTpolWeather weather;
    ACTpolWeather_default(&weather);
    if (weather_data != Py_None && 
        !unpack_weather(&weather, weather_data))
        po_raise("invalid weather specifier.");

    ACTpolArray *array = ACTpolArray_alloc(fplane.n);
    ACTpolArray_init(array, 150., 0., 0.);
    for (int i=0; i<fplane.n; i++)
        ACTpolFeedhorn_init(array->horn+i, fplane.y[i], fplane.x[i],
                            -M_PI/2 + fplane.phi[i]);
        
    ACTpolArrayCoords *coords = ACTpolArrayCoords_alloc(array);
#if LIBACTPOL_API < 20180101
    ACTpolArrayCoords_init(coords);
#else
    ACTpolArrayCoords_init(coords, ACTPOL_COORDSYS_GENERAL);
#endif    

    ACTpolState *state = ACTpolState_alloc();
    ACTpolState_init(state);

    double throw = .1;
    ACTpolScan scan;
    ACTpolScan_init(&scan, alt[0], az[0], throw);

    ACTpolArrayCoords_update_refraction(coords, &scan, &weather);

    int nfields = 2;
    npy_intp dims[1] = {n};

    PyArrayObject *output_arrays[2];
    double *output[2];
    
    for (int i=0; i<nfields; i++) {
        output_arrays[i] = (PyArrayObject*)
            PyArray_SimpleNew(1, dims, NPY_FLOAT64);
        output[i] = PyArray_DATA(output_arrays[i]);
    }

    for (int j=0; j<n; j++) {
        ACTpolState_update(state, ctime[j], alt[j], az[j]);
        //Update tracking position
        
        Quaternion *q = &state->BCRS_to_NSC_q;
        Quaternion_r3(*q, -ra[j]);
        Quaternion_r2_mul(-M_PI/2+dec[j], *q);

        ACTpolArrayCoords_update(coords, state);
        ACTpolFeedhornCoords *fc = coords->horn;
        output[0][j] = fc->m[0][2];  //dx
        output[1][j] = fc->m[1][2];  //dy
    }

    ACTpolState_free(state);
    ACTpolArrayCoords_free(coords);
    ACTpolArray_free(array);

    PyObject* tuple_out = PyTuple_New(nfields);
    for (int i=0; i<nfields; i++)
        PyTuple_SET_ITEM(tuple_out, i, (PyObject*)output_arrays[i]);

    return tuple_out;
}


/************
 * wand routines
 *
 * Pointing approximations for various projections.  These are called
 * by the main wand_get_coords routine, depending on the projection
 * specified.
 *
 * wand_pointing_single_xy is for orthographic (tangent plane).
 *
 * wand_pointing_single_ra_sindec is for cylindrical equal area
 *   projections.
 *
 * wand_pointing_single_ra_dec computes (ra,dec). Not really a projection.
 */


static void wand_pointing_single_xy(
    struct wand_data *wand,
    const double x,
    const double y,
    const double phi,
    double * restrict x_out,
    double * restrict y_out,
    double * restrict cosg,
    double * restrict sing,
    double * restrict cos2g,
    double * restrict sin2g)
{
    double cphi, sphi;
    if (cosg != NULL || cos2g != NULL) {
        cphi = cos(phi);
        sphi = sin(phi);
    }
#pragma omp parallel for
    for (int j=0; j<wand->n; j++) {
        if (x_out != NULL) {
            double dx = - x*wand->cosg[j] - y*wand->sing[j];
            double dy = + y*wand->cosg[j] - x*wand->sing[j];
            x_out[j] = wand->x[j] + dx;
            y_out[j] = wand->y[j] + dy;
        }
        if (cosg != NULL || cos2g != NULL) {
            /* double c = sphi*wand->cosg[j] - cphi*wand->sing[j]; */
            /* double s = sphi*wand->sing[j] + cphi*wand->cosg[j]; */
            double c, s;
            if (wand->cos_2hwp == NULL) {
                c = sphi*wand->cosg[j] - cphi*wand->sing[j];
                s = sphi*wand->sing[j] + cphi*wand->cosg[j];
            } else {
                /* Take gamma -> (gamma - 2*chi) and phi -> -phi. */
                double cg = wand->cosg[j]*wand->cos_2hwp[j] + wand->sing[j]*wand->sin_2hwp[j];
                double sg = wand->sing[j]*wand->cos_2hwp[j] - wand->cosg[j]*wand->sin_2hwp[j];
                c = -sphi*cg - cphi*sg;
                s = -sphi*sg + cphi*cg;
                /* This is the right idea but there's a sign flip we
                 * haven't worked out properly.*/
                /* double */
                /*     _c = c*wand->cos_2hwp[j] - s*wand->sin_2hwp[j]; */
                /* s = s*wand->cos_2hwp[j] + c*wand->sin_2hwp[j]; */
                /* c = _c; */
            }
            if (cosg != NULL) {
                cosg[j] = c;
                sing[j] = s;
            }
            if (cos2g != NULL) {
                cos2g[j] = c*c - s*s;
                sin2g[j] = 2*s*c;
            }
        }
    }
}

static void wand_pointing_single_ra_sindec(
    struct wand_data *wand,
    const double x,
    const double y,
    double * restrict ra_out,
    double * restrict sin_dec_out)
{
    double r = x*x + y*y;
    double z = sqrt(1-r);

#pragma omp parallel for
    for (int j=0; j<wand->n; j++) {
        double dx =  ( x*wand->sing[j] - y*wand->cosg[j]);
        double dy = -( y*wand->sing[j] + x*wand->cosg[j]);
        double q = dx / (z*wand->cosdec[j] - dy*wand->sindec[j]);
        ra_out[j] = wand->ra[j] + q*(1.-q*q/3);
        sin_dec_out[j] = wand->sindec[j]*(1 - 0.5*r) + wand->cosdec[j]*dy;
    }
}

static void wand_pointing_single_ra_dec(
    const struct wand_data *wand,
    const double x,
    const double y,
    double * restrict ra_out,
    double * restrict dec_out)
{
    double r = x*x + y*y;
    double z = sqrt(1-r);
#pragma omp parallel for
    for (int j=0; j<wand->n; j++) {
        double dx =  ( x*wand->sing[j] - y*wand->cosg[j]);
        double dy = -( y*wand->sing[j] + x*wand->cosg[j]);
        double q = dx / (z*wand->cosdec[j] - dy*wand->sindec[j]);
        double t = dx * wand->sindec[j] / wand->cosdec[j];
        ra_out[j] = wand->ra[j] + q*(1.-q*q/3);
        dec_out[j] = wand->dec[j] + dy - 0.5*(t*dx - dy*(dy*dy/3 - t*t));
    }
}



/* wand_project_data_to_map
 *
 * Uses pointing wand to project data into a map.
 *
 * Input:
 *   data (data_array)
 *   cuts (packed cuts info)
 *   wand_data (packed wand info)
 *   pix_data (packed pixelization info)
 *   focal_plane_data (packed focal plane coordinates)
 *   dest_map (array: 1d, float32)
 *   dest_weights (array: 1d, int32)
 *
 * Output:  None
 */
static PyObject *wand_project_data_to_map(PyObject *self, PyObject *args)
{
    PyArrayObject *data_array;
    PyObject *cuts_data;
    PyObject *pix_data;
    PyObject *wand_data;
    PyObject *fplane_data;
    PyObject *output_object;
    PyObject *outputQ_object;
    PyObject *outputU_object;
    PyObject *weights_object;

    if (!PyArg_ParseTuple(args, "O!OOOOOOOO",
                          &PyArray_Type, &data_array,
                          &cuts_data,
                          &wand_data,
                          &pix_data,
                          &fplane_data,
                          &output_object,
                          &outputQ_object,
                          &outputU_object,
                          &weights_object))
        return NULL;

    struct wand_data wand;
    if (!unpack_wand(&wand, wand_data))
        po_raise("could not unpack wand data");

    struct pixelization_data pix;
    if (!unpack_pixelization(&pix, pix_data))
        po_raise("could not unpack pixelization data");

    struct focal_plane_data fplane;
    if (!unpack_focal_plane(&fplane, fplane_data))
        po_raise("could not unpack focal plane data");

    struct cuts_set_data cuts_set = {0,NULL,NULL};
    if (cuts_data != Py_None) {
        if (!pyactpol_unpack_cuts_data(&cuts_set, cuts_data))
            po_raise("could not unpack cuts data");
    }

    float *data = PyArray_DATA(data_array);
    float *output = NULL;
    float *outputQ = NULL;
    float *outputU = NULL;
    int *weights = NULL;

    ASSERT_CARRAY_TYPE_NDIM(data_array, NPY_FLOAT32, 2);

    if (output_object != Py_None) {
        po_assert(PyArray_Check(output_object));
        PyArrayObject *output_array = (PyArrayObject*)output_object;
        ASSERT_CARRAY_TYPE_NDIM(output_array, NPY_FLOAT32, -1);
        output = PyArray_DATA(output_array);
    }
    if (weights_object != Py_None) {
        po_assert(PyArray_Check(weights_object));
        PyArrayObject *output_array = (PyArrayObject*)weights_object;
        ASSERT_CARRAY_TYPE_NDIM(output_array, NPY_INT32, -1);
        weights = PyArray_DATA(output_array);
    }

    if (outputQ_object != Py_None) {
        po_assert(PyArray_Check(outputQ_object));
        PyArrayObject *output_array = (PyArrayObject*)outputQ_object;
        ASSERT_CARRAY_TYPE_NDIM(output_array, NPY_FLOAT32, -1);
        outputQ = PyArray_DATA(output_array);
    }

    if (outputU_object != Py_None) {
        po_assert(PyArray_Check(outputU_object));
        PyArrayObject *output_array = (PyArrayObject*)outputU_object;
        ASSERT_CARRAY_TYPE_NDIM(output_array, NPY_FLOAT32, -1);
        outputU = PyArray_DATA(output_array);
    }

    /* Handle weird strides in the first dimension of the data array? Not really. */
    npy_intp *data_strides = PyArray_STRIDES(data_array);
    po_assert(data_strides[1] == PyArray_ITEMSIZE(data_array));
    int det_stride = data_strides[0] / PyArray_ITEMSIZE(data_array);
    
    /* Output arrays can be 2-d, but should be float32 and C-ordered */
    int nhorns = fplane.n;
    int nsamps = wand.n;
    int i;

    /* Enforce agreement of detectors dimension */
    if ((cuts_data != Py_None && cuts_set.n != nhorns) ||
        (PyArray_DIMS(data_array)[0] != nhorns)) {
        po_raise("Number of detectors in data, cuts, and fplane do not agree.\n");
    }
    po_assert(nsamps == PyArray_DIMS(data_array)[1]);

    // Temporary lodging for the coordinates
    double *ra = malloc(nsamps * sizeof(*ra));
    double *dec = malloc(nsamps * sizeof(*dec));
    double *cos2g = NULL;
    double *sin2g = NULL;
    if (outputQ != NULL || outputU != NULL) {
        cos2g = malloc(nsamps * sizeof(*cos2g));
        sin2g = malloc(nsamps * sizeof(*sin2g));
    }

    for (i=0; i<nhorns; i++) {
        // Get pointing
        switch(wand.proj_id) {
        case PROJECTION_TANGENT:
            wand_pointing_single_xy(&wand, fplane.x[i], fplane.y[i], fplane.phi[i],
                                    ra, dec, NULL, NULL, cos2g, sin2g);
            break;
        case PROJECTION_RADEC:
            wand_pointing_single_ra_dec(&wand, fplane.x[i], fplane.y[i],
                                        ra, dec);
            break;
        case PROJECTION_RASINDEC:
            wand_pointing_single_ra_sindec(&wand, fplane.x[i], fplane.y[i],
                                           ra, dec);
            break;
        }

        // Loop through uncut regions.
        int j = 0;
        int seg_idx = 0;
        int seg_end;
        while (j < nsamps) {
            if (cuts_data!=Py_None && seg_idx < cuts_set.ncut[i])
                seg_end = cuts_set.cuts[i][seg_idx*2];
            else
                seg_end = nsamps;

            // The actual loop.
            for (; j<seg_end; j++) {
                int i0 = (int)( (ra[j] - pix.x0) / pix.dx0
                                + pix.pix0 + 1.5) - 1;
                if (i0 < 0 || i0 >= pix.n0)
                    continue;
                int i1 = (int)( (dec[j] - pix.x1) / pix.dx1
                                + pix.pix1 + 1.5) - 1;
                if (i1 < 0 || i1 >= pix.n1)
                    continue;
                int pix_idx = i0 + i1 * pix.n0;
                if (output != NULL)
                    output[pix_idx] += data[i*det_stride + j];
                if (outputQ != NULL)
                    outputQ[pix_idx] += data[i*det_stride + j] * cos2g[i];
                if (outputU != NULL)
                    outputU[pix_idx] -= data[i*det_stride + j] * sin2g[i];
                if (weights != NULL)
                    weights[pix_idx]++;
            }

            // Cue up next region
            if (cuts_data!=Py_None && seg_idx < cuts_set.ncut[i]) {
                j = cuts_set.cuts[i][seg_idx*2+1];
                seg_idx++;
            }
        }
    }

    free(ra);
    free(dec);
    free(cos2g);
    free(sin2g);

    if (cuts_data!=Py_None)
        pyactpol_destroy_cuts_data(&cuts_set);

    Py_RETURN_NONE;
}


/* wand_project_map_to_data
 *
 * Uses pointing wand to project map into data array.
 *
 * Input:
 *   map (array: 1d, float32)
 *   cuts (packed cuts info)
 *   wand_data (packed wand info)
 *   pix_data (packed pixelization info)
 *   focal_plane_data (packed focal plane coordinates)
 *   src_data (data_array)
 *
 * Output:  None
 */
static PyObject *wand_project_map_to_data(PyObject *self, PyObject *args)
{
    PyArrayObject *data_array;
    PyObject *cuts_data;
    PyObject *pix_data;
    PyObject *wand_data;
    PyObject *fplane_data;
    PyObject *mapI_object;
    PyObject *mapQ_object;
    PyObject *mapU_object;

    if (!PyArg_ParseTuple(args, "OOOOOOOO!",
                          &mapI_object,
                          &mapQ_object,
                          &mapU_object,
                          &cuts_data,
                          &wand_data,
                          &pix_data,
                          &fplane_data,
                          &PyArray_Type, &data_array))
        return NULL;
    
    struct wand_data wand;
    if (!unpack_wand(&wand, wand_data))
        po_raise("could not unpack wand data");

    struct pixelization_data pix;
    if (!unpack_pixelization(&pix, pix_data))
        po_raise("could not unpack pixelization data");

    struct focal_plane_data fplane;
    if (!unpack_focal_plane(&fplane, fplane_data))
        po_raise("could not unpack focal plane data");

    struct cuts_set_data cuts_set = {0,NULL,NULL};
    if (cuts_data != Py_None) {
        if (!pyactpol_unpack_cuts_data(&cuts_set, cuts_data))
            po_raise("could not unpack cuts data");
    }

    ASSERT_CARRAY_TYPE_NDIM(data_array, NPY_FLOAT32, 2);
    float *dest_data = PyArray_DATA(data_array);

    float *I_map = NULL;
    float *Q_map = NULL;
    float *U_map = NULL;
    if (mapI_object != Py_None) {
        po_assert(PyArray_Check(mapI_object));
        PyArrayObject *map = (PyArrayObject*)mapI_object;
        ASSERT_CARRAY_TYPE_NDIM(map, NPY_FLOAT32, -1);
        I_map = PyArray_DATA(map);
    }
    if (mapQ_object != Py_None) {
        po_assert(PyArray_Check(mapQ_object));
        PyArrayObject *map = (PyArrayObject*)mapQ_object;
        ASSERT_CARRAY_TYPE_NDIM(map, NPY_FLOAT32, -1);
        Q_map = PyArray_DATA(map);
    }
    if (mapU_object != Py_None) {
        po_assert(PyArray_Check(mapU_object));
        PyArrayObject *map = (PyArrayObject*)mapU_object;
        ASSERT_CARRAY_TYPE_NDIM(map, NPY_FLOAT32, -1);
        U_map = PyArray_DATA(map);
    }


    /* Handle weird strides in the first dimension of the data array */
    npy_intp *data_strides = PyArray_STRIDES(data_array);
    po_assert(data_strides[1] == PyArray_ITEMSIZE(data_array));
    int det_stride = data_strides[0] / PyArray_ITEMSIZE(data_array);

    int nhorns = fplane.n;
    int nsamps = wand.n;
    int i;

    /* Enforce agreement of detectors dimension */
    if ((cuts_data != Py_None && cuts_set.n != nhorns) ||
        (PyArray_DIMS(data_array)[0] != nhorns)) {
        po_raise("Number of detectors in data, cuts, and fplane do not agree.\n");
    }

    // Temporary lodging for the coordinates
    double *ra = malloc(nsamps * sizeof(*ra));
    double *dec = malloc(nsamps * sizeof(*dec));
    double *cos2g = NULL;
    double *sin2g = NULL;
    if (Q_map != NULL || U_map != NULL) {
        cos2g = malloc(nsamps * sizeof(*cos2g));
        sin2g = malloc(nsamps * sizeof(*sin2g));
    }        
    for (i=0; i<nhorns; i++) {
        // Get pointing
        switch(wand.proj_id) {
        case PROJECTION_TANGENT:
            wand_pointing_single_xy(&wand, fplane.x[i], fplane.y[i], fplane.phi[i],
                                    ra, dec, NULL, NULL, cos2g, sin2g);
            break;
        case PROJECTION_RADEC:
            wand_pointing_single_ra_dec(&wand, fplane.x[i], fplane.y[i],
                                        ra, dec);
            break;
        case PROJECTION_RASINDEC:
            wand_pointing_single_ra_sindec(&wand, fplane.x[i], fplane.y[i],
                                           ra, dec);
            break;
        }

        // Loop through uncut regions.
        int j = 0;
        int seg_idx = 0;
        int seg_end;
        while (j < nsamps) {
            if (cuts_data!=Py_None && seg_idx < cuts_set.ncut[i])
                seg_end = cuts_set.cuts[i][seg_idx*2];
            else
                seg_end = nsamps;

            // The actual loop.
            for (; j<seg_end; j++) {
                int i0 = (int)( (ra[j] - pix.x0) / pix.dx0
                                + pix.pix0 + 1.5) - 1;
                if (i0 < 0 || i0 >= pix.n0)
                    continue;
                int i1 = (int)( (dec[j] - pix.x1) / pix.dx1
                                + pix.pix1 + 1.5) - 1;
                if (i1 < 0 || i1 >= pix.n1)
                    continue;
                int pix_idx = i0 + i1 * pix.n0;
                if (I_map != NULL)
                    dest_data[i*det_stride + j] += I_map[pix_idx];
                if (Q_map != NULL)
                    dest_data[i*det_stride + j] += Q_map[pix_idx]*cos2g[i];
                if (U_map != NULL)
                    dest_data[i*det_stride + j] -= U_map[pix_idx]*sin2g[i];
            }

            // Cue up next region
            if (cuts_data!=Py_None && seg_idx < cuts_set.ncut[i]) {
                j = cuts_set.cuts[i][seg_idx*2+1];
                seg_idx++;
            }
        }
    }

    free(ra);
    free(dec);
    free(cos2g);
    free(sin2g);
    if (cuts_data!=Py_None)
        pyactpol_destroy_cuts_data(&cuts_set);

    Py_RETURN_NONE;
}


static PyObject *wand_get_coords(PyObject *self, PyObject *args)
{
    PyObject *fplane_data;
    PyObject *wand_data;
    int get_coords = 0;
    int get_rotation = 0;
    int get_polarization = 0;

    if (!PyArg_ParseTuple(args, "OOiii",
                          &wand_data,
                          &fplane_data,
                          &get_coords,
                          &get_rotation,
                          &get_polarization
            ))
        return NULL;

    struct wand_data wand;
    if (!unpack_wand(&wand, wand_data))
        po_raise("could not unpack wand data");

    struct focal_plane_data fplane;
    if (!unpack_focal_plane(&fplane, fplane_data))
        po_raise("could not unpack focal plane data");

    int nhorns = fplane.n;
    int nsamps = wand.n;

    // Output places
    npy_intp dims_out[2];
    dims_out[0] = nhorns;
    dims_out[1] = nsamps;

    PyArrayObject *ra_out_array = NULL, *dec_out_array = NULL;
    double *ra_out = NULL, *dec_out = NULL;
    if (get_coords) {
        ra_out_array = (PyArrayObject*)
            PyArray_SimpleNew(2, dims_out, NPY_FLOAT64);
        dec_out_array = (PyArrayObject*)
            PyArray_SimpleNew(2, dims_out, NPY_FLOAT64);
        ra_out = PyArray_DATA(ra_out_array);
        dec_out = PyArray_DATA(dec_out_array);
    }

    PyArrayObject *cosg_array = NULL, *sing_array = NULL;
    double *cosg_out = NULL, *sing_out = NULL;
    if (get_rotation) {
        cosg_array = (PyArrayObject*)
            PyArray_SimpleNew(2, dims_out, NPY_FLOAT64);
        sing_array = (PyArrayObject*)
            PyArray_SimpleNew(2, dims_out, NPY_FLOAT64);
        cosg_out = PyArray_DATA(cosg_array);
        sing_out = PyArray_DATA(sing_array);
    }

    PyArrayObject *cos2g_array = NULL, *sin2g_array = NULL;
    double *cos2g_out = NULL, *sin2g_out = NULL;
    if (get_polarization) {
        cos2g_array = (PyArrayObject*)
            PyArray_SimpleNew(2, dims_out, NPY_FLOAT64);
        sin2g_array = (PyArrayObject*)
            PyArray_SimpleNew(2, dims_out, NPY_FLOAT64);
        cos2g_out = PyArray_DATA(cos2g_array);
        sin2g_out = PyArray_DATA(sin2g_array);
    }

    double *ra = NULL;
    double *dec = NULL;
    double *cosg = NULL;
    double *sing = NULL;
    double *cos2g = NULL;
    double *sin2g = NULL;
    for (int i=0; i<fplane.n; i++) {
        if (ra_out != NULL) {
            ra = ra_out + i*wand.n;
            dec = dec_out + i*wand.n;
        }
        if (cosg_out != NULL) {
            cosg = cosg_out + i*wand.n;
            sing = sing_out + i*wand.n;
        }
        if (cos2g_out != NULL) {
            cos2g = cos2g_out + i*wand.n;
            sin2g = sin2g_out + i*wand.n;
        }
        switch(wand.proj_id) {
        case PROJECTION_TANGENT:
            wand_pointing_single_xy(&wand, fplane.x[i], fplane.y[i], fplane.phi[i],
                                    ra, dec, cosg, sing, cos2g, sin2g);
            break;
        case PROJECTION_RADEC:
            wand_pointing_single_ra_dec(&wand, fplane.x[i], fplane.y[i],
                                        ra, dec);
            break;
        case PROJECTION_RASINDEC:
            wand_pointing_single_ra_sindec(&wand, fplane.x[i], fplane.y[i],
                                           ra, dec);
            break;
        }
    }
    
    // Create output tuple.
    int tuple_size = 2*((get_coords!=0) + (get_rotation!=0) +
                        (get_polarization!=0)),
        tuple_idx = 0;
    PyObject *tuple_out = PyTuple_New(tuple_size);
    if (get_coords) {
        PyTuple_SET_ITEM(tuple_out, tuple_idx++, (PyObject*)ra_out_array);
        PyTuple_SET_ITEM(tuple_out, tuple_idx++, (PyObject*)dec_out_array);
    }
    if (get_rotation) {
        PyTuple_SET_ITEM(tuple_out, tuple_idx++, (PyObject*)cosg_array);
        PyTuple_SET_ITEM(tuple_out, tuple_idx++, (PyObject*)sing_array);
    }
    if (get_polarization) {
        PyTuple_SET_ITEM(tuple_out, tuple_idx++, (PyObject*)cos2g_array);
        PyTuple_SET_ITEM(tuple_out, tuple_idx++, (PyObject*)sin2g_array);
    }

    return tuple_out;
}



/* set_bulletinA
 *
 * Pass dut1, dx, dy data to libactpol.
 *
 * Input:
 *   mjd_min (int)
 *   mjd_max (int)
 *   dut1 (double array), size mjd_max-mjd_min+1.
 *   dx (same)
 *   dy (same)
 *
 * Output:  None
 */
static PyObject *set_bulletinA(PyObject *self, PyObject *args)
{
    int mjd_min, mjd_max;
    PyArrayObject *dut1_array, *dx_array, *dy_array;
    double *dut1, *dx, *dy;

    if (!PyArg_ParseTuple(args, "iiO!O!O!",
                          &mjd_min,
                          &mjd_max,
                          &PyArray_Type, &dut1_array,
                          &PyArray_Type, &dx_array,
                          &PyArray_Type, &dy_array
            ))
        return NULL;

    int n = mjd_max - mjd_min + 1;

    ASSERT_CARRAY_TYPE_NDIM(dut1_array, NPY_FLOAT64, 1);
    ASSERT_CARRAY_TYPE_NDIM(dx_array, NPY_FLOAT64, 1);
    ASSERT_CARRAY_TYPE_NDIM(dy_array, NPY_FLOAT64, 1);

    ASSERT_AXIS_N(dut1_array, 0, n);
    ASSERT_AXIS_N(dx_array, 0, n);
    ASSERT_AXIS_N(dy_array, 0, n);

    dut1 = PyArray_DATA(dut1_array);
    dx = PyArray_DATA(dx_array);
    dy = PyArray_DATA(dy_array);

    actpol_set_iers_bulletin_a(mjd_min, mjd_max, dut1, dx, dy);
    Py_RETURN_NONE;
}




PyMethodDef pyactpol_pointing_methods[] = {
    {"get_coords", (PyCFunction)get_coords,
     METH_VARARGS | METH_KEYWORDS, ""},
    {"get_coords_inverse", (PyCFunction)get_coords_inverse,
     METH_VARARGS | METH_KEYWORDS, ""},
    {"wand_get_coords", wand_get_coords,
     METH_VARARGS, ""},
    {"wand_project_data_to_map", wand_project_data_to_map,
     METH_VARARGS, ""},
    {"wand_project_map_to_data", wand_project_map_to_data,
     METH_VARARGS, ""},
    {"set_bulletinA", set_bulletinA,
     METH_VARARGS, ""},
    {NULL, NULL, 0, NULL}        /* Sentinel */
};


