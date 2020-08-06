/* -*- mode: C; tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 4 -*-
 *      vim: sw=4 ts=8 et tw=80
 */

#include "pyactpol.h"

#include <actpol/actpol.h>

#define ASSERT_IS_CUTS_ARRAY(cuts) do {                 \
        po_assert(PyArray_NDIM(cuts) == 2 &&            \
                  PyArray_DIMS(cuts)[1] == 2 &&         \
                  PyArray_TYPE(cuts) == NPY_INT32);     \
    } while(0)

#define MAX(a, b) (a>b ? a : b)

#pragma pack(push)
#pragma pack(1)
typedef struct {
    int32_t start;
    int32_t stop;
} cut_range_t;
#pragma pack(pop)

/* Convert a boolean (char*) array into an int array of cut regions.
 * This makes no python calls and can be used from openmp contexts.
 * n_cuts is set and the malloc'd array is returned.  The array should
 * later be free()d.  When there are no cut regions, n_cuts=0 and the
 * function returns NULL.  This is declared non-static so that
 * filter.c can use it. */

int *pyactpol_cuts_vector_from_mask(char *mask, int n_mask, int *n_cuts)
{
    // Count the cuts.
    *n_cuts = 0;
    int i, on;
    on = 0;
    for (i=0; i<n_mask; i++) {
        if (on && !mask[i]) {
            on = 0;
        } else if (!on && mask[i]) {
            on = 1;
            (*n_cuts)++;
        }
    }

    // Allocate
    if (*n_cuts==0)
        return NULL;
    int *cuts = malloc(*n_cuts*2*sizeof(*cuts));
    assert(cuts != NULL);

    // Store
    int oindex = 0;
    on = 0;
    for (i=0; i<n_mask; i++) {
        if (on && !mask[i]) {
            on = 0;
            cuts[oindex++] = i;
        } else if (!on && mask[i]) {
            on = 1;
            cuts[oindex++] = i;
        }
    }
    if (on)
        cuts[oindex++] = i;
    assert(oindex==*n_cuts*2);

    return cuts;
}

/* Convert an array of bools into an CutsVector array -- start and
 * stop positions for the cuts.
 */

static PyObject *cuts_from_mask(PyObject *self, PyObject *args)
{
    PyArrayObject *mask_array;
    if (!PyArg_ParseTuple(args, "O!",
                          &PyArray_Type, &mask_array
            ))
        return NULL;
    
    ASSERT_CARRAY_TYPE_NDIM(mask_array, NPY_BOOL, 1);

    int n_mask = PyArray_SIZE(mask_array);
    char *mask = PyArray_DATA(mask_array);

    int n_cuts = 0;
    int *cuts = pyactpol_cuts_vector_from_mask(mask, n_mask, &n_cuts);
    npy_intp dims[2] = {0, 2};
    dims[0] = n_cuts;
    PyObject *arr = PyArray_SimpleNew(2, dims, NPY_INT);
    int *dest = PyArray_DATA((PyArrayObject *)arr);
    if (cuts != NULL) {
        memcpy(dest, cuts, 2*n_cuts*sizeof(*dest));
        free(cuts);
    }
    return (PyObject*)arr;
}

static PyObject *mask_from_cuts(PyObject *self, PyObject *args)
{
    PyArrayObject *mask_array;
    PyArrayObject *cuts_array;
    int n_mask;
    int invert;

    if (!PyArg_ParseTuple(args, "O!ii",
                          &PyArray_Type, &cuts_array,
                          &n_mask, &invert
            ))
        return NULL;

    ASSERT_IS_CUTS_ARRAY(cuts_array);

    int *cuts = PyArray_DATA(cuts_array);
    int n_cuts2 = 2*PyArray_DIMS(cuts_array)[0];

    npy_intp dim = n_mask;
    mask_array = (PyArrayObject*)
        PyArray_SimpleNew(1, &dim, NPY_BOOL);
    char *mask = PyArray_DATA(mask_array);

    int i = 0;
    char on = 0;
    char invertc = 0;
    if (invert) invertc = 1;
    while (i<n_mask && n_cuts2>0) {
        while (i<n_mask && i<*cuts) mask[i++] = on ^ invertc;
        on = 1-on;
        cuts++;
        n_cuts2--;
    }
    while (i<n_mask) mask[i++] = 0;
    return (PyObject*)mask_array;
}


static PyObject *fill_one_cuts(PyObject *self, PyObject *args)
{
    PyArrayObject *cuts_array;
    PyArrayObject *data_array;
    PyObject *noise_object;
    int fit_region;
    int extrapolate;

    if (!PyArg_ParseTuple(args, "O!O!Oii",
                          &PyArray_Type, &data_array,
                          &PyArray_Type, &cuts_array,
                          &noise_object,
                          &fit_region,
                          &extrapolate
            ))
        return NULL;
    
    ASSERT_IS_CUTS_ARRAY(cuts_array);
    int *cuts = PyArray_DATA(cuts_array);
    int n_cuts = PyArray_DIMS(cuts_array)[0];
    
    ASSERT_CARRAY_TYPE_NDIM(data_array, NPY_FLOAT32, 1);

    float *data = PyArray_DATA(data_array);
    int nsamps = PyArray_DIMS(data_array)[0];

    float *noise = NULL;
    int n_noise = 0;
    int noise_i = 0;
    if (noise_object != Py_None) {
        PyArrayObject *noise_array = (PyArrayObject *)noise_object;
        assert(PyArray_Check(noise_array));
        ASSERT_CARRAY_TYPE_NDIM(noise_array, NPY_FLOAT32, 1);
        noise = PyArray_DATA(noise_array);
        n_noise = PyArray_DIMS(noise_array)[0];
    }

    int last_cut = 0;

    for (int icut=0; icut<n_cuts; icut++) {
        if (cuts[icut*2] > nsamps)
            break;
        int left_edge = cuts[icut*2] - fit_region;
        if (left_edge < last_cut)
            left_edge = last_cut;
        int right_edge = cuts[icut*2+1] + fit_region;
        if (right_edge > nsamps)
            right_edge = nsamps;
        if (icut+1 < n_cuts && cuts[(icut+1)*2] < right_edge)
            right_edge = cuts[(icut+1)*2];
        /* Linear fit */
        double xbar = 0., x2bar = 0., xybar = 0., ybar = 0.;
        int n_left = 0, n_right = 0;
        int i;
        for (i=left_edge; i<cuts[icut*2]; i++) {
            n_left += 1;
            double dx = i-left_edge;
            xbar += dx;
            x2bar += dx*dx;
            xybar += dx*data[i];
            ybar += data[i];
        }
        for (i=cuts[icut*2+1]; i<right_edge; i++) {
            n_right += 1;
            double dx = i-left_edge;
            xbar += dx;
            x2bar += dx*dx;
            xybar += dx*data[i];
            ybar += data[i];
        }
        int n = n_left + n_right;
        double m, b;
        if (((n_left>0) && (n_right>0)) || (extrapolate && (n>1))) {
            float determ = (x2bar*n - xbar*xbar);
            m = (n * xybar - xbar * ybar) / determ;
            b = (x2bar * ybar - xbar * xybar) / determ;
        } else {
            m = 0;
            b = ybar / n;
        }
        // Get stdev
        double dy2 = 0;
        if (noise != NULL && n > 1) {
            for (i=left_edge; i<cuts[icut*2]; i++)
                dy2 += pow(data[i] - m*(i-left_edge) - b, 2);
            for (i=cuts[icut*2+1]; i<right_edge; i++)
                dy2 += pow(data[i] - m*(i-left_edge) - b, 2);
            dy2 = sqrt(dy2/n);
        }
        for (i=cuts[icut*2]; i<cuts[icut*2+1]; i++) {
            double n = 0;
            if (noise != NULL)  {
                if (noise_i >= n_noise)
                    po_raise("required more noise samples than were provided!\n");
                n = dy2*noise[noise_i++];
            }
            data[i] = m*(i-left_edge) + b + n;
        }
        
    }

    Py_RETURN_NONE;
}

static int cmp_cut_range(const void *v0, const void *v1)
{
    const cut_range_t *c0 = v0;
    const cut_range_t *c1 = v1;
    int delta = c0->start - c1->start;
    if (delta != 0)
        return delta;
    return c0->stop - c1->stop;
}

static PyObject *merge_cuts(PyObject *self, PyObject *args)
{
    PyObject *super_list;
    if (!PyArg_ParseTuple(args, "O!",
                          &PyList_Type, &super_list
            ))
        return NULL;
    
    // Input must be a list of CutsVector style arrays.  Check that,
    // and count the elements.
    int n_cuts = 0;
    for (int i=0; i<PyList_Size(super_list); i++) {
        PyArrayObject *cutsv = (void*)PyList_GetItem(super_list, i);
        ASSERT_IS_CUTS_ARRAY(cutsv);
        n_cuts += PyArray_DIMS(cutsv)[0];
    }

    if (n_cuts == 0) {
        npy_intp odims[2] = {0, 2};
        return PyArray_SimpleNew(2, odims, NPY_INT32);
    }

    /* Temporary storage for the cuts data, to be sorted. */
    cut_range_t *cuts = malloc(n_cuts * sizeof(*cuts));
    cut_range_t *p = cuts;
    for (int i=0; i<PyList_Size(super_list); i++) {
        PyArrayObject *cutsv = (void*)PyList_GetItem(super_list, i);
        int this_n_cuts = PyArray_DIMS(cutsv)[0];
        memcpy(p, PyArray_DATA(cutsv), this_n_cuts * sizeof(*p));
        p += this_n_cuts;
    }

    /* Sort it by starting index... */
    qsort(cuts, n_cuts, sizeof(*cuts), cmp_cut_range);
    
    /* Reduce it, by merging overlapping ranges. */
    cut_range_t *dst = cuts;
    cut_range_t *src;
    for (src = cuts + 1; src < cuts + n_cuts; src++) {
        if (src->start > dst->stop) {
            // The intervals do not overlap.  Finalize active dst;
            // initialize new dst with the src cut.
            *(++dst) = *src;
        } else {
            // Intervals overlap so update dst with the higher of the
            // stop values.
            dst->stop = MAX(src->stop, dst->stop);
        }
    }
    if (n_cuts > 0)
        dst++;
    int n_reduced = dst - cuts;

    npy_intp odims[2] = {0, 2};
    odims[0] = n_reduced;
    PyArrayObject *cutsv = (void*)PyArray_SimpleNew(2, odims, NPY_INT32);
    memcpy(PyArray_DATA(cutsv), cuts, sizeof(*cuts) * n_reduced);
    free(cuts);

    return (PyObject*) cutsv;
}


/* unpack_flags

   Inputs must be numpy(compatible) arrays corresponding to the
   stack_bounds, flag_stack, and index_stack entries in the TODFlags /
   SuperCuts format.  User must also pass the number of flag bits
   (n_flagvec) and the number of samples (n_samp) in each TOD vector.

   This returns a list containing n_flagvec elements.  Each element is
   a list of (2,n_cut) arrays suitable for conversion into CutsVector
   objects.
*/

static PyObject *unpack_flags(PyObject *self, PyObject *args)
{
    PyArrayObject *bounds_array; // uint32_t: (n_index)
    PyArrayObject *flag_array;   // uint8_t:  (n_index, n_flagbyte)
    PyArrayObject *index_array;  // uint32_t: (n_index)
    int n_flagvec, n_samp;

    if (!PyArg_ParseTuple(args, "O!O!O!ii",
                          &PyArray_Type, &bounds_array,
                          &PyArray_Type, &flag_array,
                          &PyArray_Type, &index_array,
                          &n_samp,
                          &n_flagvec
            ))
        return NULL;
    
    ASSERT_CARRAY_TYPE_NDIM(bounds_array, NPY_UINT32, 1);
    ASSERT_CARRAY_TYPE_NDIM(index_array,  NPY_UINT32, 1);
    ASSERT_CARRAY_TYPE_NDIM(flag_array,   NPY_UINT8, 2);

    int32_t *bounds = PyArray_DATA(bounds_array);
    int n_dets = PyArray_DIMS(bounds_array)[0] - 1;

    /* List for the cut lists; and the cut lists. */
    PyObject *super_list = PyList_New(n_flagvec);
    for (int i=0; i<n_flagvec; i++)
        PyList_SET_ITEM(super_list, i, PyList_New(n_dets));
    
    /* scratch space. */
    int n_scratch = 1024;
    uint32_t *scratch = malloc(n_scratch * sizeof(*scratch));

    /* For now just create trivial arrays as a time trial. */
    npy_intp odims[2];
    odims[1] = 2;
    for (int i=0; i<n_flagvec; i++) {
        PyObject *list = PyList_GetItem(super_list, i);
        int flag_stride = PyArray_DIMS(flag_array)[1];
        uint8_t flag_mask = (1 << (i % 8));
        uint32_t *index = PyArray_DATA(index_array);
        for (int j=0; j<n_dets; j++) {
            int i0 = bounds[j];
            int i1 = bounds[j+1];

            uint8_t *flag_ptr = PyArray_DATA(flag_array) + (i >> 8);
            flag_ptr += i0 * flag_stride;

            if (i1 - i0 > n_scratch) {
                n_scratch = (i1 - i0) * 2;
                scratch = realloc(scratch, n_scratch * sizeof(*scratch));
            }

            /* Record transitions into scratch vector. */
            uint8_t val = 0;
            int oidx = 0;
            for (int k = i0; k<i1; k++) {
                uint8_t new_val = (*flag_ptr) & flag_mask;
                flag_ptr += flag_stride;
                if (new_val != val) {
                    val = new_val;
                    scratch[oidx++] = index[k];
                }
            }
            /* Always record an even number of events. */
            if (val != 0)
                scratch[oidx++] = n_samp;

            /* Copy transition points into the output array. */
            odims[0] = oidx/2;
            PyArrayObject *out = (void*)PyArray_SimpleNew(2, odims, NPY_UINT32);
            uint32_t *outd = PyArray_DATA(out);
            memcpy(outd, scratch, oidx * sizeof(*outd));
            PyList_SET_ITEM(list, j, (PyObject*)out);
        }
    }        

    free(scratch);
    return super_list;
}


/* pack_flags

   Accepts a list of length n_flagvec.  Each element of this must be a
   list of CutVectors (or similar numpy-compatible objects).  Each
   list of CutVectors must be the same length.

   Returns a tuple with three elements, corresponding to the
   stack_bounds, flag_stack and index_stack arrays that are the packed
   storage blocks for the TODFlags / SuperCuts format.
*/

struct cutel {
    int32_t det;
    int32_t index;
    int32_t order;
    uint8_t bit;
    uint8_t on;
    uint8_t keep;
};

static int cmp_cutel(const void *v0, const void *v1)
{
    const struct cutel *c0 = v0;
    const struct cutel *c1 = v1;
    int delta = c0->det - c1->det;
    if (delta != 0)
        return delta;
    delta = (c0->index - c1->index);
    if (delta != 0)
        return delta;
    delta = (c0->bit - c1->bit);
    if (delta != 0)
        return delta;
    return c0->order - c1->order;
}

static PyObject *pack_flags(PyObject *self, PyObject *args)
{
    PyObject *super_list;
    int n_samp;
    if (!PyArg_ParseTuple(args, "O!i",
                          &PyList_Type, &super_list,
                          &n_samp
            ))
        return NULL;
    
    int n_flagvec = PyList_Size(super_list);

    // Get the number of cut regions for dimensioning output arrays.
    int n_det = -1;
    int n_cuts = 0;
    for (int i=0; i<n_flagvec; i++) {
        PyObject *list = PyList_GetItem(super_list, i);
        if (n_det == -1)
            n_det = PyList_Size(list);
        else
            assert(n_det == PyList_Size(list));
        for (int j=0; j<n_det; j++) {
            PyArrayObject *cutsv = (void*)PyList_GetItem(list, j);
            ASSERT_IS_CUTS_ARRAY(cutsv);
            // Two entries per cut region, plus the initialization entry.
            n_cuts += PyArray_DIMS(cutsv)[0] * 2 + 1;
        }
    }

    // Working area.
    struct cutel* cutels = calloc(n_cuts, sizeof(*cutels));

    // Copy everything in.
    struct cutel* p = cutels;
    int outi = 0;
    for (int i=0; i<n_flagvec; i++) {
        PyObject *list = PyList_GetItem(super_list, i);
        for (int j=0; j<n_det; j++) {
            PyArrayObject *cutsv = (void*)PyList_GetItem(list, j);
            int n_reg = PyArray_DIMS(cutsv)[0];
            int32_t *data = PyArray_DATA(cutsv);

            // Initialize...
            p->det = j;
            p->index = 0;
            p->on = 0;
            p->bit = i;
            p->order = outi++;
            p++;

            for (int k=0; k<n_reg; k++) {
                p->det = j;
                p->index = data[k*2];
                p->bit = i;
                p->on = 1;
                p->order = outi++;
                p++;
                p->det = j;
                p->index = data[k*2+1];
                p->bit = i;
                p->on = 0;
                p->order = outi++;
                p++;
            }
        }
    }

    /* Now sort that. */
    qsort(cutels, n_cuts, sizeof(*cutels), cmp_cutel);
    
    /* Compress overlapping entries -- rewrite the "order" member
     * represent the index into the output stack where this flag value
     * should be included. */
    int new_order = 0;
    p = cutels;

    // First point...
    p->order = new_order;
    struct cutel* last_p = p++;

    for (int i=1; i<n_cuts; i++) {
        if (p->index >= n_samp) {
            // If it's the last sample, disregard completely.
            p->order = -1;
            p++;
        } else {
            // If this flag is on a new detector, or new index... then
            // update the output pointer.
            if (p->det != last_p->det || p->index != last_p->index)
                new_order++;
            p->order = new_order;
            last_p = p++;
        }
    }
    int n_stack = new_order + 1;

    /* The bounds array is used to quickly find the flag and index for
     * a given detector. */
    npy_intp odims[2];
    odims[0] = n_det + 1;
    PyArrayObject *bounds_array = (void*)PyArray_SimpleNew(1, odims, NPY_UINT32);
    uint32_t *bounds = PyArray_DATA(bounds_array);

    /* The stack arrays. */
    int n_flagword = (n_flagvec + 7) / 8;
    odims[0] = n_stack;
    odims[1] = n_flagword;
    PyArrayObject *flag_array = (void*)PyArray_SimpleNew(2, odims, NPY_UINT8);
    PyArrayObject *index_array = (void*)PyArray_SimpleNew(1, odims, NPY_UINT32);
    uint8_t *flag = PyArray_DATA(flag_array);
    uint32_t *index = PyArray_DATA(index_array);

    const int max_flagword = 16;
    assert(n_flagword <= max_flagword);

    uint8_t f[max_flagword];
    memset(f, 0, max_flagword);

    int last_det = -1;
    for (int i=0; i<n_cuts; i++) {
        p = cutels + i;
        if (p->order < 0)
            continue;
        // If this is new det, write bounds array.
        while (p->det > last_det)
            bounds[++last_det] = p->order;
        // Update flags.
        uint8_t *dest = f + (p->bit / 8);
        if (p->on)
            *dest = (*dest) | (1 << (p->bit % 8));
        else 
            *dest = (*dest) & ~(1 << (p->bit % 8));
        index[p->order] = p->index;
        memcpy(flag + p->order * n_flagword, f, n_flagword);
    }
    while (last_det < n_det)
        bounds[++last_det] = n_stack; // Termination.

    free(cutels);

    PyObject *output = PyTuple_Pack(3, bounds_array, flag_array, index_array);
    return (void*)output;
}





PyMethodDef pyactpol_cuts_methods[] = {
    {"cuts_from_mask", cuts_from_mask, METH_VARARGS,
     ""},
    {"mask_from_cuts", mask_from_cuts, METH_VARARGS,
     ""},
    {"fill_cuts", fill_one_cuts, METH_VARARGS,
     ""},
    {"merge_cuts", merge_cuts, METH_VARARGS,
     ""},
    {"unpack_flags", unpack_flags, METH_VARARGS,
     ""},
    {"pack_flags", pack_flags, METH_VARARGS,
     ""},
    {NULL, NULL, 0, NULL}        /* Sentinel */
};


