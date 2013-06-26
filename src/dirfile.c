/* -*- mode: C; tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 4 -*-
 *      vim: sw=4 ts=8 et tw=80
 */

#include "pyactpol.h"

#include <actpol/actpol.h>
#include <actpol/dirfile.h>
#include <actpol/getdata.h>

#include "myassert.h"

/* Extraction / conversion routines, for pulling bitfields out and
 * rescaling them. */

static
int extract_float32_int32(float *dest, int32_t *src, long n,
                          long mask, float scale);

static
int extract_float32_uint32(float *dest, uint32_t *src, long n,
                          long mask, float scale);


/* ptrobj
 *
 * This is a generic object for encapsulating a C pointer, so it can
 * be passed back and forth between C and python layers.  If this is
 * abused it will certainly segfault and/or leak memory.  Maybe not in
 * that order.  Proper handling should ensure that copies of such
 * objects are associated with exactly one python object, and
 * automatically destroyed (using an appropriate destructor) when the
 * owner dies (i.e. through the owner's __del__ method).
 *
 * There are better ways to handle this but this works in simple
 * cases.  SWIG handles this automatically, if you set it up exactly
 * right, which is not easy.
 *
 */

static PyTypeObject
ptrobjType = {
    PyObject_HEAD_INIT(NULL)
#if PY_MAJOR_VERSION >= 3
#else
    0,                         /* ob_size */
#endif
    "ptrobj",                  /* tp_name */
    sizeof(pyactpol_ptrobj),   /* tp_basicsize */
    0,                         /* tp_itemsize */
    0,                         /* tp_dealloc */
    0,                         /* tp_print */
    0,                         /* tp_getattr */
    0,                         /* tp_setattr */
    0,                         /* tp_compare */
    0,                         /* tp_repr */
    0,                         /* tp_as_number */
    0,                         /* tp_as_sequence */
    0,                         /* tp_as_mapping */
    0,                         /* tp_hash */
    0,                         /* tp_call */
    0,                         /* tp_str */
    0,                         /* tp_getattro */
    0,                         /* tp_setattro */
    0,                         /* tp_as_buffer */
    Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE, /* tp_flags*/
    "ptrobj object",           /* tp_doc */
    0,                         /* tp_traverse */
    0,                         /* tp_clear */
    0,                         /* tp_richcompare */
    0,                         /* tp_weaklistoffset */
    0,                         /* tp_iter */
    0,                         /* tp_iternext */
    0,                         /* tp_methods */
    0,                         /* tp_members */
    0,                         /* tp_getset */
    0,                         /* tp_base */
    0,                         /* tp_dict */
    0,                         /* tp_descr_get */
    0,                         /* tp_descr_set */
    0,                         /* tp_dictoffset */
    0,                         /* tp_init */
    0,                         /* tp_alloc */
    PyType_GenericNew,         /* tp_new */
    0,                         /* tp_free */
    0,                         /* tp_is_gc */
    0,                         /* tp_bases */
    0,                         /* tp_mro */
    0,                         /* tp_cache */
    0,                         /* tp_subclasses */
    0,                         /* tp_weaklist */
};


/* ptrobj_decode
 *
 * For use in Py_ParseTuple -- assuming that the passed object is
 * really a ptrobj, extract the encapsulated pointer and store it in
 * dest.
 */

int pyactpol_ptrobj_decode(PyObject *o, void **dest)
{
    *dest = ((pyactpol_ptrobj*)o)->p;
    return 1;
}

/* ptrobj_new
 *
 * Create a new ptrobj, with data initialized to the "item".
 */

pyactpol_ptrobj *pyactpol_ptrobj_new(void *item)
{
    pyactpol_ptrobj* p = PyObject_New(pyactpol_ptrobj, &ptrobjType);
    p->p = item;
    return p;
}


/******
 * dirfile access
 *
 * Wraps libactpol calls.
 */


static PyObject *dirfile_open(PyObject *self, PyObject *args)
{
    char *filename;
    if (!PyArg_ParseTuple(args, "s",
                          &filename))
        return NULL;

    ACTpolDirfile *df = ACTpolDirfile_open(filename);
    if (df == NULL)
        Py_RETURN_NONE;

    return (PyObject*) pyactpol_ptrobj_new(df);
}


static PyObject *dirfile_close(PyObject *self, PyObject *args)
{
    ACTpolDirfile *df;
    if (!PyArg_ParseTuple(args, "O&",
                          pyactpol_ptrobj_decode, &df))
        return NULL;
    ACTpolDirfile_close(df);
    Py_RETURN_NONE;
}


/*
 * dirfile_get_frame_count
 *
 * Count frames and partial frames in a dirfile.
 *
 * In:
 *   dirfile (dirfile_object)
 *   channel_probe (string, possibly null) - channel to use to count partial
 *     frames
 *   
 * Out:
 *   n_full_frames (long)
 *   samples_per_frame (long)
 *   n_trailing_samples (long)
 */

static PyObject *dirfile_get_frame_count(PyObject *self, PyObject *args)
{
    ACTpolDirfile *dirfile;
    char *channel;
    if (!PyArg_ParseTuple(args, "O&s",
                          &pyactpol_ptrobj_decode, &dirfile,
                          &channel))
        Py_RETURN_NONE;

    po_assert(dirfile != NULL && channel != NULL);

    // Number of complete frames
    int status;
    long n_frames = GetNFrames(dirfile->format, &status, NULL);

    // Number of partial frames, as determined from channel.
    long spf = 0;
    long n_partial = 0;
    if (channel != NULL) {
        // Check the channel length, then ask for more.
        // You should use a real data channel, not INDEX.
        spf = GetSamplesPerFrame(dirfile->format, channel, &status);
        char *data = malloc((n_frames+1) * spf * sizeof(char));
        long nout = GetData(dirfile->format, channel,
                            0,0,
                            n_frames, spf-1,
                            'c', data,
                            &status);
        n_partial = nout - n_frames*spf;
        free(data);
    }

    return Py_BuildValue("lll",
                         n_frames,
                         spf,
                         n_partial);
}

/*
 * dirfile_get_channel_info
 *
 * Check for channel existence and get basic channel information.
 *
 * In:
 *   dirfile (dirfile_object)
 *   channel name (string)
 *
 * Out:
 *   exists (boolean)
 *   samples_per_frame (long)
 *   
 */

static PyObject *dirfile_get_channel_info(PyObject *self, PyObject *args)
{
    ACTpolDirfile *dirfile;
    char *channel;
    if (!PyArg_ParseTuple(args, "O&s",
                          &pyactpol_ptrobj_decode, &dirfile,
                          &channel))
        Py_RETURN_NONE;

    if (dirfile == NULL || channel == NULL) {
        Py_RETURN_NONE;
    }

    int status;
    bool exists = false;
    long spf = GetSamplesPerFrame(dirfile->format, channel, &status);
    exists = (status == GD_E_OK && spf > 0);
    if (!exists)
        spf = 0;

    return Py_BuildValue("Nl",
                         PyBool_FromLong((long)exists),
                         spf);
}


/*
 * dirfile_load_channel
 *
 * Load a dirfile channel.  Inputs should be sanity checked *before*
 * you get here.  E.g., to decode negative start and count parameters,
 * determine the best data type, etc.  It's easier to do that on the
 * python side.
 *
 * In:
 *   dirfile (dirfile_object)
 *   channel name (string)
 *   data_type (string)
 *   sample_start (long) - first sample to load
 *   sample_count (long) - number of samples to load
 *   dest (ndarray, possibly null) - place for the data.  Optional.
 *
 * Out:
 *   data (ndarray)
 *   
 */

static char dtype_to_typecode(int dtype) {
    switch(dtype) {
    case NPY_UINT8:
        return 'c';
    case NPY_INT16:
        return 's';
    case NPY_UINT16:
        return 'u';
    case NPY_INT32:
        return 'i';
    case NPY_UINT32:
        return 'U';
    case NPY_FLOAT32:
        return 'f';
    case NPY_FLOAT64:
        return 'd';
    }
    print_error("request for unhandled numpy type, %i\n", dtype);
    return 0;
}


static PyObject *dirfile_load_channel(PyObject *self, PyObject *args)
{
    ACTpolDirfile *dirfile;
    char *channel;
    long sample0, n_samples;
    PyArray_Descr* dtype;
    PyObject *dest_in;
    PyArrayObject *dest_array;
    if (!PyArg_ParseTuple(args, "O&sllO&O",
                          &pyactpol_ptrobj_decode, &dirfile,
                          &channel,
                          &sample0,
                          &n_samples,
                          PyArray_DescrConverter, &dtype, /* magic */
                          &dest_in))
        Py_RETURN_NONE;

    if (dirfile == NULL || channel == NULL) {
        Py_RETURN_NONE;
    }

    if (dest_in == Py_None) {
        /* Create a new array to hold the data. */
        npy_intp dims[1];
        dims[0] = n_samples;
        dest_array = (PyArrayObject*)PyArray_SimpleNewFromDescr(1, dims, dtype);
    } else if (PyArray_Check(dest_in)) {
        /* Drop data into existing array */
        dest_array = (PyArrayObject*)dest_in;
        dtype = PyArray_DESCR(dest_array);
        /* Check array dimensions and strides */
        ASSERT_CARRAY_TYPE_NDIM(dest_array, PyArray_TYPE(dest_array), 1);
        npy_intp *dims = PyArray_DIMS(dest_array);
        po_assert(dims[0] >= n_samples);
        /* Get a new reference, since we're returning this like new. */
        Py_INCREF(dest_array);
    } else {
        po_raise("unexpected object type");
    }
    void *data = PyArray_DATA(dest_array);

    int status = GD_E_OK;
    long nsamples_out = 0;
    /* Release Global Interpreter Lock during I/O */
    Py_BEGIN_ALLOW_THREADS
    nsamples_out = GetData(dirfile->format, channel,
                           0, sample0,
                           0, n_samples,
                           dtype_to_typecode(dtype->type_num),
                           data, &status);
    Py_END_ALLOW_THREADS

    if (status != GD_E_OK) {
        print_error("status = %i : %s\n", status, GD_ERROR_CODES[status]); 
        po_raise("GetData error");
    }
    if (nsamples_out != n_samples) {
        print_error("requested %i samples but only read %i\n",
                    n_samples, nsamples_out);
        po_raise("Dirfile error");
    }
    return Py_BuildValue("N",
                         dest_array);
}

struct converter {
    int raw;
    long mask;
    int is_signed;
    float scale;
};

/* decode_converter expects to encounter something like a BitField
 * object, as defined in util.mce. 
 */

int decode_converter(struct converter *conv, PyObject *convo)
{
    memset(conv, 0, sizeof(*conv));
    if (convo == NULL || convo == Py_None) {
        conv->raw = 1;
        return 0;
    }
    PyObject *start = PyObject_GetAttrString(convo, "start");
    PyObject *count = PyObject_GetAttrString(convo, "count");
    PyObject *is_signed = PyObject_GetAttrString(convo, "signed");
    PyObject *scale = PyObject_GetAttrString(convo, "scale");
    
    unsigned long bit = 1 << PyInt_AsLong(start);
    conv->scale = bit;
    for (int i=0; i < PyInt_AsLong(count); i++) {
        conv->mask |= bit;
        bit = bit << 1;
    }
    conv->scale = PyFloat_AsDouble(scale) / conv->scale;
    conv->is_signed = (int)PyInt_AsLong(is_signed);
    return 0;
}


static PyObject *dirfile_load_channels(PyObject *self, PyObject *args)
{
    ACTpolDirfile *dirfile;
    PyObject *field_list;
    long sample0, n_samples;
    PyArray_Descr* dtype;
    PyObject *dest_in;
    PyArrayObject *dest_array;
    PyObject *convo;
    struct converter conv;

    if (!PyArg_ParseTuple(args, "O&OllO&OO",
                          &pyactpol_ptrobj_decode, &dirfile,
                          &field_list,
                          &sample0,
                          &n_samples,
                          PyArray_DescrConverter, &dtype, /* magic */
                          &dest_in,
                          &convo))
        po_raise("invalid arguments.");

    if (dirfile == NULL) {
        po_raise("invalid dirfile object");
    }

    po_assert(PyList_Check(field_list));

    if (dest_in == Py_None) {
        /* Create a new array to hold the data. */
        npy_intp dims[2];
        dims[0] = PyList_Size(field_list);
        dims[1] = n_samples;
        dest_array = (PyArrayObject*)PyArray_SimpleNewFromDescr(2, dims, dtype);
    } else if (PyArray_Check(dest_in)) {
        /* Drop data into existing array */
        dest_array = (PyArrayObject*)dest_in;
        dtype = PyArray_DESCR(dest_array);
        /* Check array dimensions and strides */
        ASSERT_CARRAY_TYPE_NDIM(dest_array, PyArray_TYPE(dest_array), 2);
        npy_intp *dims = PyArray_DIMS(dest_array);
        po_assert(dims[0] == PyList_Size(field_list));
        po_assert(dims[1] == n_samples);
        /* Get a new reference, since we're returning this like new. */
        Py_INCREF(dest_array);
    } else {
        po_raise("unexpected object type");
    }
    void *data = PyArray_DATA(dest_array);

    int any_errors = 0;
    long major_stride = PyArray_STRIDES(dest_array)[0]; // Yes, it's bytes.

    char gd_type = dtype_to_typecode(dtype->type_num);

    decode_converter(&conv, convo);
    if (!conv.raw) {
        po_assert(dtype->type_num == NPY_FLOAT32);
        gd_type = 'U';
    }
    
    /* Release Global Interpreter Lock during I/O */
    Py_BEGIN_ALLOW_THREADS

#pragma omp parallel shared(any_errors)
    // Per-thread buffer is needed if we're translating the data.
    {
    void *buf = NULL;
    if (!conv.raw)
        buf = malloc(n_samples * sizeof(uint32_t));
    
#pragma omp for
    for (long i=0; i<PyList_Size(field_list); i++) {
        const char *channel = PyString_AsString(PyList_GET_ITEM(field_list, i));
        int status = GD_E_OK;
        long nsamples_out = 0;
        if (conv.raw)
            buf = data + i * major_stride;
        nsamples_out = GetData(dirfile->format, channel,
                               0, sample0, 0, n_samples,
                               gd_type, buf,
                               &status);
        if (!conv.raw) {
            if (conv.is_signed) 
                extract_float32_int32(data + i*major_stride, buf, n_samples,
                                      conv.mask, conv.scale);
            else
                extract_float32_uint32(data + i*major_stride, buf, n_samples,
                                       conv.mask, conv.scale);
        }
        if (nsamples_out != n_samples) {
            print_error("Field %s: requested %i samples but only read %i\n",
                        channel, n_samples, nsamples_out);
            any_errors = 1;
        }
        if (status != GD_E_OK) {
            print_error("Field %s: status = %i : %s\n", channel, 
                        status, GD_ERROR_CODES[status]);
            any_errors = 1;
        }
    }
    }
    
    Py_END_ALLOW_THREADS

    if (any_errors)
        po_raise("Dirfile error");

    return Py_BuildValue("N", dest_array);
}


int extract_float32_int32(float *dest, int32_t *src, long n,
                          long mask, float scale)
{
    int32_t _mask = mask;
    for (long i=0; i<n; i++)
        dest[i] = (float)(src[i] & _mask) * scale;
    return n;
}

int extract_float32_uint32(float *dest, uint32_t *src, long n,
                           long mask, float scale)
{
    uint32_t _mask = mask;
    for (long i=0; i<n; i++)
        dest[i] = (float)(src[i] & _mask) * scale;
    return n;
}



PyMethodDef pyactpol_dirfile_methods[] = {
    {"dirfile_open", dirfile_open, METH_VARARGS,
     ""},
    {"dirfile_close", dirfile_close, METH_VARARGS,
     ""},
    {"dirfile_get_channel_info", dirfile_get_channel_info, METH_VARARGS,
     "Given (dirfile, channel), returns (existence, samples_per_frame)."},
    {"dirfile_get_frame_count", dirfile_get_frame_count, METH_VARARGS,
     "Given (dirfile, channel), returns (n_samples, 0)."},
    {"dirfile_load_channel", dirfile_load_channel, METH_VARARGS,
     "Load dirfile channel."},
    {"dirfile_load_channels", dirfile_load_channels, METH_VARARGS,
     "Load multiple dirfile channels (openmp)."},
    {NULL, NULL, 0, NULL}        /* Sentinel */
};

void pyactpol_dirfile_init()
{
    if (PyType_Ready(&ptrobjType) < 0)
        return;
    Py_INCREF(&ptrobjType);
}
	
