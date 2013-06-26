/* -*- mode: C; tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 4 -*-
 *      vim: sw=4 ts=8 et tw=80
 */

/* main.c 
 *
 * The module init code.
 */


#define PYACTPOL_MAIN
#include "pyactpol.h"

#define PKGNAME "libactpol"

/* copy_methods -- assists with combining a bunch of PyMethodDef
 * arrays into a single array. */

static int copy_methods(const PyMethodDef *methods, PyMethodDef *dest)
{
    int n=0;
    while (methods[n].ml_name!=NULL) {
        if (dest != NULL)
            memcpy(dest++, methods+n, sizeof(*dest));
        n++;
    }
    return n;
}


/* Array of method arrays for different submodules.  These will be
 * combined in the module init function.
 */

static PyMethodDef *my_method_sets[] = {
    pyactpol_cuts_methods,
    pyactpol_dirfile_methods,
    pyactpol_fft_methods,
    pyactpol_filter_methods,
    pyactpol_linalg_methods,
    pyactpol_lookup_methods,
    pyactpol_mce_methods,
    pyactpol_pointing_methods,
    pyactpol_sync_methods,
    pyactpol_waterfall_methods,
    NULL
};


/* Array of init functions for different source files.  This is better
 * than exposing lots of symbols, but it's barely necessary.
 */

typedef void (*init_function)(void);

init_function my_init_functions[] = {
    pyactpol_dirfile_init,
    NULL
};

#if PY_MAJOR_VERSION >= 3
static struct PyModuleDef libactpol_module = {
    PyModuleDef_HEAD_INIT,
    "libactpol",
    NULL,
    -1,
    NULL
};
#endif

PyMODINIT_FUNC
#if PY_MAJOR_VERSION >= 3
PyInit_libactpol(void)
#else
initlibactpol(void)
#endif
{
    // Synthesize method tables.
    int n, i;
    // Pass 1, count them
    n = 0;
    for (i=0; my_method_sets[i] != NULL; i++)
        n += copy_methods(my_method_sets[i], NULL);
    n++; // terminator
    PyMethodDef *methods = malloc(n * sizeof(PyMethodDef));

    // Pass 2, copy them
    n = 0;
    for (i=0; my_method_sets[i] != NULL; i++)
        n += copy_methods(my_method_sets[i], methods+n);
    memset(methods+n, 0, sizeof(*methods)); //terminator

    // Register with python central command
#if PY_MAJOR_VERSION >= 3
    libactpol_module.m_methods = methods;
    PyObject *module = PyModule_Create(&libactpol_module);
    if (module == NULL)
        return NULL;
#else
    Py_InitModule("libactpol", methods);
#endif
    // Numpy needs this.
    import_array();

    // Call any weird per-file init functions.
    init_function* init;
    for (init = my_init_functions; (*init)!=NULL; init++)
        (*init)();

#if PY_MAJOR_VERSION >= 3
    return module;
#endif
}
