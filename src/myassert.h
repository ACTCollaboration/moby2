#pragma once

#define STRINGIZE(x) STRINGIZE2(x)
#define STRINGIZE2(x) #x
#define FILE_PLACE __FILE__ ":" STRINGIZE(__LINE__)

/* The po_assert macro is intended for use in functions that return a
 * PyObject.  By setting an error string and returning NULL, an
 * exception will be raised so user gets a stack trace as well as the
 * location of the failed assertion. */

#define po_raise(x) do {				\
		PyErr_SetString(PyExc_RuntimeError,	\
				x " at " FILE_PLACE);	\
		return NULL;				\
	} while(0)

#define po_assert(x) do {					\
		if (!(x)) po_raise("assertion failure");	\
	} while (0)

/* Error messaging */
#define print_error(fmt, args...) do {					\
		pyactpol_print_error(FILE_PLACE ": " fmt, ## args);	\
	} while (0)


/* numpy array assertions
 *
 */

/* CHECK_ARRAY_SIMPLE_1d(array)
 *
 * Fails unless array is 1-d, with stride equal to the item size.
 */

#ifndef NPY_ARRAY_CARRAY
#  define NPY_ARRAY_CARRAY 0x1
#endif 

#define IS_ARRAY_SIMPLE_1D_TYPE(x,type) (				\
		(PyArray_NDIM(x)==1) &&					\
		(PyArray_STRIDES(x)[0] == PyArray_ITEMSIZE(x)) &&	\
		(PyArray_TYPE(x) == type)				\
		)

#define XCHECK_ARRAY_SIMPLE_1D_TYPE(x,type) do {			\
		po_assert( IS_ARRAY_SIMPLE_1D_TYPE(x,type) );		\
	} while (0)

#define XCHECK_CARRAY_AND_TYPE(x,type) do {				\
		po_assert((PyArray_TYPE(x)==type) &&			\
			  ((PyArray_FLAGS(x) & NPY_ARRAY_CARRAY) ==	\
			   NPY_ARRAY_CARRAY));				\
	} while (0)

#define ASSERT_CARRAY_TYPE_NDIM(x,type,nd) do {				\
		po_assert((PyArray_TYPE(x)==type) &&			\
			  ((PyArray_FLAGS(x) & NPY_ARRAY_CARRAY) ==	\
			   NPY_ARRAY_CARRAY) &&				\
			  ((nd<0) || (PyArray_NDIM(x)==nd)));		\
	} while (0)

#define ASSERT_AXIS_N(x,axis,n) \
	po_assert(PyArray_DIMS(x)[axis] == n)
