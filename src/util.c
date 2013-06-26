#include "pyactpol.h"

#include <stdarg.h>

int pyactpol_print_error(const char *format, ... )
{
	va_list args;
	va_start(args, format);
	int i = vfprintf(stderr, format, args);
	va_end(args);
	return i;
}
