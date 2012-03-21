%module psrfits_utils
%{
#include "psrfits.h"
#include "swig_addedfunc.h"
%}
%include "cmalloc.i"
%malloc(float,floatp)
%free(float,floatp)
%malloc(unsigned char,ucharp)
%free(unsigned char,ucharp)

%include "psrfits.h"
%include "swig_addedfunc.h"

%apply double * { long double * }
