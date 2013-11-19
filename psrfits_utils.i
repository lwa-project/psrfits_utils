%module psrfits_utils
%{
#include "psrfits.h"
#include "swig_addedfunc.h"
%}
%include "cmalloc.i"
%malloc(float,floatp)
%free(float,floatp)
%malloc(double,doublep)
%free(double,doublep)
%malloc(unsigned char,ucharp)
%free(unsigned char,ucharp)

%apply double { long double }

%include "psrfits.h"
%include "swig_addedfunc.h"
