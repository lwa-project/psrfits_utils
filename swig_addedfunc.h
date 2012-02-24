#include "Python.h"

void  convert2_float_array(float *output, PyObject *input, int size);
void print_float_array(float *input, int size);
unsigned char *convert_uchar_array(unsigned char *result, PyObject *input, int size);
void print_uchar_array(unsigned char *input, int size);
long double get_ld(double in);
