#include <Python.h>
#include "swig_addedfunc.h"

void convert2_float_array(float *result, PyObject *input, int size) {
  int i;
  if (!PySequence_Check(input)) {
    printf("Expected a sequence\n");
    exit(EXIT_FAILURE);
  }
  int length=PySequence_Length(input);
  if (length > size) {
    printf("Size mismatch.\n");
    exit(EXIT_FAILURE);
  }
  for (i = 0; i < length; i++) {
    PyObject *o = PySequence_GetItem(input,i);
    if (PyNumber_Check(o)) {
      result[i] = (float) PyFloat_AsDouble(o);
    } else {
      printf("Sequence elements must be numbers\n");
      free(result);       
      exit(EXIT_FAILURE);
    }
    free(o);
  }
}

void print_float_array(float *input,int size)
{
  int i;
  for(i=0;i<size;++i)
  {
    printf("%f\n",input[i]);
  }
}

unsigned char *convert_uchar_array(unsigned char *result,PyObject *input, int size) {
  int i;
  if (!PySequence_Check(input)) {
    printf("Expected a sequence\n");
    exit(EXIT_FAILURE);
  }
  int length=PySequence_Length(input);
  if (length > size) {
    printf("Size mismatch.\n");
    exit(EXIT_FAILURE);
  }
//  result = (unsigned char *) malloc(size*sizeof(unsigned char));
  if(result==NULL)
  {
    fprintf(stderr,"Unable to allocate %d bytes.\n",(size*sizeof(unsigned char)));
    return result;
  }
  for (i = 0; i < length; i++) {
    PyObject *o = PySequence_GetItem(input,i);
    if (PyNumber_Check(o)) {
      result[i] = (unsigned char) PyFloat_AsDouble(o);
    } else {
      PyErr_SetString(PyExc_ValueError,"uchar: Sequence elements must be numbers");
      free(result);       
      return NULL;
    }
    free(o);
  }
  return result;
}

void print_uchar_array(unsigned char *input,int size)
{
  int i;
  for(i=0;i<size;++i)
  {
    printf("%u\n",input[i]);
  }
}

long double get_ld(double in)
{
  long double result=(long double) in;
  return result;
}
