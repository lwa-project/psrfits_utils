AC_INIT([psrfits_utils], [1.0])
AC_CONFIG_SRCDIR([downsample.c])
AC_CONFIG_MACRO_DIR([config])
AC_CONFIG_AUX_DIR([config])
AM_INIT_AUTOMAKE([-Wall -Werror foreign])

m4_ifdef([AM_PROG_AR], [AM_PROG_AR])
AC_PROG_CC
AC_PROG_LIBTOOL
AM_PATH_PYTHON(2.3)
AX_PKG_SWIG(1.3.17, [], [ AC_MSG_ERROR([SWIG is required to build..]) ])
AX_SWIG_MULTI_MODULE_SUPPORT
AX_SWIG_PYTHON


AX_PTHREAD

AC_CHECK_LIB([m],[sin])
AC_CHECK_LIB([cfitsio],[ffgky],[],
	     [AC_MSG_ERROR([cfitsio is required])])

AC_ARG_ENABLE([sse],
	      [AS_HELP_STRING([--enable-sse], 
			      [enable SSE folding code])],
	      [],
	      [enable_sse=no])
if test "x$enable_sse" = "xyes"; then
	AC_DEFINE([FOLD_USE_INTRINSICS], [1], [Use SSE intrinsics for folding])
fi

AC_CONFIG_FILES([Makefile])
AC_OUTPUT
