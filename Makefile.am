ACLOCAL_AMFLAGS = -I config

pkginclude_HEADERS = \
        fold.h \
        polyco.h \
        polyco_struct.h \
        psrfits.h 

noinst_HEADERS = 

dist_data_DATA = \
	guppi_PSRFITS_v3.4_fold_template.txt \
	guppi_PSRFITS_v3.4_search_template.txt

lib_LTLIBRARIES = libpsrfits_utils.la
libpsrfits_utils_la_SOURCES = \
         downsample.c \
         fold.c \
         misc_utils.c \
         polyco.c \
         read_psrfits.c \
         write_psrfits.c 

LDADD = libpsrfits_utils.la
AM_CPPFLAGS = -DPSRFITS_TEMPLATE_DIR='"$(datadir)"'

bin_PROGRAMS = fold_psrfits psrfits_singlepulse psrfits_subband combine_mocks psrfits_bandpass

fold_psrfits_CFLAGS = @PTHREAD_CFLAGS@ $(AM_CFLAGS)
fold_psrfits_LDFLAGS = @PTHREAD_CFLAGS@ $(AM_LDFLAGS)
fold_psrfits_LDADD = libpsrfits_utils.la @PTHREAD_LIBS@

psrfits_subband_SOURCES = psrfits_subband.c \
			  psrfits_subband_cmd.h psrfits_subband_cmd.c

combine_mocks_SOURCES = combine_mocks.c \
				   combine_mocks_cmd.h combine_mocks_cmd.c

psrfits_bandpass_SOURCES = psrfits_bandpass.c \
			  psrfits_bandpass_cmd.h psrfits_bandpass_cmd.c

BUILT_SOURCES = swig_addedfunc.c
SWIG_FILES = psrfits_utils.i

pkgpython_PYTHON = psrfits_utils.py
pkgpyexec_LTLIBRARIES = _psrfits_utils.la

_psrfits_utils_la_SOURCES = swig_addedfunc.c psrfits_utils_wrap.c write_psrfits.c read_psrfits.c
_psrfits_utils_la_CPPFLAGS = $(AX_SWIG_PYTHON_CPPFLAGS)
_psrfits_utils_la_LDFLAGS = -module


swig_addedfunc.c : $(SWIG_FILES)
	$(SWIG) $(SWIG_PYTHON_OPT) $(SWIG_PYTHON_CFLAGS) -o $@ $<
