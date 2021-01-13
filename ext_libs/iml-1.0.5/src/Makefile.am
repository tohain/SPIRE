# See http://www.gnu.org/software/libtool/manual/html_node/Updating-version-info.html
LIBIML_LT_CURRENT=1
LIBIML_LT_REVISION=1
LIBIML_LT_AGE=1

AM_CFLAGS = $(GMP_CFLAGS) $(CBLAS_CFLAGS) $(ATLAS_CFLAGS)

EXTERNLIB = $(GMP_LIBS) $(CBLAS_LIBS) $(ATLAS_LIBS) 

include_HEADERS = iml.h

noinst_HEADERS = basisop.h \
		 certsolve.h \
		 common.h \
		 error.h \
		 latreduce.h \
		 mtrans.h \
		 nonsysolve.h \
		 padiclift.h \
		 reconsolu.h \
		 RNSop.h \
		 nullspace.h

lib_LTLIBRARIES = libiml.la
libiml_la_SOURCES = basisop.c \
	            certsolve.c \
	            error.c \
	            latreduce.c \
	            memalloc.c \
	            mtrans.c \
	   	    nonsysolve.c \
		    padiclift.c \
		    reconsolu.c \
		    RNSop.c \
	            nullspace.c

libiml_la_CFLAGS = $(AM_CFLAGS)
libiml_la_LIBADD = $(EXTERNLIB)
libiml_la_LDFLAGS = -version-info $(LIBIML_LT_CURRENT):$(LIBIML_LT_REVISION):$(LIBIML_LT_AGE) \
                    $(LIBIML_LDFLAGS)
