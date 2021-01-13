# Check for CBLAS
# Zhuliang Chen, 2004-07
# Pascal Giorgi , 2003-03-04
# Inspired by gnome-bonobo-check.m4 by Miguel de Icaza, 99-04-12
# Stolen from Chris Lahey       99-2-5
# stolen from Manish Singh again
# stolen back from Frank Belew
# stolen from Manish Singh
# Shamelessly stolen from Owen Taylor

dnl IML_CHECK_CBLAS ([MINIMUM-VERSION [, ACTION-IF-FOUND [, ACTION-IF-NOT-FOUND]]])
dnl
dnl Test for CBLAS and define CBLAS_CFLAGS and CBLAS_LIBS

AC_DEFUN([IML_CHECK_CBLAS],
[
AC_ARG_WITH(cblas,
	     [  --with-cblas=<linker flags>
		Specify the flags to use to link to cblas.
		If argument is <empty> that means that ATLAS is installed
                and -lcblas -latlas should be used.
	      ],
	      [ 
		CBLAS_LIBS=$withval
	      ],
	      [ 
		CBLAS_LIBS="-lcblas -latlas"
	    ])


AC_ARG_WITH(cblas-include,
	     [  --with-cblas-include=<path>
		Specify the path of cblas header cblas.h
		If argument is <empty> that means the header is reachable 
            	with the default search path: /usr/include or /usr/local/include 
		or other path you add as default path using --with-default.
                If not the case a default header provided by IML will be used.
		Otherwise you give the <path> to the directory which contains
		the ATLAS header cblas.h. 
	      ],
	      [ 
		CBLAS_INPUT_HEADER=$withval
		CBLAS_HEADER_PATH="${CBLAS_INPUT_HEADER}"
		for CBLAS_HOME in ${DEFAULT_CHECKING_PATH} 
		do
			CBLAS_HEADER_PATH="${CBLAS_HEADER_PATH} ${CBLAS_HOME}/include"
		done
                CBLAS_HEADER_PATH="${CBLAS_HEADER_PATH} ${srcdir}"
	      ],
	      [ 
		for CBLAS_HOME in ${DEFAULT_CHECKING_PATH} 
		do
			CBLAS_HEADER_PATH="${CBLAS_HEADER_PATH} ${CBLAS_HOME}/include"
		done
                CBLAS_HEADER_PATH="${CBLAS_HEADER_PATH} ${srcdir}"
	    ])

AC_ARG_WITH(cblas-lib,
	    [  --with-cblas-lib=<path>
		Specify the path of CBLAS library.
	        If argument is <empty> that means the header is reachable 
                with the default search path: /usr/include or /usr/local/include
		or path you add as default path using --with-default.
		Otherwise you give the <path> to the directory which contains 
		the CBLAS library.
	     ],
	     [  
		CBLAS_INPUT_LIB=$withval
		CBLAS_LIBRARY_PATH="${CBLAS_INPUT_LIB}"
		for CBLAS_HOME in ${DEFAULT_CHECKING_PATH} 
		do
			CBLAS_LIBRARY_PATH="${CBLAS_LIBRARY_PATH} ${CBLAS_HOME}/lib"
		done
	     ],
	     [  
		for CBLAS_HOME in ${DEFAULT_CHECKING_PATH} 
		do
			CBLAS_LIBRARY_PATH="${CBLAS_LIBRARY_PATH} ${CBLAS_HOME}/lib"
		done
	    ])


dnl input data check
dnl if test  "x${CBLAS_INPUT_HEADER}" != x; then
dnl 	if test -z "${CBLAS_INPUT_LIB}"; then
dnl 		echo 'error: since you have specified the CBLAS header path, you must specify the'
dnl 		echo 'CBLAS library path using --with-cblas-lib.'
dnl 		exit 1
dnl 	fi
dnl fi
dnl 
dnl if test "x${CBLAS_INPUT_LIB}" != x; then
dnl 	if test -z "${CBLAS_INPUT_HEADER}"; then
dnl 		echo 'error: since you have specified the CBLAS library path, you must specify the'
dnl 		echo 'CBLAS header path using --with-cblas-include.'
dnl 		exit 1
dnl 	fi
dnl fi

dnl Check for existence

BACKUP_CFLAGS=${CFLAGS}
BACKUP_LIBS=${LIBS}

AC_MSG_CHECKING(for CBLAS)

cblas_found="no"

for CBLAS_HEADER in ${CBLAS_HEADER_PATH} 
do
	if test -r "${CBLAS_HEADER}/cblas.h"; then
		if test "x${CBLAS_HEADER}" != "x/usr/include" -a "x${CBLAS_HEADER}" != "x/usr/local/include"; then
			CBLAS_CFLAGS="-I${CBLAS_HEADER}"
		else
			CBLAS_CFLAGS=
		fi
		for CBLAS_LIBRARY in ${CBLAS_LIBRARY_PATH} 
		do
			if test "x${CBLAS_LIBRARY}" != "x/usr/lib" -a "x${CBLAS_LIBRARY}" != "x/usr/local/lib"; then
				CBLAS_LIBS="-L${CBLAS_LIBRARY} ${CBLAS_LIBS}"
			else
				CBLAS_LIBS="${CBLAS_LIBS}"
			fi
		done

		CFLAGS="${BACKUP_CFLAGS} ${CBLAS_CFLAGS} ${GMP_CFLAGS}" 
		LIBS="${BACKUP_LIBS} ${CBLAS_LIBS} ${GMP_LIBS}"

		AC_TRY_LINK(
		[#include <cblas.h>],
		[double a;],
		[
			   cblas_found="yes"	
			   break
		],
		[
		   unset CBLAS_CFLAGS
		   unset CBLAS_LIBS
		   ifelse([$3], , :, [$3])
		])
	fi
done

if test "x$cblas_found" = "xyes" ; then		
	AC_SUBST(CBLAS_CFLAGS)
	AC_SUBST(CBLAS_LIBS)
	AC_DEFINE(HAVE_CBLAS,1,[Define if CBLAS is installed])
	AC_DEFINE(BLAS_AVAILABLE,,[Define if BLAS routines are available])
	AC_MSG_RESULT(found, using ${CBLAS_LIBS})
	ifelse([$2], , :, [$2])
elif test "x$cblas_found" = "xno" ; then	
	AC_MSG_RESULT(not found)
	ifelse([$3], , :, [$3])
fi	

CFLAGS=${BACKUP_CFLAGS}
LIBS=${BACKUP_LIBS}
#unset LD_LIBRARY_PATH

])


