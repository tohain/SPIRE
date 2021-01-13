# Check for ATLAS
# Zhuliang Chen, 2004-07
# Pascal Giorgi , 2003-03-04
# Inspired by gnome-bonobo-check.m4 by Miguel de Icaza, 99-04-12
# Stolen from Chris Lahey       99-2-5
# stolen from Manish Singh again
# stolen back from Frank Belew
# stolen from Manish Singh
# Shamelessly stolen from Owen Taylor

dnl IML_CHECK_ATLAS ([MINIMUM-VERSION [, ACTION-IF-FOUND [, ACTION-IF-NOT-FOUND]]])
dnl
dnl Test for ATLAS and define ATLAS_CFLAGS and ATLAS_LIBS

AC_DEFUN([IML_CHECK_ATLAS],
[

AC_ARG_WITH(atlas-include,
	     [  --with-atlas-include=<path>
		Specify the path of ATLAS header cblas.h
		If argument is <empty> that means the header is reachable 
            	with the default search path: /usr/include or /usr/local/include 
		or other path you add as default path using --with-default.
		Otherwise you give the <path> to the directory which contains
		the ATLAS header cblas.h. 
	      ],
	      [ 
		ATLAS_INPUT_HEADER=$withval
		ATLAS_HEADER_PATH="${ATLAS_INPUT_HEADER}"
		for ATLAS_HOME in ${DEFAULT_CHECKING_PATH} 
		do
			ATLAS_HEADER_PATH="${ATLAS_HEADER_PATH} ${ATLAS_HOME}/include"
		done
	      ],
	      [ 
		for ATLAS_HOME in ${DEFAULT_CHECKING_PATH} 
		do
			ATLAS_HEADER_PATH="${ATLAS_HEADER_PATH} ${ATLAS_HOME}/include"
		done
	    ])

AC_ARG_WITH(atlas-lib,
	    [  --with-atlas-lib=<path>
		Specify the path of ATLAS library.
	        If argument is <empty> that means the header is reachable 
                with the default search path: /usr/include or /usr/local/include
		or path you add as default path using --with-default.
		Otherwise you give the <path> to the directory which contains 
		the ATLAS library.
	     ],
	     [  
		ATLAS_INPUT_LIB=$withval
		ATLAS_LIBRARY_PATH="${ATLAS_INPUT_LIB}"
		for ATLAS_HOME in ${DEFAULT_CHECKING_PATH} 
		do
			ATLAS_LIBRARY_PATH="${ATLAS_LIBRARY_PATH} ${ATLAS_HOME}/lib"
		done
	     ],
	     [  
		for ATLAS_HOME in ${DEFAULT_CHECKING_PATH} 
		do
			ATLAS_LIBRARY_PATH="${ATLAS_LIBRARY_PATH} ${ATLAS_HOME}/lib"
		done
	    ])


dnl input data check
if test  "x${ATLAS_INPUT_HEADER}" != x; then
	if test -z "${ATLAS_INPUT_LIB}"; then
		echo 'error: since you have specified the ATLAS header path, you must specify the'
		echo 'ATLAS library path using --with-atlas-lib.'
		exit 1
	fi
fi

if test "x${ATLAS_INPUT_LIB}" != x; then
	if test -z "${ATLAS_INPUT_HEADER}"; then
		echo 'error: since you have specified the ATLAS library path, you must specify the'
		echo 'ATLAS header path using --with-atlas-include.'
		exit 1
	fi
fi


min_atlas_version=ifelse([$1], ,3.0,$1)

dnl Check for existence

BACKUP_CFLAGS=${CFLAGS}
BACKUP_LIBS=${LIBS}

AC_MSG_CHECKING(for ATLAS >= ${min_atlas_version})

atlas_found="no"

for ATLAS_HEADER in ${ATLAS_HEADER_PATH} 
do
	if test -r "${ATLAS_HEADER}/cblas.h"; then
		if test "x${ATLAS_HEADER}" != "x/usr/include" -a "x${ATLAS_HEADER}" != "x/usr/local/include"; then
			for ATLAS_LIBRARY in ${ATLAS_LIBRARY_PATH} 
			do
				if test "x${ATLAS_LIBRARY}" != "x/usr/lib" -a "x${ATLAS_LIBRARY}" != "x/usr/local/lib"; then
					ATLAS_CFLAGS="-I${ATLAS_HEADER}"
					ATLAS_LIBS="${ATLAS_LIBS} -L${ATLAS_LIBRARY} -lcblas -latlas" 
				fi
			done
		else
			ATLAS_CFLAGS=
			ATLAS_LIBS="-lcblas -latlas" 
		fi

		CFLAGS="${BACKUP_CFLAGS} ${ATLAS_CFLAGS} ${GMP_CFLAGS}" 
		LIBS="${BACKUP_LIBS} ${ATLAS_LIBS} ${GMP_LIBS}"

		AC_TRY_LINK(
		[#include <cblas.h>],
		[double a;],
		[
			   atlas_found="yes"	
			   break
		],
		[
		   unset ATLAS_CFLAGS
		   unset ATLAS_LIBS
		   ifelse([$3], , :, [$3])
		])
	fi
done

if test "x$atlas_found" = "xyes" ; then		
	AC_SUBST(ATLAS_CFLAGS)
	AC_SUBST(ATLAS_LIBS)
	AC_DEFINE(HAVE_ATLAS,1,[Define if ATLAS is installed])
	AC_DEFINE(BLAS_AVAILABLE,,[Define if BLAS routines are available])
	HAVE_ATLAS=yes
	if test "x$atlas_cross" != "xyes"; then
		AC_MSG_RESULT(found)
	else
		AC_MSG_RESULT(unknown)
		echo 'WARNING: You appear to be cross compiling, so there is' 
		echo 'no way to determine whether your ATLAS version is new' 
		echo 'enough. I am assuming it is.'
	fi
	ifelse([$2], , :, [$2])
elif test -n "$atlas_problem"; then
	AC_MSG_RESULT(problem, your ATLAS version is too old. Disabling.)
	ifelse([$3], , :, [$3])
elif test "x$atlas_found" = "xno" ; then	
	AC_MSG_RESULT(not found)
	ifelse([$3], , :, [$3])
fi	

CFLAGS=${BACKUP_CFLAGS}
LIBS=${BACKUP_LIBS}
#unset LD_LIBRARY_PATH

])


