# Check for GMP
# Zhuliang Chen, 2004-07
# Modified by Pascal Giorgi, 2003-12-03

dnl IML_CHECK_GMP ([MINIMUM-VERSION [, ACTION-IF-FOUND [, ACTION-IF-NOT-FOUND]]])
dnl
dnl Test for the GNU Multiprecision library and define GMP_CFLAGS and GMP_LIBS

AC_DEFUN([IML_CHECK_GMP],
[

AC_ARG_WITH(gmp-include,
	     [  --with-gmp-include=<path>
		Specify the path of GMP header gmp.h
		If argument is <empty> that means the header is reachable 
		with the default search path: /usr/include or /usr/local/include
		or other path you add as default path using --with-default.
		Otherwise you give the <path> to the directory which contains 
		the header gmp.h. 
	     ],
             [  
		GMP_INPUT_HEADER=$withval
		GMP_HEADER_PATH="${GMP_INPUT_HEADER}"
		for GMP_HOME in ${DEFAULT_CHECKING_PATH} 
		do
			GMP_HEADER_PATH="${GMP_HEADER_PATH} ${GMP_HOME}/include"
		done
	     ],
	     [  
		for GMP_HOME in ${DEFAULT_CHECKING_PATH} 
		do
			GMP_HEADER_PATH="${GMP_HEADER_PATH} ${GMP_HOME}/include"
		done
	    ])

AC_ARG_WITH(gmp-lib,
	     [  --with-gmp-lib=<path>
		Specify the path of GMP library.
	        If argument is <empty> that means the header is reachable 
                with the default search path: /usr/include or /usr/local/include
		or path you add as default path using --with-default.
		Otherwise you give the <path> to the directory which contains 
		the GMP library.
	     ],
	     [  
		GMP_INPUT_LIB=$withval
		GMP_LIBRARY_PATH="${GMP_INPUT_LIB}"
		for GMP_HOME in ${DEFAULT_CHECKING_PATH} 
		do
			GMP_LIBRARY_PATH="${GMP_LIBRARY_PATH} ${GMP_HOME}/lib"
		done
	     ],
	     [  
		for GMP_HOME in ${DEFAULT_CHECKING_PATH} 
		do
			GMP_LIBRARY_PATH="${GMP_LIBRARY_PATH} ${GMP_HOME}/lib"
		done
	    ])


dnl input data check
if test  "x${GMP_INPUT_HEADER}" != x; then
	if test -z "${GMP_INPUT_LIB}"; then
		echo 'error: since you have specified the GMP header path, you must specify the'
		echo 'GMP library path using --with-gmp-lib.'
		exit 1
	fi
fi

if test "x${GMP_INPUT_LIB}" != x; then
	if test -z "${GMP_INPUT_HEADER}"; then
		echo 'error: since you have specified the GMP library path, you must specify the'
		echo 'GMP header path using --with-gmp-include.'
		exit 1
	fi
fi

min_gmp_version=ifelse([$1], ,3.1.1,$1)

dnl Check for existence
BACKUP_CFLAGS=${CFLAGS}
BACKUP_LIBS=${LIBS}

AC_MSG_CHECKING(for GMP >= $min_gmp_version)

gmp_found="no"


for GMP_HEADER in ${GMP_HEADER_PATH}
do	
	if test -r "${GMP_HEADER}/gmp.h"; then
		if test "x${GMP_HEADER}" != "x/usr/include" -a "x${GMP_HEADER}" != "x/usr/local/include"; then
			for GMP_LIBRARY in ${GMP_LIBRARY_PATH}
			do
				if test "x${GMP_LIBRARY}" != "x/usr/lib" -a "x${GMP_LIBRARY}" != "x/usr/local/lib"; then
					GMP_CFLAGS="-I${GMP_HEADER}"
					GMP_LIBS="${GMP_LIBS} -L${GMP_LIBRARY} -lgmp"	
				fi
			done
		else
			GMP_CFLAGS=
			GMP_LIBS="-lgmp"		
		fi

		CFLAGS="${CFLAGS} ${GMP_CFLAGS}"
		LIBS="${LIBS} ${GMP_LIBS}"


		AC_TRY_LINK(
		[#include <gmp.h>],
		[mpz_t a; mpz_init (a);],
		[
        		AC_TRY_RUN(
 			[#include <gmp.h>
			 int main() { if (__GNU_MP_VERSION < 3) return -1;
	 			      else return 0; }
		  	],[
				AC_MSG_RESULT(found)
				gmp_found="yes"
				AC_SUBST(GMP_CFLAGS)
		  		AC_SUBST(GMP_LIBS)
				AC_DEFINE(HAVE_GMP,1,[Define if GMP is installed])
				# See if we are running GMP 4.0
	   			AC_MSG_CHECKING(whether GMP is 4.0 or greater)
		   		AC_TRY_RUN(
		   		[#include <gmp.h>
	    			int main () { if (__GNU_MP_VERSION < 4) return -1; else return 0; }
	   			],[
					AC_MSG_RESULT(yes)
				],[
					AC_MSG_RESULT(no)
					AC_DEFINE(GMP_VERSION_3,1,[Define if GMP is version 3.xxx])
					GMP_VERSION="-DGMP_VERSION_3"
					AC_SUBST(GMP_VERSION)
				],[
					dnl This should never happen
					AC_MSG_RESULT(no)
				])
				ifelse([$2], , :, [$2])
				break
			],[			
				gmp_problem="$gmp_problem $GMP_HEADER"
				unset GMP_CFLAGS
				unset GMP_LIBS	
			],[
				AC_MSG_RESULT(unknown)
				gmp_found="yes"
				echo 'You appear to be cross compiling, so there is'
				echo 'no way to determine whether your GMP version is new enough.'
				echo 'I am assuming it is'
				AC_SUBST(GMP_CFLAGS)
				AC_SUBST(GMP_LIBS)
				AC_DEFINE(HAVE_GMP,1,[Define if GMP is installed])	
				ifelse([$2], , :, [$2])
				break
			])	
		],[
		unset GMP_CFLAGS
		unset GMP_LIBS	
		])
	fi
done

if test "x$gmp_found" = "xno"; then
	if test -n "$gmp_problem"; then
		AC_MSG_RESULT(problem: your GMP version is too old. Disabling.)
	else
		AC_MSG_RESULT(not found)
	fi
	ifelse($3, , :, $3)
fi

CFLAGS=${BACKUP_CFLAGS}
LIBS=${BACKUP_LIBS}
#unset LD_LIBRARY_PATH

])
