dnl Autoconf

AC_PREREQ(2.59)
AC_INIT(IML, 1.0.4)
AC_CONFIG_MACRO_DIR([config])
AC_CONFIG_HEADERS(config.h:config-h.in)
AC_CONFIG_SRCDIR([src/RNSop.c])
AM_INIT_AUTOMAKE

# Adds special linker flag for icc
case "x$CC" in
      x*icc) LDFLAGS="$LDFLAGS -lm -no_cpprt" 
	     ;;
      *)     LDFLAGS="$LDFLAGS -lm"
	     ;;
esac

# Checks for programs.
AC_PROG_CC
AM_PROG_CC_C_O

# Checks for libraries.


# Checks for header files.
AC_HEADER_STDC

# support windows libraries, disable shared library by default
LT_INIT([win32-dll, disable-shared])

# support shared library on Cygwin
AC_CANONICAL_HOST
case $host in
  *-*-cygwin* | *-*-mingw*)
    # By default, build only static.
    if test -z "$enable_shared"; then
      enable_shared=no
    fi
    # "-no-undefined" is required when building a DLL, see documentation on
    # AC_LIBTOOL_WIN32_DLL.
    #
    # "-Wl,--export-all-symbols" is a bit of a hack, it gets all functions
    # and variables exported. This is what libtool did in the past, and
    # it's convenient.
    if test "$enable_shared" = yes; then
      LIBIML_LDFLAGS="$LIBIML_LDFLAGS -no-undefined -Wl,--export-all-symbols"
    fi
    ;;
esac

DEFAULT_CFLAGS="-O3 -Wall"
CFLAGS="-I. ${CFLAGS:-$DEFAULT_CFLAGS}"

AC_CHECK_HEADERS([stdlib.h time.h])

AC_CHECK_HEADERS([math.h],, [echo 'error: could not find required header math.h, configuration aborted.'  exit 1])

# Checks for typedefs, structures, and compiler characteristics.
AC_C_CONST
AC_TYPE_SIZE_T

# Checks for library functions.
AC_FUNC_REALLOC

AC_CHECK_FUNCS(bzero memset, break)

AC_CHECK_FUNCS(calloc floor)

IML_VERBOSE_MODE
IML_DEFAULT_PATH
IML_CHECK_GMP(,,[AC_MSG_ERROR(
GMP not found! 
GMP version 3.1.1 or greater is required for this library to compile. 
Please make sure GMP is installed and specify the header and libraries 
location with the options --with-gmp-include=<path> and --with-gmp-lib=<path>
respectively when running configure.
)])

IML_CHECK_CBLAS(,,[AC_MSG_WARN(
CBLAS not found!
Please make sure that --with-cblas=<linker flags> and optional --with-cblas-include=<path> and --with-cblas-lib=<path> are correctly set.
Trying legacy ATLAS configuration.)
IML_CHECK_ATLAS(,,[AC_MSG_ERROR(
ATLAS not found! 
ATLAS version 3.0 or greater is required for this library to compile. Please make sure ATLAS is installed and specify the header and libraries location with the options --with-atlas-include=<path> and --with-atlas-lib=<path> respectively when running configure.)
])
])

AC_SUBST(LIBIML_LDFLAGS)
AC_SUBST(LDFLAGS)
AC_SUBST(ac_aux_dir)
	
AC_OUTPUT([
	Makefile
	config/Makefile
	src/Makefile
	doc/Makefile
	examples/Makefile
	tests/Makefile
])		

echo \
"------------------------------------------------------------------------
Configuration:
 
  Source code location:       ${srcdir}
  Compiler:                   ${CC}
  Compiler flags:             ${CFLAGS}
  Linker flags:               ${LDFLAGS}
  Host System Type:           ${host}
  Install path:               ${prefix}
 
  See config.h for further configuration information.
------------------------------------------------------------------------"
