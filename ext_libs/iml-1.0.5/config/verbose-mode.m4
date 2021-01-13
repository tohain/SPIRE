dnl IML_VERBOSE_MODE
dnl ask program to output additional information

AC_DEFUN([IML_VERBOSE_MODE],
[
AC_ARG_ENABLE(verbose-mode, 
 		[  --enable-verbose-mode
		  The program will output additional information such as
		  timings, intermidate results and so on.
		],
		[
		   AC_DEFINE([HAVE_VERBOSE_MODE], [1], [make verbose output])
		],
		[])
])
