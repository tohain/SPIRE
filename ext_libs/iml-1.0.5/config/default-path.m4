dnl setting default searching path
dnl modify from linbox

AC_DEFUN([IML_DEFAULT_PATH],
[

 AC_ARG_WITH(default,
	     [  --with-default=<path>
		Add <path> to the default path for package GMP and
		ATLAS searching. The default path without setting this 
		argument is /usr and /usr/local
	     ],
	     [  echo "Default checking path = $withval /usr /usr/local"
	        DEFAULT_CHECKING_PATH="$withval /usr /usr/local"
	     ],
	     [
	         echo "Default checking path = /usr /usr/local"
	         DEFAULT_CHECKING_PATH="/usr /usr/local"
             ])
])
