#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT(pcbr, 0.1, rmccloskey@cfenet.ubc.ca)
AC_CONFIG_SRCDIR([src/pcbr.c])
AC_CONFIG_MACRO_DIR([m4])

AM_INIT_AUTOMAKE([-Wall -Werror subdir-objects foreign 1.9.6])

# Checks for programs.
AC_PROG_CC
AC_PROG_INSTALL
AM_PROG_AR

AC_PATH_PROG([R], [R], [yes], [no])
if test "$R" = no; then 
    AC_MSG_ERROR([R is required]) 
fi

# TODO: make an AX macro out of this
AC_CACHE_CHECK([for R package mixdist], [pcbr_cv_prog_Rmixdist_have], [
    AS_IF([R --silent -e "library(mixdist)" 2>&1 > /dev/null], [pcbr_cv_prog_Rmixdist_have=yes], [
         AC_MSG_ERROR([the R package "mixdist" is required])
    ])
])

# Checks for Python and Python packages.
AM_PATH_PYTHON([3.4], [], [:])
AS_IF([test "x$PYTHON" != x:], [
    AX_PYTHON_MODULE([Bio])
    AX_PYTHON_MODULE([pyx])

    AS_IF([test "x$HAVE_PYMOD_BIO" = xno], [
        AC_MSG_WARN([BioPython not found: you will not be able to use the utility scripts])
    ])
    AS_IF([test "x$HAVE_PYMOD_PYX" = xno], [
        AC_MSG_WARN([python module PyX not found: you will not be able to use pcbr_plot])
    ])
],[
    AC_MSG_WARN([python >= 3.4 not found: you will not be able to use the utility scripts])
])
AM_CONDITIONAL([HAVE_PYMOD_BIO], [test "$HAVE_PYMOD_BIO" = yes])
AM_CONDITIONAL([HAVE_PYMOD_PYX], [test "$HAVE_PYMOD_PYX" = yes])

# Compiler and linker stuff
LT_INIT
AX_COMPILER_FLAGS()
AX_GCC_FUNC_ATTRIBUTE([noreturn])

# Environment variables
AC_DEFINE_UNQUOTED([R_HOME], ["$(R RHOME)"], [R home directory])

# Checks for libraries.
AC_CHECK_LIB([m], [cos], [], [AC_MSG_ERROR([math.h is required])])
AC_CHECK_LIB([gslcblas], [cblas_dgemm], [], [AC_MSG_ERROR([cblas is required])])
AC_CHECK_LIB([gsl], [gsl_blas_dgemm], [], [AC_MSG_ERROR([gsl is required])])
PKG_CHECK_MODULES([CHECK], [check >= 0.9.4])
PKG_CHECK_MODULES([R], [libR >= 0.26])
AC_SEARCH_LIBS(pllNewickParseFile, [pll-avx pll-sse3], [], [
    AC_MSG_ERROR([pll is required])                                      
])
AX_PTHREAD()

# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([stdlib.h string.h stddef.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_TYPE_SIZE_T

# Checks for library functions.
AC_FUNC_MALLOC
AC_CHECK_FUNCS([pow sqrt])

# Output files
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_FILES([Makefile
                 C-Thread-Pool/Makefile
                 c-cmaes/Makefile
                 src/Makefile
                 tests/Makefile])
AC_OUTPUT
