#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT(thesis, 0.1, rmccloskey@cfenet.ubc.ca, thesis, 
        https://bitbucket.org/rmcclosk/thesis)

AC_DEFUN(AC_REQUIRE_PROG, [
    AC_CHECK_PROG([$1], [$1], [yes])
    if test x"{$1}" = x"no"; then
        AC_MSG_ERROR([$1 is required])
    fi])

AC_DEFUN(AC_REQUIRE_R_PKG,
    [AC_CACHE_CHECK([for R package $1], [my_cv_have_R_pkg_$1],
                    [my_cv_have_R_pkg_$1=no
                     R -e "library($1)" &>/dev/null && my_cv_have_R_pkg_$1=yes])
     if test "x$my_cv_have_R_pkg_$1" = xyes; then
        AC_DEFINE([HAVE_R_PKG_$1], [yes], [defined if R package $1 is installed])
     else
        AC_MSG_ERROR([R package $1 is required])
     fi])

# Checks for programs.
AC_PROG_CC
AC_REQUIRE_PROG([R])
AC_REQUIRE_R_PKG([igraph])
AC_REQUIRE_R_PKG([ape])
AC_REQUIRE_R_PKG([apTreeshape])
AC_REQUIRE_R_PKG([ggplot2])
AC_REQUIRE_R_PKG([kernlab])

# Output files
#AC_CONFIG_FILES([Makefile])
#AC_OUTPUT
