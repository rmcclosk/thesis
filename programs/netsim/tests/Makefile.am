TESTS = check_util
check_PROGRAMS = check_util

check_util_SOURCES = check_util.c $(top_builddir)/src/util.h
check_util_CFLAGS = $(WARN_CFLAGS) @CHECK_CFLAGS@
check_util_LDADD = $(WARN_LDFLAGS) $(top_builddir)/src/libnetsim.la @CHECK_LIBS@
