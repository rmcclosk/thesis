TESTS = check_likelihood check_util check_tree check_thread
check_PROGRAMS = check_likelihood check_util check_tree check_thread

check_thread_SOURCES = check_thread.c $(top_builddir)/src/thread.h $(top_builddir)/src/likelihood.h
check_thread_CFLAGS = $(WARN_CFLAGS) @CHECK_CFLAGS@
check_thread_LDADD = $(WARN_LDFLAGS) $(top_builddir)/src/libpcbr.la @CHECK_LIBS@ $(top_builddir)/C-Thread-Pool/libcthreadpool.la

check_tree_SOURCES = check_tree.c $(top_builddir)/src/tree.h $(top_builddir)/src/likelihood.h
check_tree_CFLAGS = $(WARN_CFLAGS) @CHECK_CFLAGS@
check_tree_LDADD = $(WARN_LDFLAGS) $(top_builddir)/src/libpcbr.la @CHECK_LIBS@

check_util_SOURCES = check_util.c $(top_builddir)/src/util.h
check_util_CFLAGS = $(WARN_CFLAGS) @CHECK_CFLAGS@
check_util_LDADD = $(WARN_LDFLAGS) $(top_builddir)/src/libpcbr.la @CHECK_LIBS@

check_likelihood_SOURCES = check_likelihood.c $(top_builddir)/src/likelihood.h
check_likelihood_CFLAGS = $(WARN_CFLAGS) @CHECK_CFLAGS@ $(R_CFLAGS)
check_likelihood_LDADD = $(WARN_LDFLAGS) $(top_builddir)/src/libpcbr.la @CHECK_LIBS@ $(R_LIBS)
