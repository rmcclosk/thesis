TESTS = check_util check_tree check_mmpp check_treestats 
check_PROGRAMS = check_util check_tree check_mmpp check_treestats 

check_util_SOURCES = check_util.c $(top_builddir)/src/util.h
check_util_CFLAGS = $(WARN_CFLAGS) @CHECK_CFLAGS@
check_util_LDADD = $(WARN_LDFLAGS) $(top_builddir)/src/libnetsim.la @CHECK_LIBS@

check_tree_SOURCES = check_tree.c $(top_builddir)/src/tree.h
check_tree_CFLAGS = $(WARN_CFLAGS) @CHECK_CFLAGS@
check_tree_LDADD = $(WARN_LDFLAGS) $(top_builddir)/src/libnetsim.la @CHECK_LIBS@

check_treestats_SOURCES = check_treestats.c $(top_builddir)/src/treestats.h
check_treestats_CFLAGS = $(WARN_CFLAGS) @CHECK_CFLAGS@
check_treestats_LDADD = $(WARN_LDFLAGS) $(top_builddir)/src/libnetsim.la @CHECK_LIBS@

check_mmpp_SOURCES = check_mmpp.c $(top_builddir)/src/mmpp.h $(top_builddir)/src/tree.h $(top_builddir)/src/newick_parser.h
check_mmpp_CFLAGS = $(WARN_CFLAGS) @CHECK_CFLAGS@
check_mmpp_LDADD = $(WARN_LDFLAGS) $(top_builddir)/src/libnetsim.la @CHECK_LIBS@
