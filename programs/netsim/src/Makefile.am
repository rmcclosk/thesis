AM_YFLAGS = -d
BUILT_SOURCES = newick_parser.h

noinst_LTLIBRARIES = libnetsim.la
libnetsim_la_SOURCES = util.h util.c stats.c stats.h smc.c smc.h tree.h tree.c treestats.h treestats.c mmpp.h mmpp.c simulate.h simulate.c newick_parser.h newick_parser.y newick_lexer.l $(top_builddir)/c-cmaes/cmaes.h $(top_builddir)/c-cmaes/boundary_transformation.h
libnetsim_la_CFLAGS = $(PTHREAD_CFLAGS) $(GSL_CFLAGS) $(WARN_CFLAGS)
libnetsim_la_LIBADD = $(PTHREAD_LIBS) $(GSL_LIBS) $(WARN_LDFLAGS) $(top_builddir)/c-cmaes/libcmaes.la

bin_PROGRAMS = nettree treekernel pcbr netabc

netabc_SOURCES = netabc.c
netabc_CFLAGS = $(WARN_CFLAGS)
netabc_LDADD = $(WARN_LDFLAGS) libnetsim.la

pcbr_SOURCES = pcbr.c
pcbr_CFLAGS = $(WARN_CFLAGS)
pcbr_LDADD = $(WARN_LDFLAGS) libnetsim.la

nettree_SOURCES = nettree.c 
nettree_CFLAGS = $(WARN_CFLAGS) 
nettree_LDADD = $(WARN_LDFLAGS) libnetsim.la

treekernel_sources = treekernel.c
treekernel_CFLAGS = $(WARN_CFLAGS) 
treekernel_LDADD = $(WARN_LDFLAGS) libnetsim.la
