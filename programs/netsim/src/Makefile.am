AM_YFLAGS = -d
BUILT_SOURCES = newick_parser.h

noinst_LTLIBRARIES = libnetsim.la
libnetsim_la_SOURCES = util.h util.c tree.h tree.c treestats.h treestats.c mmpp.h mmpp.c simulate.h simulate.c newick_parser.h newick_parser.y newick_lexer.l
libnetsim_la_CFLAGS = $(GSL_CFLAGS) $(WARN_CFLAGS)
libnetsim_la_LIBADD = $(GSL_LIBS) $(WARN_LDFLAGS)

bin_PROGRAMS = nettree treekernel

nettree_SOURCES = nettree.c 
nettree_CFLAGS = $(WARN_CFLAGS) 
nettree_LDADD = $(WARN_LDFLAGS) libnetsim.la

treekernel_sources = treekernel.c
treekernel_CFLAGS = $(WARN_CFLAGS) 
treekernel_LDADD = $(WARN_LDFLAGS) libnetsim.la
