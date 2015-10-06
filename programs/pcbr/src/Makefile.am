noinst_LTLIBRARIES = libpcbr.la
libpcbr_la_SOURCES = likelihood.c util.c tree.c thread.c
libpcbr_la_CFLAGS = $(WARN_CFLAGS) $(GSL_CFLAGS)
libpcbr_la_LIBADD = $(WARN_LDFLAGS) $(GSL_LIBS)

bin_PROGRAMS = pcbr
pcbr_SOURCES = pcbr.c 
pcbr_CFLAGS = $(WARN_CFLAGS)
pcbr_LDADD = $(WARN_LDFLAGS) libpcbr.la $(top_builddir)/c-cmaes/libcmaes.la $(top_builddir)/C-Thread-Pool/libcthreadpool.la

if HAVE_PYMOD_BIO
if HAVE_PYMOD_PYX
dist_bin_SCRIPTS = pcbr_preproc pcbr_plot
pkgdata_DATA = ../data/colorbrewer.json
EXTRA_DIST = pcbr_preproc.in pcbr_plot.in

else
dist_bin_SCRIPTS = pcbr_preproc
EXTRA_DIST = pcbr_preproc.in
endif 

CLEANFILES = $(dist_bin_SCRIPTS)

do_subst = sed -e 's,[@]pkgdatadir[@],$(pkgdatadir),g'

pcbr_plot: pcbr_plot.in Makefile
	$(do_subst) < $(srcdir)/pcbr_plot.in > pcbr_plot
	chmod +x pcbr_plot

pcbr_preproc: pcbr_preproc.in Makefile
	cp $(srcdir)/pcbr_preproc.in pcbr_preproc
	chmod +x pcbr_preproc

endif
