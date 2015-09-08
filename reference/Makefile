CRAN_URL = "http://cran.stat.sfu.ca/"
fetch = wget -N $(1) -nd -P $(2)

.SECONDEXPANSION:

DOXYGEN_LATEST := $(basename $(shell wget ftp://ftp.stack.nl/pub/users/dimitri/ -O - \
					| grep -oP "doxygen.*?.pdf" | sort -V | tail -n 1))

AUTOCONF_ARCHIVE_LATEST := $(basename $(shell wget http://ftpmirror.gnu.org/autoconf-archive/ -O - \
							 | grep -oP "autoconf.*?.tar" | sort -V | tail -n 1))

all: bin/pandoc/README.pdf \
	 C/igraph/igraph-docs.pdf \
	 C/gsl/gsl-ref.pdf \
	 R/igraph/igraph.pdf \
	 R/apTreeshape/apTreeshape.pdf \
	 bin/doxygen/$(DOXYGEN_LATEST).pdf \
	 bin/autoconf/autoconf.pdf \
	 bin/autoconf/autoconf-archive.pdf \
	 bin/make/make.pdf \
	 spec/gml/gml-technical-report.pdf

bin/doxygen/$(DOXYGEN_LATEST).pdf: arx/doxygen/$(DOXYGEN_LATEST).pdf.zip
	unzip $^ -d $(dir $@)

arx/doxygen/$(DOXYGEN_LATEST).pdf.zip: FORCE
	$(call fetch, ftp://ftp.stack.nl/pub/users/dimitri/$(notdir $@), $(dir $@))

bin/pandoc/README.pdf: FORCE
	$(call fetch, http://pandoc.org/README.pdf, $(dir $@))

bin/autoconf/autoconf.pdf: FORCE
	$(call fetch, http://www.gnu.org/software/autoconf/manual/autoconf.pdf, $(dir $@))

bin/autoconf/autoconf-archive.pdf: arx/autoconf/$(AUTOCONF_ARCHIVE_LATEST).tar.lz
	tar xf $^ -C $(dir $^)
	cd $(basename $(basename $^)) && ./configure
	cd $(basename $(basename $^)) && make pdf
	mv $(basename $(basename $^))/doc/$(notdir $@) $@
	rm -rf $(basename $(basename $^))

bin/make/make.pdf: FORCE
	$(call fetch, https://www.gnu.org/software/make/manual/make.pdf, $(dir $@))

arx/autoconf/$(AUTOCONF_ARCHIVE_LATEST).tar.lz: FORCE
	$(call fetch, http://ftpmirror.gnu.org/autoconf-archive/$(notdir $@), $(dir $@))

C/igraph/igraph-docs.pdf: FORCE
	$(call fetch, http://igraph.org/c/doc/igraph-docs.pdf, $(dir $@))

C/gsl/gsl-ref.pdf: arx/gsl/gsl-ref.ps.gz
	gunzip -k $^
	mkdir -p $(dir $@)
	ps2pdf arx/gsl/gsl-ref.ps $@
	rm arx/gsl/gsl-ref.ps

arx/gsl/gsl-ref.ps.gz: FORCE
	$(call fetch, http://www.gnu.org/software/gsl/manual/gsl-ref.ps.gz, $(dir $@))

R/igraph/igraph.pdf: FORCE
	$(call fetch, http://igraph.org/r/doc/igraph.pdf, $(dir $@))

R/apTreeshape/apTreeshape.pdf: FORCE
	$(call fetch, $(CRAN_URL)/web/packages/apTreeshape/apTreeshape.pdf, $(dir $@))
	
spec/gml/gml-technical-report.pdf: FORCE
	$(call fetch, http://www.fim.uni-passau.de/fileadmin/files/lehrstuhl/brandenburg/projekte/gml/gml-technical-report.pdf, $(dir $@))

FORCE:
