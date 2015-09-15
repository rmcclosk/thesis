CRAN_URL = "http://cran.stat.sfu.ca/"
fetch = curl -z $(2) -o $(2) $(1)

.SECONDEXPANSION:

DOXYGEN_LATEST := $(basename $(shell wget ftp://ftp.stack.nl/pub/users/dimitri/ -O - \
					| grep -oP "doxygen.*?.pdf" | sort -V | tail -n 1))

AUTOCONF_ARCHIVE_LATEST := $(basename $(shell wget http://ftpmirror.gnu.org/autoconf-archive/ -O - \
							 | grep -oP "autoconf.*?.tar" | sort -V | tail -n 1))

.PHONY: dirs

all: dirs \
	 C/igraph-docs.pdf \
	 C/gsl-ref.pdf \
	 R/igraph.pdf \
	 R/apTreeshape.pdf \
	 R/kernlab.pdf \
	 R/R-exts.pdf \
	 bin/pandoc.pdf \
	 bin/$(DOXYGEN_LATEST).pdf \
	 bin/autoconf.pdf \
	 bin/autoconf-archive.pdf \
	 bin/make.pdf \
	 spec/gml-technical-report.pdf \
	 spec/MMformat.pdf

dirs:
	mkdir -p bin C R spec arx

bin/$(DOXYGEN_LATEST).pdf: arx/$(DOXYGEN_LATEST).pdf.zip
	unzip -f $^ -d $(dir $@)

arx/$(DOXYGEN_LATEST).pdf.zip: FORCE
	$(call fetch, ftp://ftp.stack.nl/pub/users/dimitri/$(notdir $@), $@)

bin/pandoc.pdf: FORCE
	$(call fetch, http://pandoc.org/README.pdf, $@)

bin/autoconf.pdf: FORCE
	$(call fetch, http://www.gnu.org/software/autoconf/manual/autoconf.pdf, $@)

bin/autoconf-archive.pdf: arx/$(AUTOCONF_ARCHIVE_LATEST).tar.lz
	tar xf $^ -C $(dir $^)
	cd $(basename $(basename $^)) && ./configure
	cd $(basename $(basename $^)) && make pdf
	mv $(basename $(basename $^))/doc/$(notdir $@) $@
	rm -rf $(basename $(basename $^))

bin/make.pdf: FORCE
	$(call fetch, https://www.gnu.org/software/make/manual/make.pdf, $@)

arx/$(AUTOCONF_ARCHIVE_LATEST).tar.lz: FORCE
	$(call fetch, http://gnu.mirror.iweb.com/autoconf-archive/$(notdir $@), $@)

C/igraph-docs.pdf: FORCE
	$(call fetch, http://igraph.org/c/doc/igraph-docs.pdf, $@)

C/gsl-ref.pdf: arx/gsl-ref.ps.gz
	gunzip -k $^
	ps2pdf arx/gsl-ref.ps $@
	rm arx/gsl-ref.ps

arx/gsl-ref.ps.gz: FORCE
	$(call fetch, http://www.gnu.org/software/gsl/manual/gsl-ref.ps.gz, $@)

R/R-%.pdf: FORCE
	$(call fetch, $(CRAN_URL)/doc/manuals/r-release/R-$*.pdf, $@)

R/%.pdf: FORCE
	$(call fetch, $(CRAN_URL)/web/packages/$*/$*.pdf, $@)

spec/gml-technical-report.pdf: FORCE
	$(call fetch, http://www.fim.uni-passau.de/fileadmin/files/lehrstuhl/brandenburg/projekte/gml/gml-technical-report.pdf, $@)

spec/MMformat.pdf: FORCE
	if test `curl -I -z $@ http://math.nist.gov/MatrixMarket/reports/MMformat.ps -w %{http_code} -o /dev/null` -eq 200; then \
		curl http://math.nist.gov/MatrixMarket/reports/MMformat.ps | ps2pdf - $@; \
	fi

FORCE:
