CRAN_URL = http://cran.stat.sfu.ca/
fetch = curl -L -z $(2) -o $(2) $(1)

.SECONDEXPANSION:

DOXYGEN_LATEST := $(basename $(shell wget ftp://ftp.stack.nl/pub/users/dimitri/ -O - \
					| grep -oP "doxygen.*?.pdf" | sort -V | tail -n 1))

AUTOCONF_ARCHIVE_LATEST := $(basename $(shell wget http://ftpmirror.gnu.org/autoconf-archive/ -O - \
							 | grep -oP "autoconf.*?.tar" | sort -V | tail -n 1))

JUDY_LATEST := $(basename $(shell wget http://sourceforge.net/projects/judy/files/ -O - \
				 | grep -oP -m 1 "Judy-[0-9.]*?.tar"))

PYTHON_LATEST := $(basename $(shell wget https://docs.python.org/3/download.html -O - \
				 | grep -oP -m 1 "python-.*?letter.tar"))

.PHONY: dirs

all: dirs \
	 C/igraph-docs.pdf \
	 C/gsl-ref.pdf \
	 C/judy.pdf \
	 R/igraph.pdf \
	 R/apTreeshape.pdf \
	 R/kernlab.pdf \
	 R/R-exts.pdf \
	 R/parallel.pdf \
	 R/ape.pdf \
	 R/yaml.pdf \
	 R/ergm.pdf \
	 R/network.pdf \
	 R/data.table.pdf \
	 R/reshape2.pdf \
	 R/roxygen2.pdf \
	 R/ggplot2.pdf \
	 R/ks.pdf \
	 R/vignettes/ergm/ergm.pdf \
	 python/python-library.pdf \
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

bin/autoconf-archive.pdf: arx/$(AUTOCONF_ARCHIVE_LATEST).tar.xz
	tar xf $^ -C $(dir $^)
	cd $(basename $(basename $^)) && ./configure
	cd $(basename $(basename $^)) && make pdf
	mv $(basename $(basename $^))/doc/$(notdir $@) $@
	rm -rf $(basename $(basename $^))

arx/$(AUTOCONF_ARCHIVE_LATEST).tar.xz: FORCE
	$(call fetch, http://gnu.mirror.iweb.com/autoconf-archive/$(notdir $@), $@)

bin/make.pdf: FORCE
	$(call fetch, https://www.gnu.org/software/make/manual/make.pdf, $@)

C/igraph-docs.pdf: FORCE
	$(call fetch, http://igraph.org/c/doc/igraph-docs.pdf, $@)

C/gsl-ref.pdf: arx/gsl-ref.ps.gz
	gunzip -k $^
	ps2pdf arx/gsl-ref.ps $@
	rm arx/gsl-ref.ps

arx/gsl-ref.ps.gz: FORCE
	$(call fetch, http://www.gnu.org/software/gsl/manual/gsl-ref.ps.gz, $@)

C/judy.pdf: arx/$(JUDY_LATEST).tar.gz
	tar xf $^ -C $(dir $^)
	mv arx/$(shell echo $(JUDY_LATEST) | tr J j) arx/$(JUDY_LATEST)
	ls -1d $(basename $(basename $^))/doc/ext/*.htm | sort -t_ -k1,1 | xargs pandoc -o$@
	rm -rf $(basename $(basename $^))

arx/$(JUDY_LATEST).tar.gz: FORCE
	$(call fetch, http://sourceforge.net/projects/judy/files/latest/download, $@)

python/python-library.pdf: arx/$(PYTHON_LATEST).tar.bz2
	tar xf $^ -C $(dir $^)
	mv arx/docs-pdf/library.pdf $@
	rm -rf arx/docs-pdf

arx/$(PYTHON_LATEST).tar.bz2: FORCE
	$(call fetch, https://docs.python.org/3/archives/$(PYTHON_LATEST).tar.bz2, $@)

R/R-%.pdf: FORCE
	$(call fetch, $(CRAN_URL)/doc/manuals/r-release/R-$*.pdf, $@)

R/vignettes/%.pdf: FORCE
	mkdir -p $(dir $@)
	$(call fetch, $(CRAN_URL)/web/packages/$(subst /,/vignettes/,$*), $@)

R/%.pdf: FORCE
	$(call fetch, $(CRAN_URL)/web/packages/$*/$*.pdf, $@)

spec/gml-technical-report.pdf: FORCE
	$(call fetch, http://www.fim.uni-passau.de/fileadmin/files/lehrstuhl/brandenburg/projekte/gml/gml-technical-report.pdf, $@)

spec/MMformat.pdf: FORCE
	if test `curl -I -z $@ http://math.nist.gov/MatrixMarket/reports/MMformat.ps -w %{http_code} -o /dev/null` -eq 200; then \
		curl http://math.nist.gov/MatrixMarket/reports/MMformat.ps | ps2pdf - $@; \
	fi

FORCE:
