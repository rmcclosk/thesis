fetch = wget -N $(1) -nd -P $(2)

DOXYGEN_LATEST := $(shell wget ftp://ftp.stack.nl/pub/users/dimitri/ -O - | \
	grep -oP 'doxygen.*?.pdf' | sort -V | tail -n 1)

all: bin/pandoc/README.pdf \
	 C/igraph/igraph-docs.pdf \
	 R/igraph/igraph.pdf \
	 bin/doxygen/$(DOXYGEN_LATEST)

bin/pandoc/README.pdf: FORCE
	$(call fetch, http://pandoc.org/README.pdf, $(dir $@))

bin/doxygen/$(DOXYGEN_LATEST):
	$(call fetch, ftp://ftp.stack.nl/pub/users/dimitri/$(DOXYGEN_LATEST).zip, $(dir $@))
	unzip $@.zip -d $(dir $@)
	rm $@.zip

C/igraph/igraph-docs.pdf: FORCE
	$(call fetch, http://igraph.org/c/doc/igraph-docs.pdf, $(dir $@))

R/igraph/igraph.pdf: FORCE
	$(call fetch, http://igraph.org/r/doc/igraph.pdf, $(dir $@))

FORCE:
