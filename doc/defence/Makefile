LATEX=$(shell ls *.tex)


all: test.pdf presentation.pdf 

%.pdf: %.tex $(LATEX)
	$(eval TMP := $(shell mktemp))
	pdflatex $*
	R -e "extrafont::embed_fonts('$@')"
	pdf2ps $@ - | ps2pdf - $(TMP)
	mv $(TMP) $@

papers.bib:
	ln -s ../../reading/$@ .

clean:
	rm -f *.pdf *.aux *.log *.bcf *.nav *.out *.xml *.snm *.blg *.toc *.bbl
