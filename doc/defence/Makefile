LATEX=$(shell ls *.tex)

all: test.pdf presentation.pdf

#pdflatex $* && \
#biber $* && \
#pdflatex $* && \
#pdflatex $* && \
#R -e "extrafont::embed_fonts('$@')"
%.pdf: %.tex $(LATEX)
	pdflatex $*

papers.bib:
	ln -s ../../reading/$@ .

clean:
	rm -f *.pdf *.aux *.log *.bcf *.nav *.out *.xml *.snm *.blg *.toc *.bbl
