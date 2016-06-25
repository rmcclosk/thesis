all: presentation.pdf

#R -e "extrafont::embed_fonts('$@')"
%.pdf: %.tex papers.bib
	pdflatex $* && \
	biber $* && \
	pdflatex $* && \
	pdflatex $* 

papers.bib:
	ln -s ../../reading/$@ .

clean:
	rm -f *.pdf *.aux *.log *.bcf *.nav *.out *.xml *.snm *.blg *.toc
