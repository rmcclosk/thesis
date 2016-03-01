all: poster.pdf

clean:
	rm -f *.aux *.bcf *.log *.nav *.out *.snm *.toc poster.pdf

poster.pdf: poster.tex papers.bib
	pdflatex $^
	biber poster
	pdflatex $^
	pdflatex $^

papers.bib:
	ln -s ../../reading/$@ $@
