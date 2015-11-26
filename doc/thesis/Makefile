all: main.pdf

main.pdf: main.tex intro.tex papers.bib
	pdflatex main.tex

papers.bib:
	ln -s ../../reading/$@ .
