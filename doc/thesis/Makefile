all: main.pdf

main.pdf: main.tex intro.tex methods.tex papers.bib
	pdflatex main && biber main && pdflatex main

papers.bib:
	ln -s ../../reading/$@ .
