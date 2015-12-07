all: main.pdf

main.pdf: main.tex intro.tex methods.tex papers.bib figures
	pdflatex main && biber main && pdflatex main

papers.bib:
	ln -s ../../reading/$@ .

figures: figures/nltt.pdf

figures/nltt.pdf: figures/nltt.R
	cd figures && ./nltt.R

