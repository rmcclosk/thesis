all: summary.pdf presentation.pdf

summary.pdf: summary.tex papers.bib
	pdflatex summary && biber summary && pdflatex summary

presentation.pdf: presentation.tex
	pdflatex presentation

papers.bib:
	ln -sf ../../reading/$@ .
