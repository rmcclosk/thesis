.INTERMEDIATE: results.tex

all: main.pdf

clean:
	rm -f *.acn *.acr *.alg *.aux *.bbl *.bcf *.blg *.glg *.glo *.gls *.ist \
		  *.lof *.log *.out *.pdf *.xml *.toc results.tex

outline.pdf: outline.md outline.sty
	pandoc -V geometry:margin=2cm --number-sections -H $(word 2, $^) $(word 1, $^) -o $@

main.pdf: main.tex intro.tex methods.tex objective.tex results.tex papers.bib figures
	pdflatex main && \
	biber main && \
	makeglossaries main && \
	pdflatex main && \
	pdflatex main && \
	R -e "extrafont::embed_fonts('main.pdf')"

results.tex: results.Rnw
	R -e "knitr::knit('$^')"

papers.bib:
	ln -sf ../../reading/$@ .
