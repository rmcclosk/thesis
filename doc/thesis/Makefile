#.INTERMEDIATE: results.tex

all: main.pdf

clean:
	rm -f *.acn *.acr *.alg *.aux *.bbl *.bcf *.blg *.glg *.glo *.gls *.ist \
		  *.lof *.log *.out *.pdf *.xml *.toc results.tex

outline.pdf: outline.md outline.sty
	pandoc -V geometry:margin=2cm --number-sections -H $(word 2, $^) $(word 1, $^) -o $@

#pdflatex $1 && biber $1 && makeglossaries $1 && pdflatex $1 && pdflatex $1 && R -e "extrafont::embed_fonts('$1.pdf')"

define compile_latex
pdflatex $1
endef

main.pdf: main.tex preamble.tex preface.tex intro.tex objective.tex methods.tex results.tex appendix.tex papers.bib
	$(call compile_latex,$(patsubst %.pdf,%,$@))

results.tex: results.Rnw
	R -e "knitr::knit('$^')"

papers.bib:
	ln -sf ../../reading/$@ .
