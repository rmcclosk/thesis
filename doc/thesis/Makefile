#.INTERMEDIATE: results.tex

all: main.pdf

clean:
	rm -f *.acn *.acr *.alg *.aux *.bbl *.bcf *.blg *.glg *.glo *.gls *.ist \
		  *.lof *.log *.out *.pdf *.xml *.toc 
	rm -f results.tex discussion.tex abstract.tex

check: main.tex preamble.tex abstract.Rnw preface.tex intro.tex objective.tex methods_results.Rnw appendix.Rnw discussion.Rnw conclusion.tex
	for FILE in $^; do \
		aspell check $$FILE --add-tex-command="autocite op" --add-tex-command="textcite op"; \
	done

outline.pdf: outline.md outline.sty
	pandoc -V geometry:margin=2cm --number-sections -H $(word 2, $^) $(word 1, $^) -o $@

# && R -e "extrafont::embed_fonts('$1.pdf')"
define compile_latex
pdflatex $1 && biber $1 && makeglossaries $1 && pdflatex $1 && pdflatex $1
endef

main.pdf: main.tex preamble.tex abstract.tex preface.tex intro.tex objective.tex methods_results.tex appendix.tex discussion.tex conclusion.tex papers.bib
	$(call compile_latex,$(patsubst %.pdf,%,$@))

%.tex: %.Rnw
	R -e "knitr::knit('$^')"

papers.bib:
	ln -sf ../../reading/$@ .
