#.INTERMEDIATE: results.tex

all: main.pdf

clean:
	rm -f *.acn *.acr *.alg *.aux *.bbl *.bcf *.blg *.glg *.glo *.gls *.ist \
		  *.lof *.log *.out *.pdf *.xml *.toc results.tex

check: main.tex preamble.tex preface.tex intro.tex objective.tex methods_results.tex appendix.tex discussion.tex conclusion.tex
	for FILE in $^; do \
		aspell check $$FILE --add-tex-command="autocite op" --add-tex-command="textcite op" \
	done

outline.pdf: outline.md outline.sty
	pandoc -V geometry:margin=2cm --number-sections -H $(word 2, $^) $(word 1, $^) -o $@

define compile_latex
pdflatex $1 && biber $1 && makeglossaries $1 && pdflatex $1 && pdflatex $1 && R -e "extrafont::embed_fonts('$1.pdf')"
endef

main.pdf: main.tex preamble.tex preface.tex intro.tex objective.tex methods_results.tex appendix.tex discussion.tex conclusion.tex papers.bib
	$(call compile_latex,$(patsubst %.pdf,%,$@))

methods_results.tex: methods_results.Rnw
	R -e "knitr::knit('$^')"

discussions.tex: discussion.Rnw
	R -e "knitr::knit('$^')"

appendix.tex: appendix.Rnw
	R -e "knitr::knit('$^')"

papers.bib:
	ln -sf ../../reading/$@ .
