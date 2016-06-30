.INTERMEDIATE: paper.tex

#all: main.pdf supp.pdf coverletter.pdf
#all: paper.pdf

all: changes.pdf

clean:
	rm -f *.aux *.pdf *.log *.acn *.bbl *.bcf *.blg *.glo *.idx *.ist *.xml *.bak *.dvi

check: paper.Rnw abstract.tex
	aspell check --mode=tex --add-tex-command="citep op" --add-tex-command="citet op" $(word 1,$^)
	aspell check --mode=tex --add-tex-command="citep op" --add-tex-command="citet op" $(word 2,$^)
	aspell check --mode=tex --add-tex-command="citep op" --add-tex-command="citet op" $(word 3,$^)

main.pdf: paper.pdf
	$(eval last := $(shell echo `pdfgrep -n Supplementary $^ | cut -d ':' -f 1` - 1 | bc))
	pdftk $^ cat 1-$(last) output $@

supp.pdf: paper.pdf
	$(eval first := $(shell pdfgrep -n Supplementary $^ | cut -d ':' -f 1))
	pdftk $^ cat $(first)-end output $@

paper.pdf: paper.tex papers.bib abstract.tex
	pdflatex $(word 1, $^)
	biber $(patsubst %.tex,%,$(word 1, $^))
	pdflatex $(word 1, $^)
	pdflatex $(word 1, $^)

coverletter.pdf: coverletter.tex
	pdflatex $^

paper.tex: paper.Rnw
	R -e "library(knitr); knit('$^')"

papers.bib:
	ln -s ../../reading/papers.bib .

changes.pdf: changes.tex
	pdflatex $^
