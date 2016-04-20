.INTERMEDIATE: paper.tex

all: paper.pdf supp-f1-9.pdf supp-f10.pdf mbe_submission.pdf

clean:
	rm -f *.aux *.pdf *.log *.acn *.bbl *.bcf *.blg *.glo *.idx *.ist *.xml *.bak *.dvi

check: paper.Rnw mbe_submission.tex abstract.tex
	aspell check --mode=tex --add-tex-command="citep op" --add-tex-command="citet op" $(word 1,$^)
	aspell check --mode=tex --add-tex-command="citep op" --add-tex-command="citet op" $(word 2,$^)
	aspell check --mode=tex --add-tex-command="citep op" --add-tex-command="citet op" $(word 3,$^)

mbe_submission.pdf: mbe_submission.tex abstract.tex
	pdflatex --enable-write18 $^

paper.pdf: paper.tex papers.bib MBE.cls figures abstract.tex
	latex $(word 1, $^)
	bibtex $(patsubst %.tex,%,$(word 1, $^))
	latex $(word 1, $^)
	latex $(word 1, $^)
	dvipdf $(patsubst %.tex,%.dvi,$(word 1, $^))

supp-f1-9.pdf: supp-f1-9.tex supp-f1-9
	pdflatex $^

paper.tex: paper.Rnw
	R -e "library(knitr); knit('$^')"

papers.bib:
	ln -s ../../reading/papers.bib .

figures: kernel-rsquared.eps \
	abc-point-estimate-m2.eps \
	abc-posterior-example.eps \
	realdata-hpd.eps

supp-f1-9: leventhal2012fig1.pdf \
	smc-test.pdf \
	kernel-alpha-crossv.pdf \
	kernel-m-crossv.pdf \
	kernel-I-crossv.pdf \
	kernel-N-crossv.pdf \
	abc-point-estimate-m3.pdf \
	abc-point-estimate-m4.pdf \
	alpha-gamma.pdf

supp-f10.pdf:
	ln -s ../figures/abc-posterior.pdf $@

%.eps:
	ln -s ../figures/$@ .

%.pdf:
	ln -s ../figures/$@ .
