all: paper.pdf

paper.pdf: paper.tex papers.bib figures 
	pdflatex $(word 1, $^)
	biber $(patsubst %.tex,%,$(word 1, $^))
	pdflatex $(word 1, $^)
	pdflatex $(word 1, $^)

paper.tex: paper.Rnw
	R -e "knitr::knit('$^')"

papers.bib:
	ln -s ../../reading/papers.bib .

figures: figures/leventhal2012fig1.pdf \
	figures/kernel-alpha-crossv.pdf \
	figures/smc_test.pdf

figures/leventhal2012fig1.pdf:
	cp ../../simulations/leventhal2012-fig1/sackin-plot/* $@ || cp figures/missing.pdf $@

figures/kernel-alpha-crossv.pdf:
	cp ../../simulations/kernel-alpha/crossv-plot/* $@ || cp figures/missing.pdf $@

figures/smc_test.pdf:
	cp ../../programs/netsim/build/tests/check_smc_hist.pdf $@ || cp figures/missing.pdf $@
