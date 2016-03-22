.INTERMEDIATE: paper.tex

all: paper.pdf

paper.pdf: paper.tex papers.bib figures 
	pdflatex $(word 1, $^)
	biber $(patsubst %.tex,%,$(word 1, $^))
	pdflatex $(word 1, $^)
	pdflatex $(word 1, $^)

paper.tex: paper.Rnw
	R -e "library(knitr); knit('$^')"

papers.bib:
	ln -s ../../reading/papers.bib .

figures: figures/leventhal2012fig1.pdf \
	figures/kernel-alpha-crossv.pdf \
	figures/smc_test.pdf \
	figures/abc-point-est.pdf \
	figures/abc-example.pdf

figures/leventhal2012fig1.pdf:
	cp ../../simulations/leventhal2012-fig1/sackin-plot/* $@

figures/kernel-alpha-crossv.pdf:
	cp ../../simulations/kernel-alpha/crossv-plot/* $@

figures/smc_test.pdf:
	cp ../../programs/netsim/build/tests/check_smc_hist.pdf $@

figures/abc-point-est.pdf:
	cp ../../simulations/`echo "SELECT path FROM data WHERE step = 'point-estimate-plot' AND parameter = 'true_m' AND value = 2;" | sqlite3 ../../simulations/abc-pa-free-m.sqlite` $@

figures/abc-example.pdf:
	cp ../../simulations/`echo "SELECT path FROM data WHERE step = 'posterior-plot' AND parameter = 'true_m' AND value = 2 UNION SELECT path FROM data WHERE step = 'posterior-plot' AND parameter = 'true_alpha' AND value = 1.0;" | sqlite3 ../../simulations/abc-pa-free-m.sqlite | head -n 1` $@
