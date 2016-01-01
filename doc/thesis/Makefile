all: main.pdf

outline.pdf: outline.md outline.sty
	pandoc -V geometry:margin=2cm --number-sections -H $(word 2, $^) $(word 1, $^) -o $@

main.pdf: main.tex intro.tex methods.tex papers.bib figures
	pdflatex main && biber main && pdflatex main

papers.bib:
	ln -s ../../reading/$@ .

figures: figures/speciestree.pdf figures/contactnet.pdf

figures/nltt.pdf: figures/nltt.R
	cd figures && ./nltt.R

figures/speciestree.pdf: figures/speciestree.tex
	cd figures && pdflatex $(notdir $^)
	pdfcrop $@ $@

figures/contactnet.pdf: figures/contactnet.tex
	cd figures && pdflatex $(notdir $^)
	pdfcrop $@ $@

figures/kernel-pa.pdf:
	$(eval FIG_PATH := $(shell echo 'select path from data where parameter = "decay_factor" and value = 0.3 intersect select path from data where parameter = "rbf_variance" and value = 2 intersect select path from data where parameter = "nltt" and value = "no" intersect select path from data where parameter = "ntip" and value = 1000 intersect select path from data where step = "kpca-plot";' | sqlite3 ../../simulations/kernel-pa-power.sqlite))
	cp ../../simulations/$(FIG_PATH) $@

figures/crossv-pa.pdf:
	$(eval FIG_PATH := $(shell echo 'select path from data where parameter = "ntip" and value = 1000 intersect select path from data where step = "crossv-plot";' | sqlite3 ../../simulations/kernel-pa-power.sqlite))
	cp ../../simulations/$(FIG_PATH) $@

figures/gridsearch-pa.pdf:
	$(eval FIG_PATH := $(shell ls ../../simulations/gridsearch-pa-power/kernel-plot | head -n 1))
	cp ../../simulations/gridsearch-pa-power/kernel-plot/$(FIG_PATH) $@

figures/abc-pa.pdf:
	$(eval FIG_PATH := $(shell echo 'select path from data where parameter = "ntip" and value = 500 intersect select path from data where step = "accuracy-plot" limit 1;' | sqlite3 ../../simulations/abc-pa.sqlite))
	cp ../../simulations/$(FIG_PATH) $@

figures/sackin.pdf:
	cp ../../simulations/leventhal2012-fig1/sackin-plot/* $@

figures/smc.pdf:
	cp ../../programs/netsim/build/tests/check_smc_hist.pdf $@
