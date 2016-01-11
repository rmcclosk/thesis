all: main.pdf

outline.pdf: outline.md outline.sty
	pandoc -V geometry:margin=2cm --number-sections -H $(word 2, $^) $(word 1, $^) -o $@

main.pdf: main.tex intro.tex methods.tex papers.bib figures
	pdflatex main && biber main && pdflatex main

papers.bib:
	ln -s ../../reading/$@ .

figures: figures/speciestree.pdf \
	figures/contactnet.pdf \
	figures/pa_power_bounds.pdf \
	figures/kernel_pa_kpca.pdf \
	figures/kernel_pa_crossv.pdf \
	figures/nltt.pdf

figures/nltt.pdf: figures/nltt.R
	cd figures && ./nltt.R

figures/speciestree.pdf: figures/speciestree.tex
	cd figures && pdflatex $(notdir $^)
	pdfcrop $@ $@

figures/contactnet.pdf: figures/contactnet.tex
	cd figures && pdflatex $(notdir $^)
	pdfcrop $@ $@

figures/kernel_pa_kpca.pdf:
	$(eval FIG_PATH := $(shell echo 'SELECT DISTINCT path FROM data WHERE step = "kpca-plot";' \
		| sqlite3 ../../simulations/kernel-pa-power.sqlite))
	cp ../../simulations/$(FIG_PATH) $@

figures/kernel_pa_crossv.pdf:
	$(eval FIG_PATH := $(shell echo 'SELECT path FROM data WHERE parameter = "mean_degree" AND value = "4" INTERSECT SELECT path FROM data WHERE step = "crossv-plot";' | sqlite3 ../../simulations/kernel-pa-power.sqlite))
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
