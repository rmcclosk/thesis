all: main.pdf

clean:
	rm -f *.acn *.acr *.alg *.aux *.bbl *.bcf *.blg *.glg *.glo *.gls *.ist \
		  *.lof *.log *.out *.pdf *.xml *.toc 

outline.pdf: outline.md outline.sty
	pandoc -V geometry:margin=2cm --number-sections -H $(word 2, $^) $(word 1, $^) -o $@

main.pdf: main.tex intro.tex methods.tex objective.tex results.tex papers.bib figures
	pdflatex main && \
	biber main && \
	pdflatex main && \
	makeglossaries main && \
	pdflatex main

papers.bib:
	ln -sf ../../reading/$@ .

figures: figures/speciestree.pdf \
	figures/virustrees.pdf \
	figures/contactnet.pdf \
	figures/pa_power_bounds.pdf \
	figures/kernel_pa_kpca.pdf \
	figures/kernel_pa_crossv.pdf \
	figures/kernel_pa_trees.pdf \
	figures/kernel_nsimnode_kpca.pdf \
	figures/kernel_nsimnode_crossv.pdf \
	figures/kernel_nsimnode_trees.pdf \
	figures/gridsearch_pa_kernel.pdf \
	figures/gridsearch_pa_estimates.pdf \
	figures/nltt.pdf \
	figures/kernel_expt.pdf \
	figures/gridsearch_expt.pdf \
	figures/acc_prec.pdf \
	figures/separable.pdf

figures/nltt.pdf: figures/nltt.R
	cd figures && ./nltt.R

figures/speciestree.pdf: figures/speciestree.tex
figures/contactnet.pdf: figures/contactnet.tex
figures/missing.pdf: figures/missing.tex
figures/virustrees.pdf: figures/virustrees.tex

figures/speciestree.pdf figures/contactnet.pdf figures/missing.pdf figures/virustrees.pdf:
	cd figures && pdflatex $(notdir $^)
	pdfcrop $@ $@

figures/kernel_pa_kpca.pdf:
	cp ../../simulations/kernel-pa-power/kpca-plot/* $@ || cp figures/missing.pdf $@

figures/kernel_pa_crossv.pdf:
	cp ../../simulations/kernel-pa-power/crossv-plot/* $@ || cp figures/missing.pdf $@

figures/kernel_pa_trees.pdf:
	cp ../../simulations/kernel-pa-power/tree-plot/* $@ || cp figures/missing.pdf $@
	pdfcrop $@ $@

figures/kernel_nsimnode_kpca.pdf:
	cp ../../simulations/kernel-pa-nsimnode/kpca-plot/* $@ || cp figures/missing.pdf $@

figures/kernel_nsimnode_crossv.pdf:
	cp ../../simulations/kernel-pa-nsimnode/crossv-plot/* $@ || cp figures/missing.pdf $@

figures/kernel_nsimnode_trees.pdf:
	cp ../../simulations/kernel-pa-nsimnode/tree-plot/* $@ || cp figures/missing.pdf $@
	pdfcrop $@ $@

figures/gridsearch_pa_kernel.pdf:
	cp ../../simulations/gridsearch-pa-power/kernel-plot/* $@ || cp figures/missing.pdf $@

figures/gridsearch_pa_estimates.pdf:
	cp ../../simulations/gridsearch-pa-power/error-plot/* $@ || cp figures/missing.pdf $@

figures/kernel_expt.pdf: figures/kernel_expt.R figures/kernel_expt.tex
	cd figures && ./$(notdir $(word 1, $^))
	cd figures && pdflatex $(notdir $(word 2, $^))
	cd figures && pdfcrop $(notdir $@) $(notdir $@)

figures/gridsearch_expt.pdf: figures/gridsearch_expt.R figures/gridsearch_expt.tex
	cd figures && ./$(notdir $(word 1, $^))
	cd figures && pdflatex $(notdir $(word 2, $^))
	cd figures && pdfcrop $(notdir $@) $(notdir $@)

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

figures/%.pdf: figures/%.R
	cd figures && ./$(notdir $^)
