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
	pdflatex main && \
	pdflatex main && \
	R -e "extrafont::embed_fonts('main.pdf')"

papers.bib:
	ln -sf ../../reading/$@ .

figures: figures/speciestree.pdf \
	figures/virustrees.pdf \
	figures/contactnet.pdf \
	figures/pa_power_bounds.pdf \
	figures/kernel_alpha_kpca.pdf \
	figures/kernel_alpha_crossv.pdf \
	figures/kernel_alpha_tree.pdf \
	figures/kernel_I_kpca.pdf \
	figures/kernel_I_crossv.pdf \
	figures/kernel_I_tree.pdf \
	figures/kernel_m_kpca.pdf \
	figures/kernel_m_crossv.pdf \
	figures/kernel_m_tree.pdf \
	figures/kernel_N_kpca.pdf \
	figures/kernel_N_crossv.pdf \
	figures/kernel_N_tree.pdf \
	figures/gridsearch_pa_kernel.pdf \
	figures/gridsearch_pa_estimates.pdf \
	figures/nltt.pdf \
	figures/kernel_expt.pdf \
	figures/gridsearch_expt.pdf \
	figures/acc_prec.pdf \
	figures/separable.pdf \
	figures/kernel_rsquared.pdf

figures/nltt.pdf: figures/nltt.R
	cd figures && ./nltt.R

figures/speciestree.pdf: figures/speciestree.tex
figures/contactnet.pdf: figures/contactnet.tex
figures/missing.pdf: figures/missing.tex
figures/virustrees.pdf: figures/virustrees.tex

figures/speciestree.pdf figures/contactnet.pdf figures/missing.pdf figures/virustrees.pdf:
	cd figures && pdflatex $(notdir $^)
	pdfcrop $@ $@

figures/kernel_rsquared.pdf:
	cp ../../simulations/aggregates/classifier-plot/* $@ || cp figures/missing.pdf $@
	pdfcrop $@ $@

figures/kernel_%.pdf:
	cp ../../simulations/kernel-$(word 1, $(subst _, ,$*))/$(word 2, $(subst _, ,$*))-plot/* $@ || cp figures/missing.pdf $@

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
