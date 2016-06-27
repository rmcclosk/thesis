.SECONDARY:

LATEX_FILES = $(basename $(wildcard *.tex))
R_FILES = $(basename $(wildcard *.R))
ABC_POSTERIOR_TARGETS = $(foreach alpha,0.0 0.5 1.0 1.5,\
	$(foreach I,1000 2000,\
		$(foreach m,2 3 4,\
			$(foreach N,5000,\
				$(foreach replicate,0 1 2,\
					abc-posterior/$(alpha)_$(I)_$(m)_$(N)_$(replicate).pdf)))))
ABC_POSTERIOR_2D_TARGETS = $(foreach alpha,0.0 0.5 1.0 1.5,\
	$(foreach I,1000 2000,\
		$(foreach m,2 3 4,\
			$(foreach N,5000,\
				$(foreach replicate,0 1 2,\
					abc-posterior-2d/$(alpha)_$(I)_$(m)_$(N)_$(replicate).pdf)))))

define FINALIZE_PDF
	pdfcrop --bbox "`pdfcrop --verbose $(1) /dev/null | grep HiResBoundingBox | cut -d ' ' -f 2-5 | sort -n | head -n 1`" $(1) $(1)
endef

define MAKE_LATEX
$(1).pdf: $(1).tex
	pdflatex --enable-write18 $$^
	$$(call FINALIZE_PDF,$$@)
endef

define MAKE_R
$(1).pdf: $(1).R
	./$$^
	$$(call FINALIZE_PDF,$$@)
endef

define ABC_POSTERIOR
abc-posterior/$(1)_$(2)_$(3)_$(4)_$(5).pdf:
	cp ../../simulations/`echo "SELECT path FROM data WHERE step == 'posterior-plot' AND parameter == 'true_alpha' AND value == $(1) INTERSECT \
	      SELECT path FROM data WHERE step == 'posterior-plot' AND parameter == 'true_I' AND value == $(2) INTERSECT \
	      SELECT path FROM data WHERE step == 'posterior-plot' AND parameter == 'true_m' AND value == $(3) INTERSECT \
	      SELECT path FROM data WHERE step == 'posterior-plot' AND parameter == 'true_N' AND value == $(4) INTERSECT \
		  SELECT path FROM data WHERE step == 'posterior-plot' AND parameter == 'replicate' AND value == $(5);" | \
		sqlite3 ../../simulations/abc-pa-free-m.sqlite` $$@
	$$(call FINALIZE_PDF,$$@)
endef

define ABC_POSTERIOR_2D
abc-posterior-2d/$(1)_$(2)_$(3)_$(4)_$(5).pdf:
	cp ../../simulations/`echo "SELECT path FROM data WHERE step == 'posterior-plot-2d' AND parameter == 'true_alpha' AND value == $(1) INTERSECT \
	      SELECT path FROM data WHERE step == 'posterior-plot-2d' AND parameter == 'true_I' AND value == $(2) INTERSECT \
	      SELECT path FROM data WHERE step == 'posterior-plot-2d' AND parameter == 'true_m' AND value == $(3) INTERSECT \
	      SELECT path FROM data WHERE step == 'posterior-plot-2d' AND parameter == 'true_N' AND value == $(4) INTERSECT \
		  SELECT path FROM data WHERE step == 'posterior-plot-2d' AND parameter == 'replicate' AND value == $(5);" | \
		sqlite3 ../../simulations/abc-pa-free-m.sqlite` $$@
	$$(call FINALIZE_PDF,$$@)
endef

# http://stackoverflow.com/questions/9551416/gnu-make-how-to-join-list-and-separate-it-with-separator
space :=
space +=
join-with = $(subst $(space),$1,$(strip $2))

all: dirs \
	contactnet.pdf \
	transtree.pdf \
    contactnet-empty.pdf \
	bctree-posterior.pdf \
	cuevas2009-posterior.pdf \
	li2015-posterior.pdf \
	little2014-posterior.pdf \
	kao2011-posterior.pdf \
	grabowski2014-posterior.pdf \
	mccormack2002-posterior.pdf \
	zetterberg2004-posterior.pdf \
	niculescu2015-posterior.pdf \
	novitsky2014-posterior.pdf \
	wang2015-posterior.pdf \
	epitrees.pdf \
	speciestree.pdf \
	objects.pdf \
	nltt.pdf \
	pa.pdf \
	kernel-idea.pdf \
	pa-example.pdf \
	sir-trajectories.pdf \
	sir-trajectories-vertical.pdf \
	compartments.pdf \
	separable.pdf \
	alpha-bounds.pdf \
	smc-test.pdf \
	smc-test-bimodal.pdf \
	acc-prec.pdf \
	abc-idea.pdf \
	kernel-idea.pdf \
	abc-smc.pdf \
	kernel-expt.pdf \
	kernel-kpca.pdf \
	kernel-rsquared.pdf \
	kernel-alpha-crossv.pdf \
	kernel-alpha-tree.pdf \
	kernel-alpha-kpca.pdf \
	kernel-m-crossv.pdf \
	kernel-m-kpca.pdf \
	kernel-m-tree.pdf \
	kernel-I-crossv.pdf \
	kernel-I-kpca.pdf \
	kernel-I-tree.pdf \
	kernel-N-crossv.pdf \
	kernel-N-kpca.pdf \
	kernel-N-tree.pdf \
	gridsearch-expt.pdf \
	gridsearch-example.pdf \
	gridsearch-alpha-kernel.pdf \
	gridsearch-alpha-point-estimate.pdf \
	gridsearch-I-kernel.pdf \
	gridsearch-I-point-estimate.pdf \
	gridsearch-m-kernel.pdf \
	gridsearch-m-point-estimate.pdf \
	gridsearch-N-kernel.pdf \
	gridsearch-N-point-estimate.pdf \
	pipeline.pdf \
	realdata-hpd-bc.pdf \
	realdata-hpd-bc-m2.pdf \
	realdata-hpd.pdf \
	realdata-hpd-m2.pdf \
	realdata-genes.pdf \
	abc-point-estimate-m2.pdf \
	abc-point-estimate-m3.pdf \
	abc-point-estimate-m4.pdf \
	abc-posterior-example.pdf \
	abc-posterior.pdf \
	leventhal2012fig1.pdf \
	alpha-gamma.pdf \
	mixed-posterior.pdf \
	peerdriven-posterior.pdf \
	derivatives.pdf \
	powerlaw.pdf \
	abc-marginal.pdf \
	homophily.pdf \
	vaccinate.pdf \
	nettypes.pdf \
	$(ABC_POSTERIOR_TARGETS) \
	$(ABC_POSTERIOR_2D_TARGETS)

clean:
	/bin/rm -f *.png *.pdf *.log *.aux *.eps *.tiff abc-posterior/*.pdf

dirs:
	mkdir -p abc-posterior abc-posterior-2d

# http://stackoverflow.com/questions/20128656/how-to-convert-a-pdf-to-grayscale-from-command-line-avoiding-to-be-rasterized
%.eps: %.pdf
	gs -sDEVICE=pdfwrite \
   		-sProcessColorModel=DeviceGray \
   		-sColorConversionStrategy=Gray \
   		-dOverrideICC \
   		-o temp.pdf \
   		-f $^
	pdftops temp.pdf $@

%.tiff: %.pdf
	convert -density 1200 $(word 1,$^) $@

$(foreach f,$(LATEX_FILES),$(eval $(call MAKE_LATEX,$(f))))
$(foreach f,$(R_FILES),$(eval $(call MAKE_R,$(f))))
$(foreach alpha,0.0 0.5 1.0 1.5,\
	$(foreach I,1000 2000,\
		$(foreach m,2 3 4,\
			$(foreach N,5000,\
				$(foreach replicate,0 1 2,\
					$(eval $(call ABC_POSTERIOR,$(alpha),$(I),$(m),$(N),$(replicate))))))))
$(foreach alpha,0.0 0.5 1.0 1.5,\
	$(foreach I,1000 2000,\
		$(foreach m,2 3 4,\
			$(foreach N,5000,\
				$(foreach replicate,0 1 2,\
					$(eval $(call ABC_POSTERIOR_2D,$(alpha),$(I),$(m),$(N),$(replicate))))))))

kernel-rsquared.pdf:
	cp ../../simulations/aggregates/classifier-plot/* $@
	$(call FINALIZE_PDF,$@)

kernel-kpca.pdf:
	cp ../../simulations/aggregates/kpca-plot/* $@
	$(call FINALIZE_PDF,$@)

kernel-%.pdf:
	cp `ls ../../simulations/kernel-$(word 1, $(subst -, ,$*))/$(word 2, $(subst -, ,$*))-plot/* | head -n 1` $@ || cp stock/missing.pdf $@
	$(call FINALIZE_PDF,$@)

gridsearch-example.pdf:
	cp ../../simulations/aggregates/gridsearch-plot/* $@ || cp stock/missing.pdf $@
	$(call FINALIZE_PDF,$@)

gridsearch-%.pdf:
	cp ../../simulations/gridsearch-$(word 1, $(subst -, ,$*))/$(call join-with,-,$(wordlist 2,100,$(subst -, ,$*)))-plot/* $@ || cp stock/missing.pdf $@
	$(call FINALIZE_PDF,$@)

realdata-hpd-bc.pdf:
	cp ../../simulations/`echo "SELECT path FROM data WHERE step == 'hpd-plot' AND parameter == 'bctree' AND value == 1 INTERSECT SELECT path FROM data WHERE step == 'hpd-plot' AND parameter == 'm_min' AND value == 1;" | sqlite3 ../../simulations/aggregates.sqlite` $@
	$(call FINALIZE_PDF,$@)

realdata-hpd-bc-m2.pdf:
	cp ../../simulations/`echo "SELECT path FROM data WHERE step == 'hpd-plot' AND parameter == 'bctree' AND value == 1 INTERSECT SELECT path FROM data WHERE step == 'hpd-plot' AND parameter == 'm_min' AND value == 2;" | sqlite3 ../../simulations/aggregates.sqlite` $@
	$(call FINALIZE_PDF,$@)

realdata-hpd.pdf:
	cp ../../simulations/`echo "SELECT path FROM data WHERE step == 'hpd-plot' AND parameter == 'bctree' AND value == 0 INTERSECT SELECT path FROM data WHERE step == 'hpd-plot' AND parameter == 'm_min' AND value == 1;" | sqlite3 ../../simulations/aggregates.sqlite` $@
	$(call FINALIZE_PDF,$@)

realdata-genes.pdf:
	cp ../../simulations/aggregates/gene-plot/* $@
	$(call FINALIZE_PDF,$@)

realdata-hpd-m2.pdf:
	cp ../../simulations/`echo "SELECT path FROM data WHERE step == 'hpd-plot' AND parameter == 'bctree' AND value == 0 INTERSECT SELECT path FROM data WHERE step == 'hpd-plot' AND parameter == 'm_min' AND value == 2;" | sqlite3 ../../simulations/aggregates.sqlite` $@
	$(call FINALIZE_PDF,$@)

abc-point-estimate-m%.pdf:
	cp ../../simulations/`echo "SELECT path FROM data WHERE step == 'point-estimate-plot' AND parameter == 'true_m' AND value == $*;" | sqlite3 ../../simulations/abc-pa-free-m.sqlite` $@
	$(call FINALIZE_PDF,$@)

abc-posterior-example.pdf:
	cp ../../simulations/`echo "SELECT path FROM data WHERE step == 'posterior-plot' AND parameter == 'true_m' AND value == 2 INTERSECT SELECT path FROM data WHERE parameter == 'true_alpha' AND value == 1.0 INTERSECT SELECT path FROM data WHERE parameter == 'true_I' AND value == 1000 LIMIT 1;" | sqlite3 ../../simulations/abc-pa-free-m.sqlite` $@
	$(call FINALIZE_PDF,$@)

abc-posterior.pdf:
	pdfnup --nup 1x1 ../../simulations/abc-pa-free-m/posterior-plot/* -o $@
	$(call FINALIZE_PDF,$@)

leventhal2012fig1.pdf:
	cp ../../simulations/leventhal2012-fig1/sackin-plot/* $@
	$(call FINALIZE_PDF,$@)

alpha-gamma.pdf:
	cp ../../simulations/pa-gamma/alpha-gamma-plot/* $@
	$(call FINALIZE_PDF,$@)

sir-trajectories.pdf:
	cp ../../simulations/`echo "SELECT path FROM data WHERE step == 'trajectory-plot' AND parameter == 'model' AND value == 'sir' INTERSECT SELECT path FROM data WHERE step == 'trajectory-plot' AND parameter == 'orientation' AND value == 'horizontal';" | sqlite3 ../../simulations/sir.sqlite` $@
	$(call FINALIZE_PDF,$@)

sir-trajectories-vertical.pdf:
	cp ../../simulations/`echo "SELECT path FROM data WHERE step == 'trajectory-plot' AND parameter == 'model' AND value == 'sir' INTERSECT SELECT path FROM data WHERE step == 'trajectory-plot' AND parameter == 'orientation' AND value == 'vertical';" | sqlite3 ../../simulations/sir.sqlite` $@
	$(call FINALIZE_PDF,$@)

mixed-posterior.pdf:
	cp ../../simulations/abc-pa-mixed-alpha/posterior-plot/* $@
	$(call FINALIZE_PDF,$@)

peerdriven-posterior.pdf:
	cp ../../simulations/`echo "SELECT path FROM data WHERE step == 'posterior-plot' AND parameter == 'sample_peer' AND value == 1;" | sqlite3 ../../simulations/abc-pa-peerdriven.sqlite` $@
	$(call FINALIZE_PDF,$@)

grabowski2014-posterior.pdf:
	cp ../../simulations/`echo "SELECT path FROM data WHERE step == 'posterior-plot' AND parameter == 'm_min' AND value == 1 INTERSECT SELECT path FROM data WHERE step == 'posterior-plot' AND parameter == 'nparticle' AND value == 10000 INTERSECT SELECT path FROM data WHERE step == 'posterior-plot' AND parameter == 'gene' AND value == 'gag';" | sqlite3 ../../simulations/grabowski2014.sqlite` $@

zetterberg2004-posterior.pdf:
	cp ../../simulations/`echo "SELECT path FROM data WHERE step == 'posterior-plot' AND parameter == 'm_min' AND value == 1 INTERSECT SELECT path FROM data WHERE step == 'posterior-plot' AND parameter == 'nparticle' AND value == 10000 INTERSECT SELECT path FROM data WHERE step == 'posterior-plot' AND parameter == 'gene' AND value == 'gag';" | sqlite3 ../../simulations/zetterberg2004.sqlite` $@

mccormack2002-posterior.pdf:
	cp ../../simulations/`echo "SELECT path FROM data WHERE step == 'posterior-plot' AND parameter == 'm_min' AND value == 1 INTERSECT SELECT path FROM data WHERE step == 'posterior-plot' AND parameter == 'nparticle' AND value == 10000 INTERSECT SELECT path FROM data WHERE step == 'posterior-plot' AND parameter == 'gene' AND value == 'gag';" | sqlite3 ../../simulations/mccormack2002.sqlite` $@

bctree-posterior.pdf:
	cp ../../simulations/`echo "SELECT path FROM data WHERE step == 'posterior-plot' AND parameter == 'm_min' AND value == 1 INTERSECT SELECT path FROM data WHERE step == 'posterior-plot' AND parameter == 'nparticle' AND value == 10000 INTERSECT SELECT path FROM data WHERE step == 'posterior-plot' AND parameter == 'cluster' AND value == 0;" | sqlite3 ../../simulations/bctree.sqlite` $@

%-posterior.pdf:
	cp ../../simulations/`echo "SELECT path FROM data WHERE step == 'posterior-plot' AND parameter == 'm_min' AND value == 1 INTERSECT SELECT path FROM data WHERE step == 'posterior-plot' AND parameter == 'nparticle' AND value == 10000;" | sqlite3 ../../simulations/$*.sqlite` $@
	$(call FINALIZE_PDF,$@)

derivatives.pdf:
	cp ../../simulations/I-vs-N/derivative-plot/* $@
	$(call FINALIZE_PDF,$@)

powerlaw.pdf:
	cp ../../simulations/powerlaw/fit-plot/* $@
	$(call FINALIZE_PDF,$@)

abc-marginal.pdf:
	cp ../../simulations/abc-marginal/point-estimate-plot/* $@
	$(call FINALIZE_PDF,$@)
