.SECONDARY:

LATEX_FILES = $(basename $(wildcard *.tex))
R_FILES = $(basename $(wildcard *.R))

define MAKE_LATEX
$(1).pdf: $(1).tex
	pdflatex --enable-write18 $$^
endef

define MAKE_R
$(1).pdf: $(1).R
	./$$^
endef

# http://stackoverflow.com/questions/9551416/gnu-make-how-to-join-list-and-separate-it-with-separator
space :=
space +=
join-with = $(subst $(space),$1,$(strip $2))

all: contactnet.png \
	epitrees.png \
	speciestree.png \
	objects.png \
	nltt.png \
	pa.png \
	separable.png \
	alpha-bounds.png \
	smc-test.png \
	acc-prec.png \
	abc-idea.png \
	abc-smc.png \
	kernel-expt.png \
	kernel-kpca.png \
	kernel-rsquared.png \
	kernel-alpha-crossv.png \
	kernel-alpha-tree.png \
	kernel-alpha-kpca.png \
	kernel-m-crossv.png \
	kernel-m-kpca.png \
	kernel-m-tree.png \
	kernel-I-crossv.png \
	kernel-I-kpca.png \
	kernel-I-tree.png \
	kernel-N-crossv.png \
	kernel-N-kpca.png \
	kernel-N-tree.png \
	gridsearch-expt.png \
	gridsearch-alpha-kernel.png \
	gridsearch-alpha-point-estimate.png \
	gridsearch-I-kernel.png \
	gridsearch-I-point-estimate.png \
	gridsearch-m-kernel.png \
	gridsearch-m-point-estimate.png \
	gridsearch-N-kernel.png \
	gridsearch-N-point-estimate.png \
	pipeline.png \
	realdata-posterior-bc.png \
	realdata-hpd-bc.png \
	realdata-hpd.png \
	abc-point-estimate-m2.png \
	abc-point-estimate-m3.png \
	abc-point-estimate-m4.png \
	abc-posterior-example.png \
	abc-posterior.zip \
	leventhal2012fig1.png \
	alpha-gamma.png

clean:
	/bin/rm -f *.png *.pdf *.log *.aux

%.png: %.pdf
	R --silent --vanilla -e "library(extrafont); embed_fonts('$^')"
	pdfcrop --bbox "`pdfcrop --verbose $^ /dev/null | grep HiResBoundingBox | cut -d ' ' -f 2-5 | sort -n | head -n 1`" $^ $^
	convert -density 500 -transparent white $^ $@

$(foreach f,$(LATEX_FILES),$(eval $(call MAKE_LATEX,$(f))))
$(foreach f,$(R_FILES),$(eval $(call MAKE_R,$(f))))

kernel-rsquared.pdf:
	cp ../../simulations/aggregates/classifier-plot/* $@

kernel-kpca.pdf:
	cp ../../simulations/aggregates/kpca-plot/* $@

kernel-%.pdf:
	cp ../../simulations/kernel-$(word 1, $(subst -, ,$*))/$(word 2, $(subst -, ,$*))-plot/* $@ || cp stock/missing.pdf $@

gridsearch-%.pdf:
	cp ../../simulations/gridsearch-$(word 1, $(subst -, ,$*))/$(call join-with,-,$(wordlist 2,100,$(subst -, ,$*)))-plot/* $@ || cp stock/missing.pdf $@

realdata-posterior-bc.pdf:
	cp ../../simulations/`echo "SELECT path FROM data WHERE step == 'posterior-plot' AND parameter == 'bctree' AND value == 1;" | sqlite3 ../../simulations/aggregates.sqlite` $@

realdata-hpd-bc.pdf:
	cp ../../simulations/`echo "SELECT path FROM data WHERE step == 'hpd-plot' AND parameter == 'bctree' AND value == 1;" | sqlite3 ../../simulations/aggregates.sqlite` $@

realdata-hpd.pdf:
	cp ../../simulations/`echo "SELECT path FROM data WHERE step == 'hpd-plot' AND parameter == 'bctree' AND value == 0;" | sqlite3 ../../simulations/aggregates.sqlite` $@

abc-point-estimate-m%.pdf:
	cp ../../simulations/`echo "SELECT path FROM data WHERE step == 'point-estimate-plot' AND parameter == 'true_m' AND value == $*;" | sqlite3 ../../simulations/abc-pa-free-m.sqlite` $@

abc-posterior-example.pdf:
	cp ../../simulations/`echo "SELECT path FROM data WHERE step == 'posterior-plot' AND parameter == 'true_m' AND value == 2 INTERSECT SELECT path FROM data WHERE parameter == 'true_alpha' AND value == 1.0 INTERSECT SELECT path FROM data WHERE parameter == 'true_I' AND value == 1000;" | sqlite3 ../../simulations/abc-pa-free-m.sqlite` $@

abc-posterior.zip:
	mkdir -p $(basename $@)
	cp ../../simulations/abc-pa-free-m/posterior-plot/* $(basename $@)
	zip -r $(basename $@) $(basename $@)
	rm -rf $(basename $@)

leventhal2012fig1.pdf:
	cp ../../simulations/leventhal2012-fig1/sackin-plot/* $@

alpha-gamma.pdf:
	cp ../../simulations/pa-gamma/alpha-gamma-plot/* $@
