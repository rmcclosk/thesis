.SECONDARY:

all: contactnet.png \
	speciestree.png \
	kernel-kpca.png \
	kernel-rsquared.png \
	kernel-alpha-crossv.png \
	kernel-m-crossv.png \
	kernel-I-crossv.png \
	kernel-N-crossv.png \
	realdata-posterior-bc.png \
	realdata-hpd-bc.png \
	realdata-hpd.png \
	abc-point-estimate.png \
	abc-posterior-example.png \
	pipeline.png \
	alpha-bounds.png \
	abc-smc.png \
	leventhal2012fig1.png \
	abc-idea.png

%.png: %.pdf
	R --silent --vanilla -e "library(extrafont); embed_fonts('$^')"
	pdfcrop $^ $^
	convert -density 500 -transparent white $^ $@

contactnet.pdf: contactnet.tex
	pdflatex $^

pipeline.pdf: pipeline.tex
	pdflatex $^

speciestree.pdf: speciestree.tex
	pdflatex $^

alpha-bounds.pdf: alpha-bounds.R
	./$^

abc-idea.pdf: abc-idea.R
	./$^

abc-smc.pdf: abc-smc.R
	./$^

kernel-rsquared.pdf:
	cp ../../simulations/aggregates/classifier-plot/* $@

kernel-kpca.pdf:
	cp ../../simulations/aggregates/kpca-plot/* $@

kernel-%.pdf:
	cp ../../simulations/kernel-$(word 1, $(subst -, ,$*))/$(word 2, $(subst -, ,$*))-plot/* $@ || cp figures/missing.pdf $@

realdata-posterior-bc.pdf:
	cp ../../simulations/`echo "SELECT path FROM data WHERE step == 'posterior-plot' AND parameter == 'bctree' AND value == 1;" | sqlite3 ../../simulations/aggregates.sqlite` $@

realdata-hpd-bc.pdf:
	cp ../../simulations/`echo "SELECT path FROM data WHERE step == 'hpd-plot' AND parameter == 'bctree' AND value == 1;" | sqlite3 ../../simulations/aggregates.sqlite` $@

realdata-hpd.pdf:
	cp ../../simulations/`echo "SELECT path FROM data WHERE step == 'hpd-plot' AND parameter == 'bctree' AND value == 0;" | sqlite3 ../../simulations/aggregates.sqlite` $@

abc-point-estimate.pdf:
	cp ../../simulations/`echo "SELECT path FROM data WHERE step == 'point-estimate-plot' AND parameter == 'true_m' AND value == 2;" | sqlite3 ../../simulations/abc-pa-free-m.sqlite` $@

abc-posterior-example.pdf:
	cp ../../simulations/`echo "SELECT path FROM data WHERE step == 'posterior-plot' AND parameter == 'true_m' AND value == 2 INTERSECT SELECT path FROM data WHERE parameter == 'true_alpha' AND value == 1.0 INTERSECT SELECT path FROM data WHERE parameter == 'true_I' AND value == 1000;" | sqlite3 ../../simulations/abc-pa-free-m.sqlite` $@

leventhal2012fig1.pdf:
	cp ../../simulations/leventhal2012-fig1/sackin-plot/* $@
