all: contactnet.png \
	kernel-kpca.png \
	realdata-posterior-bc.png \
	realdata-hpd-bc.png \
	abc-point-estimate.png \
	abc-posterior-example.png \
	pipeline.png \
	alpha-bounds.png \
	abc-smc.png

%.png: %.pdf
	pdfcrop $^ $^
	R --silent --vanilla -e "library(extrafont); embed_fonts('$^')"
	convert -density 500 -transparent white $^ $@

contactnet.pdf: contactnet.tex
	pdflatex $^

pipeline.pdf: pipeline.tex
	pdflatex $^

alpha-bounds.pdf: alpha-bounds.R
	./$^

abc-smc.pdf: abc-smc.R
	./$^

kernel-kpca.pdf:
	cp ../../simulations/aggregates/kpca-plot/* $@

realdata-posterior-bc.pdf:
	cp ../../simulations/`echo "SELECT path FROM data WHERE step == 'posterior-plot' AND parameter == 'bctree' AND value == 0;" | sqlite3 ../../simulations/aggregates.sqlite` $@

realdata-hpd-bc.pdf:
	cp ../../simulations/`echo "SELECT path FROM data WHERE step == 'hpd-plot' AND parameter == 'bctree' AND value == 0;" | sqlite3 ../../simulations/aggregates.sqlite` $@

abc-point-estimate.pdf:
	cp ../../simulations/`echo "SELECT path FROM data WHERE step == 'point-estimate-plot' AND parameter == 'true_m' AND value == 2;" | sqlite3 ../../simulations/abc-pa-free-m.sqlite` $@

abc-posterior-example.pdf:
	cp ../../simulations/`echo "SELECT path FROM data WHERE step == 'posterior-plot' AND parameter == 'true_m' AND value == 2 INTERSECT SELECT path FROM data WHERE parameter == 'true_alpha' AND value == 1.0 INTERSECT SELECT path FROM data WHERE parameter == 'true_I' AND value == 1000;" | sqlite3 ../../simulations/abc-pa-free-m.sqlite` $@
