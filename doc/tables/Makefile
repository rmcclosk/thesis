all: papers.bib \
	 abc-hpd.pdf \
	 pa-gamma-glm.pdf \
	 kernel-expt.pdf \
	 gridsearch-expt.pdf \
	 abc-expt.pdf \
	 gamma-expt.pdf \
	 realdata.pdf \
	 realdata-bc.pdf \
	 realdata-gamma.pdf 

clean:
	rm -f *.pdf abc-hpd.tex pa-gamma-glm.tex realdata-gamma.tex realdata.tex

%.pdf: %.tex
	$(eval TMPFILE := $(shell mktemp --suffix=.tex))
	echo '\\documentclass[12pt]{article}' > $(TMPFILE)
	echo '\\usepackage[osf]{garamondx}' >> $(TMPFILE)
	echo '\\usepackage[garamondx,cmbraces]{newtxmath}' >> $(TMPFILE)
	echo '\\usepackage{array}' >> $(TMPFILE)
	echo '\\usepackage{pdflscape}' >> $(TMPFILE)
	echo '\\usepackage{nicefrac}' >> $(TMPFILE)
	echo '\\usepackage[backend=biber]{biblatex}' >> $(TMPFILE)
	echo '\\usepackage{multirow}' >> $(TMPFILE)
	echo '\\addbibresource{papers.bib}' >> $(TMPFILE)
	echo '\\begin{document}' >> $(TMPFILE)
	echo '\\pagestyle{empty}' >> $(TMPFILE)
	echo '\\begin{landscape}' >> $(TMPFILE)
	echo '\\input{$^}' >> $(TMPFILE)
	echo '\\end{landscape}' >> $(TMPFILE)
	echo '\\end{document}' >> $(TMPFILE)
	pdflatex -output-directory $(dir $(TMPFILE)) $(TMPFILE)
	mv `echo $(TMPFILE) | sed s/.tex/.pdf/` $@
	pdfcrop $@ $@

abc-hpd.tex:
	cp ../../simulations/abc-pa-free-m/hpd-table/* $@

pa-gamma-glm.tex:
	cp ../../simulations/pa-gamma/glm/* $@

realdata-gamma.tex:
	cp ../../simulations/aggregates/gamma-table/* $@

papers.bib:
	ln -s ../../reading/$@ .

realdata.tex: realdata-bc.tex
	grep -v unpublished $^ > $@
