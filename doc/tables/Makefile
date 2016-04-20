all: abc-hpd.pdf \
	 pa-gamma-glm.pdf

%.pdf: %.tex
	$(eval TMPFILE := $(shell mktemp --suffix=.tex))
	echo '\\documentclass[12pt]{article}' > $(TMPFILE)
	echo '\\usepackage[osf]{garamondx}' >> $(TMPFILE)
	echo '\\usepackage[garamondx,cmbraces]{newtxmath}' >> $(TMPFILE)
	echo '\\usepackage{array}' >> $(TMPFILE)
	echo '\\begin{document}' >> $(TMPFILE)
	echo '\\pagestyle{empty}' >> $(TMPFILE)
	echo '\\input{$^}' >> $(TMPFILE)
	echo '\\end{document}' >> $(TMPFILE)
	pdflatex -output-directory $(dir $(TMPFILE)) $(TMPFILE)
	mv `echo $(TMPFILE) | sed s/.tex/.pdf/` $@
	pdfcrop $@ $@

abc-hpd.tex:
	cp ../../simulations/abc-pa-free-m/hpd-table/* $@

pa-gamma-glm.tex:
	cp ../../simulations/pa-gamma/glm/* $@
