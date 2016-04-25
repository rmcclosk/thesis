#!/bin/bash

echo '\documentclass[12pt]{article}' > /tmp/tmp.OxKysAvb3p.tex
echo '\usepackage[osf]{garamondx}' >> /tmp/tmp.OxKysAvb3p.tex
echo '\usepackage[garamondx,cmbraces]{newtxmath}' >> /tmp/tmp.OxKysAvb3p.tex
echo '\usepackage{array}' >> /tmp/tmp.OxKysAvb3p.tex
echo '\begin{document}' >> /tmp/tmp.OxKysAvb3p.tex
echo '\pagestyle{empty}' >> /tmp/tmp.OxKysAvb3p.tex
echo '\input{abc-hpd.tex}' >> /tmp/tmp.OxKysAvb3p.tex
echo '\end{document}' >> /tmp/tmp.OxKysAvb3p.tex
pdflatex -output-directory /tmp/ /tmp/tmp.OxKysAvb3p.tex
mv `echo /tmp/tmp.OxKysAvb3p.tex | sed s/.tex/.pdf/` abc-hpd.pdf
pdfcrop abc-hpd.pdf abc-hpd.pdf
