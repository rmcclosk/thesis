all: committee-meeting-1-summary.pdf presentation.pdf

clean:
	rm -f *.aux *.bbl *.bcf *.blg *.log *.xml

presentation.pdf: presentation.tex
	pdflatex $^

%.pdf: %.tex
	pdflatex $^
	biber $(basename $^).bcf
	pdflatex $^
	pdflatex $^
