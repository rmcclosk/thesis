PDF_FOLDER=/home/rmcclosk/gdrive

all: notes.pdf check

.PHONY: check

check:
	@./check.sh

notes.pdf: notes.md papers.bib
	pandoc -f markdown $(word 1, $^) --bibliography=$(word 2, $^) -o $@
