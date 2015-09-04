PDF_FOLDER=/home/rmcclosk/gdrive

all: check 

.PHONY: check push pull

check:
	@./check.sh

notes.pdf: notes.md papers.bib
	pandoc -f markdown $(word 1, $^) --bibliography=$(word 2, $^) -o $@

pull:
	drive pull pdf todo

push:
	drive push pdf todo

