#!/bin/bash

# This script checks that the bibliography file, PDF files, and notes are all
# in sync.

# First we check the completed readings.

BIBFILE=papers.bib
NOTEFILE=notes.md
PDFDIR=pdf
TODODIR=todo
TODOFILE=todo/list.md

BIBKEYS=`grep @ papers.bib | cut -d '{' -f 2 | sed s/,$//`
NOTEKEYS=`grep '\*\*\[@' notes.md | sed s/[[:punct:]]//g`
PDFKEYS=`ls pdf -1 | cut -d '.' -f 1`

ERROR=0

for KEY in $BIBKEYS; do
    if [[ ! -f $PDFDIR/$KEY.pdf ]]; then
        echo "PDF file for $KEY is missing"
        ERROR=1
    fi

    grep "\*\*\[@$KEY" $NOTEFILE &> /dev/null
    if [[ $? -ne 0 ]]; then
        echo "Notes for $KEY are missing"
        ERROR=1
    fi
done

for KEY in $NOTEKEYS; do
    grep "@.*{$KEY" $BIBFILE &>/dev/null
    if [[ $? -ne 0 ]]; then
        echo "Bibliography entry for $KEY is missing"
        ERROR=1
    fi

    if [[ ! -f pdf/$KEY.pdf ]]; then
        echo "PDF file for $KEY is missing"
        ERROR=1
    fi
done

for KEY in $PDFKEYS; do
    grep "@.*{$KEY" $BIBFILE &>/dev/null
    if [[ $? -ne 0 ]]; then
        echo "Bibliography entry for $KEY is missing"
        ERROR=1
    fi

    grep "\*\*\[@$KEY" notes.md &> /dev/null
    if [[ $? -ne 0 ]]; then
        echo "Notes for $KEY are missing"
        ERROR=1
    fi
done

# Then we check the todo readings.

TODOPDFKEYS=`ls $TODODIR -1 | grep pdf | cut -d '.' -f 1`
TODOLISTKEYS=`grep '##' $TODOFILE | cut -d ' ' -f 2`

for KEY in $TODOPDFKEYS; do
    grep "## $KEY" $TODOFILE &> /dev/null
    if [[ $? -ne 0 ]]; then
        echo "$KEY is missing from todo list"
        ERROR=1
    fi
done

for KEY in $TODOLISTKEYS; do
    if [[ ! -f $TODODIR/$KEY.pdf ]]; then
        echo "$TODODIR/$KEY.pdf is missing"
        ERROR=1
    fi
done

exit $ERROR
