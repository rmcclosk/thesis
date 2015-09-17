SHELL := /bin/bash

.PHONY: push pull

push:
	find . -name *.pdf | xargs drive push
	OLD_FILES=`drive diff -ignore-name-clashes 2>&1 | grep remote | grep -oP .*.pdf`; \
	if [[ $$? -eq 0 ]]; then \
		echo $$OLD_FILES | tr ' ' $$'\n' | cut -d '/' -f 3- | xargs drive delete; \
	fi

pull:
	dirname `drive list -r | grep '\.pdf *$$' | cut -d '/' -f 3-` | xargs mkdir -p
	drive list -r | grep '\.pdf *$$' | cut -d '/' -f 3- | xargs drive pull
	OLD_FILES=`drive diff -ignore-name-clashes 2>&1 | grep local | grep -oP '.*.pdf '`; \
	if [[ $$? -eq 0 ]]; then \
		echo $$OLD_FILES | tr ' ' $$'\n' | cut -d '/' -f 3- | xargs rm; \
	fi
