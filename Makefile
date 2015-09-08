.PHONY: push pull

push:
	find . -name *.pdf | xargs drive push
	drive diff -ignore-name-clashes 2>&1 | grep remote | grep -oP .*.pdf | cut -d '/' -f 3- | xargs drive delete

pull:
	drive list -r | grep '\.pdf *$' | xargs drive pull
