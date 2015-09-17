SHELL := /bin/bash
SEED = 0

.SECONDARY:
.SECONDEXPANSION:

%.png: %.pdf
	convert -density 300 -trim -resize @250000 +repage $^ $@

%.bz2: %
	bzip2 $^

include kernel-shapes.mk
include kernel-mean-degree.mk

.Rprofile:
	echo "set.seed($(SEED))" > $@
	echo "suppressPackageStartupMessages(library(netabc))" >> $@
	echo "suppressPackageStartupMessages(library(igraph))" >> $@
