SHELL := /bin/bash

NNODE = 5000
MEAN_DEGREE = 8
WS_REWIRE_PROB = 0.01
T_VALUES=$(shell seq 0.1 0.1 0.9)
REPLICATES=$(shell seq -w 0 99)
NET_TYPES=WS ER BA
REMOVE_RATE = 2e-3

define \n


endef

.SECONDARY:

TREES = $(foreach T, $(T_VALUES), $(foreach i, $(REPLICATES), $(foreach nettype, $(NET_TYPES),$(nettype)/$(T)/$(i).nwk)))

all: sackin.png

%.png: %.pdf
	convert -density 300 $^ $@

sackin.pdf: src/sackin_plot.R sackin.tsv
	src/sackin_plot.R sackin.tsv sackin.pdf

sackin.tsv: src/sackin.R $(TREES)
	src/sackin.R $(NET_TYPES) sackin.tsv

kernel.tsv: $(foreach i, $(REPLICATES), $(foreach nettype, $(NET_TYPES),$(nettype)/0.9/$(i).nwk))
	echo tree1\\ttree2\\tK > $@
	for T1 in $^; do \
		for T2 in $^; do \
			K=`treekernel -d -b mean $$T1 $$T2`; \
			if [[ $$? -eq 0 ]]; then \
				echo -e $$T1\\t$$T2\\t$$K >> $@; \
			else \
				echo -e $$T1\\t$$T2\\t0.0 >> $@; \
			fi \
		done \
	done

%.nwk: %.gml
	nettree $*.gml --seed $(notdir $*) $*.nwk

# folder structure is NET_TYPE/T/SEED.gml
%.gml: src/make_net.R
	mkdir -p $(dir $@)
	$(eval NET_TYPE := $(word 1, $(subst /, ,$@)))
	$(eval T := $(word 2, $(subst /, ,$*)))
	$(eval SEED := $(notdir $*))
	$^ $(SEED) $(NET_TYPE) $(NNODE) $(MEAN_DEGREE) $(REMOVE_RATE) $(T) $*.gml $(WS_REWIRE_PROB)
