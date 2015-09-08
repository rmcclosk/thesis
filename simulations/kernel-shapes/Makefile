SHELL := /bin/bash

SEED = 0
NNODE = 5000
MEAN_DEGREE = 8
REPLICATES = $(shell seq -w 0 999)

.SECONDARY:
net_boiler = -e "g <- as.directed(g, 'mutual')" \
	  -e "V(g)\$$remove <- 0" \
	  -e "E(g)\$$transmit <- 0.01" \
	  -e "write_graph(g, '$1', 'gml')"

all: kernel.tsv

kernel.tsv: $(foreach nettype,BA ER WS,$(foreach i,$(REPLICATES),$(nettype)/$(i).kernel.tsv))
	cat $^ > $@

%.kernel.tsv: $(foreach nettype,BA ER WS,$(foreach i,$(REPLICATES),$(nettype)/$(i).nwk))
	for T in $^; do \
		if [[ "$*.nwk" < "$$T" || "$*.nwk" == "$$T" ]]; then \
			K=`treekernel -d -b mean $*.nwk $$T`; \
			if [[ $$? -eq 0 ]]; then \
				echo -e $*.nwk\\t$$T\\t$$K >> $@; \
			else \
				echo -e $*.nwk\\t$$T\\t0.0 >> $@; \
			fi \
		fi \
	done

.SECONDEXPANSION:
%.nwk: $$(subst /,,$$(dir $$@)).gml
	nettree --tree-tips 100 --seed $(notdir $*) $^ $@

BA.gml: .Rprofile
	R -e "g <- barabasi.game($(NNODE), m=$(MEAN_DEGREE)/2, directed=FALSE)" \
		$(call net_boiler,$@)

ER.gml: .Rprofile
	R -e "g <- erdos.renyi.game($(NNODE), $(MEAN_DEGREE)/$(NNODE))" \
		$(call net_boiler,$@)

WS.gml: .Rprofile
	R -e "g <- watts.strogatz.game(1, $(NNODE), $(MEAN_DEGREE)/2, 0.01)" \
		$(call net_boiler,$@)
