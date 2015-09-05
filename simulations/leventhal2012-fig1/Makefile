all: $(foreach T, $(shell seq 0.1 0.1 1), \
		$(foreach i, $(shell seq -w 0 99), ws/T$(T)/$(i).nwk))

ws/%.gml:
	mkdir -p $(dir $@)
	$(eval T := $(patsubst T%,%,$(word 2, $(subst /, ,$@))))
	Rscript --default-packages=igraph \
	-e "set.seed($(basename $(notdir $@)))" \
	-e "g <- watts.strogatz.game(1, 5000, 8, 0.01)" \
	-e "V(g)\$$remove <- 2e-3" \
	-e "E(g)\$$transmit <- $(T) / (1-$(T)) * 2e-3" \
	-e "write.graph(g, '$@', format='gml')"

%.nwk: %.gml
	nettree $^ --seed $(basename $(notdir $@)) $@
