define FETCH
	rsync --delete -avz rmcclosk@$(1):/home/rmcclosk/thesis/simulations/$(2).yaml `pwd`
	rsync --delete -avz rmcclosk@$(1):/home/rmcclosk/thesis/simulations/$(2).sqlite `pwd`
	rsync --delete -avz rmcclosk@$(1):/home/rmcclosk/thesis/simulations/$(2) `pwd`
endef

define RUN
	./expt.py -vv $(1).yaml
endef

.PHONY: bulbasaur breezy.westgrid.ca jasper.westgrid.ca local

all: bctree novitsky2014 cuevas2009 li2015 wang2015 niculescu2015 \
	abc-pa-free-m abc-pa-mixed-alpha abc-pa-peerdriven \
	kernel-alpha kernel-I kernel-m kernel-N \
	gridsearch-alpha gridsearch-I gridsearch-m gridsearch-n \
	aggregates

# kernel classifier experiments

kernel-%: local
	$(call RUN,$@)

# grid search experiments

gridsearch-%: local
	$(call RUN,$@)

# ABC simulation experiments

abc-pa-free-m: breezy.westgrid.ca
	$(call FETCH,$^,$@)
	$(call RUN,$@)

abc-pa-mixed-alpha: jasper.westgrid.ca
	$(call FETCH,$^,$@)
	$(call RUN,$@)

abc-pa-peerdriven: jasper.westgrid.ca
	$(call FETCH,$^,$@)
	$(call RUN,$@)

# real data

bctree: bulbasaur
	$(call FETCH,$^,$@)
	$(call RUN,$@)

novitsky2014: jasper.westgrid.ca
	$(call FETCH,$^,$@)
	$(call RUN,$@)

cuevas2009: jasper.westgrid.ca
	$(call FETCH,$^,$@)
	$(call RUN,$@)

li2015: jasper.westgrid.ca
	$(call FETCH,$^,$@)
	$(call RUN,$@)

wang2015: jasper.westgrid.ca
	$(call FETCH,$^,$@)
	$(call RUN,$@)

niculescu2015: jasper.westgrid.ca
	$(call FETCH,$^,$@)
	$(call RUN,$@)

# misc

aggregates: local
	$(call RUN,$@)

leventhal2012-fig1: local
	$(call RUN,$@)

pa-gamma: local
	$(call RUN,$@)

sir: local
	$(call RUN,$@)
