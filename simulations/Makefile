all: kernel-edge-density \
	kernel-mean-degree \
	kernel-nbhd-size \
	kernel-pa-power \
	kernel-rewire \
	kernel-nnode \
	gridsearch-edge-density \
	gridsearch-pa-power \
	abc-pa \
	abc-pa-bctree

%-bctree:
	rsync --delete -avz rmcclosk@bulbasaur:/home/rmcclosk/thesis/simulations/$@.sqlite `pwd`
	rsync --delete -avz rmcclosk@bulbasaur:/home/rmcclosk/thesis/simulations/$@ `pwd`

pcbr-goodreau pcbr-cliques:
	rsync --delete -avz rmcclosk@hermes:/home/rmcclosk/thesis/simulations/$@.sqlite `pwd`
	rsync --delete -avz rmcclosk@hermes:/home/rmcclosk/thesis/simulations/$@ `pwd`

pcbr-%:
	rsync --delete -avz rmcclosk@bulbasaur:/home/rmcclosk/thesis/simulations/$@.sqlite `pwd`
	rsync --delete -avz rmcclosk@bulbasaur:/home/rmcclosk/thesis/simulations/$@ `pwd`

kernel-nnode:
	rsync --delete -avz rmcclosk@bulbasaur:/home/rmcclosk/thesis/simulations/$@.sqlite `pwd`
	rsync --delete -avz rmcclosk@bulbasaur:/home/rmcclosk/thesis/simulations/$@ `pwd`

kernel-pa-nsimnode:
	rsync --delete -avz rmcclosk@bulbasaur:/home/rmcclosk/thesis/simulations/$@.sqlite `pwd`
	rsync --delete -avz rmcclosk@bulbasaur:/home/rmcclosk/thesis/simulations/$@ `pwd`

kernel-pa-power:
	rsync --delete -avz rmcclosk@bulbasaur:/home/rmcclosk/thesis/simulations/$@.sqlite `pwd`
	rsync --delete -avz rmcclosk@bulbasaur:/home/rmcclosk/thesis/simulations/$@ `pwd`

kernel-edge-density:
	rsync --delete -avz rmcclosk@bulbasaur:/home/rmcclosk/thesis/simulations/$@.sqlite `pwd`
	rsync --delete -avz rmcclosk@bulbasaur:/home/rmcclosk/thesis/simulations/$@ `pwd`

kernel-nbhd-size:
	rsync --delete -avz rmcclosk@bulbasaur:/home/rmcclosk/thesis/simulations/$@.sqlite `pwd`
	rsync --delete -avz rmcclosk@bulbasaur:/home/rmcclosk/thesis/simulations/$@ `pwd`

kernel-%:
	rsync --delete -avz rmcclosk@hermes.westgrid.ca:/home/rmcclosk/thesis/simulations/$@.sqlite `pwd`
	rsync --delete -avz rmcclosk@hermes.westgrid.ca:/home/rmcclosk/thesis/simulations/$@ `pwd`

gridsearch-pa-power:
	rsync --delete -avz rmcclosk@bulbasaur:/home/rmcclosk/thesis/simulations/$@.sqlite `pwd`
	rsync --delete -avz rmcclosk@bulbasaur:/home/rmcclosk/thesis/simulations/$@ `pwd`

gridsearch-%:
	rsync --delete -avz rmcclosk@hermes.westgrid.ca:/home/rmcclosk/thesis/simulations/$@.sqlite `pwd`
	rsync --delete -avz rmcclosk@hermes.westgrid.ca:/home/rmcclosk/thesis/simulations/$@ `pwd`

abc-pa:
	rsync --delete -avz rmcclosk@breezy:/home/rmcclosk/thesis/simulations/$@.sqlite `pwd`
	rsync --delete -avz rmcclosk@breezy:/home/rmcclosk/thesis/simulations/$@ `pwd`

abc-pa-%:
	rsync --delete -avz rmcclosk@bulbasaur:/home/rmcclosk/thesis/simulations/$@.sqlite `pwd`
	rsync --delete -avz rmcclosk@bulbasaur:/home/rmcclosk/thesis/simulations/$@ `pwd`

abc-genbank:
	rsync --delete -avz rmcclosk@hermes:/home/rmcclosk/thesis/simulations/$@.sqlite `pwd`
	rsync --delete -avz rmcclosk@hermes:/home/rmcclosk/thesis/simulations/$@ `pwd`

abc-%:
	rsync --delete -avz rmcclosk@breezy.westgrid.ca:/home/rmcclosk/thesis/simulations/$@.sqlite `pwd`
	rsync --delete -avz rmcclosk@breezy.westgrid.ca:/home/rmcclosk/thesis/simulations/$@ `pwd`

clean:
	rm -f Rplots.pdf
	rm -f .expttmp/*
	rm -f *.[eo][0-9][0-9][0-9][0-9][0-9][0-9]
	rm -f *.[eo][0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9]
