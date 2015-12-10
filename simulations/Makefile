all: kernel-pa-power kernel-mean-degree

kernel-%:
	rsync --delete -avz rmcclosk@hermes.westgrid.ca:/home/rmcclosk/thesis/simulations/$@.sqlite `pwd`
	rsync --delete -avz rmcclosk@hermes.westgrid.ca:/home/rmcclosk/thesis/simulations/$@ `pwd`

clean:
	rm -f Rplots.pdf
	rm -f *.[eo][0-9][0-9][0-9][0-9][0-9][0-9]
	rm -f *.[eo][0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9]
