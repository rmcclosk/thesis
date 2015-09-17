ANALYSIS = kernel-mean-degree

SEED = 0
NNODE = 5000
MEAN_DEGREES = 2 4 8 16
TREE_TIPS = 100
BRANCH_SCALING = mean
REPLICATES = $(shell seq -w 0 99)
LAMBDA_VALUES = 0.1 0.2 0.3
SIGMA_VALUES = 0.5 1 2 5 10 20 50 100
N_CROSSV = 1000

.PHONY: $(ANALYSIS)-dirs

ALL_TREES = $(sort \
	$(foreach deg, $(MEAN_DEGREES), \
		$(foreach i, $(REPLICATES), \
			$(ANALYSIS)/data/$(deg)/$(i).nwk)))

$(ANALYSIS): $(ANALYSIS)-dirs \
	$(foreach s,$(SIGMA_VALUES),$(foreach l,$(LAMBDA_VALUES),$(ANALYSIS)/results/sigma-$(s)_lambda-$(l).tsv))

$(ANALYSIS)-dirs:
	mkdir -p $(ANALYSIS)/{plots,results}
	mkdir -p $(ANALYSIS)/data/{1,2,4,8}

$(ANALYSIS)/plots/cv.pdf: .Rprofile $(foreach s,$(SIGMA_VALUES),$(foreach l,$(LAMBDA_VALUES),$(ANALYSIS)/results/sigma-$(s)_lambda-$(l).tsv)) 
	Rscript -e "files <- strsplit('$(wordlist 2, $(words $^), $^)', ' ')[[1]]" \
			-e "summary.plot(collect.data(files), 'sigma', 'accuracy', 'lambda', fun='mean')" \
			-e "ggsave('$@')"

$(ANALYSIS)/results/sigma-%.pdf: $(ANALYSIS)/results/sigma-%.mtx.bz2 $(ANALYSIS)/data/degree.tsv .Rprofile
	Rscript -e "k <- read.mm('$<')" \
		    -e "kpca.plot(k, color=list(net.type=read.table('$(word 2, $^)')[,1]))" \
			-e "ggsave('$@')"

$(ANALYSIS)/results/%.tsv: $(ANALYSIS)/results/%.mtx.bz2 $(ANALYSIS)/data/degree.tsv .Rprofile
	Rscript -e "k <- read.mm('$<')" \
			-e "y <- read.table('$(word 2, $^)')[,1]" \
			-e "write.tsv(ksvm.cv(k, y, n.cv=$(N_CROSSV), stats=c('rsquared')), '$@')"

# results/sigma-0.5_lambda-0.1.mtx
$(ANALYSIS)/results/%.mtx: $(ALL_TREES)
	# parse sigma and lambda from file name
	$(eval SIGMA=$(patsubst sigma-%,%,$(word 1, $(subst _, ,$*))))
	$(eval LAMBDA=$(patsubst lambda-%,%,$(word 2, $(subst _, ,$*))))
	
	# compare each pair of trees
	echo "%%MatrixMarket matrix array real symmetric" > $@
	echo $(words $^) $(words $^) >> $@
	for T1 in $(sort $^); do \
		for T2 in $(sort $^); do \
			if [[ "$$T1" < "$$T2" || "$$T1" == "$$T2" ]]; then \
				treekernel -d -l $(LAMBDA) -g $(SIGMA) -b $(BRANCH_SCALING) $$T1 $$T2 >> $@; \
			fi \
		done \
	done

$(ANALYSIS)/data/degree.tsv: $(ALL_TREES)
	echo $^ | tr ' ' $$'\n' | cut -d '/' -f 3 > $@

$(ANALYSIS)/data/%.nwk: $$(@D).gml
	nettree --tree-tips $(TREE_TIPS) --seed $(basename $(@F)) $^ $@

$(ANALYSIS)/data/%.gml: .Rprofile
	Rscript -e "g <- barabasi.game($(NNODE), m=$*/2, directed=F)" \
			-e "write.graph(SI.net(g), '$@', format='gml')"
