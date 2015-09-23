#!/usr/bin/env python3

from pprint import pprint
import sys
from makefile import *

data_params = {"mean_degree": [2, 4, 8, 16],
               "nnode": [5000],
               "ntip": [100],
               "seed": list(range(100))}

data_targets = {
    ("mean_degree", "nnode"): [
        ("{file_basename}.gml",
         "",
         """
         Rscript -e "g <- barabasi.game({nnode}, m={mean_degree}/2, directed=F)" \
                 -e 'graph_attr(g, "comment") <- "{yaml}"' \
                 -e "write.graph(SI.net(g), '{file_basename}.gml', format='gml')"
         """)
    ],
    ("mean_degree", "nnode", "ntip", "seed"): [
        ("{file_basename}.nwk",
         "test/data/{mean_degree}_{nnode}.gml",
         """
         echo "#{yaml}" > $@
         nettree --tree-tips {ntip} --seed {seed} $^ >> $@
         """)
    ]
}

data_makefile, data_files = create_makefile(data_params, data_targets, "test/data")

results_params = {"branch_scaling": ["mean"],
                  "decay_factor": [0.1, 0.2, 0.3],
                  "rbf_variance": [0.0625, 0.125, 0.25, 0.5, 1, 2, 5],
                  "sst_control": [0, 0.5, 1],
                  "coal_power": [0, 1, 2],
                  "ncrossv": [1000]}

results_targets = {
    ("branch_scaling", "decay_factor", "rbf_variance", "sst_control", "coal_power"): [
        ("{file_basename}.mtx.bz2",
         " ".join(f for f in data_files if f.endswith("nwk")),
         """
         echo "%%MatrixMarket matrix array real symmetric" > {file_basename}.mtx
         echo "%{yaml}" >> {file_basename}.mtx
         echo $(words $^) $(words $^) >> {file_basename}.mtx
         for T1 in $(sort $^); do \
             for T2 in $(sort $^); do \
                 if [[ "$$T1" < "$$T2" || "$$T1" == "$$T2" ]]; then \
                     treekernel -d -l {decay_factor} -g {rbf_variance} -b {branch_scaling} -s {sst_control} -c {coal_power} $$T1 $$T2 >> {file_basename}.mtx; \
                 fi \
             done \
         done
         bzip2 {file_basename}.mtx
         """)],

    ("branch_scaling", "decay_factor", "rbf_variance", "sst_control", "coal_power", "ncrossv"): [
        ("{file_basename}.tsv.bz2",
         "test/results/{branch_scaling}_{coal_power}_{decay_factor}_{rbf_variance}_{sst_control}.mtx.bz2",
         """
         echo "#{yaml}" > {file_basename}.tsv
         Rscript -e "k <- read.mm('$^')" \
                 -e "y <- collect.metadata('test/data/*.nwk')[,'mean_degree']" \
                 -e "write.tsv(ksvm.cv(k, y, n.cv={ncrossv}, stats=c('rsquared')), '{file_basename}.tsv', append=TRUE)"
         bzip2 {file_basename}.tsv
         """)
    ]
}

results_makefile, results_files = create_makefile(results_params, results_targets, "test/results")

all_depends = [f for f in results_files if f.endswith("tsv.bz2")]
generic_targets = {
    (): [("all", "dirs " + " ".join(all_depends), ""),
         ("dirs", "", "mkdir -p test/data test/results"),
         ("%.png", "%.pdf", "convert -density 300 -trim -transparent white -resize @250000 +repage $^ $@")]
}

generic_makefile, _ = create_makefile({}, generic_targets, "test")
print("SHELL := /bin/bash\n")
print(".SECONDARY:\n")
sys.stdout.write(generic_makefile)
sys.stdout.write(data_makefile)
sys.stdout.write(results_makefile)
