#!/usr/bin/env python3

from pprint import pprint
import sys
from makefile import *

data_params = {"mean_degree": [2, 4, 8, 16],
               "nnode": [5000],
               "nsimnode": [1000],
               "ntip": [200],
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
    ("mean_degree", "nnode", "ntip", "nsimnode", "seed"): [
        ("{file_basename}.nwk",
         "kernel-mean-degree/data/{mean_degree}_{nnode}.gml",
         """
         echo "#{yaml}" > $@
         nettree --sim-nodes {nsimnode} --tree-tips {ntip} --seed {seed} $^ >> $@
         """)
    ]
}

data_makefile, data_files = create_makefile(data_params, data_targets, "kernel-mean-degree/data")

results_params = {"branch_scaling": ["mean"],
                  "decay_factor": [0.1, 0.2, 0.4, 0.5],
                  "rbf_variance": [0.0625, 0.125, 1, 2],
                  "sst_control": [0, 1],
                  "coal_variance": ["INFINITY", 2],
                  "ncrossv": [1000]}

results_targets = {
    ("branch_scaling", "decay_factor", "rbf_variance", "sst_control", "coal_variance"): [
        ("{file_basename}.mtx.bz2",
         " ".join(f for f in data_files if f.endswith("nwk")),
         """
         echo "%%MatrixMarket matrix array real symmetric" > {file_basename}.mtx
         echo "%{yaml}" >> {file_basename}.mtx
         echo $(words $^) $(words $^) >> {file_basename}.mtx
         for T1 in kernel-mean-degree/data/*.nwk; do \
             for T2 in kernel-mean-degree/data/*.nwk; do \
                 if [[ "$$T1" < "$$T2" || "$$T1" == "$$T2" ]]; then \
                     treekernel -d -l {decay_factor} -g {rbf_variance} -b {branch_scaling} -s {sst_control} -c {coal_variance} $$T1 $$T2 >> {file_basename}.mtx; \
                 fi; \
             done; \
         done
         bzip2 {file_basename}.mtx
         """)],

    ("branch_scaling", "decay_factor", "rbf_variance", "sst_control", "coal_variance", "ncrossv"): [
        ("{file_basename}.tsv.bz2",
         "kernel-mean-degree/results/{branch_scaling}_{coal_variance}_{decay_factor}_{rbf_variance}_{sst_control}.mtx.bz2",
         """
         echo "#{yaml}" > {file_basename}.tsv
         Rscript -e "k <- read.mm('$^')" \
                 -e "y <- collect.metadata('kernel-mean-degree/data/*.nwk')[,'mean_degree']" \
                 -e "write.tsv(ksvm.cv(k, y, n.cv={ncrossv}, stats=c('rsquared')), '{file_basename}.tsv', append=TRUE)"
         bzip2 {file_basename}.tsv
         """)
    ]
}

results_makefile, results_files = create_makefile(results_params, results_targets, "kernel-mean-degree/results")

plots_targets = {
    ("branch_scaling", "decay_factor", "rbf_variance", "sst_control", "coal_variance"): [
        ("{file_basename}.pdf",
         "kernel-mean-degree/results/{branch_scaling}_{coal_variance}_{decay_factor}_{rbf_variance}_{sst_control}.mtx.bz2",
         """
         Rscript -e "k <- read.mm('$^')" \
                 -e 'kpca.plot(k, color=list(mean.degree=collect.metadata("kernel-mean-degree/data/*.nwk")[,"mean_degree"]), yaml="{yaml}")' \
                 -e "ggsave('$@')"
         """)
    ],
    ("ncrossv",): [
        ("{file_basename}.pdf",
         " ".join(x for x in results_files if x.endswith(".tsv.bz2")),
         """
	 Rscript -e "summary.plot(collect.data('kernel-mean-degree/results/*.tsv.bz2'), x='rbf_variance', y='rsquared', group='decay_factor', facet.x='sst_control', facet.y='coal_variance', fun='mean')" \
		 -e "ggsave('$@', width=10)"
         """)
    ]
}

plots_makefile, plots_files = create_makefile(results_params, plots_targets, "kernel-mean-degree/plots")

all_depends = [f.replace("pdf", "png") for f in plots_files if f.endswith("pdf")]

generic_targets = {
    (): [("all", "dirs " + " ".join(all_depends), ""),
         ("dirs", "", "mkdir -p kernel-mean-degree/data kernel-mean-degree/results kernel-mean-degree/plots"),
         ("%.png", "%.pdf", "convert -density 300 -trim -transparent white -resize @250000 +repage $^ $@")]
}

generic_makefile, _ = create_makefile({}, generic_targets, "kernel-mean-degree")
print("SHELL := /bin/bash\n")
print(".SECONDARY:\n")
sys.stdout.write(generic_makefile)
sys.stdout.write(plots_makefile)
sys.stdout.write(results_makefile)
sys.stdout.write(data_makefile)
