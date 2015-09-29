#!/usr/bin/env python3

from pprint import pprint
import sys
from makefile import *

data_params = {"mean_degree": [8],
               "nnode": [5000],
               "nsimnode": [1000],
               "net_type": ["BA"],
               "seed": list(range(100)),
               "num_clusters": [1, 2, 5, 10],
               "cluster_size": [5, 10, 20, 50, 100],
               "cluster_rate": [1.5, 2, 4],
               "connect": ["TRUE", "FALSE"],
               "tree_tips": [200]}

data_targets = {
    ("mean_degree", "nnode", "net_type", "seed", "num_clusters", "cluster_size", "cluster_rate", "connect"): [
        ("{file_basename}.gml",
         "",
         """
         Rscript -e "set.seed({seed})" \
                 -e "g <- barabasi.game({nnode}, m={mean_degree}/2, directed=F)" \
                 -e 'graph_attr(g, "comment") <- "{yaml}"' \
                 -e 'g <- SI.net(add.transmission.clusters(g, {cluster_size}, {num_clusters}, {cluster_rate}, {connect}))' \
                 -e "write.graph(g, '{file_basename}_BA.gml', format='gml')"
         """)
    ],
    ("mean_degree", "nnode", "net_type", "seed", "num_clusters", "cluster_size", "cluster_rate", "connect", "tree_tips", "nsimnode"): [
        ("{file_basename}.nwk",
         "pcbr-validation/data/{cluster_rate}_{cluster_size}_{connect}_{mean_degree}_{net_type}_{nnode}_{num_clusters}_{seed}.gml",
         """
         echo "#{yaml}" > $@
         nettree --sim-nodes {nsimnode} --tree-tips {tree_tips} --seed {seed} $^ >> $@
         """)
    ]
}

data_makefile, data_files = create_makefile(data_params, data_targets, "pcbr-validation/data")

results_params = {"nrates": [2],
                  "use_tips": ["yes", "trans", "no"],
                  "trans_per_branch": ["one", "any"],
                  "tip_pdf": ["yes", "no"],
                  "cheating": ["yes", "no"]}

results_targets = {
    ("nrates", "use_tips", "trans_per_branch", "tip_pdf", "cheating"): [
        ("{file_basename}/%.tsv.bz2",
         "{pcbr-validation/%.nwk",
         """
         mkdir -p {file_basename}
         pcbr --rates {nrates} --use-tips {use_tips} --trans-per-branch {trans_per_branch} \
              --tip-pdf {tip_pdf} --cheating {cheating} --newick-outfile /dev/null \
              --annot-outfile {file_basename}/$*.tsv
         bzip $*.tsv
         """)
    ]
}

#results_targets = {
#    ("branch_scaling", "decay_factor", "rbf_variance", "sst_control", "coal_variance"): [
#        ("{file_basename}.mtx.bz2",
#         " ".join(f for f in data_files if f.endswith("nwk")),
#         """
#         echo "%%MatrixMarket matrix array real symmetric" > {file_basename}.mtx
#         echo "%{yaml}" >> {file_basename}.mtx
#         echo $(words $^) $(words $^) >> {file_basename}.mtx
#         for T1 in $(sort $^); do \
#             for T2 in $(sort $^); do \
#                 if [[ "$$T1" < "$$T2" || "$$T1" == "$$T2" ]]; then \
#                     treekernel -d -l {decay_factor} -g {rbf_variance} -b {branch_scaling} -s {sst_control} -c {coal_variance} $$T1 $$T2 >> {file_basename}.mtx; \
#                 fi \
#             done \
#         done
#         bzip2 {file_basename}.mtx
#         """)],
#
#    ("branch_scaling", "decay_factor", "rbf_variance", "sst_control", "coal_variance", "ncrossv"): [
#        ("{file_basename}.tsv.bz2",
#         "pcbr-validation/results/{branch_scaling}_{coal_variance}_{decay_factor}_{rbf_variance}_{sst_control}.mtx.bz2",
#         """
#         echo "#{yaml}" > {file_basename}.tsv
#         Rscript -e "k <- read.mm('$^')" \
#                 -e "y <- collect.metadata('pcbr-validation/data/*.nwk')[,'mean_degree']" \
#                 -e "write.tsv(ksvm.cv(k, y, n.cv={ncrossv}, stats=c('rsquared')), '{file_basename}.tsv', append=TRUE)"
#         bzip2 {file_basename}.tsv
#         """)
#    ]
#}
#
#results_makefile, results_files = create_makefile(results_params, results_targets, "pcbr-validation/results")
#
#plots_targets = {
#    ("branch_scaling", "decay_factor", "rbf_variance", "sst_control", "coal_variance"): [
#        ("{file_basename}.pdf",
#         "pcbr-validation/results/{branch_scaling}_{coal_variance}_{decay_factor}_{rbf_variance}_{sst_control}.mtx.bz2",
#         """
#         Rscript -e "k <- read.mm('$^')" \
#                 -e 'kpca.plot(k, color=list(mean.degree=collect.metadata("pcbr-validation/data/*.nwk")[,"mean_degree"]), yaml="{yaml}")' \
#                 -e "ggsave('$@')"
#         """)
#    ],
#    ("ncrossv",): [
#        ("{file_basename}.pdf",
#         " ".join(x for x in results_files if x.endswith(".tsv.bz2")),
#         """
#	 Rscript -e "summary.plot(collect.data('pcbr-validation/results/*.tsv.bz2'), x='rbf_variance', y='rsquared', group='decay_factor', facet.x='sst_control', facet.y='coal_variance', fun='mean')" \
#		 -e "ggsave('$@', width=10)"
#         """)
#    ]
#}
#
#plots_makefile, plots_files = create_makefile(results_params, plots_targets, "pcbr-validation/plots")

all_depends = [f for f in data_files if f.endswith("nwk")]

generic_targets = {
    (): [("all", "dirs " + " ".join(all_depends), ""),
         ("dirs", "", "mkdir -p pcbr-validation/data pcbr-validation/results pcbr-validation/plots"),
         ("%.png", "%.pdf", "convert -density 300 -trim -transparent white -resize @250000 +repage $^ $@")]
}

generic_makefile, _ = create_makefile({}, generic_targets, "pcbr-validation")
print("SHELL := /bin/bash\n")
print(".SECONDARY:\n")
sys.stdout.write(generic_makefile)
#sys.stdout.write(plots_makefile)
#sys.stdout.write(results_makefile)
sys.stdout.write(data_makefile)
