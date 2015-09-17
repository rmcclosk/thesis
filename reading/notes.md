# Tree Kernel

**[@moschitti2006making]**

* Fast algorithm for computing the tree kernel by ordering nodes by production
rule (Table 1).

**[@poon2013mapping]**

* Introduces the tree kernel.

# ABC

**[@sisson2007sequential]**

* Introduces the sequential Monte Carlo approach to ABC.
* I don't understand the backwards kernel thing, but the correction replaces it
with a sum over the forward kernel.
* If the forward and backward kernels are equal, and the prior is uniform, then
all the weights are equal and the algorithm becomes quite trivial (end of
Methods).

**[@marin2012approximate]**

* Overview of developed ABC algorithms.
* Three main types: rejection, MCMC, and SMC.
* There has been no observed advantage to simulating multiple data points from
the same parameter vector.
* Reference to Ratmann et al. 2009, who describe a method to use ABC for model
selection.
* Evidence from multiple groups that SMC seems to outperform MCMC.
* Different regression-based approaches to estimating true parameters from the
output.
* Gibbs random fields? Something to look into.
* It would be easier to directly simulate the summary statistics, rather than 
simulate the data and then calculate them.

# Kernel-ABC

**[@poon2015phylodynamic]**

* Kernel-ABC was able to distinguish differing contact rates in a compartmental
model much more effectively than Sackin's index (Figure 1).

# Phylodynamics

## Contact Networks

**[@leventhal2012inferring]**

* Compare Sackin's index of trees generated from Erdős–Rényi (ER) networks,
Barabási–Albert (BA) networks, Watts-Strogatz (WS) networks, and full graphs
(FG).
* Sackin's index has a local maximum when the transmissibility is set to the
critical value between the epidemic dying out and taking off (Figure 1). This
means Sackin's index is not useful in determining whether an epidemic will
emerge or not - like Figure 1A in [@poon2015phylodynamic].
* "it is not always possible to distinguish between the BA and WS models when
considering the Sackin index as a measure of tree topology"
* When a small number of tips (5%) are sampled from a larger tree, the network
types become indistinguishable by Sackin's index (Figure 4).
* Sampling only alive/infectious individuals, as opposed to including
dead/removed samples in the tree, increases the difference in Sackin's index
between the network types (Figure 6).
* Also examined a real HIV tree, but were only able to say that the contact
structure was non-random.

**[@colijn2014phylogenetic]**

* Examine effect of transmission patterns on tree topology only (not branch
lengths).
* Use a "duration of infection" parameter rather than a recovery rate.
* Define eleven tree summary metrics (Figure 3). Useful to supplement the
kernel?
* Reference 30 has information about degree distributions in real sexual
contact networks.
* Reduced sampling density has a strong effect on the ability to distinguish 
tree shapes from different epidemic types.


# Other Statistics

**[@duchene2015evaluating]**

* Traditional model selection can only indicate relative goodness of models,
not tell us if a particular model fits well or not. This is a problem if all
the models being tested explain the data poorly.
* Compare simulated data sets from a model to real data to assess model fit.
* Could be an effective method of model selection with ABC.
* In a Bayesian setting, underparametrization is more damaging than
overparametrization, due to integration over uncertainty in unnecessary
parameters.
* The accuracy statistic proposed by the authors doesn't give a hard "yes or
no" answer to whether a model is a good fit. It's more like a tool to assist in
model selection by eye.

# HIV Epidemiology

**[@beyrer2012global]**

* General information about prevalence and practices among MSM worldwide.
* Different infection rates for insertive vs. receptive intercourse - edge
attributes?
* Fast spreading of HIV among MSM is mostly explained by high per-act infection
rates, rather than network structure. Indicates broad applicability of method?

**[@keener2014connectivity]**

* Describes a recently developed transmission network score assigned to each
person based on their viral DNA (presumably from a phylogeny-derived network).

**[@centre2011mtrack]**

* Surveillance data about MSM in Canada.
* Tables 10, 11, 14, 18 give numbers of contacts.
