# Tree Kernel

**[@moschitti2006making]**

* Fast algorithm for computing the tree kernel by ordering nodes by production
rule (Table 1).

**[@poon2013mapping]**

* Introduces the tree kernel.

# Other Kernels

**[@lodhi2002text]**

* Kernel methods provide a way to avoid explicit feature extraction on complex
data, by operating on dot-products of data points in a high-dimensional feature
space.
* This approach won't work on sequences of real numbers, since it considers the
feature space to be all possible k-tuples, which is only countable when the
alphabet itself is countable. I think I'm looking for something more similar to
n-grams.

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

**[@wang2014approximate]**

* An exponential random graph (ERG) model assigns a probability to each network
as the exponential of a weighted sum of several network statistics.
* Assumes the set of possible edge weights is finite, so that the number of
possible networks is also finite.
* There seems to be no methods to exactly sample from ERG probability
distributions?
* There is too much math in here.
* Propose a M-H move on graphs, which either adds an edge, deletes an edge, or
completes a triangle.

**[@beaumont2009adaptive]**

* Develops an unbiased importance sampling version of SMC.
* This is mostly a theoretical paper demonstrating and correcting a bias in the
version developed by Sisson et al. 2007.

**[@beaumont2010approximate]**

* Pop-gen oveview of ABC.
* See Joyce & Marjoram 2008 for a method of systematically choosing among
summary statistics.

**[@del2012adaptive]**

* Develop a version of ABC-SMC with an adaptive tolerance schedule.
* This is the version I ended up implementing.
* Particle weights are proportional to the improvement at each step (Equation
11).
* Applies an MCMC kernel at each step. That is, particle perturbations are
accepted with probability equal to a Metropolis-Hastings ratio.
* Tolerances at each step are chosen so that the estimated sample size of the
particle population decreases by a small amount (between 1% and 10%) each
iteration (Equation 12).

# Kernel-ABC

**[@poon2015phylodynamic]**

* Kernel-ABC was able to distinguish differing contact rates in a compartmental
model much more effectively than Sackin's index (Figure 1).

# Phylodynamics

## Stochastic Mapping

**[@nielsen2002mapping]**

* Analysing reconstructed characters along a phylogeny is statistically
unsound, because it doesn't take the uncertainty in the character estimation
into account.
* Simulate N mappings of mutations along the tree in proportion to their
probability under some model, then perform inference on the set of mappings.

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

**[@o2010contact]**

* Consider biases in epimediological parameter estimates caused by contact
heterogeneity.
* Can use this as the citation for how epidemics were simulated over networks.
* Network heterogeneity tends to result in shorter transmission times near the
root compared to near the tips.
* Transmission bottlenecks increase the likelihood that the transmission tree
is concordant with the virus phylogeny.

**[@worby2014within]**

* Transmission tree doesn't correspond to virus tree - check refs 27-29 for
more.
* See also refs 4, 6, 7, and 31 for other attempts to infer transmission
networks.
* This just says that genetic distance is a poor proxy for evolutionary
distance - I think we knew this already.

**[@little2014using]**

* Develop a network-based score (TNS) to determine an individual's transmission
risk.
* This paper is a circlejerk. High network degree leads to high network degree
in the future, says nothing about actual transmission.
* Specifically targeting ART to high-risk people isn't any better than
targeting people at random.
* They don't say whether TNS was correlated with risk factors. I'm guessing
not.

**[@robinson2013dynamics]**

* Reference 5 has some information about the relationship between a long-term
aggregated contact network and an instantaneous contact network.
* 5 also has a method to create networks based on a real survey of sexual
contacts.
* Apparently the method of choosing infection times from an exponential
distribution, and then choosing infected nodes at random, is called a Gillespie
simulation.
* Network structure influences the relationship between an LTT plot and the
underlying prevalence plot: could affect using LTT plots in the kernel.
* Distribution of phylogenetic cluster sizes does not line up with degree
distribution in the contact network.
* Sackin's index is pretty bad at distinguishing network types, especially if
the networks are dynamic (Figure 8).
* It's still possible to distinguish trees from realistic vs. random networks
if the incidence nand prevalence curves are very closely matched.
* Imbalance, or rather the relationship between epidemic dynamics and
imbalance, is significantly affected by sampling schemes (Figure 11).

## Other Phylodynamics

**[@hughes2009molecular]**

* Phylogenetic clustering of 11,000 UK HIV sequences by genetic distance.
* Were able to fit a power-law distribution to cluster sizes by varying the
distance cutoff.
* Found few transmissions in acute infection, concluding that acute
infectiousness doesn't play a role in the heterosexual UK epidemic (Figure 4).
* Estimated the shape parameter of a Pareto degree distribution as 2.1.

# ERGMs

**[@robins2007introduction]**

* General overview of ERGMs.
* See Snijders et al: New specifications for exponential random graph models.
* ERGMs offer a way to model social networks while incorporating realistic
social dynamics.
* Models are fitted to observed networks (we don't have this).
* Network "ties" (ie. edges) depend on local structure.
* We assign each configuration of interest (eg. a reciprocal tie) a parameter,
and then equate or constrain a bunch of the parameters.
* Markov random graph model: includes eg. edge, 2-star, 3-star, and triangle
effects (Figure 1 and Equation 4).

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
