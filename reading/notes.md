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

**[@fearnhead2012constructing]**

* Search for summary statistics which maximise the accuracy of ABC parameter 
estimates, not necessarily the accuracy of the approximation to the posterior.
* Introduce the term "calibrated" to indicate the property that (equation 7)
the integral of the true posterior over a subset of parameter space agrees with
the integral of the ABC posterior over the same subset.
* Calibration of the ABC posterior implies that we can correctly infer 
eg. confidence intervals.
* It's actually a fairly straightforward modification to the standard importance
sampling (not sure how it would fit in adaptive ABC-SMC though).
* Summary statistics are linear models fitted to observed data, based on a 
preliminary grid search. Not too complicated.
* There's a note under theorem 5 that "ABC gives estimates that are at least as
accurate or more accurate than any other estimators based on the same summary
statistics."

# Kernel-ABC

**[@poon2015phylodynamic]**

* Kernel-ABC was able to distinguish differing contact rates in a compartmental
model much more effectively than Sackin's index (Figure 1).

# Phylodynamics

**[@grenfell2004unifying]**

* This is the reference to use for the "definition" of phylodynamics.
* Sound byte: "epidemiological and population genetic processes occur on a
similar time scale".
* Identifies four broad classes of pathogen: short infections with strong
cross-immunity, short infections with partial cross-immunity, infections with
immune enhancement, and persistent infections. These each generate
qualitatively different epidemic trajectories and phylogenies.
* Only persistent infections have both intra- and inter-host phylogenies to
consider.

**[@ypma2013relating]**

* Talks about how transmission trees are not phylogenies.
* Lots of good references, including some starting points for the parallel
problem in classical pop gen (ie. gene trees vs. species trees).
* Discrepancy between gene tree and phylogeny is worst when sampling is
highest.
* Develop a method for joint estimation of the transmission tree and phylogeny.
* Reference to "transmission tree reconstruction methods" of Wallinga and
Teunis 2004. Seems worth checking out, except the cited paper has nothing to do
with trees. Wrong reference maybe?
* Evaluate the effect of assuming at coalescent events coincide with
transmission events on the accuracy of reconstructing the transmission tree (!).
Turns out it's not that bad.
* There are some other references in the discussion about methods to
simultaneously estimate both trees.

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

**[@robinson2012dynamics]**

* Develops a method for modelling dynamic sexual contact networks.
* There is a link to a survey of sexual partnering information.

**[@welch2011statistical]**

* Says in the intro that "network-based studies to date have largely focused on
the impact of network structure on disease dynamics". Probably should mine the
references.
* The last paragraph of the section "Defining contact networks" gives some
justification for using static contact networks instead of dynamic, as a first
pass.
* Discusses the difference between network models which explicitly define the
probability of a network (eg. random graphs) vs. those which define it
implicitly via construction rules (eg. preferential attachment).
* Subsection "Inference about contact networks from disease data" is a goldmine
of references, several past papers which have tried to infer network parameters
from epidemic data.

**[@britton2002bayesian]**

* Develops a Bayesian method to infer beta and gamma in an SIR model, given
that the epidemic is spreading over a random (GNP, they call it Bernoulli)
contact network. The network itself is a nuisance parameter.
* They also optionally want to infer the infection times of each individual.
It's assumed that the removal times are known.
* Aha, they also infer the edge probability p. This is a good one for previous
work.
* They weren't able to get much information about p in their application. 

## Other Phylodynamics

**[@hughes2009molecular]**

* Phylogenetic clustering of 11,000 UK HIV sequences by genetic distance.
* Were able to fit a power-law distribution to cluster sizes by varying the
distance cutoff.
* Found few transmissions in acute infection, concluding that acute
infectiousness doesn't play a role in the heterosexual UK epidemic (Figure 4).
* Estimated the shape parameter of a Pareto degree distribution as 2.1.

**[@volz2012simple]**

* Modelling shows that viruses from acute/early infections cluster at higher
rates even when there is no elevated transmission risk during acute infection.
* This might just be because of the use of a distance cutoff for clustering.
Acute infections = shorter branches = more clustering.

**[@volz2013viral]**

* General overview of viral phylodynamics.
* Some good tidbits about the effects of various processes on phylogeny shape.
Strong directional selection in flu causes a ladder-like shape. An epidemic
which underwend rapid expansion, like HIV, tends to have a star-like phylogeny.
* Of course a bunch of unnecessarily detailed coalescent stuff is included.
* The coalescence rate should be correlated with the rate of infection of new
hosts. But is this "I" or "beta", or some combination of those? Do we have to
vary the transmission rate to get any signal?
* Coalescence rates can be calculated at equilibrium, but what about before the
model gets to equilibrium? In our case, when the epidemic continues until
saturation, there is no equilibrium except the end state.
* The ladder-like structure of flu phylogenies isn't because they all go
extinct except for one, it's because the old strains are outcompeted by the new
strain.

**[@pybus2009evolutionary]**

* Another review of viral phylodynamics. Many many references.
* Box 2 gives some alternative explanations for the within-host evolutionary
rate being higher than the between-host rate.

**[@mooers1997inferring]**

* This is a way older source for the idea of phylodynamics than Grenfell et al.,
although it doesn't use that word.
* Can use this is a reference for "nodes represent speciation events".
* Tree balance has been used to accept or reject (the latter in almost all
cases) the Markov (better known as Yule) model in macro-organism phylogenies.
* Heard 1996 talks about diversification rates depending on an evolving trait.
This is like the MMPP, but simpler (a Cox process).
* Phylogenetic estimation (admittedly they only discuss parsimony and UPGMA)
tends to be biased towards producing more unbalanced trees.

**[@kirkpatrick1993searching]**

* Another early phylodynamics paper, without the word "phylodynamics".
* Derives expectations and confidence intervals under the Yule model for
various tree balance statistics, including Sackin's (here called N-bar), and
Colless'.
* The B1 statistic is the most powerful in the particular scenario they
evaluated, but they caution against generalizing that and advise that all six
statistics be used.

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

# Other Graph Theory

**[@erdos1960evolution]**

* Didn't read this one; it describes random graphs.

**[@barabasi1999emergence]**

* Didn't read. Describes preferential attachment/scale-free networks.

**[@watts1998collective]**

* Didn't read. Describes small world networks.

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

**[@gillespie1977exact]**

* I didn't actually read this one, but it's the reference for performing
simulations with concurrent Poisson processes.

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

# HIV Biology

**[@mcmichael2001cellular]**

* This is a well-known reference for the fact that RNA viruses (in this case,
HIV) adapt to selection pressure exerted by the host immune response.

**[@domingo2012viral]**

* Overview of the concept of treating a within-host viral population as a
quasispecies.
* The quasispecies idea was apparently developed by Manfred Eigen, although not
in a virological context. Refs 240-242 possibly have some "original" uses of
viral quasispecies, good for references. Also ref 206.
* Also some good references for the fact that quasispecies theory is really the
same as species theory, with selection etc.
* On pg. 167 some good references for inter- and intra-host bottlenecking.
* We should make a distinction between the "mutation rate", which refers to
changes in individual templates, and the "rate of accumulation of mutations",
which refers to changes in the consensus sequence.
* There is some evidence that a viral quasispecies itself, rather than the
strains of which it is composed, may be a target of selection, although this is
debated (p. 175).
* Refs 230, 544 are the original work showing escape mutations in response to
ART.

# General Phylogenetics

**[@nee1992tempo]**

* This is as reference for the idea that branching times correspond to
speciation events in a time-scaled phylogeny.
* Might also be the first instance of a lineages-through-time plot (Figure 1)?
Yes, according to mooers1997inferring.
* Reference 20 sounds like it's talking about root-to-tip regression, might be
worth checking out as a new reference preceeding Rambaut and Pybus.

**[@drummond2003measurably]**

* Defines heterochronous vs. isochronous sampling, and measurably evolving
populations (MEPs).
* MEPs are "operationally defined", meaning that what exactly constitutes an
MEP depends on the amount of data available and the range of sampling times.
* Some useful references for "seminal" papers in here.

**[@janzen2015approximate]**

* Develop the normalized lineages-through-time statistic and show that it's
better than imbalance.

**[@coyne2004speciation]**

* I'm using this as a reference for what allopatric speciation is, and the fact
that it's common. See chapter 3, particularly page 84.

# Software

**[@csardi2006igraph]**

* Describes the igraph library. Just needed the reference, I didn't read it in
detail.

**[@baskins2004judy]**

* Judy arrays. Didn't read, just used.

**[@gough2009gnu]**

* The GSL manual, which is apparently the right way to cite GSL. Check on this,
it says on the GSL page that I should be citing one by a different guy?
