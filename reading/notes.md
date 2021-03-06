# Tree Kernel

**[@collins2002new]**
* Original definition of the subset tree kernel for NLP.

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

**[@scholkopf1998nonlinear]**

* The paper which introduces kernel PCA. Only skimmed.

**[@burges1998tutorial]**

* A commonly used reference for the "kernel trick" as applied to SVM.

# ABC

**[@rubin1984bayesianly]**

* Introduces the concept of approximate Bayesian computation, although not by
that name. Didn't read.

**[@tavare1997inferring]**

* First "real" use of ABC.

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

**[@sunnaker2013approximate]**

* Another ABC overview.
* TODO

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

**[@del2006sequential]**

* This is an important paper in the development of SMC. It goes over some of
the problems with MCMC in the introduction.
* In importance sampling, the variance in estimates of the expected value (or
any statistic) is proportional to the variance in importance weights. For this
reason, it is important to select importance distributions which are "close" to
the true distribution of interest.
* Ah this is beautiful.
* In its original formulation, SMC is meant for sequentially sampling from
distributions defined on spaces of increasing dimension. To sample from
distributions on the *same* dimension, the authors introduce $n$-dimensional
distributions which admit the $n$th distribution of interest as a marginal.

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

**[@doucet2001introduction]**

* One of the primary motivations of SMC was that MCMC is not designed for
recursive problems, like sequential Bayes.
* Develops SMC by first describing IS, then SIS, and then finally the
"bootstrap filter" which combines SIS with resampling according to the
importance weights.

**[@liu2001theoretical]**

* Gives some history of the development of SIS.
* Nice definition of Monte Carlo: "any computation of expectations... can be
replaced to an acceptable degree of accuracy by using the empirical
distribution resulting from the discrete sample."
* A more intuitive explanation of SIS: instead of sampling from a series of
distributions on increasing dimensions, we only want to sample from a single
distribution on n dimensions, but we break this up into n terms using the
"telescope" law of probability.

**[@liu2008monte]**

* TODO

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

**[@lintusaari2016identifiability]**

* Previous ABC experiments with simulated data (they list several) have only
allowed one parameter to vary.
* Precedent for stopping simulations once a pre-determined number of
individuals has become infected.
* Also precedent for doing grid search.
* Transformation of variables can change uninformative priors into informative
ones: a uniform prior on the SI model parameters does not translate to a
uniform prior on the basic reproductive number, because of how these are
related.
* Different priors give wildly different posteriors.

# Kernel-ABC

**[@nakagome2013kernel]**

* First use of the kernel-ABC concept.

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

**[@ypma2012unravelling]**

* Combines genetic, geographic, and epidemiological data to estimate
transmission trees.
* Assume that the tree data types are independent and take their product to get
a total likelihood.
* Apply to a flu outbreak on farms.

**[@cottam2008integrating]**

* "Integrates" genetic data into transmission tree construction by eliminating
all trees which are inconsistent with known epidemiological data.
* Uses something like UPGMA to reconstruct the trees.

**[@leitner1996accurate]**

* Use phylogenetic methods to reconstruct a known HIV transmission history by
eliminating transmission trees not consistent with genetic data.
* Most of the estimated phylogenies were quite accurate, even though the
methods were crude.
* Choice of gene has more of an impact than choice of method or model; more
data was better.

**[@bernard2007hiv]**

* Not possible to phylogenetically ascertain who transmitted to whom (didn't
read).

**[@ou1992molecular]**

* First use of genetic data to establish transmission links.

**[@jombart2011reconstructing]**

* Phylogenetic methods will fail to reconstruct the transmission tree if both
ancestors and descendants are present in a sample, which the authors claim is
likely during the early stages of an epidemic.
* Develop a method to find the optimal phylogeny consistent with known
infection times, where "optimal" is defined as optimizing a weight function
such as genetic distance.
* Worked much better than phylogenetics when there were sampled ancestors.
Phylogenetics did very poorly in this scenario.
* Branch lengths are not estimated, only topology.
* They suggest to use their method in densely sampled early epidemics. 

**[@didelot2014bayesian]**

* Another method to combine genomic and epidemiological data to get a
transmission tree, this time a Bayesian one.
* Considers within-host diversity, in contrast to some of the other developed
methods.
* Points out that coalescences in the viral phylogeny coincide with
transmissions when within-host diversity is zero.
* They also infer the SIR parameters and effective population size.

**[@stadler2011estimating]**

* A good example of the cutting-edge use of phylodynamics to infer
epidemiological parameters. 
* Develops the birth-death model for serially sampled data.
* The likeilhood can be explicitly calculated. This is very pretty math.

**[@minin2008smooth]**

* Introduces the Bayesian skyride (just used for background, didn't read).

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
if the incidence and prevalence curves are very closely matched.
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

**[@keeling2005networks]**

* Overview of the interface between contact networks and epidemiology.
* There are some published networks obtained through contact tracing, where
infected people listed their contacts for the purposes of identifying new
infections. However these don't ask anything from the uninfected people.
* Some other groups seems to have independently discovered the BA model in
different contexts, should read the references.
* This is a really useful review with a lot of good references, although it's
old.

**[@pastor2001epidemic]**

* Describes how an epidemic spreads on a scale-free network, in the context of
a computer virus spreading over the internet.
* References for epidemic threshold, should mention in intro.
* No epidemic threshold for scale-free networks with alpha = 1.

**[@britton2002bayesian]**

* Develops a Bayesian method to infer beta and gamma in an SIR model, given
that the epidemic is spreading over a random (GNP, they call it Bernoulli)
contact network. The network itself is a nuisance parameter.
* They also optionally want to infer the infection times of each individual.
It's assumed that the removal times are known.
* Aha, they also infer the edge probability p. This is a good one for previous
work.
* They weren't able to get much information about p in their application. 

**[@groendyke2011bayesian]**

* Extension of Britton 2002 towards an "inferential methodology of practical
use."
* Assume an ER network and explicitly calculate the likelihood of an observed
epidemic. This is possible only because it's ER I think.
* The MCMC algorithm is very complicated. But they have an R package! It's
called epinet.
* They marginalize over the transmission tree. This is troublesome, it seems to
be a superset of what we are doing.
* Extensible to ERGMs.

**[@goodreau2006assessing]**

* Uses ERGMs to model HIV contact networks.
* Check Grassly et al. for some arguments against what we're doing (they
argue in favour of the panmixis / effective population size assumptions).
* This has actual values for the parameters! Very useful for validating
clustering.
* Also provides a significantly more sophisticated method of simulating
epidemics over a network: instead of transmission trees, they actually simulate
viral evolution including within-host. This actually results in a simulated
phylogeny, not a transmission tree. Could be very useful!
* Contact network structure can affect the estimated effective population size
by multiple orders of magnitude, so it's important to take it into account.
* I should definitely do this one for journal club.

**[@welch2011network]**

* Simulates epidemics over networks with a constant degree sequence but varying
clustering coefficients.
* Finds that the trees cannot be distinguished from one another. I wonder
though, if this would be the same with the tree kernel.
* Lots of good references in here for how contact network structure influences
epidemic spread.
* Suggests that inference of contact network properties should be focused on
"degree distribution or broader notions of population structure", which is what
I'm doing.
* Use some odd summary statistics to show there isn't much difference between
trees simulated with different levels of clustering.

**[@bansal2007individual]**

* Contrasts contact networks with homogeneous mixing models, explaining
situations where each type is useful.
* Exponential degree distribution, now power-law, was best fit to six published
contact networks.
* Should reread this one.

**[@volz2008sir]**

* Ah Erik. Of course you have done this already. Le sigh.
* Uses probability generating functions to derive explicit formulae for SIR
epidemic trajectories on random networks.
* The word "random" here is like in random variable, it doesn't mean ER
networks. They are defined in terms of their degree distributions.
* He does this by considering what happens in the limit of effective population
size.
* Final DEs are written in terms of probability generating functions of the
degree distribution. Ah it's pretty.

**[@volz2007susceptible]**

* Argh. Apparently there are tonnes of papers which fit contact network models
to epi data. Check the references.
* Extends the above to dynamic networks (but was published sooner? dunno).

**[@klovdahl1985social]**

* Canonical reference for the (modern) use of social networks in an
epidimiological context.
* This is a very politically incorrect article.

**[@barthelemy2005dynamical]**

* This is basically about how different contact network structures can
influence the prevalence curve.
* Turns out R0 is proportional to the second moment of the degree distribution.
That is E(k^2) / E(k).
* The speed of epidemic spread in a structured network can be dramatically
different than in a homogeneously mixing population.
* Seems gamma is bounded between 2 and 3?
* Most of this isn't relevant to me and was skimmed.
* Basically the epidemic makes its way to the superspreaders, and from there is
rapidly disseminated to everybody else.

**[@kemper1980identification]**

* Canonical reference for the "superspreader" concept.

**[@brown2011transmission]**

* The authors previously showed that the degree distribution of the
phylogenetic network of HIV+ MSM in the UK followed a power law.
* Here, they show that a preferential attachment model fit the data better than
a random attachment model.
* The Waring distribution, which results from a somewhat more complex
preferential attachment model, was the best fit.

**[@yirrell1998molecular]**

* Finds very little concordance between the phylogenetic and actual (based on
contact tracing) networks in a Ugandan town.

**[@resik2007limitations]**

* The phylogenetic network doesn't always agree with the contact network,
especially in cases where there is a long delay between transmission and
sampling.

**[@schneeberger2004scale]**

* Networks with a power-law degree distribution are only "scale free" if 2 <
gamma < 3.
* References 12 and 13 for the fact that gamma in "real systems" lies between 1
and 4.
* References 24-27 for arguments for and against real networks being scale-free.
* This has real data values for gamma.

**[@colgate1989risk]**

* Growth of HIV in early 1980's was cubic, which is not consistent with a
homogeneously mixed population, even taking increasing education and behaviour
change into account.
* Find that a biased mixing model can produce the kind of cubic growth seen.
* A power law distribution with exponent between 3 and 4 gives a good fit to
surveys of number of sex partners among MSM.

**[@liljeros2001web]**

* First study to indicate that sexual networks follow a power law. Has real
data values for gamma.

**[@jones2003assessment]**

* Disagrees with the other literature on PA networks for real data - it says
that PA provides a "very poor fit" to sexual networks in the USA.
* Fixing the minimum degree above 1 yields "wildly increasing confidince
intervals" apparently. Dammit.
* Estimated power-law exponents varied from 3.03 to 17.04 (!).

**[@handcock2004likelihood]**

* Persistence of STDs in the population is because of contact heterogeneity.
Without it, they would die out.
* Good historical references for some basic results (STDs only survive because
of heterogeneity; R0 increases with variance of degree distribution).
* Most people have 1 or fewer sex partners, so spread is driven primarily by
the tail of the degree distribution.
* Consider three types of network model: non-homogeneous Poisson, preferential
attachment, and vetting.
* Fit all three models to six survey datasets.
* Most networks were best fit by the non-homogeneous Poisson model over
preferential attachment. However the difference in fit was very small for all
but one dataset (American men).

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

**[@volz2012complex]**

* Develops the likelihood of a phylogeny under a general compartmental model.

**[@rasmussen2014phylodynamic]**

* The coalescent process can be used to fit arbitrarily complex demographic
models to phylogenies.
* Uses a particle filter to sample population state trajectories given
compartmental model parameters at each iteration of MCMC.
* Based on the framework of [@volz2012complex].
* Too complicated. I have no idea.

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

**[@yule1925mathematical]**

* Origin of the Yule model.

**[@kirkpatrick1993searching]**

* Another early phylodynamics paper, without the word "phylodynamics".
* Derives expectations and confidence intervals under the Yule model for
various tree balance statistics, including Sackin's (here called N-bar), and
Colless'.
* The B1 statistic is the most powerful in the particular scenario they
evaluated, but they caution against generalizing that and advise that all six
statistics be used.

**[@holmes1995revealing]**

* First reference to using the LTT to estimate epidemic prevalence.
* This is nostalgic to read. 72 sequences from all over the world. Alignment
with clustal, then tree building with phylip.
* HIV underwent exponential growth in the past but more constant in the
present, whereas HCV was endemic in the past but became epidemic recently.

**[@poon2014impact]**

* Reference for clustering by patristic distance cutoff.

**[@kenah2015algorithms]**

* Overview of relationship between transmission trees and viral phylogenies.
* This would make a good journal club paper, although it's really long.
* Great source of references for all the "standard" phylodynamics stuff.
* Data from uninfected individuals is not used by any phylodynamic methods.
* This develops a framework for the formation of transmission trees from
contact networks, including epidemiological data.

**[@stadler2013uncovering]**

* Develops the multi-type birth-death branching process.
* Detects evidence for superspreaders in MSM in Latvia.

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

**[@simon1955class]**

* This is a much older description of the preferential attachment mechanism.

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

**[@gillespie1976general]**

* I didn't actually read this one, but it's the reference for performing
simulations with concurrent Poisson processes.

**[@smola1997support]**

* The paper which introduced support vector regression. Didn't read.

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

**[@wang2015targeting]**

* Pretty solid analysis of acutely infected MSM from Beijing.
* Authors claim that they have "deeply sampled" the MSM subnetwork.
* "Similar to previously inferred HIV networks among MSM [3,23], the inferred
transmission network was best described by a preferential attachment model (ie. 
Waring distribution, rho = 3.19) [24]."
* See refs 3, 5, 33-37 for other deeply sampled MSM cohorts.

**[@novitsky2013phylogenetic]**

* Tested ~6000 residents of Mochudi, Botswana, and found ~1300 HIV+. 75\% of
those who tested positive were female.
* They hypothesize that there were "multiple introductions of HIV into the
community", which would explain the low PA power.
* Mother-to-child transmission pairs can be used as a control for clustering
analysis.
* The town is also close to Gabarone (the capital of Botswana), so it's likely
that the sexual network extends into the city.
* However they do say that the epidemic is "dominated by locally circulating
viral variants."

**[@novitsky2014impact]**

* The concept of an HIV cluster is poorly defined in the current literature.
* Very densely sampled a village in Botswana.

**[@cuevas2009hiv]**

* Analysed pol sequences from 261 newly diagnosed individuals from six
hospitals in the Basque Country, Spain.
* About half of the sequences (47%) grouped in transmission clusters.

**[@shiino2014phylodynamic]**

* Large drug resistance surveillance network in Japan.
* 3618 individuals. I assumed the second part of the Genbank identifier (eg.
"2069136" in "C12-2069136-1") was the patient ID, because there are
3563 unique IDs in the pol sequences.

**[@li2015hiv]**

* 1265 newly infected MSM from Shangai, China. There are exactly 1265 linked
sequences in Genbank.

**[@niculescu2015recent]**

* Outbreak of HIV among IDU in Romainia.

**[@morris1993epidemiology]**

* Early reference arguing for the use of social networks in epidemic modelling.
* "The diffusion of disease through a human population traces the structure of
social networks."
* Nice historical overview.
* Talks about the introduction of compartmental models and preferential mixing.
* Basically this is extending the compartmental model idea all the way down to
one differential equation per person.

# Other Contact Networks

**[@wasserman1994social]**

* Book about social network analysis, primarily in the social science context.

**[@barnes1954class]**

* Apparently the first use of the term "social network", according to
[@wasserman1994social].

**[@moreno1953shall]**

* Invention of the graphical form of social network (he called it a
"sociogram").

**[@jeong2000large]**

* Metabolic networks are scale-free (didn't read).

**[@shen2004superspreading]**

* Contribution of super-spreaders to the 2003 SARS epidemic in Beijing (didn't
read).

# HIV Biology

**[@mcmichael2001cellular]**

* This is a well-known reference for the fact that RNA viruses (in this case,
HIV) adapt to selection pressure exerted by the host immune response.

**[@domingo2012viral]**

* Overview of the concept of treating a within-host viral population as a
quasispecies.
* Refs 240-242 possibly have some "original" uses of viral quasispecies, good
for references. Also ref 206.
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

**[@drummond2003inference]**

* Review of dated-tips methods for estimating viral evolutionary rates (also a
good reference for root-to-tip regression).
* Three main methods of estimating mutation rates from serially sample data:
linear regression, ML, and Bayesian inference.
* Substitution rate is not the same as mutation rate - mutation rate is much
higher and depends only on replication errors, while substitution rate is an
expression of how many polymorphisms become fixed in a population per unit time.

**[@janzen2015approximate]**

* Develop the normalized lineages-through-time statistic and show that it's
better than imbalance.

**[@coyne2004speciation]**

* I'm using this as a reference for what allopatric speciation is, and the fact
that it's common. See chapter 3, particularly page 84.

**[@fitzpatrick2008if]**

* Speciation processes are a continuum from completely allopatric to completely
sympatric.
* Points out the implausability of completely sympatric speciation, which
requires complete panmixia. Unstructured host populations are similarly
unlikely.

**[@haeckel1866generelle]**

* The first instance of the word "phylogeny".

**[@cavalli1967phylogenetic]**

* Seems to be the first instance of the idea of the "topology" of a tree.
Possibly also of the "root" concept, although that might be the below.

**[@harding1971probabilities]**

* Might be the first instance of the term "root" for phylogenetics. At any rate
it gives a nice definition.

**[@buneman1974note]**

* Introduces the 4-point condition for ultrametric trees.

**[@nei2000molecular]**

* Good general reference book for a lot of the basic phylogenetics concepts.

**[@sokal1958statistical]**

* First instance of UPGMA, and distance-based phylogenetic methods in general.

**[@shao1990tree]**

* Defines Sackin's index.

**[@korber2000timing]**

* One of the earliest uses of a root-to-tip regression.

**[@shankarappa1999consistent]**

* The other early use of root-to-tip regression, also one of the most famous
HIV datasets.

**[@li1988rates]**

* First use of rooting with an outgroup (I think).

# Population Genetics

**[@kendall1948generalized]**

* Birth-death process.

**[@kingman1982coalescent]**

* Coalescent process.

**[@kermack1927contribution]**

* The SIR model.

# Software

**[@csardi2006igraph]**

* Describes the igraph library. Just needed the reference, I didn't read it in
detail.

**[@baskins2004judy]**

* Judy arrays. Didn't read, just used.

**[@gough2009gnu]**

* The GSL manual, which is apparently the right way to cite GSL. Check on this,
it says on the GSL page that I should be citing one by a different guy?

**[@karatzoglou2004kernlab]**

* The reference for kernlab. Didn't read.

**[@cock2009biopython]**

* Biopython.

**[@katoh2013mafft]**

* MAFFT.

**[@edgar2004muscle]**

* MUSCLE.

**[@gouy2010seaview]**

* Seaview.

**[@price2010fasttree]**

* FastTree2.

**[@capella2009trimal]**

* trimAl.

**[@drummond2007beast]**

* BEAST.

**[@plummer2006coda]**

* coda.

**[@to2015fast]**

* Least-square dating.

**[@bouckaert2014beast]**

* BEAST2.

**[@ronquist2012mrbayes]**

* MrBayes.
