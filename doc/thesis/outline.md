# Background

## Phylodynamics

### Species trees, transmission trees

- What a species tree is, and some of the forces which act to shape it, like
selection and speciation. 
- What transmission trees are, and how they are analogous to species trees. 
- The parallels between some of the forces shaping species trees and those
shaping transmission trees, such as between transmission and speciation.

### Gene trees, viral phylogenies

- What a gene tree is, how it's made, and how it differs from a species tree.

### "Converting" a gene tree into a species tree

- The problem of incomplete lineage sorting, and how this applies to virus
trees.
- Rooting and time-scaling, molecular clock models, and root-to-tip regression,
including a reference that root-to-tip works okay in practice although there
are criticisms.

### Tree shapes

- Overview of known causes of particular tree shapes, such as seasonal
bottlenecking causing a ladder shape, or exponential growth causing a star.
- Also shapes of particular parts of the tree, such as bursts of
speciation/transmission.
- Ways of quantifying tree shape, including the tree kernel and the nLTT.

## Contact networks

- What they are, and how they capture more information about structured
populations.
- Simple generative models for contact networks: random graphs, preferential
attachment, small world, full graphs (which represent homogeneous mixing).
- Review some of the prior work on investigating epidemic spread over contact
networks.

## Approximate Bayesian computation

### What ABC is used for

- What does it mean to "fit" a model - either find coefficients which optimize
an objective function (such as likelihood), or draw a sample of coefficients
from a distribution (such as the posterior).
- An example of a model which can be fit explicitly with calculus, such as
least squares or SI.
- An example of a model which can be fit through numerical optimization,
because the objective function can be calculated but not differentiated, such
as the SIR model.
- An example of a model which falls into none of the above categories, because
the likelihood cannot be calculated explicitly, like the chess player example.
- Explain that ABC is intended to address the last type of situation, and give
the basic idea (try parameters until something looks sort of right).

### Algorithms for ABC

- Rejection, MCMC, and SMC.
- Statistical points of importance, such as the sensitivity to appropriate
choices of summary statistics.
- Some recent advances in SMC - automatic summary statistics and tolerance
schedules.

### Applications of ABC in phylodynamics

- Overview of some studies which have used ABC, including some of the
population genetics stuff.

# Objective

**Infer parameters of contact networks from viral phylogenies.**

## Prior work

- Studies that show different epidemic parameters produce differently shaped
trees, eg. flu stuff trying to explain the ladder shape, without trying to
estimate the parameters.
- Estimation of actual values of epidemic parameters from trees, eg.
kernel-ABC.
- Studies that show different *network* parameters produce differently shaped
trees, eg. Leventhal's and Gardy's stuff.
- At least one paper tried to infer a network parameter (edge density) from
infection and removal times; this is actually one step further because they
include the transmission tree as a nuisance parameter (Britton & O'Niell).

## Our approach

- Method for inferring parameters of an arbitrary contact network model.
- Use ABC: for a set of parameters, simulate a network, then simulate a tree
over that network, then compare it to the true observed tree (illustrate this
with a schematic figure).
- Validate by generating trees from a known model, and recovering the
parameters.
- Apply to a tree of HIV-infected MSM in BC.

# Methods

## Programs

### Simulating epidemics over networks

- Gillespie simulation - everything is a Poisson process.
- Point out a few recent papers which have used this, but not provided public
implementations.
- Validate by replicating Figure 1A from Leventhal 2012.

### Adaptive ABC-SMC

- Implement the version from Del Moral et al. 2012 with an adaptive tolerance
schedule.
- Validate by replicating Figure 1 from Del Moral 2012.
- Some components: fast implementation of the tree kernel, slightly modified
nLTT.

## Simulation experiments

- Common parameters throughout:
    - number of nodes in the network = 5000
    - number of infected nodes = 1000
    - number of tips in tree = 500 (50% sampling, unless otherwise specified)
    - transmission rate = 1 (we scale the trees by their mean branch length, so
            it doesn't matter)
    - removal rate = 0
    - mean degree of nodes = 4 (unless otherwise specified)

### Identifying separable parameters and kernel meta-parameters

- First step, to tell us which parameters are possibly identifiable by the tree
kernel.
- Five network types, with one or two parameters each: random, preferential
attachment, small world, Pareto, full.
- Kernel-SVM cross-validation to identify optimal meta-parameters ($\lambda$,
$\sigma$, and nLTT).

### Marginal parameter estimation with grid search

- The previous step told us which parameters might be identified and what
kernel meta-parameters to use.
- This step tells us what kind of precision and accuracy we can expect.
- Generate many simulated trees at various parameter values, compare to a
single tree with the tree kernel.

### Multi-variate parameter estimation with ABC

- Simulate trees on preferential attachment networks with known parameters.
- Allow all remaining parameters to be estimated with ABC-SMC, except for the
number of edges per vertex in preferential attachment networks (this is fixed).

## Application to real data

- Describe the HIV BC tree, how it was made, how we pruned it to include only
MSM, and display it in a figure.
- ABC-SMC on this tree, with the preferential attachment network model, and
priors on $N$ and $I$ informed by epidemiological estimates.

# Results

## Simulation experiments

### Separable parameters and kernel meta-parameters

- Basically every parameter we look at is separable in kernel space, and has a
reasonably accurate classifier with appropriate choices of meta-parameters.
- Figure: trees with different preferential attachment power separated in
kernel space. Remaining kernel space separation plots go in supplements.
- Figure: cross-validation accuracy of a kernel-SVM classifier on preferential
attachment power, under several combinations of kernel meta-parameters.
Remaining plots of this type go in supplements.

### Marginal parameter estimation with grid search

- Some comments about which parameters, and which values of those parameters, 
it's easy or hard to estimate.
- For example, we can be quite precise when the true preferential attachment
power is equal to 1, but less precise when it is less than or greater than 1.
- Figure: plot of "distribution" of kernel scores for three different
preferential attachment powers, with marginal boxplots.

### Multi-variate parameter estimation with ABC

- Again, discussion about which parameters are and are not estimable.
- For example, for networks with two edges per vertex, we can estimate the
preferential attachment power quite well, but the number of nodes in the
network poorly.
- Figure: boxplots showing marginal posterior distributions of estimated
attachment power, number of total nodes, and number of infected nodes, for
various numbers of edges per vertex. It should hopefully show that it gets
easier to estimate number of infected nodes, but harder to estimate attachment
power, as the edge density increases.

## Application to real data

- Discuss the estimated parameters, their precision, and how they line up (or
don't) with other estimates. 
- Figure: density plots of marginal posterior distributions of preferential
attachment power, number of total nodes, and number of infected nodes.

# Conclusions

## Summary of results

- The tree kernel is pretty good at differentiating most network parameters, at
least into broad classes.
- We provide a general method to estimate contact network parameters in a
phylodynamic setting.
- We were able to estimate some parameters with ABC-SMC with fairly good
accuracy, but not others.
- What we saw in the BC tree.

## Future work

- We have only examined very simple network models here, but there exist many
more realistic ones, such as ERGMs, to which this method might be applied.
- Evaluate the robustness of the parameter estimates to amount of sampling.
- There exist methods to jointly estimate the transmission tree and viral
phylogeny (eg. Ypma 2013), perhaps these could be used instead of "converting"
the phylogeny into a transmission tree beforehand.

# Figures

## Background 

1. Schematic of a species tree, to illustrate concepts (extant taxa, branching
times, etc.).

## Methods

2. Schematic of our ABC approach to estimating network parameters.
3. Phylogeny of MSM HIV infections in BC.

## Results

4. Kernel PCA projection of trees with different preferential attachment powers.
5. Cross-validation accuracy of kSVM classifier on preferential attachment power.
6. "Distribution" of grid search kernel scores for preferential attachment
powers on simulated trees.
7. Boxplots of marginal posterior distributions of preferential attachment
network parameters, on simulated trees, for various (fixed) numbers of edges
per vertex.
8. Density plots of marginal posterior distributions for network parameters
estimated for the BC tree.

# Supplemental Figures

- Our version of Figure 1A from Leventhal 2012, to demonstrate that the
epidemic simulation works as expected.
- Our version of Figure 1 from Del Moral 2012, to demonstrate that the adaptive
ABC-SMC works as expected.
- Kernel PCA projections for all parameters except preferential attachment
power: edge density (random networks), edges per vertex (preferential
attachment), neighbourhood size and rewiring probability (small world), number
of infected nodes (full graph), Pareto distribution parameter, and network
type.
- Kernel-SVM cross-validation accuracy for all parameters except preferential
attachment power.
- "Distributions" of grid search kernel scores for all parameters except
preferential attachment power.

\newpage

# Figure 1 {-}
![](figures/speciestree.pdf)

\newpage

# Figure 2 {-}

\vspace{0.5\textheight}

# Figure 3 {-}

\newpage

# Figure 4 {-}
![](figures/kernel-pa.pdf)

\newpage

# Figure 5 {-}

![](../../simulations/kernel-edge-density/crossv-plot/aa55adf9e9775187b347d11aab237acf.pdf)

\newpage

# Figure 6 {-}
![](figures/gridsearch-pa.pdf)

\newpage

# Figure 7 {-}
![](figures/abc-pa.pdf)

\newpage

# Figure 8 {-}

\newpage

# Supplemental Figure 1 {-}
![](figures/sackin.pdf)

\newpage

# Supplemental Figure 2 {-}
![](figures/smc.pdf)
