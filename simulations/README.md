# Simulations

This folder contains all simulations and computational experiments.

The experiments are described in the YAML files and run with the script
`expt.py`. For example, to run the "kernel-shapes" experiment, run
`./expt.py kernel-shapes.yaml`. Each experiment has its own folder, with
subfolders for each of its steps.

The `expt.py` script is hacked together and not very robust. If a step fails or
is interrupted, it will not know to delete the file, and you will have to delete
it manually in order to run the step properly. 

Each experiment has an accompanying SQLite database, which is created and
updated by `expt.py`. It contains the location of each file and all the
parameters which went into creating the file. The names of the files in the
experiment directories are not themselves descriptive, although most of them
have their parameters embedded in some way.

Below I describe what each experiment is on a high level. In the "YAML files"
section, I detail the structure of the YAML files specifying the experiments.

## Experiments

All of the kernel-* experiments are testing the tree kernel's ability to
estimate a particular network parameter. Three or four different values of the
parameter are chosen, and 100 different networks are generated for each value.
A tree is simulated from each network. A kernel matrix is computed for the
trees for several combinations of its two parameters (the decay factor and
radial basis function variance). Each of these matrices is used to evaluate the
performance of a kSVM classifier by performing 1000 two-fold cross-validations.
We also calculate and plot some tree summary statistics.

### kernel-edge-density

Edge density parameter of a random (also known as Erdos-Renyi or GNP) network.
It is the probability that a particular edge exists. 

### kernel-mean-degree

The mean degree of a preferential attachment network, which is twice the number
of edges added for each vertex.

### kernel-pa-power

The preferential attachment power in a preferential attachment network. The
probability of joining with a vertex of degree d is proportional to d^(1+power).

### kernel-pareto

The shape parameter of a Pareto distribution from which degrees are drawn.

### kernel-shapes

The overall network type (random vs. preferential attachment vs. small world).

### kernel-smallworld

The rewiring probability for a small world network.

## YAML files

The YAML files should be fairly self explanatory. Experiments are organized into
steps, which have parameters and recipes, and depend on other steps. Below I
describe each of the fields in the YAML files. To create a new experiment, it is
probably easiest to just copy and modify one of the existing files.

### Global fields

These fields apply to the entire experiment.

#### Name

The name of the experiment. This is used to create the data folder and database.

#### Description

Just for documentation, it's not used for anything.

#### Processes

The number of processes to spawn, if you're using a multi-threaded machine. Note
that these are actual processes, not threads. We let the operating system handle
which process is distributed to which core. Also, tasks are distributed to
processes in a very unintelligent way - we just give each one an equal number of
tasks and wait for them all to finish. So it's possible that at the end of the
step, only one process will be left running, because it happened to be assigned
tasks that took longer. It was more trouble than it was worth to implement
process monitoring and so on.

### Step-specific fields

These parameters apply to individual steps.

#### Extension

The created files will all be named with this extension. Note that a step can
only create one file per rule. I realize that this is a pretty severe
limitation, but I had no need to overcome it, so I did not.

#### Processes

If you want to override the global number of processes, it can be set per-step.
This is useful if your rule is a multi-threaded command run within a single
interpreter.

#### Parameters

Parameters for a step are expanded combinatorially. For example, consider a step
called `network` with three parameters as follows.

    Parameters:
        nnode: 5000
        mean_degree: [4, 8]
        net_type: ["BA", "ER"]

This means that the step will be run 6 times, once for each possible combination
of parameters: `{nnode: 5000, mean_degree: 4, net_type: "BA"}`, 
`{nnode: 5000, mean_degree: 4, net_type: "ER"}`, and so on. The parameters
are substituted into the rule using placeholders in curly braces (see Rule,
below).

Parameters also affect which files from previous steps will be used to run the
step's rule. Continuing the example, suppose we have another step called `tree`
with three parameters, two of which overlap with the step we're looking at now.

    Parameters:
        nnode: 5000
        mean_degree: [4, 8]
        replicate: [1, 2, 3, 4, 5]

This step will be run 10 times, for each combination of parameters. The first
combination is `{nnode: 5000, mean_degree: 4, replicate: 1}'. This step depends
on `network`, so we pull in all the matching files from that step. There are two
of them, with parameters `{nnode: 5000, mean_degree: 4, net_type: "BA"}` and 
`{nnode: 5000, mean_degree: 4, net_type: "ER"}`. 

Notice that the `replicate` parameter has no effect on which files are pulled
in, because the previous `tree` step doesn't have that parameter. Likewise,
since the `tree` step doesn't have a `net_type` parameter, we don't filter out
files from the `network` step based on that parameter. Importantly, this means
that _all_ files in the prerequisite step are used by default, and the
parameters act as filters.

Currently, `expt.py` does not handle steps with no parameters. If the step has
no parameters, you should just specify one and not use it, like so:

    Parameters:
        placeholder: 0

#### Exclusions

Combinations of parameters which aren't to be used. In the `tree` step above,
the below exclusions will result in only the first three replicates being
created for `mean_degree = 4`.

    Exclusions:
        -
            mean_degree: 4
            replicate: [4, 5]

#### Interpreter

The program which is used to run the rule. Some examples are `bash` (for an
ordinary shell command) or `R --silent --vanilla` for an R script.

#### Startup

Any input to feed to the interpreter right after it's opened, that only needs to
be done once at the start. For example, loading packages.

#### Rule

The actual command, or sequence of commands, used to carry out the step. It will
be fed to the open interpreter. The parameters are available by placeholders in
curly braces. For example, for the `network` step above, the `nnode` parameter
can be accessed like so.

    x <- rep(0, {nnode})

To include actual curly braces in the rule, you have to double them up.

    if ({nnode} == 1000) {{
        ...
    }}

The name of the output file is accessible by the name of the step. Here, the
step is called `tree`.

    echo "this is not a tree" >> {tree}

The prerequisite files are also available by their step names. If there are
multiple prerequisites for a step, the placeholder will be replaced by a
_space-delimited_ list of the files. For example, if `tree` depends on
`network`, we can do this.

    for N in {network}; do
        echo "I am a network!"
    done

Note that, because the file names are separated by a space, you may have to
split them up depending on the language.

    network_list = "{network}".split(" ")

There are currently two special placeholders. `{yaml}` will give you a
YAML-formatted string containing all the parameters for the current step.
And `{$#}` will be replaced with the total number of prequisite files.
