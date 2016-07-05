# Simulations

This folder contains all simulations and computational experiments.

The experiments are described in the YAML files and run with the script
`expt.py`. For example, to run the "kernel-shapes" experiment, run
`./expt.py kernel-shapes.yaml`. Each experiment has its own folder, with
subfolders for each of its steps.

Each experiment has an accompanying SQLite database, which is created and
updated by `expt.py`. It contains the location of each file and all the
parameters which went into creating the file. The names of the files in the
experiment directories are not themselves descriptive, although most of the
files have their parameters recorded in a comment inside the file.
