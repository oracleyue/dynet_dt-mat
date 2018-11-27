# Toolbox for Dynamic Network Reconstruction 

This is a MATLAB toolbox of *discrete-time methods* for dynamic network reconstruction. 

(Note that this is not the original git repository for development. It is only used for released versions on github. The original one is maintained in a private gitLab repository.)

The core idea is to model causal interactions between measurable variables by a
dynamical network model, called *dynamical structure function* (DSF), and to
infer network structures (with feedback loops allowed) by identifying DSF
models. Moreover, it takes into account particular properties of biological time series:

-  limited lengths of time series,
-  low sampling frequencies, and
-  heterogeneity in experiment repetitions.

This is an interdisciplinary study, involving control engineering (system
identification, dynamical system theory), Bayesian statistics
(empirical/variational Bayes, EM, sampling methods, etc.)  and theoretical tools
from machine learning (RVM, Gaussian processes, etc.).

## File organisation

The folders in the project root are mainly two types:

- folders with name starting with capital letters are used for storing application/simualation data, outputs/results, which you may skip if you want to use this toolbox;
- folders with lower-case names (incl. starting with `@`) are part of toolbox.

The toolbox consists of the following essential directories:

- `@dynet` and `@idnet`: classes to define network models for identification, verification and analysis;
- `simulators`: functions to simulate DSF models or ARX networks to generate data;
- `solvers`: essential solvers for sparse parameter estimation;
- `auxiliary`: miscellany of the third-party MATLAB functions;
- `demo`: a collection of demos, which help you to quick start using the toolbox.

And the data/build directories:

- `Applications`: real biological/medical data and analysis scripts; it should
  be empty if released in public (like github);
- `Data`: the default folder to store simulation data of network models for further inference;
- `Results`: the default folder to save inference results;
- `Backup` and `Test`: the folders for backup or temporary scripts, which shouldn't appear in released versions.

If you didn't see the folders `Data` and `Results`, you may create them by `mkdir Data Results`.

## Optional dependency

If you like to use CVX solver in the toolbox, you need to setup CVX in your MATLAB, which can be download from the [CVX website](http://cvxr.com/cvx/download/).

## Examples

You may go to `PROJECTROOT/demos/` to find demo examples to use this toolbox.

**More simple demos** (to be completed)

## References

1. [Dynamic Network Reconstruction from Heterogeneous Datasets](https://arxiv.org/abs/1612.01963)
2. [Dynamic Network Reconstruction in Systems Biology: Methods and Algorithms](http://publications.uni.lu/handle/10993/35580)