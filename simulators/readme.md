# Guideline for Generation and Simulation of Random Sparse Dynamic Network Model s


## Dynamical Structure Functions

1. Run the script `sSysSimDemo.m` （in `../demos/Basics/`） to generate random network models (DSF) and simulate them to generate time series (principle: generating A-matrix and then simulating the state-space system). 
    -  The data are saved in the folders `dataCSV`, `dataMAT` and the figures are in `figures`. 
    -  The path to save data/figures can be set in the call `dynet_sim.m`. See the demo for more examples.
    -  The writing and plotting behaviors of data/figures can be controlled via `export` argument in `dynet_sim.m`.

2. More about the function call of `dynet_sim.m`:
    -  The available types of input signals are *step signals* or *random Gaussian white noise*, controlled in `dynet_sim.m` via `Input` argument.
    -  The hidden states can be chosen randomly or fixed, controlled by the argument `sys_info.IdxNodes`.
    -  The input signals is fed into the system at randomly chosen positions (*B*: randomly choose *m* columns of the identity matrix). 
    -  The process noise randomly chooses *p* states to perturb in the *innovation form of state-space representations* (*K*: randomly choose *p* columns of the identity matrix).
    -  The initial state value are the standard Gaussian random vector multiplied by the magnitude given by users.

3. `dynet_sim.m` uses the following function to generate stable sparse A-matrix "randomly":
```
SparsityDensity = .001;
A  = dynet_sprandstab(n, 'nicolo-overlap', SparsityDensity, 5);
```

4. Modify `dynet_sim.m` if requiring more freedom, e.g. other input signals.


## ARX Models

1. The demo `sARXSimDemo.m` is placed in `ROOT/demos/Basics/`. One may run it to quickly generate models and obtain signals from simulations.

2. More about the supporting functions in `ROOT/simulators/`:
     - The stable sparse ARX models are generated by `arx_gen.m`. If you don't need to perform simulations, you may just call this function.
     - The function `arx_sim.m` generates models and performs simulations, and save input/output signals, by default, in `ARXSimData.mat`, if you didn't specify a new file name. In order to modify the simulation options, you have to provide three arguments: `dataOpt`, `netOpt`, `signalOpt`. To ease its usage, we have a function `arxsimOptions.m` to quickly generate these three structures. Then you may modify specific items, referring to the help document of `arx_sim.m` to see their member variables.
     - Regarding the implementation of `arx_gen.m`, it requires two supportive functions: `pathfind.m` and `stabpoly.m`, which are placed under `ROOT/simulators/members/`.