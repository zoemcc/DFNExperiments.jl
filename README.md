# DFNExperiments

Tighter integration between PyBaMM for model spec generation and simulation and NeuralPDE.jl for PINN training.

Requires specific branches of both PyBaMM and NeuralPDE.jl:

PyBaMM branch:
https://github.com/zoemcc/PyBaMM/tree/issue-1129-julia

NeuralPDE branch:
https://github.com/zoemcc/NeuralPDE.jl/tree/hyperparam

To get PyBaMM to work I used Anaconda python3 on python version 3.7, with a specific conda env just for this project.  
The scripts in the install folder are designed to help with this process.
Git clone the above branches, install Anaconda and then set the following environment variables before you run the install script
(where the actual paths you used to install the various things are used instead of the placeholder strings):

```sh
export PYBAMM_LOCATION="/path/to/pybamm_repo"
export DFNEXPERIMENTS_LOCATION="/path/to/this_repo"
export NEURALPDE_LOCATION="/path/to/neuralpde_repo"
export ANACONDA_LOCATION="/path/to/anaconda3"
```

and then run:

```sh
sh $DFNEXPERIMENTS_LOCATION/install/install_script.sh
```

and hopefully everything should get setup. 
The steps are basically :

create conda env with the right python version
activate that python env
install the right pybamm version to that python env 
and then in a julia environment:
set NeuralPDE to develop with the right version
set the PYTHON environment variable to the new python env
build PyCall in julia using the right python that has the right pybamm installed
