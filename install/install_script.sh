#!/bin/bash
cd $PYBAMM_LOCATION
pip install -e .
cd $DFNEXPERIMENTS_LOCATION
julia ./install/initialize_julia.jl
