conda create -n pybamm_dev python=3.7
conda activate pybamm_dev
cd $PYBAMM_LOCATION
pip install -e .
cd $DFNEXPERIMENTS_LOCATION
julia ./install/initialize_julia.jl
