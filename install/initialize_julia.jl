using Pkg
Pkg.activate(abspath(joinpath(@__DIR__, "..")))
Pkg.develop(ENV["NEURALPDE_LOCATION"]) # make sure the hyperparams branch is checked out and then develop using that branch
ENV["PYTHON"] = joinpath(ENV["ANACONDA_LOCATION"], "envs", "pybamm_dev", "bin", "python") # this is the env var PyCall will use to find python
Pkg.build("PyCall") # initialize PyCall to use the above python
