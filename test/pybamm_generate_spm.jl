#ENV["PYTHON"] = joinpath(homedir(), "anaconda3", "envs", "pybamm_dev", "bin", "python")
#using Pkg
#Pkg.build("PyCall")
using DFNExperiments
using ModelingToolkit, Symbolics, DomainSets
using PyCall
using JSON
pybamm = pyimport("pybamm")

current_input = false
model = SPMModel()
models_dir = abspath(joinpath(@__DIR__, "..", "models"))
output_dir = joinpath(models_dir, "$(pybamm_func_str(model))")
sim_data, pde_system = generate_sim_model(model; current_input=current_input, output_dir=output_dir)

sim_json_data = read_sim_data(output_dir)

nothing

# Plot pybamm solution
"""
include("unpickle.jl")
pybamm_sols = myunpickle("pybamm/hardcoded_models/MTK_format/pybamm_solutions/SPM.pickle")

t_pb = pybamm_sols["Time"]
rn_pb = pybamm_sols["r_n"][:,1]
rp_pb = pybamm_sols["r_p"][:,1]
c_s_n_pb = pybamm_sols["X-averaged negative particle concentration"]
c_s_p_pb = pybamm_sols["X-averaged positive particle concentration"]
Q_pb = pybamm_sols["Discharge capacity [A.h]"]

using Plots
anim = @animate for i in 1:length(t_pb)
    p1 = plot(rn_pb,c_s_n_pb[:,i];ylims=(0,1),ylabel="c_s_n",xlabel="r_n",legend=false)
    p2 = plot(rp_pb,c_s_p_pb[:,i];ylims=(0,1),ylabel="c_s_p",xlabel="r_p",legend=false)
    p3 = plot(t_pb,Q_pb;ylims=(0,1),ylabel="Q",xlabel="t",legend=false)
    plot(p1,p2,p3)
end
gif(anim, "pybamm/hardcoded_models/MTK_format/animations/SPM.gif",fps=30)
"""

