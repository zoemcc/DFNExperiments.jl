#ENV["PYTHON"] = joinpath(homedir(), "anaconda3", "envs", "pybamm_dev", "bin", "python")
#using Pkg
#Pkg.build("PyCall")
using Distributed
addprocs(1)
@everywhere begin 
    using Pkg
    Pkg.activate(abspath(joinpath(@__DIR__, "..")))
end
@everywhere begin 
    using DFNExperiments
    using ModelingToolkit, Symbolics, DomainSets
    import ModelingToolkit: Interval, infimum, supremum
    using LabelledArrays
    using NeuralPDE
    using PyCall
    using JSON
end

model = SPMModel()
models_dir = abspath(joinpath(@__DIR__, "..", "models"))
output_dir = joinpath(models_dir, "$(pybamm_func_str(model))")
regenerate_sim = true
if regenerate_sim
    pybamm = pyimport("pybamm")
    current_input = false
    sim_data, pde_system = generate_sim_model(model; current_input=current_input, output_dir=output_dir)
else 
    include(joinpath(output_dir, "model.jl"))
end

sim_data = read_sim_data(output_dir)
@everywhere sim_data = $sim_data
@everywhere pde_system = $pde_system

sg = StructGenerator(
    :CompositeHyperParameter,
    RandomChoice(1:2^10), # seed
    StructGenerator( # nn
        :SimpleFeedForwardNetwork, # type/constructor name
        RandomChoice(1:2),
        RandomChoice(10, 20, 30),
        RandomChoice(:GELUNonLin, :SigmoidNonLin),
        :GlorotUniformParams
    ),
    StructGenerator( # training
        :GridTraining,
        RandomChoice(0.1, 0.2, 0.06)
    ),
    StructGenerator( # adaptive loss
        :MiniMaxAdaptiveLoss,
        40
    ),
    RandomChoice( # optimizer
        StructGenerator(:ADAMOptimiser, 10000, 1e-2),
        StructGenerator(:ADAMOptimiser, 10000, 1e-3)
    )
)


hyperparametersweep = StructGeneratorHyperParameterSweep(1, 16, sg)
hyperparameters = generate_hyperparameters(hyperparametersweep)

@everywhere function get_cb()
    cb = function (p,l)
        return false
    end
    return cb
end

ts = sim_data[:ivs][:t]
r_ns = sim_data[:ivs][:r_n]
r_ps = sim_data[:ivs][:r_p]
qs = sim_data[:dvs][:Q_Ah]
csns = sim_data[:dvs][:c_s_n_xav]
csps = sim_data[:dvs][:c_s_p_xav]
@everywhere function get_plot_function()
    #xs,ys = [infimum(d.domain):0.01:supremum(d.domain) for d in domains]
    
    analytic_sol_func(x,y) = (sin(pi*x)*sin(pi*y))/(2pi^2)
    u_real = reshape([analytic_sol_func(x,y) for x in xs for y in ys], (length(xs),length(ys)))
    p1 = plot(xs, ys, u_real, linetype=:contourf,title = "analytic");
    function plot_function(logger, step, phi, θ, adaloss)
        u_predict = reshape([first(phi[1]([x,y],θ)) for x in xs for y in ys],(length(xs),length(ys)))
        diff_u = abs.(u_predict .- u_real)


        p2 = plot(xs, ys, u_predict, linetype=:contourf,title = "predict");
        p3 = plot(xs, ys, diff_u,linetype=:contourf,title = "error");
        [(name="analytic", image=p1), (name="predict", image=p2), (name="error", image=p3)]
    end
    return plot_function
end

log_options = NeuralPDE.LogOptions(;plot_function=get_plot_function())

neuralpde_workers = map(NeuralPDE.NeuralPDEWorker, workers())
cb_func = get_cb()
experiment_manager = NeuralPDE.ExperimentManager(pde_system, hyperparameters, cb_func, log_options, neuralpde_workers)


NeuralPDE.run_experiment_queue(experiment_manager)

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

