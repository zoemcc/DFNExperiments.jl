begin 
    using DFNExperiments, NeuralPDE, DiffEqFlux
    using ModelingToolkit, Symbolics, DomainSets
    import ModelingToolkit: Interval, infimum, supremum
    using IfElse
    import Term as Terminal
end

model = ReducedCPhiJModel()
model_file = abspath(joinpath(@__DIR__, "..", "models", "reduced_c_phi_j", "model_subdomain_design.jl"))
include(model_file)


pde_system
subdomain_relations
eqs_integration_domains
ics_bcs_integration_domains

nn_spec = SimpleFeedForwardNetwork(
    1,
    8,
    SigmoidNonLin(),
    GlorotUniformParams()
)

nn_fc, init_params = NeuralPDE.getfunction(nn_spec, fill(2,4))
flat_init_params = vcat(init_params...)
strategy = StochasticTraining(16, 16)
discretization = PhysicsInformedNN(nn_fc, strategy; init_params = init_params)
sym_prob = NeuralPDE.symbolic_discretize(pde_system, discretization)
prob = NeuralPDE.discretize(pde_system, discretization)
sym_prob = NeuralPDE.symbolic_discretize(pde_system, discretization)
nothing
# current error:
# ERROR: BoundsError: attempt to access 2×16 Matrix{Float32} at index [[3], 1:16]
prob.f(flat_init_params, Float32[]) # currently doesn't work but will work once I fix this
nothing
