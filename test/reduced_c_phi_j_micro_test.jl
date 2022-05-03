begin 
    using DFNExperiments, NeuralPDE, DiffEqFlux
    using ModelingToolkit, Symbolics, DomainSets
    import ModelingToolkit: Interval, infimum, supremum
    using IfElse
end

model = ReducedCPhiJModel()
model_file = abspath(joinpath(@__DIR__, "..", "models", "reduced_c_phi_j", "model_subdomain_design.jl"))
include(model_file)
nothing





