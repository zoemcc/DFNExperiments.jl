begin 
    using Term
    using DFNExperiments
    using ModelingToolkit, NeuralPDE, DiffEqFlux
    import ModelingToolkit: Interval, infimum, supremum
    using IfElse
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

dummy_eq = [u1(t,x) ~ u1(t,x)]
dummy_ic = [u1(0.0,x) ~ u1(0.0,x)]
strategy = StochasticTraining(16, 16)
discretization = PhysicsInformedNN(nn_fc, strategy; init_params = init_params)
pde_system_main = pde_system
pde_system = pde_system_main
pde_system = PDESystem(pde_system_main.eqs[[1]], pde_system_main.bcs[[1]], pde_system_main.domain, pde_system_main.ivs, pde_system_main.dvs; name=:pde_system_small)
pde_system = PDESystem(pde_system_main.eqs[[1]], dummy_ic, pde_system_main.domain, pde_system_main.ivs, pde_system_main.dvs; name=:pde_system_small)
pde_system = PDESystem(pde_system_main.eqs[[1]], dummy_eq, pde_system_main.domain, pde_system_main.ivs, pde_system_main.dvs; name=:pde_system_small)
pde_system = PDESystem(dummy_eq, pde_system_main.bcs[[1]], pde_system_main.domain, pde_system_main.ivs, pde_system_main.dvs; name=:pde_system_small)
pde_system = PDESystem(dummy_eq, dummy_ic, pde_system_main.domain, pde_system_main.ivs, pde_system_main.dvs; name=:pde_system_small)
begin 
    sym_prob = NeuralPDE.symbolic_discretize(pde_system, discretization; 
        subdomain_relations=subdomain_relations, 
        eqs_integration_domains=eqs_integration_domains,
        ics_bcs_integration_domains=ics_bcs_integration_domains,
        #ics_bcs_integration_domains=[[:t, :x]],
    )
    println(" ")
    println(sym_prob)
    #println(sym_prob.args[2].args[3].args[3])
end
prob = NeuralPDE.discretize(pde_system, discretization; 
    subdomain_relations=subdomain_relations, 
    eqs_integration_domains=eqs_integration_domains,
    ics_bcs_integration_domains=ics_bcs_integration_domains,
    #ics_bcs_integration_domains=[[:t, :x]],
)
#inspect(sym_prob)
#expressiontree(sym_prob)
nothing
# current error:
# ERROR: BoundsError: attempt to access 2×16 Matrix{Float32} at index [[3], 1:16]
prob.f(flat_init_params, Float32[]) # currently doesn't work but will work once I fix this
nothing
