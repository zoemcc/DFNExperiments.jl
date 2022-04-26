using DFNExperiments
using NeuralPDE, DiffEqFlux
using IfElse
using ModelingToolkit, Symbolics, SymbolicUtils, DomainSets

@parameters t 
@variables u(..)
Dt = Differential(t)

eqs = [
   Dt(u(t)) ~ IfElse.ifelse(t > 1.0, IfElse.ifelse(t > 2.0, 3.0, 2.0), 1.0),
]


ics_bcs = [
   u(0) ~ 0.0,
]

t_domain = Interval(0.0,3.0)

domains = [
   t in t_domain,
]

ind_vars = [t]
dep_vars = [u(t)]

@named ifelse_pde_system = PDESystem(eqs, ics_bcs, domains, ind_vars, dep_vars)
pde_system = ifelse_pde_system

analytic_sol(i) = IfElse.ifelse(i > 1, IfElse.ifelse(i > 2, 3 * (i - 2) + 3, 2 * (i - 1) + 1), i)

struct IfElseSolAnalytic
end

function (::IfElseSolAnalytic)(ts, θ) 
    @show ts
    @show size(ts)
    @show θ
    sol = similar(ts)
    for i in eachindex(ts)
        @show i
        @show ts[i]
        sol[i] = analytic_sol(ts[i])
    end
    sol
end

DiffEqFlux.initial_params(::IfElseSolAnalytic) = Float64[]

analytic_fc = [FastChain(IfElseSolAnalytic())]

strategy =  NeuralPDE.StochasticTraining(8, 8)
#derivative =  NeuralPDE.get_numeric_derivative(epsvar=1e-6)
discretization = NeuralPDE.PhysicsInformedNN(analytic_fc, strategy)
prob = NeuralPDE.discretize(pde_system, discretization)
sym_prob = NeuralPDE.symbolic_discretize(pde_system, discretization)
prob.f(Float64[], Float64[])

