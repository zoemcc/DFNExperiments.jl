begin
using IfElse
# ('negative electrode', 'separator', 'positive electrode') -> x
@parameters t x
independent_variables_to_pybamm_names = Dict(
  :t => "Time",
  :x => "negative electrode",
)
# 'Electrolyte concentration' -> u1
# 'Electrolyte potential' -> u2
@variables u1(..) u2(..)
dependent_variables_to_pybamm_names = Dict(
  :u1 => "Electrolyte concentration",
  :u2 => "Electrolyte potential",
)
dependent_variables_to_dependencies = Dict(
  :u1 => (:t, :x),
  :u2 => (:t, :x),
)
Dt = Differential(t)
Dx = Differential(x)

# 'Electrolyte concentration' equation

function concatenation(x, n, s, p)
   # A concatenation in the electrolyte domain
   IfElse.ifelse(
      x < 0.4444444444444445, n, IfElse.ifelse(
         x < 0.5555555555555556, s, p
      )
   )
end

cache_4457050546309239314 = concatenation(x, 2.25, 0.0, -2.25)
cache_m4416472977406629345 = (Dx(Dx(u1(t, x)))) + cache_4457050546309239314

# 'Electrolyte potential' equation
cache_4457050546309239314 = concatenation(x, 2.25, 0.0, -2.25)
cache_m1466633705168541415 = (Dx(Dx(u1(t, x)) / u1(t, x) - Dx(u2(t, x)))) - cache_4457050546309239314


eqs = [
   Dt(u1(t, x)) ~ cache_m4416472977406629345,
   0 ~ cache_m1466633705168541415,
]


ics_bcs = [
   # initial conditions
   u1(0, x) ~ 1.0,
   u2(0, x) ~ 0.0,
   # boundary conditions
   Dx(u1(t, 0.0)) ~ 0.0,
   Dx(u1(t, 1.0)) ~ 0.0,
   u2(t, 0.0) ~ 0.0,
   Dx(u2(t, 1.0)) ~ 0.0,
]

t_domain = Interval(0.000,1.000)
x_domain = Interval(0.0, 1.0)

domains = [
   t in t_domain,
   x in x_domain,
]

ind_vars = [t, x]
dep_vars = [u1(t, x), u2(t, x)]

@named reduced_c_phi_pde_system = PDESystem(eqs, ics_bcs, domains, ind_vars, dep_vars)
pde_system = reduced_c_phi_pde_system

nothing
end
