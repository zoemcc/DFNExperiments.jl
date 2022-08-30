begin
using IfElse
# ('negative electrode', 'separator', 'positive electrode') -> x
@parameters t x
independent_variables_to_pybamm_names = Dict(
  :t => "Time",
  :x => "negative electrode",
)
# 'Electrolyte concentration' -> c_e
@variables c_e(..)
dependent_variables_to_pybamm_names = Dict(
  :c_e => "Electrolyte concentration",
)
dependent_variables_to_dependencies = Dict(
  :c_e => (:t, :x),
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

cache_m8587173813837221631 = concatenation(x, 2.25, 0.0, -2.25)
cache_m5374597158501241961 = (Dx(Dx(c_e(t, x)))) + cache_m8587173813837221631


eqs = [
   Dt(c_e(t, x)) ~ cache_m5374597158501241961,
]


ics_bcs = [
   # initial conditions
   c_e(0, x) ~ 1.0,
   # boundary conditions
   Dx(c_e(t, 0.0)) ~ 0.0,
   Dx(c_e(t, 1.0)) ~ 0.0,
]

t_domain = Interval(0.000,1.000)
x_domain = Interval(0.0, 1.0)

domains = [
   t in t_domain,
   x in x_domain,
]

ind_vars = [t, x]
dep_vars = [c_e(t, x)]

@named reduced_c_pde_system = PDESystem(eqs, ics_bcs, domains, ind_vars, dep_vars)
pde_system = reduced_c_pde_system

nothing
end
