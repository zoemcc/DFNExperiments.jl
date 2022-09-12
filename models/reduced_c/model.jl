begin
using IfElse
# ('negative electrode',) -> x_n
# ('separator',) -> x_s
# ('positive electrode',) -> x_p
@parameters t x_n x_s x_p
independent_variables_to_pybamm_names = Dict(
  :t => "Time",
  :x_n => "negative electrode",
  :x_s => "separator",
  :x_p => "positive electrode",
)
# 'Electrolyte concentration neg' -> c_e_n
# 'Electrolyte concentration sep' -> c_e_s
# 'Electrolyte concentration pos' -> c_e_p
@variables c_e_n(..) c_e_s(..) c_e_p(..)
dependent_variables_to_pybamm_names = Dict(
  :c_e_n => "Electrolyte concentration neg",
  :c_e_s => "Electrolyte concentration sep",
  :c_e_p => "Electrolyte concentration pos",
)
dependent_variables_to_dependencies = Dict(
  :c_e_n => (:t, :x_n),
  :c_e_s => (:t, :x_s),
  :c_e_p => (:t, :x_p),
)
Dt = Differential(t)
Dx_n = Differential(x_n)
Dx_s = Differential(x_s)
Dx_p = Differential(x_p)

# 'Electrolyte concentration' equation
cache_m5729566733887177216_n = 2.25
cache_m5729566733887177216_s = 0.0
cache_m5729566733887177216_p = -2.25
cache_8327581641087914053_n = (Dx_n(Dx_n(c_e_n(t, x_n)))) + cache_m5729566733887177216_n
cache_8327581641087914053_s = (Dx_s(Dx_s(c_e_s(t, x_s)))) + cache_m5729566733887177216_s
cache_8327581641087914053_p = (Dx_p(Dx_p(c_e_p(t, x_p)))) + cache_m5729566733887177216_p
eqs = [
   Dt(c_e_n(t, x_n)) ~ cache_8327581641087914053_n,
   Dt(c_e_s(t, x_s)) ~ cache_8327581641087914053_s,
   Dt(c_e_p(t, x_p)) ~ cache_8327581641087914053_p,
]


ics_bcs = [
   # initial conditions
   c_e_n(0, x_n) ~ 1.0,
   c_e_s(0, x_s) ~ 1.0,
   c_e_p(0, x_p) ~ 1.0,
   # boundary conditions
   Dx_n(c_e_n(t, 0.0)) ~ 0.0,
   c_e_n(t, 0.4444444444444445) ~ c_e_s(t, 0.4444444444444445),
   c_e_p(t, 0.5555555555555556) ~ c_e_s(t, 0.5555555555555556),
   Dx_p(c_e_p(t, 1.0)) ~ 0.0,
]

t_domain = Interval(0.000,1.000)
x_n_domain = Interval(0.0, 0.4444444444444445)
x_s_domain = Interval(0.4444444444444445, 0.5555555555555556)
x_p_domain = Interval(0.5555555555555556, 1.0)

domains = [
   t in t_domain,
   x_n in x_n_domain,
   x_s in x_s_domain,
   x_p in x_p_domain,
]

ind_vars = [t, x_n, x_s, x_p]
dep_vars = [c_e_n(t, x_n), c_e_s(t, x_s), c_e_p(t, x_p)]

@named reduced_c_pde_system = PDESystem(eqs, ics_bcs, domains, ind_vars, dep_vars)
pde_system = reduced_c_pde_system

nothing
end
