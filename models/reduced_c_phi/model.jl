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
# 'Negative electrolyte potential' -> phi_e_n
# 'Separator electrolyte potential' -> phi_e_s
# 'Positive electrolyte potential' -> phi_e_p
@variables c_e_n(..) c_e_s(..) c_e_p(..) phi_e_n(..) phi_e_s(..) phi_e_p(..)
dependent_variables_to_pybamm_names = Dict(
  :c_e_n => "Electrolyte concentration neg",
  :c_e_s => "Electrolyte concentration sep",
  :c_e_p => "Electrolyte concentration pos",
  :phi_e_n => "Negative electrolyte potential",
  :phi_e_s => "Separator electrolyte potential",
  :phi_e_p => "Positive electrolyte potential",
)
dependent_variables_to_dependencies = Dict(
  :c_e_n => (:t, :x_n),
  :c_e_s => (:t, :x_s),
  :c_e_p => (:t, :x_p),
  :phi_e_n => (:t, :x_n),
  :phi_e_s => (:t, :x_s),
  :phi_e_p => (:t, :x_p),
)
Dt = Differential(t)
Dx_n = Differential(x_n)
Dx_s = Differential(x_s)
Dx_p = Differential(x_p)

# 'Electrolyte concentration' equation
cache_2342231486479082649_n = 2.25
cache_2342231486479082649_s = 0.0
cache_2342231486479082649_p = -2.25
cache_6108101611517005270_n = (Dx_n(Dx_n(c_e_n(t, x_n)))) + cache_2342231486479082649_n
cache_6108101611517005270_s = (Dx_s(Dx_s(c_e_s(t, x_s)))) + cache_2342231486479082649_s
cache_6108101611517005270_p = (Dx_p(Dx_p(c_e_p(t, x_p)))) + cache_2342231486479082649_p

# 'Electrolyte potential' equation
cache_2342231486479082649_n = 2.25
cache_2342231486479082649_s = 0.0
cache_2342231486479082649_p = -2.25
cache_8603744982470494737_n = (Dx_n(Dx_n(c_e_n(t, x_n)) / c_e_n(t, x_n) - Dx_n(phi_e_n(t, x_n)))) - cache_2342231486479082649_n
cache_8603744982470494737_s = (Dx_s(Dx_s(c_e_s(t, x_s)) / c_e_s(t, x_s) - Dx_s(phi_e_s(t, x_s)))) - cache_2342231486479082649_s
cache_8603744982470494737_p = (Dx_p(Dx_p(c_e_p(t, x_p)) / c_e_p(t, x_p) - Dx_p(phi_e_p(t, x_p)))) - cache_2342231486479082649_p
eqs = [
   Dt(c_e_n(t, x_n)) ~ cache_6108101611517005270_n,
   Dt(c_e_s(t, x_s)) ~ cache_6108101611517005270_s,
   Dt(c_e_p(t, x_p)) ~ cache_6108101611517005270_p,
   0 ~ cache_8603744982470494737_n,
   0 ~ cache_8603744982470494737_s,
   0 ~ cache_8603744982470494737_p,
]


ics_bcs = [
   # initial conditions
   c_e_n(0, x_n) ~ 1.0,
   c_e_s(0, x_s) ~ 1.0,
   c_e_p(0, x_p) ~ 1.0,
   phi_e_n(0, x_n) ~ 0.0,
   phi_e_s(0, x_s) ~ 0.0,
   phi_e_p(0, x_p) ~ 0.0,
   # boundary conditions
   Dx_n(c_e_n(t, 0.0)) ~ 0.0,
   c_e_n(t, 0.4444444444444445) ~ c_e_s(t, 0.4444444444444445),
   c_e_p(t, 0.5555555555555556) ~ c_e_s(t, 0.5555555555555556),
   Dx_p(c_e_p(t, 1.0)) ~ 0.0,
   phi_e_n(t, 0.0) ~ 0.0,
   phi_e_n(t, 0.4444444444444445) ~ phi_e_s(t, 0.4444444444444445),
   phi_e_p(t, 0.5555555555555556) ~ phi_e_s(t, 0.5555555555555556),
   Dx_p(phi_e_p(t, 1.0)) ~ 0.0,
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
dep_vars = [c_e_n(t, x_n), c_e_s(t, x_s), c_e_p(t, x_p), phi_e_n(t, x_n), phi_e_s(t, x_s), phi_e_p(t, x_p)]

@named reduced_c_phi_pde_system = PDESystem(eqs, ics_bcs, domains, ind_vars, dep_vars)
pde_system = reduced_c_phi_pde_system

nothing
end
