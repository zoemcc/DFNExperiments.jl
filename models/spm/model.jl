begin
using IfElse
# ('negative particle',) -> r_n
# ('positive particle',) -> r_p
@parameters t r_n r_p
independent_variables_to_pybamm_names = Dict(
  :t => "Time",
  :r_n => "negative particle",
  :r_p => "positive particle",
)
# 'X-averaged negative particle concentration' -> c_s_n_xav
# 'X-averaged positive particle concentration' -> c_s_p_xav
@variables c_s_n_xav(..) c_s_p_xav(..)
dependent_variables_to_pybamm_names = Dict(
  :c_s_n_xav => "X-averaged negative particle concentration",
  :c_s_p_xav => "X-averaged positive particle concentration",
)
dependent_variables_to_dependencies = Dict(
  :c_s_n_xav => (:t, :r_n),
  :c_s_p_xav => (:t, :r_p),
)
Dt = Differential(t)
Dr_n = Differential(r_n)
Dr_p = Differential(r_p)

# 'X-averaged negative particle concentration' equation
cache_m1744780769296494788 = 8.813457647415216 * (1 / r_n^2 * Dr_n(r_n^2 * Dr_n(c_s_n_xav(t, r_n))))

# 'X-averaged positive particle concentration' equation
cache_m2497120790061496988 = 22.598609352346717 * (1 / r_p^2 * Dr_p(r_p^2 * Dr_p(c_s_p_xav(t, r_p))))
eqs = [
   Dt(c_s_n_xav(t, r_n)) ~ cache_m1744780769296494788,
   Dt(c_s_p_xav(t, r_p)) ~ cache_m2497120790061496988,
]


ics_bcs = [
   # initial conditions
   c_s_n_xav(0, r_n) ~ 0.8000000000000016,
   c_s_p_xav(0, r_p) ~ 0.6000000000000001,
   # boundary conditions
   Dr_n(c_s_n_xav(t, 0.0)) ~ 0.0,
   Dr_n(c_s_n_xav(t, 1.0)) ~ -0.14182855923368468,
   Dr_p(c_s_p_xav(t, 0.0)) ~ 0.0,
   Dr_p(c_s_p_xav(t, 1.0)) ~ 0.03237700710041634,
]

t_domain = Interval(0.000,0.159)
r_n_domain = Interval(0.0, 1.0)
r_p_domain = Interval(0.0, 1.0)

domains = [
   t in t_domain,
   r_n in r_n_domain,
   r_p in r_p_domain,
]

ind_vars = [t, r_n, r_p]
dep_vars = [c_s_n_xav(t, r_n), c_s_p_xav(t, r_p)]

@named SPM_pde_system = PDESystem(eqs, ics_bcs, domains, ind_vars, dep_vars)
pde_system = SPM_pde_system

nothing
end
