begin
using IfElse
@parameters t
independent_variables_to_pybamm_names = Dict(
  :t => "Time",
)
# 'Discharge capacity [A.h]' -> Q_Ah
# 'Average negative particle concentration' -> c_s_n_av
# 'Average positive particle concentration' -> c_s_p_av
@variables Q_Ah(..) c_s_n_av(..) c_s_p_av(..)
dependent_variables_to_pybamm_names = Dict(
  :Q_Ah => "Discharge capacity [A.h]",
  :c_s_n_av => "Average negative particle concentration",
  :c_s_p_av => "Average positive particle concentration",
)
dependent_variables_to_dependencies = Dict(
  :Q_Ah => (:t,),
  :c_s_n_av => (:t,),
  :c_s_p_av => (:t,),
)
Dt = Differential(t)

# 'Average negative particle concentration' equation
cache_143377372744660759 = massc_s_n_av(t) @ -3.7500000000000004

# 'Average positive particle concentration' equation
cache_5101339567594643364 = massc_s_p_av(t) @ 2.195026006381394


eqs = [
   Dt(Q_Ah(t)) ~ 4.27249308415467,
   Dt(c_s_n_av(t)) ~ cache_143377372744660759,
   Dt(c_s_p_av(t)) ~ cache_5101339567594643364,
]


ics_bcs = [
   # initial conditions
   Q_Ah(0) ~ 0.0,
   c_s_n_av(0, z) ~ 0.8000000000000016,
   c_s_p_av(0, z) ~ 0.6000000000000001,
   # boundary conditions
]

t_domain = Interval(0.000,0.159)

domains = [
   t in t_domain,
]

ind_vars = [t]
dep_vars = [Q_Ah(t), c_s_n_av(t), c_s_p_av(t)]

@named SPM_no_r_pde_system = PDESystem(eqs, ics_bcs, domains, ind_vars, dep_vars)
pde_system = SPM_no_r_pde_system

nothing
end
