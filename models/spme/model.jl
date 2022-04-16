begin
using IfElse
# ('negative electrode', 'separator', 'positive electrode') -> x
# ('negative particle',) -> r_n
# ('positive particle',) -> r_p
@parameters t x r_n r_p
independent_variables_to_pybamm_names = Dict(
  :t => "Time",
  :x => "negative electrode",
  :r_n => "negative particle",
  :r_p => "positive particle",
)
# 'Discharge capacity [A.h]' -> Q_Ah
# 'X-averaged negative particle concentration' -> c_s_n_xav
# 'X-averaged positive particle concentration' -> c_s_p_xav
# 'Porosity times concentration' -> eps_c_e
@variables Q_Ah(..) c_s_n_xav(..) c_s_p_xav(..) eps_c_e(..)
dependent_variables_to_pybamm_names = Dict(
  :Q_Ah => "Discharge capacity [A.h]",
  :c_s_n_xav => "X-averaged negative particle concentration",
  :c_s_p_xav => "X-averaged positive particle concentration",
  :eps_c_e => "Porosity times concentration",
)
dependent_variables_to_dependencies = Dict(
  :Q_Ah => (:t,),
  :c_s_n_xav => (:t, :r_n),
  :c_s_p_xav => (:t, :r_p),
  :eps_c_e => (:t, :x),
)
Dt = Differential(t)
Dx = Differential(x)
Dr_n = Differential(r_n)
Dr_p = Differential(r_p)

# 'X-averaged negative particle concentration' equation
cache_m8822410815049285743 = 8.813457647415216 * (1 / r_n^2 * Dr_n(r_n^2 * Dr_n(c_s_n_xav(t, r_n))))

# 'X-averaged positive particle concentration' equation
cache_6957192350346960075 = 22.598609352346717 * (1 / r_p^2 * Dr_p(r_p^2 * Dr_p(c_s_p_xav(t, r_p))))

# 'Porosity times concentration' equation

function concatenation(x, n, s, p)
   # A concatenation in the electrolyte domain
   IfElse.ifelse(
      x < 0.4444444444444445, n, IfElse.ifelse(
         x < 0.5555555555555556, s, p
      )
   )
end

function D_e(c_e, T)
   5.34e-10 * exp(-0.00065 * c_e) * exp(14.941767670181717 - (4454.8880308646785 * 1.0 / T))
end

cache_461017418220814889 = concatenation(x, -589429731.3893579, -3587155110.512914, -589429731.3893579)
cache_m1984000167394293568 = concatenation(x, eps_c_e(t, x) * 3333.3333333333335, 
   eps_c_e(t, x) * 1000.0
, 
   eps_c_e(t, x) * 3333.3333333333335
)
cache_6602808265225132919 = concatenation(x, eps_c_e(t, x) / 0.3, eps_c_e(t, x), 
   eps_c_e(t, x) / 0.3
)
cache_1221441731471020359 = concatenation(x, 0.1806862603261852 * x, 0.08030500458941565, 
   (0.007232292579356106 - (0.007232292579356106 * x)) / 0.04002679875215723
)
cache_9041100221364284082 = concatenation(x, 56.212339486148316, 0.0, -56.212339486148316)
cache_m5489219686130520664 = ((Dx(((cache_461017418220814889 * (D_e(cache_m1984000167394293568, 298.15))) * Dx(cache_6602808265225132919)) + cache_1221441731471020359)) / -0.008035880643729006) + cache_9041100221364284082


eqs = [
   Dt(Q_Ah(t)) ~ 4.27249308415467,
   Dt(c_s_n_xav(t, r_n)) ~ cache_m8822410815049285743,
   Dt(c_s_p_xav(t, r_p)) ~ cache_6957192350346960075,
   Dt(eps_c_e(t, x)) ~ cache_m5489219686130520664,
]

# 'Porosity times concentration' initial condition
cache_5346377623449515321 = concatenation(x, 0.3, 1.0, 0.3)


ics_bcs = [
   # initial conditions
   Q_Ah(0) ~ 0.0,
   c_s_n_xav(0, r_n) ~ 0.8000000000000016,
   c_s_p_xav(0, r_p) ~ 0.6000000000000001,
   eps_c_e(0, x) ~ cache_5346377623449515321,
   # boundary conditions
   Dr_n(c_s_n_xav(t, 0.0)) ~ 0.0,
   Dr_n(c_s_n_xav(t, 1.0)) ~ -0.14182855923368468,
   Dr_p(c_s_p_xav(t, 0.0)) ~ 0.0,
   Dr_p(c_s_p_xav(t, 1.0)) ~ 0.03237700710041634,
]

t_domain = Interval(0.000,0.159)
x_domain = Interval(0.0, 1.0)
r_n_domain = Interval(0.0, 1.0)
r_p_domain = Interval(0.0, 1.0)

domains = [
   t in t_domain,
   x in x_domain,
   r_n in r_n_domain,
   r_p in r_p_domain,
]

ind_vars = [t, x, r_n, r_p]
dep_vars = [Q_Ah(t), c_s_n_xav(t, r_n), c_s_p_xav(t, r_p), eps_c_e(t, x)]

@named SPMe_pde_system = PDESystem(eqs, ics_bcs, domains, ind_vars, dep_vars)
pde_system = SPMe_pde_system

nothing
end
