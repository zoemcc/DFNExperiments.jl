begin
using IfElse
# ('negative particle',) -> r_n
# ('positive particle',) -> r_p
# ('negative electrode',) -> x_n
# ('separator',) -> x_s
# ('positive electrode',) -> x_p
@parameters t r_n r_p x_n x_s x_p
independent_variables_to_pybamm_names = Dict(
  :t => "Time",
  :r_n => "negative particle",
  :r_p => "positive particle",
  :x_n => "negative electrode",
  :x_s => "separator",
  :x_p => "positive electrode",
)
# 'X-averaged negative particle concentration' -> c_s_n_xav
# 'X-averaged positive particle concentration' -> c_s_p_xav
# 'Negative electrolyte concentration' -> c_e_n
# 'Separator electrolyte concentration' -> c_e_s
# 'Positive electrolyte concentration' -> c_e_p
@variables c_s_n_xav(..) c_s_p_xav(..) c_e_n(..) c_e_s(..) c_e_p(..)
dependent_variables_to_pybamm_names = Dict(
  :c_s_n_xav => "X-averaged negative particle concentration",
  :c_s_p_xav => "X-averaged positive particle concentration",
  :c_e_n => "Negative electrolyte concentration",
  :c_e_s => "Separator electrolyte concentration",
  :c_e_p => "Positive electrolyte concentration",
)
dependent_variables_to_dependencies = Dict(
  :c_s_n_xav => (:t, :r_n),
  :c_s_p_xav => (:t, :r_p),
  :c_e_n => (:t, :x_n),
  :c_e_s => (:t, :x_s),
  :c_e_p => (:t, :x_p),
)
Dt = Differential(t)
Dr_n = Differential(r_n)
Dr_p = Differential(r_p)
Dx_n = Differential(x_n)
Dx_s = Differential(x_s)
Dx_p = Differential(x_p)

# 'X-averaged negative particle concentration' equation
cache_m4119150365269624523 = 8.813457647415216 * (1 / r_n^2 * Dr_n(r_n^2 * Dr_n(c_s_n_xav(t, r_n))))

# 'X-averaged positive particle concentration' equation
cache_5866437115522166420 = 22.598609352346717 * (1 / r_p^2 * Dr_p(r_p^2 * Dr_p(c_s_p_xav(t, r_p))))

# 'Electrolyte concentration' equation

function D_e(c_e, T)
   5.34e-10 * exp(-0.00065 * c_e) * exp(14.941767670181717 - (4454.8880308646785 * 1.0 / T))
end

cache_5380421293040872985_n = -589429731.3893579
cache_5380421293040872985_s = -3587155110.512914
cache_5380421293040872985_p = -589429731.3893579
cache_2434126157835956296_n = 
   0.1806862603261852 * x_n

cache_2434126157835956296_s = 
   0.08030500458941565

cache_2434126157835956296_p = 
   (0.007232292579356106 - (0.007232292579356106 * x_p)) / 0.04002679875215723

cache_2742393921021778050_n = 
   56.212339486148316

cache_2742393921021778050_s = 0.0
cache_2742393921021778050_p = 
   -56.212339486148316

cache_7037718863933443715_n = 0.3
cache_7037718863933443715_s = 1.0
cache_7037718863933443715_p = 0.3
cache_m9197988089807064229_n = (((Dx_n(((cache_5380421293040872985_n * (D_e(c_e_n(t, x_n) * 1000.0, 298.15))) * Dx_n(c_e_n(t, x_n))) + cache_2434126157835956296_n)) / -0.008035880643729006) + cache_2742393921021778050_n) / cache_7037718863933443715_n
cache_m9197988089807064229_s = (((Dx_s(((cache_5380421293040872985_s * (D_e(c_e_s(t, x_s) * 1000.0, 298.15))) * Dx_s(c_e_s(t, x_s))) + cache_2434126157835956296_s)) / -0.008035880643729006) + cache_2742393921021778050_s) / cache_7037718863933443715_s
cache_m9197988089807064229_p = (((Dx_p(((cache_5380421293040872985_p * (D_e(c_e_p(t, x_p) * 1000.0, 298.15))) * Dx_p(c_e_p(t, x_p))) + cache_2434126157835956296_p)) / -0.008035880643729006) + cache_2742393921021778050_p) / cache_7037718863933443715_p
eqs = [
   Dt(c_s_n_xav(t, r_n)) ~ cache_m4119150365269624523,
   Dt(c_s_p_xav(t, r_p)) ~ cache_5866437115522166420,
   Dt(c_e_n(t, x_n)) ~ cache_m9197988089807064229_n,
   Dt(c_e_s(t, x_s)) ~ cache_m9197988089807064229_s,
   Dt(c_e_p(t, x_p)) ~ cache_m9197988089807064229_p,
]


ics_bcs = [
   # initial conditions
   c_s_n_xav(0, r_n) ~ 0.8000000000000016,
   c_s_p_xav(0, r_p) ~ 0.6000000000000001,
   c_e_n(0, x_n) ~ 1.0,
   c_e_s(0, x_s) ~ 1.0,
   c_e_p(0, x_p) ~ 1.0,
   # boundary conditions
   Dr_n(c_s_n_xav(t, 0.0)) ~ 0.0,
   Dr_n(c_s_n_xav(t, 1.0)) ~ -0.14182855923368468,
   Dr_p(c_s_p_xav(t, 0.0)) ~ 0.0,
   Dr_p(c_s_p_xav(t, 1.0)) ~ 0.03237700710041634,
   Dx_n(c_e_n(t, 0.0)) ~ 0.0,
   c_e_n(t, 0.4444444444444445) ~ c_e_s(t, 0.4444444444444445),
   c_e_p(t, 0.5555555555555556) ~ c_e_s(t, 0.5555555555555556),
   Dx_p(c_e_p(t, 1.0)) ~ 0.0,
]

t_domain = Interval(0.000,0.159)
r_n_domain = Interval(0.0, 1.0)
r_p_domain = Interval(0.0, 1.0)
x_n_domain = Interval(0.0, 0.4444444444444445)
x_s_domain = Interval(0.4444444444444445, 0.5555555555555556)
x_p_domain = Interval(0.5555555555555556, 1.0)

domains = [
   t in t_domain,
   r_n in r_n_domain,
   r_p in r_p_domain,
   x_n in x_n_domain,
   x_s in x_s_domain,
   x_p in x_p_domain,
]

ind_vars = [t, r_n, r_p, x_n, x_s, x_p]
dep_vars = [c_s_n_xav(t, r_n), c_s_p_xav(t, r_p), c_e_n(t, x_n), c_e_s(t, x_s), c_e_p(t, x_p)]

@named SPMe_pde_system = PDESystem(eqs, ics_bcs, domains, ind_vars, dep_vars)
pde_system = SPMe_pde_system

nothing
end
