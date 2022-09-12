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
# 'Negative particle concentration' -> c_s_n
# 'Positive particle concentration' -> c_s_p
# 'Negative electrolyte concentration' -> c_e_n
# 'Separator electrolyte concentration' -> c_e_s
# 'Positive electrolyte concentration' -> c_e_p
# 'Negative electrode potential' -> phi_s_n
# 'Positive electrode potential' -> phi_s_p
# 'Negative electrolyte potential' -> phi_e_n
# 'Separator electrolyte potential' -> phi_e_s
# 'Positive electrolyte potential' -> phi_e_p
@variables c_s_n(..) c_s_p(..) c_e_n(..) c_e_s(..) c_e_p(..) phi_s_n(..) phi_s_p(..) phi_e_n(..) phi_e_s(..) phi_e_p(..)
dependent_variables_to_pybamm_names = Dict(
  :c_s_n => "Negative particle concentration",
  :c_s_p => "Positive particle concentration",
  :c_e_n => "Negative electrolyte concentration",
  :c_e_s => "Separator electrolyte concentration",
  :c_e_p => "Positive electrolyte concentration",
  :phi_s_n => "Negative electrode potential",
  :phi_s_p => "Positive electrode potential",
  :phi_e_n => "Negative electrolyte potential",
  :phi_e_s => "Separator electrolyte potential",
  :phi_e_p => "Positive electrolyte potential",
)
dependent_variables_to_dependencies = Dict(
  :c_s_n => (:t, :r_n, :x_n),
  :c_s_p => (:t, :r_p, :x_p),
  :c_e_n => (:t, :x_n),
  :c_e_s => (:t, :x_s),
  :c_e_p => (:t, :x_p),
  :phi_s_n => (:t, :x_n),
  :phi_s_p => (:t, :x_p),
  :phi_e_n => (:t, :x_n),
  :phi_e_s => (:t, :x_s),
  :phi_e_p => (:t, :x_p),
)
Dt = Differential(t)
Dr_n = Differential(r_n)
Dr_p = Differential(r_p)
Dx_n = Differential(x_n)
Dx_s = Differential(x_s)
Dx_p = Differential(x_p)

# 'Negative particle concentration' equation
cache_m3772969026142279169 = 8.813457647415216 * (1 / r_n^2 * Dr_n(r_n^2 * Dr_n(c_s_n(t, r_n, x_n))))

# 'Positive particle concentration' equation
cache_m4733341015708365011 = 22.598609352346717 * (1 / r_p^2 * Dr_p(r_p^2 * Dr_p(c_s_p(t, r_p, x_p))))

# 'Electrolyte concentration' equation

function D_e(c_e, T)
   5.34e-10 * exp(-0.00065 * c_e) * exp(14.941767670181717 - (4454.8880308646785 * 1.0 / T))
end

function kappa_e(c_e, T)
   (0.0911 + 0.0019100999999999999 * c_e - (1.052e-06 * (c_e ^ 2.0)) + 1.554e-10 * (c_e ^ 3.0)) * exp(13.9978223044089 - (4173.450720059513 * 1.0 / T))
end

function j0_n(c_e, c_s_surf, T)
   2e-05 * exp(15.119261670583445 - (4507.807867084453 * 1.0 / T)) * (c_e ^ 0.5) * (c_s_surf ^ 0.5) * ((24983.2619938437 - c_s_surf) ^ 0.5)
end

function U_n(sto)
   0.194 + 1.5 * exp(-120.0 * sto) + 0.0351 * tanh((sto - 0.286) / 0.083) - (0.0045 * tanh((sto - 0.849) / 0.119)) - (0.035 * tanh((sto - 0.9233) / 0.05)) - (0.0147 * tanh((sto - 0.5) / 0.034)) - (0.102 * tanh((sto - 0.194) / 0.142)) - (0.022 * tanh((sto - 0.9) / 0.0164)) - (0.011 * tanh((sto - 0.124) / 0.0226)) + 0.0155 * tanh((sto - 0.105) / 0.029)
end

function j0_p(c_e, c_s_surf, T)
   6e-07 * exp(15.962358172491644 - (4759.177089128383 * 1.0 / T)) * (c_e ^ 0.5) * (c_s_surf ^ 0.5) * ((51217.9257309275 - c_s_surf) ^ 0.5)
end

function U_p(sto)
   2.16216 + 0.07645 * tanh(30.834 - (57.858397200000006 * sto)) + 2.1581 * tanh(52.294 - (53.412228 * sto)) - (0.14169 * tanh(11.0923 - (21.0852666 * sto))) + 0.2051 * tanh(1.4684 - (5.829105600000001 * sto)) + 0.2531 * tanh((-1.062 * sto + 0.56478) / 0.1316) - (0.02167 * tanh(((1.062 * sto) - 0.525) / 0.006))
end

cache_m7615344741526561470_n = -589429731.3893579
cache_m7615344741526561470_s = -3587155110.512914
cache_m7615344741526561470_p = -589429731.3893579
cache_m8640978618337745505_n = 0.7818002858515763
cache_m8640978618337745505_s = 4.757885022498838
cache_m8640978618337745505_p = 0.7818002858515763
cache_8053574672103279471_n = 
   ((((2.0 * max(1.0,1e-08)) / 0.5925925925925926) * (j0_n(max(1e-08,c_e_n(t, x_n)) * 1000.0, (max(1e-08,min(c_s_n(t, 1.0, x_n),0.99999999)) * 24983.2619938437), 298.15))) * (sinh(0.5 * ((phi_s_n(t, x_n) - phi_e_n(t, x_n)) - (((U_n(c_s_n(t, 1.0, x_n)) + ((1e-06 * (1.0 / c_s_n(t, 1.0, x_n))) + (1e-06 * (1.0 / (c_s_n(t, 1.0, x_n) - 1.0))))) - 0.175189434028335) / 0.025692579121493725))))) / 0.04002679875215723

cache_8053574672103279471_s = 0.0
cache_8053574672103279471_p = 
   ((((2.0 * max(1.0,1e-08)) * 1.4062500000000002) * (j0_p(max(1e-08,c_e_p(t, x_p)) * 1000.0, (max(1e-08,min(c_s_p(t, 1.0, x_p),0.99999999)) * 51217.9257309275), 298.15))) * (sinh(0.5 * ((phi_s_p(t, x_p) - phi_e_p(t, x_p)) - (((U_p(c_s_p(t, 1.0, x_p)) + ((1e-06 * (1.0 / c_s_p(t, 1.0, x_p))) + (1e-06 * (1.0 / (c_s_p(t, 1.0, x_p) - 1.0))))) - 4.027013014008729) / 0.025692579121493725))))) / 0.04002679875215723

cache_m6677540497496108090_n = 0.3
cache_m6677540497496108090_s = 1.0
cache_m6677540497496108090_p = 0.3
cache_m8025144623980950_n = (((Dx_n(((cache_m7615344741526561470_n * (D_e(c_e_n(t, x_n) * 1000.0, 298.15))) * Dx_n(c_e_n(t, x_n))) + (0.08030500458941565 * (((kappa_e(c_e_n(t, x_n) * 1000.0, 298.15)) * cache_m8640978618337745505_n) * (((1.2 * Dx_n(c_e_n(t, x_n))) / c_e_n(t, x_n)) - Dx_n(phi_e_n(t, x_n))))))) / -0.008035880643729006) + cache_8053574672103279471_n) / cache_m6677540497496108090_n
cache_m8025144623980950_s = (((Dx_s(((cache_m7615344741526561470_s * (D_e(c_e_s(t, x_s) * 1000.0, 298.15))) * Dx_s(c_e_s(t, x_s))) + (0.08030500458941565 * (((kappa_e(c_e_s(t, x_s) * 1000.0, 298.15)) * cache_m8640978618337745505_s) * (((1.2 * Dx_s(c_e_s(t, x_s))) / c_e_s(t, x_s)) - Dx_s(phi_e_s(t, x_s))))))) / -0.008035880643729006) + cache_8053574672103279471_s) / cache_m6677540497496108090_s
cache_m8025144623980950_p = (((Dx_p(((cache_m7615344741526561470_p * (D_e(c_e_p(t, x_p) * 1000.0, 298.15))) * Dx_p(c_e_p(t, x_p))) + (0.08030500458941565 * (((kappa_e(c_e_p(t, x_p) * 1000.0, 298.15)) * cache_m8640978618337745505_p) * (((1.2 * Dx_p(c_e_p(t, x_p))) / c_e_p(t, x_p)) - Dx_p(phi_e_p(t, x_p))))))) / -0.008035880643729006) + cache_8053574672103279471_p) / cache_m6677540497496108090_p

# 'Negative electrode potential' equation
cache_m6533399985655703878 = (-Dx_n(221.12651346369245 * Dx_n(phi_s_n(t, x_n)))) + ((((2.0 * max(1.0,1e-08)) / 0.5925925925925926) * (j0_n(max(1e-08,c_e_n(t, x_n)) * 1000.0, (max(1e-08,min(c_s_n(t, 1.0, x_n),0.99999999)) * 24983.2619938437), 298.15))) * (sinh(0.5 * ((phi_s_n(t, x_n) - phi_e_n(t, x_n)) - (((U_n(c_s_n(t, 1.0, x_n)) + ((1e-06 * (1.0 / c_s_n(t, 1.0, x_n))) + (1e-06 * (1.0 / (c_s_n(t, 1.0, x_n) - 1.0))))) - 0.175189434028335) / 0.025692579121493725)))))

# 'Positive electrode potential' equation
cache_7285051224364462583 = (-Dx_p(16.821663817574187 * Dx_p(phi_s_p(t, x_p)))) + ((((2.0 * max(1.0,1e-08)) * 1.4062500000000002) * (j0_p(max(1e-08,c_e_p(t, x_p)) * 1000.0, (max(1e-08,min(c_s_p(t, 1.0, x_p),0.99999999)) * 51217.9257309275), 298.15))) * (sinh(0.5 * ((phi_s_p(t, x_p) - phi_e_p(t, x_p)) - (((U_p(c_s_p(t, 1.0, x_p)) + ((1e-06 * (1.0 / c_s_p(t, 1.0, x_p))) + (1e-06 * (1.0 / (c_s_p(t, 1.0, x_p) - 1.0))))) - 4.027013014008729) / 0.025692579121493725)))))

# 'Electrolyte potential' equation
cache_m8640978618337745505_n = 0.7818002858515763
cache_m8640978618337745505_s = 4.757885022498838
cache_m8640978618337745505_p = 0.7818002858515763
cache_m3518860117650756109_n = 
   (((2.0 * max(1.0,1e-08)) / 0.5925925925925926) * (j0_n(max(1e-08,c_e_n(t, x_n)) * 1000.0, (max(1e-08,min(c_s_n(t, 1.0, x_n),0.99999999)) * 24983.2619938437), 298.15))) * (sinh(0.5 * ((phi_s_n(t, x_n) - phi_e_n(t, x_n)) - (((U_n(c_s_n(t, 1.0, x_n)) + ((1e-06 * (1.0 / c_s_n(t, 1.0, x_n))) + (1e-06 * (1.0 / (c_s_n(t, 1.0, x_n) - 1.0))))) - 0.175189434028335) / 0.025692579121493725))))

cache_m3518860117650756109_s = 0.0
cache_m3518860117650756109_p = 
   (((2.0 * max(1.0,1e-08)) * 1.4062500000000002) * (j0_p(max(1e-08,c_e_p(t, x_p)) * 1000.0, (max(1e-08,min(c_s_p(t, 1.0, x_p),0.99999999)) * 51217.9257309275), 298.15))) * (sinh(0.5 * ((phi_s_p(t, x_p) - phi_e_p(t, x_p)) - (((U_p(c_s_p(t, 1.0, x_p)) + ((1e-06 * (1.0 / c_s_p(t, 1.0, x_p))) + (1e-06 * (1.0 / (c_s_p(t, 1.0, x_p) - 1.0))))) - 4.027013014008729) / 0.025692579121493725))))

cache_5541103136077085898_n = (Dx_n(((kappa_e(c_e_n(t, x_n) * 1000.0, 298.15)) * cache_m8640978618337745505_n) * (((1.2 * Dx_n(c_e_n(t, x_n))) / c_e_n(t, x_n)) - Dx_n(phi_e_n(t, x_n))))) - cache_m3518860117650756109_n
cache_5541103136077085898_s = (Dx_s(((kappa_e(c_e_s(t, x_s) * 1000.0, 298.15)) * cache_m8640978618337745505_s) * (((1.2 * Dx_s(c_e_s(t, x_s))) / c_e_s(t, x_s)) - Dx_s(phi_e_s(t, x_s))))) - cache_m3518860117650756109_s
cache_5541103136077085898_p = (Dx_p(((kappa_e(c_e_p(t, x_p) * 1000.0, 298.15)) * cache_m8640978618337745505_p) * (((1.2 * Dx_p(c_e_p(t, x_p))) / c_e_p(t, x_p)) - Dx_p(phi_e_p(t, x_p))))) - cache_m3518860117650756109_p

eq_fake = c_s_n(t, r_n, x_n) ~ c_s_n(t, r_n, x_n)
eqs = [
   Dt(c_s_n(t, r_n, x_n)) ~ cache_m3772969026142279169,
   Dt(c_s_p(t, r_p, x_p)) ~ cache_m4733341015708365011,
   Dt(c_e_n(t, x_n)) ~ cache_m8025144623980950_n,
   Dt(c_e_s(t, x_s)) ~ cache_m8025144623980950_s,
   Dt(c_e_p(t, x_p)) ~ cache_m8025144623980950_p,
   0 ~ cache_m6533399985655703878,
   0 ~ cache_7285051224364462583,
   0 ~ cache_5541103136077085898_n,
   0 ~ cache_5541103136077085898_s,
   0 ~ cache_5541103136077085898_p,
]

# 'Negative particle concentration' boundary condition
cache_4232793693923302623 = -0.06303491521497095 * ((((2.0 * max(1.0,1e-08)) / 0.5925925925925926) * (j0_n(max(1e-08,c_e_n(t, x_n)) * 1000.0, (max(1e-08,min(c_s_n(t, 1.0, x_n),0.99999999)) * 24983.2619938437), 298.15))) * (sinh(0.5 * ((phi_s_n(t, x_n) - phi_e_n(t, x_n)) - (((U_n(c_s_n(t, 1.0, x_n)) + ((1e-06 * (1.0 / c_s_n(t, 1.0, x_n))) + (1e-06 * (1.0 / (c_s_n(t, 1.0, x_n) - 1.0))))) - 0.175189434028335) / 0.025692579121493725)))))

# 'Positive particle concentration' boundary condition
cache_2857815746009979266 = -0.014389780933518372 * ((((2.0 * max(1.0,1e-08)) * 1.4062500000000002) * (j0_p(max(1e-08,c_e_p(t, x_p)) * 1000.0, (max(1e-08,min(c_s_p(t, 1.0, x_p),0.99999999)) * 51217.9257309275), 298.15))) * (sinh(0.5 * ((phi_s_p(t, x_p) - phi_e_p(t, x_p)) - (((U_p(c_s_p(t, 1.0, x_p)) + ((1e-06 * (1.0 / c_s_p(t, 1.0, x_p))) + (1e-06 * (1.0 / (c_s_p(t, 1.0, x_p) - 1.0))))) - 4.027013014008729) / 0.025692579121493725)))))

ic_bc_fake = c_s_n(0, r_n, x_n) ~ c_s_n(0, r_n, x_n)

ics_bcs = [
   # initial conditions
   c_s_n(0, r_n, x_n) ~ 0.8000000000000016,  # ~ 0
   c_s_p(0, r_p, x_p) ~ 0.6000000000000001,
   phi_s_n(0, x_n) ~ 0.0,  # 9.77 for this and below
   phi_s_p(0, x_p) ~ 0.0,
   c_e_n(0, x_n) ~ 1.0,  # ~ 0
   c_e_s(0, x_s) ~ 1.0,
   c_e_p(0, x_p) ~ 1.0,
   phi_e_n(0, x_n) ~ -0.0,
   phi_e_s(0, x_s) ~ -0.0,
   phi_e_p(0, x_p) ~ -0.0,
   # boundary conditions
   Dr_n(c_s_n(t, 0.0, x_n)) ~ 0.0,
   Dr_n(c_s_n(t, 1.0, x_n)) ~ cache_4232793693923302623,
   Dr_p(c_s_p(t, 0.0, x_p)) ~ 0.0,
   Dr_p(c_s_p(t, 1.0, x_p)) ~ cache_2857815746009979266,
   phi_s_n(t, 0.0) ~ 0.0,
   Dx_n(phi_s_n(t, 0.4444444444444445)) ~ 0.0,
   Dx_p(phi_s_p(t, 0.5555555555555556)) ~ 0.0,
   Dx_p(phi_s_p(t, 1.0)) ~ -0.05944715165186363,
   Dx_n(c_e_n(t, 0.0)) ~ 0.0,
   c_e_n(t, 0.4444444444444445) ~ c_e_s(t, 0.4444444444444445),
   c_e_p(t, 0.5555555555555556) ~ c_e_s(t, 0.5555555555555556),
   Dx_p(c_e_p(t, 1.0)) ~ 0.0,
   Dx_n(phi_e_n(t, 0.0)) ~ 0.0,
   phi_e_n(t, 0.4444444444444445) ~ phi_e_s(t, 0.4444444444444445),
   phi_e_p(t, 0.5555555555555556) ~ phi_e_s(t, 0.5555555555555556),
   Dx_p(phi_e_p(t, 1.0)) ~ 0.0,
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
dep_vars = [c_s_n(t, r_n, x_n), c_s_p(t, r_p, x_p), c_e_n(t, x_n), c_e_s(t, x_s), c_e_p(t, x_p), phi_s_n(t, x_n), phi_s_p(t, x_p), phi_e_n(t, x_n), phi_e_s(t, x_s), phi_e_p(t, x_p)]

@named DFN_pde_system = PDESystem(eqs, ics_bcs, domains, ind_vars, dep_vars)
pde_system = DFN_pde_system

nothing
end
