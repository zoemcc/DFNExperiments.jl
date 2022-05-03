begin
using IfElse
# ('negative electrode', 'separator', 'positive electrode') -> x
# ('negative electrode',) -> x_n
# ('positive electrode',) -> x_p
@parameters t x x_n x_p
independent_variables_to_pybamm_names = Dict(
  :t => "Time",
  :x => "negative electrode",
  :x_n => "negative electrode",
  :x_p => "positive electrode",
)
# 'Discharge capacity [A.h]' -> Q_Ah
# 'R-averaged negative particle concentration' -> c_s_n_rav
# 'R-averaged positive particle concentration' -> c_s_p_rav
# 'Electrolyte concentration' -> c_e
# 'Negative electrode potential' -> phi_s_n
# 'Positive electrode potential' -> phi_s_p
# 'Electrolyte potential' -> phi_e
@variables Q_Ah(..) c_s_n_rav(..) c_s_p_rav(..) c_e(..) phi_s_n(..) phi_s_p(..) phi_e(..)
dependent_variables_to_pybamm_names = Dict(
  :Q_Ah => "Discharge capacity [A.h]",
  :c_s_n_rav => "R-averaged negative particle concentration",
  :c_s_p_rav => "R-averaged positive particle concentration",
  :c_e => "Electrolyte concentration",
  :phi_s_n => "Negative electrode potential",
  :phi_s_p => "Positive electrode potential",
  :phi_e => "Electrolyte potential",
)
dependent_variables_to_dependencies = Dict(
  :Q_Ah => (:t,),
  :c_s_n_rav => (:t, :x_n),
  :c_s_p_rav => (:t, :x_p),
  :c_e => (:t, :x),
  :phi_s_n => (:t, :x_n),
  :phi_s_p => (:t, :x_p),
  :phi_e => (:t, :x),
)
Dt = Differential(t)
Dx = Differential(x)
Dx_n = Differential(x_n)
Dx_p = Differential(x_p)

# 'R-averaged negative particle concentration' equation

function j0_n(c_e, c_s_surf, T)
   2e-05 * exp(15.119261670583445 - (4507.807867084453 * 1.0 / T)) * (c_e ^ 0.5) * (c_s_surf ^ 0.5) * ((24983.2619938437 - c_s_surf) ^ 0.5)
end

function U_n(sto)
   0.194 + 1.5 * exp(-120.0 * sto) + 0.0351 * tanh((sto - 0.286) / 0.083) - (0.0045 * tanh((sto - 0.849) / 0.119)) - (0.035 * tanh((sto - 0.9233) / 0.05)) - (0.0147 * tanh((sto - 0.5) / 0.034)) - (0.102 * tanh((sto - 0.194) / 0.142)) - (0.022 * tanh((sto - 0.9) / 0.0164)) - (0.011 * tanh((sto - 0.124) / 0.0226)) + 0.0155 * tanh((sto - 0.105) / 0.029)
end

cache_m5586716010616745556 = -1.6666666666666667 * ((((2.0 * max(1.0,1e-08)) / 0.5925925925925926) * (j0_n(max(1e-08,c_e(t, x)) * 1000.0, (max(1e-08,min(c_s_n_rav(t, x_n),0.99999999)) * 24983.2619938437), 298.15))) * (sinh(0.5 * ((phi_s_n(t, x_n) - phi_e(t, x)) - (((U_n(c_s_n_rav(t, x_n)) + ((1e-06 * (1.0 / c_s_n_rav(t, x_n))) + (1e-06 * (1.0 / (c_s_n_rav(t, x_n) - 1.0))))) - 0.175189434028335) / 0.025692579121493725)))))

# 'R-averaged positive particle concentration' equation

function j0_p(c_e, c_s_surf, T)
   6e-07 * exp(15.962358172491644 - (4759.177089128383 * 1.0 / T)) * (c_e ^ 0.5) * (c_s_surf ^ 0.5) * ((51217.9257309275 - c_s_surf) ^ 0.5)
end

function U_p(sto)
   2.16216 + 0.07645 * tanh(30.834 - (57.858397200000006 * sto)) + 2.1581 * tanh(52.294 - (53.412228 * sto)) - (0.14169 * tanh(11.0923 - (21.0852666 * sto))) + 0.2051 * tanh(1.4684 - (5.829105600000001 * sto)) + 0.2531 * tanh((-1.062 * sto + 0.56478) / 0.1316) - (0.02167 * tanh(((1.062 * sto) - 0.525) / 0.006))
end

cache_m3435911063506985360 = -0.9755671139472863 * ((((2.0 * max(1.0,1e-08)) * 1.4062500000000002) * (j0_p(max(1e-08,c_e(t, x)) * 1000.0, (max(1e-08,min(c_s_p_rav(t, x_p),0.99999999)) * 51217.9257309275), 298.15))) * (sinh(0.5 * ((phi_s_p(t, x_p) - phi_e(t, x)) - (((U_p(c_s_p_rav(t, x_p)) + ((1e-06 * (1.0 / c_s_p_rav(t, x_p))) + (1e-06 * (1.0 / (c_s_p_rav(t, x_p) - 1.0))))) - 4.027013014008729) / 0.025692579121493725)))))

# 'Electrolyte concentration' equation

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

function kappa_e(c_e, T)
   (0.0911 + 0.0019100999999999999 * c_e - (1.052e-06 * (c_e ^ 2.0)) + 1.554e-10 * (c_e ^ 3.0)) * exp(13.9978223044089 - (4173.450720059513 * 1.0 / T))
end

cache_6137284500702632209 = concatenation(x, -589429731.3893579, -3587155110.512914, -589429731.3893579)
cache_7441480552996732851 = concatenation(x, 0.7818002858515763, 4.757885022498838, 0.7818002858515763)
cache_7972285907166862482 = concatenation(x, ((((2.0 * max(1.0,1e-08)) / 0.5925925925925926) * (j0_n(max(1e-08,c_e(t, x)) * 1000.0, (max(1e-08,min(c_s_n_rav(t, x_n),0.99999999)) * 24983.2619938437), 298.15))) * (sinh(0.5 * ((phi_s_n(t, x_n) - phi_e(t, x)) - (((U_n(c_s_n_rav(t, x_n)) + ((1e-06 * (1.0 / c_s_n_rav(t, x_n))) + (1e-06 * (1.0 / (c_s_n_rav(t, x_n) - 1.0))))) - 0.175189434028335) / 0.025692579121493725))))) / 0.04002679875215723, 0.0, 
   ((((2.0 * max(1.0,1e-08)) * 1.4062500000000002) * (j0_p(max(1e-08,c_e(t, x)) * 1000.0, (max(1e-08,min(c_s_p_rav(t, x_p),0.99999999)) * 51217.9257309275), 298.15))) * (sinh(0.5 * ((phi_s_p(t, x_p) - phi_e(t, x)) - (((U_p(c_s_p_rav(t, x_p)) + ((1e-06 * (1.0 / c_s_p_rav(t, x_p))) + (1e-06 * (1.0 / (c_s_p_rav(t, x_p) - 1.0))))) - 4.027013014008729) / 0.025692579121493725))))) / 0.04002679875215723
)
cache_m5058333278220336370 = concatenation(x, 0.3, 1.0, 0.3)
cache_m5020652664574118404 = (((Dx(((cache_6137284500702632209 * (D_e(c_e(t, x) * 1000.0, 298.15))) * Dx(c_e(t, x))) + (0.08030500458941565 * (((kappa_e(c_e(t, x) * 1000.0, 298.15)) * cache_7441480552996732851) * (((1.2 * Dx(c_e(t, x))) / c_e(t, x)) - Dx(phi_e(t, x))))))) / -0.008035880643729006) + cache_7972285907166862482) / cache_m5058333278220336370

# 'Negative electrode potential' equation
cache_m2639882052141794349 = (-Dx_n(221.12651346369245 * Dx_n(phi_s_n(t, x_n)))) + ((((2.0 * max(1.0,1e-08)) / 0.5925925925925926) * (j0_n(max(1e-08,c_e(t, x)) * 1000.0, (max(1e-08,min(c_s_n_rav(t, x_n),0.99999999)) * 24983.2619938437), 298.15))) * (sinh(0.5 * ((phi_s_n(t, x_n) - phi_e(t, x)) - (((U_n(c_s_n_rav(t, x_n)) + ((1e-06 * (1.0 / c_s_n_rav(t, x_n))) + (1e-06 * (1.0 / (c_s_n_rav(t, x_n) - 1.0))))) - 0.175189434028335) / 0.025692579121493725)))))

# 'Positive electrode potential' equation
cache_m7743774010292684701 = (-Dx_p(16.821663817574187 * Dx_p(phi_s_p(t, x_p)))) + ((((2.0 * max(1.0,1e-08)) * 1.4062500000000002) * (j0_p(max(1e-08,c_e(t, x)) * 1000.0, (max(1e-08,min(c_s_p_rav(t, x_p),0.99999999)) * 51217.9257309275), 298.15))) * (sinh(0.5 * ((phi_s_p(t, x_p) - phi_e(t, x)) - (((U_p(c_s_p_rav(t, x_p)) + ((1e-06 * (1.0 / c_s_p_rav(t, x_p))) + (1e-06 * (1.0 / (c_s_p_rav(t, x_p) - 1.0))))) - 4.027013014008729) / 0.025692579121493725)))))

# 'Electrolyte potential' equation
cache_7441480552996732851 = concatenation(x, 0.7818002858515763, 4.757885022498838, 0.7818002858515763)
cache_88799198024298508 = concatenation(x, (((2.0 * max(1.0,1e-08)) / 0.5925925925925926) * (j0_n(max(1e-08,c_e(t, x)) * 1000.0, (max(1e-08,min(c_s_n_rav(t, x_n),0.99999999)) * 24983.2619938437), 298.15))) * (sinh(0.5 * ((phi_s_n(t, x_n) - phi_e(t, x)) - (((U_n(c_s_n_rav(t, x_n)) + ((1e-06 * (1.0 / c_s_n_rav(t, x_n))) + (1e-06 * (1.0 / (c_s_n_rav(t, x_n) - 1.0))))) - 0.175189434028335) / 0.025692579121493725)))), 0.0, 
   (((2.0 * max(1.0,1e-08)) * 1.4062500000000002) * (j0_p(max(1e-08,c_e(t, x)) * 1000.0, (max(1e-08,min(c_s_p_rav(t, x_p),0.99999999)) * 51217.9257309275), 298.15))) * (sinh(0.5 * ((phi_s_p(t, x_p) - phi_e(t, x)) - (((U_p(c_s_p_rav(t, x_p)) + ((1e-06 * (1.0 / c_s_p_rav(t, x_p))) + (1e-06 * (1.0 / (c_s_p_rav(t, x_p) - 1.0))))) - 4.027013014008729) / 0.025692579121493725))))
)
cache_m7704575124354296774 = (Dx(((kappa_e(c_e(t, x) * 1000.0, 298.15)) * cache_7441480552996732851) * (((1.2 * Dx(c_e(t, x))) / c_e(t, x)) - Dx(phi_e(t, x))))) - cache_88799198024298508


eqs = [
   Dt(Q_Ah(t)) ~ 4.27249308415467,
   Dt(c_s_n_rav(t, x_n)) ~ cache_m5586716010616745556,
   Dt(c_s_p_rav(t, x_p)) ~ cache_m3435911063506985360,
   Dt(c_e(t, x)) ~ cache_m5020652664574118404,
   0 ~ cache_m2639882052141794349,
   0 ~ cache_m7743774010292684701,
   0 ~ cache_m7704575124354296774,
]


ics_bcs = [
   # initial conditions
   Q_Ah(0) ~ 0.0,
   c_s_n_rav(0, x_n) ~ 0.8000000000000016,
   c_s_p_rav(0, x_p) ~ 0.6000000000000001,
   phi_s_n(0, x_n) ~ 0.0,
   phi_s_p(0, x_p) ~ 0.0,
   c_e(0, x) ~ 1.0,
   phi_e(0, x) ~ -0.0,
   # boundary conditions
   phi_s_n(t, 0.0) ~ 0.0,
   Dx_n(phi_s_n(t, 0.4444444444444445)) ~ 0.0,
   Dx_p(phi_s_p(t, 0.5555555555555556)) ~ 0.0,
   Dx_p(phi_s_p(t, 1.0)) ~ -0.05944715165186363,
   Dx(c_e(t, 0.0)) ~ 0.0,
   Dx(c_e(t, 1.0)) ~ 0.0,
   Dx(phi_e(t, 0.0)) ~ 0.0,
   Dx(phi_e(t, 1.0)) ~ 0.0,
]

t_domain = Interval(0.000,0.159)
x_domain = Interval(0.0, 1.0)
x_n_domain = Interval(0.0, 0.4444444444444445)
x_p_domain = Interval(0.5555555555555556, 1.0)

domains = [
   t in t_domain,
   x in x_domain,
   x_n in x_n_domain,
   x_p in x_p_domain,
]

subdomain_relations = Dict(
   :x_n => :x,
   :x_p => :x,
)

eqs_integration_domains = [
   [:t,],
   [:t, :x_n,],
   [:t, :x_p,],
   [:t, :x,],
   [:t, :x_n,],
   [:t, :x_p,],
   [:t, :x,],
]

ics_bcs_integration_domains = [
   [:t => 0.0,],
   [:t => 0.0, :x_n,],
   [:t => 0.0, :x_p,],
   [:t => 0.0, :x_n,],
   [:t => 0.0, :x_p,],
   [:t => 0.0, :x,],
   [:t => 0.0, :x,],
   [:t, :x_n => 0.0,],
   [:t, :x_n => 0.4444444444444445,],
   [:t, :x_p => 0.5555555555555556,],
   [:t, :x_p => 1.0,],
   [:t, :x => 0.0,],
   [:t, :x => 1.0,],
   [:t, :x => 0.0,],
   [:t, :x => 1.0,],
]


ind_vars = [t, x, x_n, x_p]
dep_vars = [Q_Ah(t), c_s_n_rav(t, x_n), c_s_p_rav(t, x_p), c_e(t, x), phi_s_n(t, x_n), phi_s_p(t, x_p), phi_e(t, x)]

@named DFN_no_r_pde_system = PDESystem(eqs, ics_bcs, domains, ind_vars, dep_vars)
pde_system = DFN_no_r_pde_system

nothing
end
