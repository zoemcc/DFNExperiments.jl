#begin
#using DrWatson
#import DrWatson.savename
#end
#@quickactivate "NeuralPDEWatson"
begin
push!(LOAD_PATH, "/home/zobot/.julia/dev/NeuralPDE.jl/src")
using Revise
using LinearAlgebra
using IfElse
using PyCall
using Flux
println("NNPDE_tests_heterogeneous")
using DiffEqFlux
println("Starting Soon!")
using ModelingToolkit
using DiffEqBase
using Test, NeuralPDE
using GalacticOptim
using Optim
using Quadrature,Cubature, Cuba
using QuasiMonteCarlo
using SciMLBase
using DelimitedFiles
using CSV
using OrdinaryDiffEq
using Plots
using DataFrames
using LineSearches
using Zygote

using Random
end

begin
@pyimport pickle
@pyimport numpy

function myunpickle(filename)
    r = nothing
    @pywith pybuiltin("open")(filename,"rb") as f begin
        r = pickle.load(f)
    end
    return r
end
const rampupr = 0.05
const righttargetrn = 0.8000000000000016
const righttargetrp = 0.6000000000000001
const rightfluxrn = -0.14182855923368468
const rightfluxrp = 0.03237700710041634
function rampupcubic(r, targetflux, targetval)
    IfElse.ifelse(r < 1.0 .- rampupr, targetval, targetval .+ (r .- (1.0 .- rampupr)).^3 .* targetflux ./ (3 .* rampupr.^2) )
end
@register rampupcubic(r,f,v)
function rampupquadratic(r, targetflux, targetval)
    IfElse.ifelse(r < 1.0 .- rampupr, targetval, targetval .+ (r .- (1.0 .- rampupr)).^2 .* targetflux ./ (2rampupr) )
end
@register rampupquadratic(r,f,v)
function concatenation(x, n, s, p)
# A concatenation in the electrolyte domain
IfElse.ifelse(
    x < 0.4444444444444445, n, IfElse.ifelse(
        x < 0.5555555555555556, s, p
    )
)
end
@register concatenation(x, n, s, p)
end

#function main_train()
begin
    # PyBaMM solution
    pybamm_sols = myunpickle("/home/zobot/.julia/dev/DFN.jl/pybamm/hardcoded_models/MTK_format/pybamm_solutions/DFN.pickle")

    t_pb = pybamm_sols["Time"]
    x_pb = pybamm_sols["x"][:,1]
    x_n_pb = pybamm_sols["x_n"][:,1]
    x_p_pb = pybamm_sols["x_p"][:,1]
    r_n_pb = pybamm_sols["r_n"][:, 1, 1]
    r_p_pb = pybamm_sols["r_p"][:, 1, 1]
    c_s_n_pb = pybamm_sols["Negative particle concentration"]
    c_s_p_pb = pybamm_sols["Positive particle concentration"]
    c_e_pb = pybamm_sols["Electrolyte concentration"]
    phi_e_pb = pybamm_sols["Electrolyte potential"]
    phi_s_n_pb = pybamm_sols["Negative electrode potential"]
    phi_s_p_pb = pybamm_sols["Positive electrode potential"]

    begin
    # ('negative electrode',) -> xn
    # ('positive particle',) -> rp
    # ('positive electrode',) -> xp
    # ('negative electrode', 'separator', 'positive electrode') -> x
    # ('negative particle',) -> rn
    @parameters t xn rp xp x rn
    # 'Discharge capacity [A.h]' -> Q
    # 'Negative particle concentration' -> c_s_n
    # 'Positive particle concentration' -> c_s_p
    # 'Porosity times concentration' -> eps_c_e
    # 'Negative electrode potential' -> phi_s_n
    # 'Positive electrode potential' -> phi_s_p
    # 'Electrolyte potential' -> phi_e
    @variables Q(..) c_s_n(..) c_s_p(..) eps_c_e(..) phi_s_n(..) phi_s_p(..) phi_e(..)
    Dt = Differential(t)
    Dxn = Differential(xn)
    Drp = Differential(rp)
    Dxp = Differential(xp)
    Dx = Differential(x)
    Drn = Differential(rn)

    # 'Negative particle concentration' equation
    #cache_m2473104508719688576 = 8.813457647415216 * (1 / rn^2 * Drn(rn^2 * Drn(c_s_n(t, rn, xn))))
    cache_m2473104508719688576 = 8.813457647415216 * (2 / rn * Drn(c_s_n(t, rn, xn)) + Drn(Drn(c_s_n(t, rn, xn))))  #simplified

    # 'Positive particle concentration' equation
    #cache_1772252443521738373 = 22.598609352346717 * (1 / rp^2 * Drp(rp^2 * Drp(c_s_p(t, rp, xp))))
    cache_1772252443521738373 = 22.598609352346717 * (2 / rp * Drp(c_s_p(t, rp, xp)) + Drp(Drp(c_s_p(t, rp, xp))))  #simplified

    # 'Porosity times concentration' equation


    cache_2496679329654508972_f = concatenation(x, -0.31475547656191716, -1.9155408290138962, -0.31475547656191716)
    cache_2496679329654508972_n = -0.31475547656191716
    cache_2496679329654508972_s = -1.9155408290138962
    cache_2496679329654508972_p = -0.31475547656191716
    cache_m7835415840857063672_f = concatenation(x, -2.166666666666667 * eps_c_e(t, x), 
    -0.65 * eps_c_e(t, x)
    , 
    -2.166666666666667 * eps_c_e(t, x)
    )
    cache_m7835415840857063672_n = -2.166666666666667 * eps_c_e(t, x)
    cache_m7835415840857063672_s = -0.65 * eps_c_e(t, x)
    cache_m7835415840857063672_p = -2.166666666666667 * eps_c_e(t, x)
    cache_187754896691324381_f = concatenation(x, eps_c_e(t, x) / 0.3, eps_c_e(t, x), 
    eps_c_e(t, x) / 0.3
    )
    cache_187754896691324381_n = eps_c_e(t, x) / 0.3
    cache_187754896691324381_s = eps_c_e(t, x) 
    cache_187754896691324381_p = eps_c_e(t, x) / 0.3
    # 'Electrolyte potential' equation
    cache_m3686922780758088698_f = concatenation(x, 0.0911 + (6.367 * eps_c_e(t, x)), 
    0.0911 + (1.9101 * eps_c_e(t, x))
    , 
    0.0911 + (6.367 * eps_c_e(t, x))
    )
    cache_m3686922780758088698_n = 0.0911 + (6.367 * eps_c_e(t, x)) 
    cache_m3686922780758088698_s = 0.0911 + (1.9101 * eps_c_e(t, x))
    cache_m3686922780758088698_p = 0.0911 + (6.367 * eps_c_e(t, x))
    cache_m351415169828981746_f = concatenation(x, 11.68888888888889 * (eps_c_e(t, x) ^ 2.0), 
    1.052 * (eps_c_e(t, x) ^ 2.0)
    , 
    11.68888888888889 * (eps_c_e(t, x) ^ 2.0)
    )
    cache_m351415169828981746_n = 11.68888888888889 * (eps_c_e(t, x) ^ 2.0)
    cache_m351415169828981746_s = 1.052 * (eps_c_e(t, x) ^ 2.0)
    cache_m351415169828981746_p = 11.68888888888889 * (eps_c_e(t, x) ^ 2.0)
    cache_4937909445635284222_f = concatenation(x, 5.755555555555557 * (eps_c_e(t, x) ^ 3.0), 
    0.1554 * (eps_c_e(t, x) ^ 3.0)
    , 
    5.755555555555557 * (eps_c_e(t, x) ^ 3.0)
    )
    cache_4937909445635284222_n = 5.755555555555557 * (eps_c_e(t, x) ^ 3.0) 
    cache_4937909445635284222_s = 0.1554 * (eps_c_e(t, x) ^ 3.0)
    cache_4937909445635284222_p = 5.755555555555557 * (eps_c_e(t, x) ^ 3.0)
    cache_m6635399386264426189_f = concatenation(x, 0.7818002858515763, 4.757885022498838, 0.7818002858515763)
    cache_m6635399386264426189_n = 0.7818002858515763
    cache_m6635399386264426189_s = 4.757885022498838
    cache_m6635399386264426189_p = 0.7818002858515763
    cache_m891205910910975400_f = concatenation(x, ((3.3750000000000004 * (((0.0006324555320336759 * (max(1e-08,eps_c_e(t, x) / 0.3) ^ 0.5)) * ((max(1e-08,min(c_s_n(t, 1.0, x),0.99999999)) ^ 0.5) * 158.06094392304414)) * (24983.2619938437 - (max(1e-08,min(c_s_n(t, 1.0, x),0.99999999)) * 24983.2619938437) ^ 0.5))) * (sinh(0.5 * ((phi_s_n(t, x) - phi_e(t, x)) - ((((((((((((0.194 + (1.5 * (exp(-120.0 * c_s_n(t, 1.0, x))))) + (0.0351 * (tanh((c_s_n(t, 1.0, x) - 0.286) / 0.083)))) - (0.0045 * (tanh((c_s_n(t, 1.0, x) - 0.849) / 0.119)))) - (0.035 * (tanh((c_s_n(t, 1.0, x) - 0.9233) / 0.05)))) - (0.0147 * (tanh((c_s_n(t, 1.0, x) - 0.5) / 0.034)))) - (0.102 * (tanh((c_s_n(t, 1.0, x) - 0.194) / 0.142)))) - (0.022 * (tanh((c_s_n(t, 1.0, x) - 0.9) / 0.0164)))) - (0.011 * (tanh((c_s_n(t, 1.0, x) - 0.124) / 0.0226)))) + (0.0155 * (tanh((c_s_n(t, 1.0, x) - 0.105) / 0.029)))) + ((1e-06 * (1.0 / c_s_n(t, 1.0, x))) + (1e-06 * (1.0 / (c_s_n(t, 1.0, x) - 1.0))))) - 0.175189434028335) / 0.025692579121493725))))) / 0.04002679875215723, 0.0, 
    ((2.8125000000000004 * (((1.8973665961010275e-05 * (max(1e-08,eps_c_e(t, x) / 0.3) ^ 0.5)) * ((max(1e-08,min(c_s_p(t, 1.0, x),0.99999999)) ^ 0.5) * 226.313777156689)) * (51217.9257309275 - (max(1e-08,min(c_s_p(t, 1.0, x),0.99999999)) * 51217.9257309275) ^ 0.5))) * (sinh(0.5 * ((phi_s_p(t, x) - phi_e(t, x)) - (((((((((2.16216 + (0.07645 * (tanh(30.834 - (57.858397200000006 * c_s_p(t, 1.0, x)))))) + (2.1581 * (tanh(52.294 - (53.412228 * c_s_p(t, 1.0, x)))))) - (0.14169 * (tanh(11.0923 - (21.0852666 * c_s_p(t, 1.0, x)))))) + (0.2051 * (tanh(1.4684 - (5.829105600000001 * c_s_p(t, 1.0, x)))))) + (0.2531 * (tanh(((-1.062 * c_s_p(t, 1.0, x)) + 0.56478) / 0.1316)))) - (0.02167 * (tanh(((1.062 * c_s_p(t, 1.0, x)) - 0.525) / 0.006)))) + ((1e-06 * (1.0 / c_s_p(t, 1.0, x))) + (1e-06 * (1.0 / (c_s_p(t, 1.0, x) - 1.0))))) - 4.027013014008729) / 0.025692579121493725))))) / 0.04002679875215723
    )
    cache_m891205910910975400_n = ((3.3750000000000004 * (((0.0006324555320336759 * (max(1e-08,eps_c_e(t, x) / 0.3) ^ 0.5)) * ((max(1e-08,min(c_s_n(t, 1.0, x),0.99999999)) ^ 0.5) * 158.06094392304414)) * (24983.2619938437 - (max(1e-08,min(c_s_n(t, 1.0, x),0.99999999)) * 24983.2619938437) ^ 0.5))) * (sinh(0.5 * ((phi_s_n(t, x) - phi_e(t, x)) - ((((((((((((0.194 + (1.5 * (exp(-120.0 * c_s_n(t, 1.0, x))))) + (0.0351 * (tanh((c_s_n(t, 1.0, x) - 0.286) / 0.083)))) - (0.0045 * (tanh((c_s_n(t, 1.0, x) - 0.849) / 0.119)))) - (0.035 * (tanh((c_s_n(t, 1.0, x) - 0.9233) / 0.05)))) - (0.0147 * (tanh((c_s_n(t, 1.0, x) - 0.5) / 0.034)))) - (0.102 * (tanh((c_s_n(t, 1.0, x) - 0.194) / 0.142)))) - (0.022 * (tanh((c_s_n(t, 1.0, x) - 0.9) / 0.0164)))) - (0.011 * (tanh((c_s_n(t, 1.0, x) - 0.124) / 0.0226)))) + (0.0155 * (tanh((c_s_n(t, 1.0, x) - 0.105) / 0.029)))) + ((1e-06 * (1.0 / c_s_n(t, 1.0, x))) + (1e-06 * (1.0 / (c_s_n(t, 1.0, x) - 1.0))))) - 0.175189434028335) / 0.025692579121493725))))) / 0.04002679875215723
    cache_m891205910910975400_s = 0.0
    cache_m891205910910975400_p = ((2.8125000000000004 * (((1.8973665961010275e-05 * (max(1e-08,eps_c_e(t, x) / 0.3) ^ 0.5)) * ((max(1e-08,min(c_s_p(t, 1.0, x),0.99999999)) ^ 0.5) * 226.313777156689)) * (51217.9257309275 - (max(1e-08,min(c_s_p(t, 1.0, x),0.99999999)) * 51217.9257309275) ^ 0.5))) * (sinh(0.5 * ((phi_s_p(t, x) - phi_e(t, x)) - (((((((((2.16216 + (0.07645 * (tanh(30.834 - (57.858397200000006 * c_s_p(t, 1.0, x)))))) + (2.1581 * (tanh(52.294 - (53.412228 * c_s_p(t, 1.0, x)))))) - (0.14169 * (tanh(11.0923 - (21.0852666 * c_s_p(t, 1.0, x)))))) + (0.2051 * (tanh(1.4684 - (5.829105600000001 * c_s_p(t, 1.0, x)))))) + (0.2531 * (tanh(((-1.062 * c_s_p(t, 1.0, x)) + 0.56478) / 0.1316)))) - (0.02167 * (tanh(((1.062 * c_s_p(t, 1.0, x)) - 0.525) / 0.006)))) + ((1e-06 * (1.0 / c_s_p(t, 1.0, x))) + (1e-06 * (1.0 / (c_s_p(t, 1.0, x) - 1.0))))) - 4.027013014008729) / 0.025692579121493725))))) / 0.04002679875215723
    cache_3233777004814567773_f = ((Dx(((cache_2496679329654508972_f * exp(cache_m7835415840857063672_f)) * Dx(cache_187754896691324381_f)) + (0.08030500458941565 * ((((cache_m3686922780758088698_f - cache_m351415169828981746_f) + cache_4937909445635284222_f) * cache_m6635399386264426189_f) * (((1.2 * Dx(cache_187754896691324381_f)) / cache_187754896691324381_f) - Dx(phi_e(t, x))))))) / -0.008035880643729006) + cache_m891205910910975400_f
    cache_3233777004814567773_n = ((Dx(((cache_2496679329654508972_n * exp(cache_m7835415840857063672_n)) * Dx(cache_187754896691324381_n)) + (0.08030500458941565 * ((((cache_m3686922780758088698_n - cache_m351415169828981746_n) + cache_4937909445635284222_n) * cache_m6635399386264426189_n) * (((1.2 * Dx(cache_187754896691324381_n)) / cache_187754896691324381_n) - Dx(phi_e(t, x))))))) / -0.008035880643729006) + cache_m891205910910975400_n
    cache_3233777004814567773_s = ((Dx(((cache_2496679329654508972_s * exp(cache_m7835415840857063672_s)) * Dx(cache_187754896691324381_s)) + (0.08030500458941565 * ((((cache_m3686922780758088698_s - cache_m351415169828981746_s) + cache_4937909445635284222_s) * cache_m6635399386264426189_s) * (((1.2 * Dx(cache_187754896691324381_s)) / cache_187754896691324381_s) - Dx(phi_e(t, x))))))) / -0.008035880643729006) + cache_m891205910910975400_s
    cache_3233777004814567773_p = ((Dx(((cache_2496679329654508972_p * exp(cache_m7835415840857063672_p)) * Dx(cache_187754896691324381_p)) + (0.08030500458941565 * ((((cache_m3686922780758088698_p - cache_m351415169828981746_p) + cache_4937909445635284222_p) * cache_m6635399386264426189_p) * (((1.2 * Dx(cache_187754896691324381_p)) / cache_187754896691324381_p) - Dx(phi_e(t, x))))))) / -0.008035880643729006) + cache_m891205910910975400_p
    cache_3233777004814567773_concat = concatenation(x, cache_3233777004814567773_n, cache_3233777004814567773_s, cache_3233777004814567773_p)

    # 'Negative electrode potential' equation
    cache_m3101737402143440636 = (Dxn(-221.12651346369245 * Dxn(phi_s_n(t, xn)))) + ((3.3750000000000004 * (((0.0006324555320336759 * (max(1e-08,eps_c_e(t, xn) / 0.3) ^ 0.5)) * ((max(1e-08,min(c_s_n(t, 1.0, xn),0.99999999)) ^ 0.5) * 158.06094392304414)) * (24983.2619938437 - (max(1e-08,min(c_s_n(t, 1.0, xn),0.99999999)) * 24983.2619938437) ^ 0.5))) * (sinh(0.5 * ((phi_s_n(t, xn) - phi_e(t, xn)) - ((((((((((((0.194 + (1.5 * (exp(-120.0 * c_s_n(t, 1.0, xn))))) + (0.0351 * (tanh((c_s_n(t, 1.0, xn) - 0.286) / 0.083)))) - (0.0045 * (tanh((c_s_n(t, 1.0, xn) - 0.849) / 0.119)))) - (0.035 * (tanh((c_s_n(t, 1.0, xn) - 0.9233) / 0.05)))) - (0.0147 * (tanh((c_s_n(t, 1.0, xn) - 0.5) / 0.034)))) - (0.102 * (tanh((c_s_n(t, 1.0, xn) - 0.194) / 0.142)))) - (0.022 * (tanh((c_s_n(t, 1.0, xn) - 0.9) / 0.0164)))) - (0.011 * (tanh((c_s_n(t, 1.0, xn) - 0.124) / 0.0226)))) + (0.0155 * (tanh((c_s_n(t, 1.0, xn) - 0.105) / 0.029)))) + ((1e-06 * (1.0 / c_s_n(t, 1.0, xn))) + (1e-06 * (1.0 / (c_s_n(t, 1.0, xn) - 1.0))))) - 0.175189434028335) / 0.025692579121493725)))))

    # 'Positive electrode potential' equation
    cache_2225606857123748425 = (Dxp(-16.82166381757419 * Dxp(phi_s_p(t, xp)))) + ((2.8125000000000004 * (((1.8973665961010275e-05 * (max(1e-08,eps_c_e(t, xp) / 0.3) ^ 0.5)) * ((max(1e-08,min(c_s_p(t, 1.0, xp),0.99999999)) ^ 0.5) * 226.313777156689)) * (51217.9257309275 - (max(1e-08,min(c_s_p(t, 1.0, xp),0.99999999)) * 51217.9257309275) ^ 0.5))) * (sinh(0.5 * ((phi_s_p(t, xp) - phi_e(t, xp)) - (((((((((2.16216 + (0.07645 * (tanh(30.834 - (57.858397200000006 * c_s_p(t, 1.0, xp)))))) + (2.1581 * (tanh(52.294 - (53.412228 * c_s_p(t, 1.0, xp)))))) - (0.14169 * (tanh(11.0923 - (21.0852666 * c_s_p(t, 1.0, xp)))))) + (0.2051 * (tanh(1.4684 - (5.829105600000001 * c_s_p(t, 1.0, xp)))))) + (0.2531 * (tanh(((-1.062 * c_s_p(t, 1.0, xp)) + 0.56478) / 0.1316)))) - (0.02167 * (tanh(((1.062 * c_s_p(t, 1.0, xp)) - 0.525) / 0.006)))) + ((1e-06 * (1.0 / c_s_p(t, 1.0, xp))) + (1e-06 * (1.0 / (c_s_p(t, 1.0, xp) - 1.0))))) - 4.027013014008729) / 0.025692579121493725)))))

    cache_187754896691324381_f = concatenation(x, eps_c_e(t, x) / 0.3, eps_c_e(t, x), 
    eps_c_e(t, x) / 0.3
    )
    cache_187754896691324381_n = eps_c_e(t, x) / 0.3
    cache_187754896691324381_s = eps_c_e(t, x) 
    cache_187754896691324381_p = eps_c_e(t, x) / 0.3
    cache_7504237941213128410_f = concatenation(x, (3.3750000000000004 * (((0.0006324555320336759 * (max(1e-08,eps_c_e(t, x) / 0.3) ^ 0.5)) * ((max(1e-08,min(c_s_n(t, 1.0, x),0.99999999)) ^ 0.5) * 158.06094392304414)) * (24983.2619938437 - (max(1e-08,min(c_s_n(t, 1.0, x),0.99999999)) * 24983.2619938437) ^ 0.5))) * (sinh(0.5 * ((phi_s_n(t, x) - phi_e(t, x)) - ((((((((((((0.194 + (1.5 * (exp(-120.0 * c_s_n(t, 1.0, x))))) + (0.0351 * (tanh((c_s_n(t, 1.0, x) - 0.286) / 0.083)))) - (0.0045 * (tanh((c_s_n(t, 1.0, x) - 0.849) / 0.119)))) - (0.035 * (tanh((c_s_n(t, 1.0, x) - 0.9233) / 0.05)))) - (0.0147 * (tanh((c_s_n(t, 1.0, x) - 0.5) / 0.034)))) - (0.102 * (tanh((c_s_n(t, 1.0, x) - 0.194) / 0.142)))) - (0.022 * (tanh((c_s_n(t, 1.0, x) - 0.9) / 0.0164)))) - (0.011 * (tanh((c_s_n(t, 1.0, x) - 0.124) / 0.0226)))) + (0.0155 * (tanh((c_s_n(t, 1.0, x) - 0.105) / 0.029)))) + ((1e-06 * (1.0 / c_s_n(t, 1.0, x))) + (1e-06 * (1.0 / (c_s_n(t, 1.0, x) - 1.0))))) - 0.175189434028335) / 0.025692579121493725)))), 0.0, 
    (2.8125000000000004 * (((1.8973665961010275e-05 * (max(1e-08,eps_c_e(t, x) / 0.3) ^ 0.5)) * ((max(1e-08,min(c_s_p(t, 1.0, x),0.99999999)) ^ 0.5) * 226.313777156689)) * (51217.9257309275 - (max(1e-08,min(c_s_p(t, 1.0, x),0.99999999)) * 51217.9257309275) ^ 0.5))) * (sinh(0.5 * ((phi_s_p(t, x) - phi_e(t, x)) - (((((((((2.16216 + (0.07645 * (tanh(30.834 - (57.858397200000006 * c_s_p(t, 1.0, x)))))) + (2.1581 * (tanh(52.294 - (53.412228 * c_s_p(t, 1.0, x)))))) - (0.14169 * (tanh(11.0923 - (21.0852666 * c_s_p(t, 1.0, x)))))) + (0.2051 * (tanh(1.4684 - (5.829105600000001 * c_s_p(t, 1.0, x)))))) + (0.2531 * (tanh(((-1.062 * c_s_p(t, 1.0, x)) + 0.56478) / 0.1316)))) - (0.02167 * (tanh(((1.062 * c_s_p(t, 1.0, x)) - 0.525) / 0.006)))) + ((1e-06 * (1.0 / c_s_p(t, 1.0, x))) + (1e-06 * (1.0 / (c_s_p(t, 1.0, x) - 1.0))))) - 4.027013014008729) / 0.025692579121493725))))
    )
    cache_7504237941213128410_n = (3.3750000000000004 * (((0.0006324555320336759 * (max(1e-08,eps_c_e(t, x) / 0.3) ^ 0.5)) * ((max(1e-08,min(c_s_n(t, 1.0, x),0.99999999)) ^ 0.5) * 158.06094392304414)) * (24983.2619938437 - (max(1e-08,min(c_s_n(t, 1.0, x),0.99999999)) * 24983.2619938437) ^ 0.5))) * (sinh(0.5 * ((phi_s_n(t, x) - phi_e(t, x)) - ((((((((((((0.194 + (1.5 * (exp(-120.0 * c_s_n(t, 1.0, x))))) + (0.0351 * (tanh((c_s_n(t, 1.0, x) - 0.286) / 0.083)))) - (0.0045 * (tanh((c_s_n(t, 1.0, x) - 0.849) / 0.119)))) - (0.035 * (tanh((c_s_n(t, 1.0, x) - 0.9233) / 0.05)))) - (0.0147 * (tanh((c_s_n(t, 1.0, x) - 0.5) / 0.034)))) - (0.102 * (tanh((c_s_n(t, 1.0, x) - 0.194) / 0.142)))) - (0.022 * (tanh((c_s_n(t, 1.0, x) - 0.9) / 0.0164)))) - (0.011 * (tanh((c_s_n(t, 1.0, x) - 0.124) / 0.0226)))) + (0.0155 * (tanh((c_s_n(t, 1.0, x) - 0.105) / 0.029)))) + ((1e-06 * (1.0 / c_s_n(t, 1.0, x))) + (1e-06 * (1.0 / (c_s_n(t, 1.0, x) - 1.0))))) - 0.175189434028335) / 0.025692579121493725))))
    cache_7504237941213128410_s = 0.0 
    cache_7504237941213128410_p = (2.8125000000000004 * (((1.8973665961010275e-05 * (max(1e-08,eps_c_e(t, x) / 0.3) ^ 0.5)) * ((max(1e-08,min(c_s_p(t, 1.0, x),0.99999999)) ^ 0.5) * 226.313777156689)) * (51217.9257309275 - (max(1e-08,min(c_s_p(t, 1.0, x),0.99999999)) * 51217.9257309275) ^ 0.5))) * (sinh(0.5 * ((phi_s_p(t, x) - phi_e(t, x)) - (((((((((2.16216 + (0.07645 * (tanh(30.834 - (57.858397200000006 * c_s_p(t, 1.0, x)))))) + (2.1581 * (tanh(52.294 - (53.412228 * c_s_p(t, 1.0, x)))))) - (0.14169 * (tanh(11.0923 - (21.0852666 * c_s_p(t, 1.0, x)))))) + (0.2051 * (tanh(1.4684 - (5.829105600000001 * c_s_p(t, 1.0, x)))))) + (0.2531 * (tanh(((-1.062 * c_s_p(t, 1.0, x)) + 0.56478) / 0.1316)))) - (0.02167 * (tanh(((1.062 * c_s_p(t, 1.0, x)) - 0.525) / 0.006)))) + ((1e-06 * (1.0 / c_s_p(t, 1.0, x))) + (1e-06 * (1.0 / (c_s_p(t, 1.0, x) - 1.0))))) - 4.027013014008729) / 0.025692579121493725))))
    cache_5792085376191370772_f = (Dx((((cache_m3686922780758088698_f - cache_m351415169828981746_f) + cache_4937909445635284222_f) * cache_m6635399386264426189_f) * (((1.2 * Dx(cache_187754896691324381_f)) / cache_187754896691324381_f) - Dx(phi_e(t, x))))) - cache_7504237941213128410_f
    cache_5792085376191370772_n = (Dx((((cache_m3686922780758088698_n - cache_m351415169828981746_n) + cache_4937909445635284222_n) * cache_m6635399386264426189_n) * (((1.2 * Dx(cache_187754896691324381_n)) / cache_187754896691324381_n) - Dx(phi_e(t, x))))) - cache_7504237941213128410_n
    cache_5792085376191370772_s = (Dx((((cache_m3686922780758088698_s - cache_m351415169828981746_s) + cache_4937909445635284222_s) * cache_m6635399386264426189_s) * (((1.2 * Dx(cache_187754896691324381_s)) / cache_187754896691324381_s) - Dx(phi_e(t, x))))) - cache_7504237941213128410_s
    cache_5792085376191370772_p = (Dx((((cache_m3686922780758088698_p - cache_m351415169828981746_p) + cache_4937909445635284222_p) * cache_m6635399386264426189_p) * (((1.2 * Dx(cache_187754896691324381_p)) / cache_187754896691324381_p) - Dx(phi_e(t, x))))) - cache_7504237941213128410_p
    cache_5792085376191370772_concat = concatenation(x, cache_m3686922780758088698_n, cache_m3686922780758088698_s, cache_m3686922780758088698_p)


    eqs = [
    Dt(Q(t)) ~ 4.27249308415467,
    Dt(c_s_n(t, rn, xn)) ~ cache_m2473104508719688576, #only integrate xn
    Dt(c_s_p(t, rp, xp)) ~ cache_1772252443521738373, #only integrate xp
    Dt(eps_c_e(t, x)) ~ cache_3233777004814567773_concat, # all x
    0 ~ cache_m3101737402143440636, # should only be integrated along xn
    0 ~ cache_2225606857123748425, # should only be integrated along xp
    0 ~ cache_5792085376191370772_concat, # should be integrated all x
    ]
    # 'Porosity times concentration' initial condition

    cache_3549472971836658861 = concatenation(x, 0.3, 1.0, 0.3)

    # 'Negative particle concentration' boundary condition
    cache_m4573998876320545433 = -0.06303491521497095 * ((3.3750000000000004 * (((0.0006324555320336759 * (max(1e-08,eps_c_e(t, xn) / 0.3) ^ 0.5)) * ((max(1e-08,min(c_s_n(t, 1.0, xn),0.99999999)) ^ 0.5) * 158.06094392304414)) * (24983.2619938437 - (max(1e-08,min(c_s_n(t, 1.0, xn),0.99999999)) * 24983.2619938437) ^ 0.5))) * (sinh(0.5 * ((phi_s_n(t, xn) - phi_e(t, xn)) - ((((((((((((0.194 + (1.5 * (exp(-120.0 * c_s_n(t, 1.0, xn))))) + (0.0351 * (tanh((c_s_n(t, 1.0, xn) - 0.286) / 0.083)))) - (0.0045 * (tanh((c_s_n(t, 1.0, xn) - 0.849) / 0.119)))) - (0.035 * (tanh((c_s_n(t, 1.0, xn) - 0.9233) / 0.05)))) - (0.0147 * (tanh((c_s_n(t, 1.0, xn) - 0.5) / 0.034)))) - (0.102 * (tanh((c_s_n(t, 1.0, xn) - 0.194) / 0.142)))) - (0.022 * (tanh((c_s_n(t, 1.0, xn) - 0.9) / 0.0164)))) - (0.011 * (tanh((c_s_n(t, 1.0, xn) - 0.124) / 0.0226)))) + (0.0155 * (tanh((c_s_n(t, 1.0, xn) - 0.105) / 0.029)))) + ((1e-06 * (1.0 / c_s_n(t, 1.0, xn))) + (1e-06 * (1.0 / (c_s_n(t, 1.0, xn) - 1.0))))) - 0.175189434028335) / 0.025692579121493725)))))

    # 'Positive particle concentration' boundary condition
    cache_2508578313141055385 = -0.014389780933518372 * ((2.8125000000000004 * (((1.8973665961010275e-05 * (max(1e-08,eps_c_e(t, xp) / 0.3) ^ 0.5)) * ((max(1e-08,min(c_s_p(t, 1.0, xp),0.99999999)) ^ 0.5) * 226.313777156689)) * (51217.9257309275 - (max(1e-08,min(c_s_p(t, 1.0, xp),0.99999999)) * 51217.9257309275) ^ 0.5))) * (sinh(0.5 * ((phi_s_p(t, xp) - phi_e(t, xp)) - (((((((((2.16216 + (0.07645 * (tanh(30.834 - (57.858397200000006 * c_s_p(t, 1.0, xp)))))) + (2.1581 * (tanh(52.294 - (53.412228 * c_s_p(t, 1.0, xp)))))) - (0.14169 * (tanh(11.0923 - (21.0852666 * c_s_p(t, 1.0, xp)))))) + (0.2051 * (tanh(1.4684 - (5.829105600000001 * c_s_p(t, 1.0, xp)))))) + (0.2531 * (tanh(((-1.062 * c_s_p(t, 1.0, xp)) + 0.56478) / 0.1316)))) - (0.02167 * (tanh(((1.062 * c_s_p(t, 1.0, xp)) - 0.525) / 0.006)))) + ((1e-06 * (1.0 / c_s_p(t, 1.0, xp))) + (1e-06 * (1.0 / (c_s_p(t, 1.0, xp) - 1.0))))) - 4.027013014008729) / 0.025692579121493725)))))

    left_r_domain = 0.001

    ics_bcs = [
    # initial conditions
    Q(0) ~ 0.0,
    c_s_n(0, rn, xn) ~ 0.8000000000000016, # I changed these to have the xn/ xp
    c_s_p(0, rp, xp) ~ 0.6000000000000001,
    #phi_s_n(0, xn) ~ 0.0,  # phi eqs not helpful, actually hurtful later on to optimize to, since no initial guess needed for this
    #phi_s_p(0, xp) ~ -0.0,
    eps_c_e(0, x) ~ cache_3549472971836658861,
    #phi_e(0, x) ~ -0.0,
    # boundary conditions
    Drn(c_s_n(t, left_r_domain, xn)) ~ 0.0,
    Drn(c_s_n(t, 1.0, xn)) ~ cache_m4573998876320545433,
    Drp(c_s_p(t, left_r_domain, xp)) ~ 0.0,
    Drp(c_s_p(t, 1.0, xp)) ~ cache_2508578313141055385,
    phi_s_n(t, 0.0) ~ 0.0, 
    Dxn(phi_s_n(t, 0.4444444444444445)) ~ 0.0,
    Dxp(phi_s_p(t, 0.5555555555555556)) ~ 0.0,
    Dxp(phi_s_p(t, 1.0)) ~ -0.059447151651863615,
    Dx(phi_e(t, 0.0)) ~ 0.0,
    Dx(phi_e(t, 1.0)) ~ 0.0,
    ]

    t_domain = IntervalDomain(0.000,0.159)
    xn_domain = IntervalDomain(0.0, 0.4444444444444445)
    rp_domain = IntervalDomain(left_r_domain, 1.0)
    xp_domain = IntervalDomain(0.5555555555555556, 1.0)
    x_domain = IntervalDomain(0.0, 1.0)
    rn_domain = IntervalDomain(left_r_domain, 1.0)

    #something about xp ⊂ x and xn ⊂ x

    domains = [
    t in t_domain,
    xn in xn_domain,
    rp in rp_domain,
    xp in xp_domain,
    x in x_domain,
    rn in rn_domain,
    ]
    ind_vars = [t, xn, rp, xp, x, rn]
    dep_vars = [Q(t), c_s_n(t, rn, xn), c_s_p(t, rp, xp), eps_c_e(t, x), phi_s_n(t, xn), phi_s_p(t, xp), phi_e(t, x)]
    chain_input_dimensions = [[1], [1,6,2], [1,3,4], [1,5], [1,2], [1,4], [1,5]]

    DFN_pde_system = PDESystem(eqs, ics_bcs, domains, ind_vars, dep_vars)

    rampupr = 0.05
end
end
begin


    # solve for initial derivative at 0, end derivative at targetflux


    righttargetrn = 0.8000000000000016
    righttargetrp = 0.6000000000000001
    rightfluxrn = -0.14182855923368468
    rightfluxrp = 0.03237700710041634

    leftrdomain = 0.001




    ## PINN Part


    begin
    num_dim = 50
    nonlin = Flux.gelu
    #in_dim = 3
    #out_dim = 3
    strategy = NeuralPDE.QuadratureTraining(;quadrature_alg=HCubatureJL(),abstol=1e-3, reltol=1, maxiters=1000, batch=0)
    #strategy_ = NeuralPDE.QuadratureTraining(;abstol=1e-6, reltol=1e-8, maxiters=2000)
    #strategy = NeuralPDE.QuadratureTraining(;quadrature_alg=HCubatureJL(), batch=0)
    #strategy = NeuralPDE.StochasticTraining(128)
    #strategy = NeuralPDE.QuadratureTraining()
    #strategy = NeuralPDE.QuadratureTraining(;quadrature_alg=HCubatureJL(),
                                                        #reltol=1e-3,abstol=1e-5,
                                                        #maxiters =2000, batch=0)

    #in_dims = [3, 3, 3]
    #in_dims = [1, 2, 2]
    in_dims = [1, 3, 3, 2, 2, 2, 2]
    num_hid = 2
    chains_ = [FastChain(FastDense(in_dim,num_dim,nonlin),
                        [FastDense(num_dim,num_dim,nonlin) for i in 1:num_hid]...,
                        FastDense(num_dim,1)) for in_dim in in_dims]
    #adalosspoisson = NeuralPDE.LossGradientsAdaptiveLoss(20; α=0.9f0)
    adaloss = NeuralPDE.MiniMaxAdaptiveLoss(20; α_pde=1e-3, α_bc=1e-1, pde_weights_start=1e-6, bc_weights_start=1e3)
    data_dir = "/home/zobot/.julia/dev/NeuralPDEWatson/data"
    exp_name = "dfn_first_testexperimental"
    exp_folder = joinpath(data_dir, exp_name)
    wipe_logs = false
    add_loss = generate_supervised_loss()[1]
    discretization = NeuralPDE.PhysicsInformedNN(chains_,
                                                    strategy; adaptive_loss=adaloss, exp_name=exp_name, data_dir=data_dir, wipe=wipe_logs, additional_loss=add_loss)
end
end
begin
    sym_prob = NeuralPDE.symbolic_discretize(DFN_pde_system,discretization)
end
    prob = NeuralPDE.discretize(DFN_pde_system,discretization)
begin
                                                    
    end

    #@run sym_prob = NeuralPDE.symbolic_discretize(SPM_pde_system,discretization)
    #@run sym_prob = NeuralPDE.symbolic_discretize(SPM_pde_system,discretization)
    begin
    initθ = vcat(discretization.init_params...)
    #opt = Flux.Optimiser(ClipValue(1e-3), ExpDecay(1, 0.5, 25_000), ADAM(3e-4))
    opt = ADAM(3e-4)
    saveevery = 100
    loss = zeros(Float64, saveevery)

    losssavefile = joinpath(exp_folder, "loss.csv")

    if wipe_logs
        rm(losssavefile; force=true)
    end
    iteration_count_arr = [1]
    cb = function (p,l)
        iteration_count = iteration_count_arr[1]


        println("Current loss is: $l, iteration is: $(iteration_count)")
        loss[((iteration_count - 1) % saveevery) + 1] = l
        if iteration_count % saveevery == 0
            cursavefile = joinpath(exp_folder, string(iteration_count, base=10, pad=5) * ".csv")
            writedlm(cursavefile, p, ",")
            df = DataFrame(loss=loss)
            CSV.write(losssavefile, df, writeheader=false, append=true)
        end
        iteration_count_arr[1] += 1
        if isnan(l)
            println("Aborting training, NaN found")
            return true
        else
            return false 
        end
    end
    #prob.f(initθ, [])
    #prob_pretrained.f(supervised_res, [])
    prob_pretrained = remake(prob; u0=supervised_res)
    res = GalacticOptim.solve(prob_pretrained, ADAM(3e-4); cb = cb, maxiters=20_000)
end

function generate_supervised_loss()
    pybamm_sols = myunpickle("/home/zobot/.julia/dev/DFN.jl/pybamm/hardcoded_models/MTK_format/pybamm_solutions/DFN.pickle")
    t_pb = pybamm_sols["Time"]
    x_pb = pybamm_sols["x"][:,1]
    x_n_pb = pybamm_sols["x_n"][:,1]
    x_p_pb = pybamm_sols["x_p"][:,1]
    r_n_pb = pybamm_sols["r_n"][:, 1, 1]
    r_p_pb = pybamm_sols["r_p"][:, 1, 1]
    c_s_n_pb = pybamm_sols["Negative particle concentration"]
    c_s_p_pb = pybamm_sols["Positive particle concentration"]
    c_e_pb = pybamm_sols["Electrolyte concentration"]
    phi_e_pb = pybamm_sols["Electrolyte potential"]
    phi_s_n_pb = pybamm_sols["Negative electrode potential"]
    phi_s_p_pb = pybamm_sols["Positive electrode potential"]

    c_s_n_inputs = hcat(reshape([[t_i, r_n_i, x_n_i] for x_n_i in x_n_pb, r_n_i in r_n_pb, t_i in t_pb], (length(c_s_n_pb),))...)
    c_s_n_sup = reshape(c_s_n_pb, (1, length(c_s_n_pb),)) 
    c_s_p_inputs = hcat(reshape([[t_i, r_p_i, x_p_i] for x_p_i in x_p_pb, r_p_i in r_p_pb, t_i in t_pb], (length(c_s_p_pb),))...)
    c_s_p_sup = reshape(c_s_p_pb, (1, length(c_s_p_pb),)) 
    eps_c_e_inputs = hcat(reshape([[t_i, x_i] for x_i in x_pb, t_i in t_pb], (length(c_e_pb),))...)
    eps_c_e_sup = reshape(c_e_pb, (1, length(c_e_pb),)) 
    phi_s_n_inputs = hcat(reshape([[t_i, x_n_i] for x_n_i in x_n_pb, t_i in t_pb], (length(phi_s_n_pb),))...)
    phi_s_n_sup = reshape(phi_s_n_pb, (1, length(phi_s_n_pb),)) 
    phi_s_p_inputs = hcat(reshape([[t_i, x_p_i] for x_p_i in x_p_pb, t_i in t_pb], (length(phi_s_p_pb),))...)
    phi_s_p_sup = reshape(phi_s_p_pb, (1, length(phi_s_p_pb),)) 
    phi_e_inputs = hcat(reshape([[t_i, x_i] for x_i in x_pb, t_i in t_pb], (length(phi_e_pb),))...)
    phi_e_sup = reshape(phi_e_pb, (1, length(phi_e_pb),)) 

    param_lengths = (length ∘ initial_params).(chains_)
    indices_in_params = map(zip(param_lengths, cumsum(param_lengths))) do (param_length, cumsum_param)
            cumsum_param - (param_length - 1) : cumsum_param
    end

    #ind_vars = [t, xn, rp, xp, x, rn]
    #dep_vars = [Q(t), c_s_n(t, rn, xn), c_s_p(t, rp, xp), eps_c_e(t, x), phi_s_n(t, xn), phi_s_p(t, xp), phi_e(t, x)]
    #chain_input_dimensions = [[1], [1,6,2], [1,3,4], [1,5], [1,2], [1,4], [1,5]]
    function supervised_loss(θ)
        c_s_n_evals = chains_[2](c_s_n_inputs, θ[indices_in_params[2]])
        c_s_p_evals = chains_[3](c_s_p_inputs, θ[indices_in_params[3]])
        eps_c_e_evals = chains_[4](eps_c_e_inputs, θ[indices_in_params[4]])
        phi_s_n_evals = chains_[5](phi_s_n_inputs, θ[indices_in_params[5]])
        phi_s_p_evals = chains_[6](phi_s_p_inputs, θ[indices_in_params[6]])
        phi_e_evals = chains_[7](phi_e_inputs, θ[indices_in_params[7]])

        c_s_n_l2 = sum(abs2, c_s_n_evals .- c_s_n_sup) / length(c_s_n_pb)
        c_s_p_l2 = sum(abs2, c_s_p_evals .- c_s_p_sup) / length(c_s_p_pb)
        eps_c_e_l2 = sum(abs2, eps_c_e_evals .- eps_c_e_sup) / length(c_e_pb)
        phi_s_n_l2 = sum(abs2, phi_s_n_evals .- phi_s_n_sup) / length(phi_s_n_pb)
        phi_s_p_l2 = sum(abs2, phi_s_p_evals .- phi_s_p_sup) / length(phi_s_p_pb)
        phi_e_l2 = sum(abs2, phi_e_evals .- phi_e_sup) / length(phi_e_pb)
        c_s_n_l2 + c_s_p_l2 + eps_c_e_l2 + phi_s_n_l2 + phi_s_p_l2 + phi_e_l2
    end
    function supervised_additional_loss(phi, θ)
        supervised_loss(θ)
    end
    function supervised_loss_only(θ, p)
        supervised_loss(θ)
    end

    (supervised_additional_loss, supervised_loss_only)
end


f = OptimizationFunction(generate_supervised_loss()[2], GalacticOptim.AutoZygote())
supervised_opt_prob = GalacticOptim.OptimizationProblem(f, initθ)
supervisedcb = function (p,l)
    println("Current loss is: $l")
    return false
end
supervised_res = GalacticOptim.solve(supervised_opt_prob, ADAM(3e-4); cb = supervisedcb, maxiters=10_000)
#=
begin
    pretrained_params_file = joinpath(exp_folder, "20000.csv")
    pretrained_pinn_params = Array(CSV.read(pretrained_params_file, DataFrame; header=false))[:,1]
    params = pretrained_pinn_params

    dts = sols_t
    drns = sols_r_n
    drps = sols_r_p

    param_lengths = (length ∘ initial_params).(chains_)
    indices_in_params = map(zip(param_lengths, cumsum(param_lengths))) do (param_length, cumsum_param)
            cumsum_param - (param_length - 1) : cumsum_param
    end

    #phi_i_params = [params[indices_in_params_i] for indices_in_params_i in indices_in_params]

    #Q_evals = [chains_[1]([dt], phi_i_params[1])[1] for dt in sols_t]

    #c_s_n_evals_pretrained_pinn = [chains_[2]([dt, drn], phi_i_params[2])[1] for dt in dts, drn in drns]
    c_s_p_evals_pretrained_pinn = [chains_[1]([dt, drp], params)[1] for dt in dts, drp in drps]

    anim = @animate for i in 1:length(dts)
        #p1 = plot(drns,sols_c_s_n[:, i];ylims=(0,1.3),ylabel="c_s_n",xlabel="r_n",legend=true, label="FDM")
        #p1 = plot(drns,c_s_n_evals_pretrained_pinn[i, :];ylims=(0,1.3),ylabel="c_s_n",xlabel="r_n",legend=true, label="FDM")
        #p1 = plot!(p1,drns,c_s_n_evals_pretrained_pinn[i,:];ylims=(0,1.3),ylabel="c_s_n",xlabel="r_n",legend=true, label="Fulltrained pinn", linestyle=:dashdot)
        p2 = plot(drps,sols_c_s_p[:, i];ylims=(0,1.3),ylabel="c_s_p",xlabel="r_p",legend=true, label="FDM")
        #p2 = plot(drps,c_s_p_evals_pretrained_pinn[i, :];ylims=(0,1.3),ylabel="c_s_p",xlabel="r_p",legend=true, label="FDM")
        p2 = plot!(p2,drps,c_s_p_evals_pretrained_pinn[i,:];ylims=(0,1.3),ylabel="c_s_p",xlabel="r_p",legend=true, label="Fulltrained pinn", linestyle=:dashdot)
        #plot(p1,p2)
    end
    gif(anim, joinpath(exp_folder, "SPM_fulltrain_only_c_s_p_bigdomain.gif"),fps=30)
end
=#
main_train()

