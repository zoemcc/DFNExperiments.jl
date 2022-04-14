begin 
    #using DFNExperiments
    #using ModelingToolkit, DomainSets
    using Symbolics
    #import ModelingToolkit: Interval, infimum, supremum
    #using LabelledArrays
    #using NeuralPDE
    #using JSON
    #using Random
    using Flux, DiffEqFlux
    using LinearAlgebra
    #using Plots
    #using GalacticOptim
end

abstract type AbstractFuncType end

struct MatrixApplyFuncType <: AbstractFuncType end

# Neural network
chain = [FastChain(FastDense(2,8,Flux.tanh),FastDense(8,1))]
initθ = map(c -> Float64.(c), DiffEqFlux.initial_params.(chain))


function gen_range_slices(vardomains::AbstractVector{Symbolics.VarDomainPairing}, num_points_per_dim=128)
    map(vardomain -> range(vardomain.domain, num_points_per_dim), vardomains)
end

num_points_per_dim = 4


function eval_func_slice(func, vardomains::AbstractVector{Symbolics.VarDomainPairing}, symbols::AbstractVector{Num}, num_points_per_dim)
    sym_indices = map(sym->findmax(map(vardomain->isequal(vardomain.variables, sym.val), vardomains))[2], symbols)
    
    linear_indices = axes(ones(repeat([num_points_per_dim], 4)...))
    cart_indices = CartesianIndices(linear_indices)
    for I in eachindex(cart_indices)
        @show I
    end

    range_slices = gen_range_slices(vardomains, num_points_per_dim)
    @show range_slices

    value_from_index(index) = []

end

function gen_subslices(constant_values, nonconstant_ranges)
    # 
end



slice_vars = [x2, x1, x3]
func_sliced = eval_func_slice(sol_func2, domains, slice_vars, num_points_per_dim)

@generated function eval_pde_func(f, x::Union{Real, AbstractVector{<:Real}}...) 
    broadcast_dims = map(x_i->x_i <: AbstractVector, x)
    broadcast_indices = Set(map(first, filter(i_b->i_b[2], collect(enumerate(broadcast_dims)))))
    N = length(x)
    eltypes = map(x_i -> x_i <: AbstractVector ? eltype(x_i) : x_i, x)
    promotion_type = promote_type(eltypes...)
    quote
        # need to make a 2D array with the cartesian product of all the broadcasted 
        broadcast_lengths = length.(x)
        num_points = prod(broadcast_lengths)
        point_array = Array{$promotion_type, 2}(undef, $N, num_points) # TODO: make the AbstractArray type be more general
        iter_indices = CartesianIndices(tuple(broadcast_lengths...))
        for (i, index) in enumerate(iter_indices)
            for j in 1:$N
                if j in $(broadcast_indices)
                    point_array[j, i] = x[j][index.I[j]]
                else
                    point_array[j, i] = x[j]
                end
            end
        end

        # apply the function to the point array
        output_pre_reshape = f(point_array)
        out_dim = size(output_pre_reshape, 1)
        output = reshape(output_pre_reshape, (out_dim, broadcast_lengths...))
        return (N=$N, broadcast_dims=$(broadcast_dims), broadcast_indices=$(broadcast_indices), broadcast_lengths=broadcast_lengths, 
            num_points=num_points, point_array=point_array, shape=size(point_array), iter_indices=iter_indices,
            output_pre_reshape=output_pre_reshape, out_dim=out_dim, output=output,
            eltypes=$eltypes, promotion_type=$promotion_type)
    end

end

@generated function (f::MultiDimensionalFunction{Func, MatrixApplyFuncType})(x::Union{<:Real, AbstractVector{<:Real}}...; flat=Val{false}) where {Func}
    quote
        # need to make a 2D array with the cartesian product of all the broadcasted 
        # and then apply the function to the point array
        # and reshape into the correct shape if flat is nothing
        if flat == Val{false}
            point_array, resize_info = cartesian_product(x...; flat=Val{true}, resize_info=Val{true}, debug=Val{false}) 
            transformed_array = f.f(point_array)
            reshaped_transformed_array = reshape(transformed_array, :, resize_info...)
            return reshaped_transformed_array
        elseif flat == Val{true}
            point_array = cartesian_product(x...; flat=Val{true}, resize_info=Val{false}, debug=Val{false}) 
            transformed_array = f.f(point_array)
            return transformed_array
        else
            throw("Flat must be either Val{false} or Val{true}")
        end
    end
end


@generated function cartesian_product(x::Union{<:Real, AbstractVector{<:Real}}...; flat=Val{false}, debug=Val{false}, resize_info=Val{false}) 
    broadcast_dims = map(x_i->x_i <: AbstractVector, x)
    broadcast_indices = Set(map(first, filter(i_b->i_b[2], collect(enumerate(broadcast_dims)))))
    N = length(x)
    eltypes = map(x_i -> x_i <: AbstractVector ? eltype(x_i) : x_i, x)
    promotion_type = promote_type(eltypes...)
    quote
        # need to make a 2D array with the cartesian product of all the broadcasted 
        broadcast_lengths = length.(x)
        num_points = prod(broadcast_lengths)
        iter_indices = CartesianIndices(tuple(broadcast_lengths...))

        if flat == Val{false} 
            point_array = Array{$promotion_type, ($N + 1)}(undef, $N, broadcast_lengths...) # TODO: make the AbstractArray type be more general
            for (i, index) in enumerate(iter_indices)
                for j in 1:$N
                    if j in $(broadcast_indices)
                        point_array[j, index] = x[j][index.I[j]]
                    else
                        point_array[j, index] = x[j]
                    end
                end
            end
        elseif flat == Val{true}
            point_array = Array{$promotion_type, 2}(undef, $N, num_points) # TODO: make the AbstractArray type be more general
            for (i, index) in enumerate(iter_indices)
                for j in 1:$N
                    if j in $(broadcast_indices)
                        point_array[j, i] = x[j][index.I[j]]
                    else
                        point_array[j, i] = x[j]
                    end
                end
            end
        else
            throw("Flat must be either Val{false} or Val{true}")
        end

        if debug == Val{false}
            if resize_info == Val{false}
                return point_array
            elseif resize_info == Val{true}
                return (point_array=point_array, resize_shape=broadcast_lengths)
            else
                throw("resize_info must be either Val{false} or Val{true}")
            end
        elseif debug == Val{true}
            return (N=$N, broadcast_dims=$(broadcast_dims), broadcast_indices=$(broadcast_indices), broadcast_lengths=broadcast_lengths, 
                num_points=num_points, point_array=point_array, shape=size(point_array), iter_indices=iter_indices)
        else
            throw("Debug must be either Val{false} or Val{true}")
        end
    end

end

struct MultiDimensionalFunction{Func, FType <: AbstractFuncType, IVsNum, DVsNum, DVDepsTuple} 
    f::Func
    ftype::FType
    ivs::NTuple{IVsNum, Symbol}
    dvs::NTuple{DVsNum, Symbol}
    dv_deps::DVDepsTuple 

    function MultiDimensionalFunction(f, ftype::FType, ivs::NTuple{IVsNum, Symbol}, dvs::NTuple{DVsNum, Symbol}, dv_deps::Tuple) where
        {FType <: AbstractFuncType, IVsNum, DVsNum}

        @assert length(dv_deps) == DVsNum
        for dv_dep in dv_deps
            @assert dv_dep isa Tuple
            @assert length(dv_dep) <= IVsNum
            for iv in dv_dep
                @assert iv isa Symbol
            end
        end

        new{typeof(f), FType, IVsNum, DVsNum, typeof(dv_deps)}(f, ftype, ivs, dvs, dv_deps)
    end
end

ivs = (:x1, :x2)
dvs = (:u1,)
dv_deps = ((:x1, :x2),)

chain_apply(x) = chain[1](x, initθ[1])
mdf = MultiDimensionalFunction(chain_apply, MatrixApplyFuncType(), ivs, dvs, dv_deps)

x1s = 0:0.1:1
x1i = 0.1
x2s = 0:0.2:1
x2i = 0.4
ys = 0:0.1:0.2
da = mdf(ys, ys; flat=Val{false})
grids = cartesian_product(x1s, x2i, ys; flat=nothing)
grids = cartesian_product(x1s, x2i, ys; flat=true)
grids = cartesian_product(ys, ys; flat=true)
f_eval(x) = phi[1](x, res.u)
allys = eval_pde_func(f_eval, ys, ys, ys, Float32(x1i))
eltypes = allys.eltypes
promote_type(eltypes...)
#eval_pde_func(x1i, x2i)
#ful2 = eval_pde_func(Float32.(x1s), x2s)
#ful3 = eval_pde_func(Float32.(x1s), x2s, x2s)
#eval_pde_func(x1i, x2s)
#bd = eval_pde_func(x1i)
#bd = eval_pde_func(x1s)

lengths2 =  (length).([x1s, x2s])
indices = CartesianIndices(tuple(lengths2...))


"""
things I want to plot:

pde values and slices of them (heterogenous)



"""


"""
desired interface:

AbstractPDESolution <: AbstractNoTimeSolution
FunctionPDESolution <: AbstractPDESolution
PointwisePDESolution <: AbstractPDESolution

would like to treat both as a function and both as a grid

res::FunctionPDESolution = solve(prob, opt)
sample_slice()


"""



"""
const TYPE_LENGHT_LIMIT = Ref(20)
function Base.print_type_stacktrace(io, type; color=:normal)
str = first(sprint(show, type, context=io), TYPE_LENGHT_LIMIT[])
i = findfirst('{', str)
if isnothing(i) || !get(io, :backtrace, false)::Bool
printstyled(io, str; color=color)
else
printstyled(io, str[1:prevind(str,i)]; color=color)
printstyled(io, str[i:end]; color=:light_black)
end
end """
