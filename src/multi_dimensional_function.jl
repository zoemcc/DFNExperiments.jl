
struct ParameterizedMatrixApplyFuncType <: AbstractApplyFuncType end
struct VectorOfParameterizedMDFApplyFuncType <: AbstractApplyFuncType end

struct MultiDimensionalFunction{Func, FType <: AbstractApplyFuncType, Extra, IVsNT <: NamedTuple, DVsNT <: NamedTuple, DVDepsNT <: NamedTuple} 
    f::Func
    ftype::FType
    extra_data::Extra
    ivs::IVsNT # NT of IV_Symbol -> IV_LocalIndex
    dvs::DVsNT # NT of DV_Symbol -> DV_LocalIndex
    dv_deps::DVDepsNT # NT of DV_Symbol -> (NT of IV_Symbol -> IV_LocalIndex)

    function MultiDimensionalFunction(f, ftype::FType, ivs::IVsNT, dvs::DVsNT, dv_deps::DVDepsNT) where 
        {FType <: AbstractApplyFuncType, IVsNT <: NamedTuple, DVsNT <: NamedTuple, DVDepsNT <: NamedTuple}

        fntyivs = keys(ivs)
        num_ivs = length(ivs)
        for i in 1:num_ivs
            @assert ivs[i] isa Int
            @assert ivs[i] == i
        end

        fntydvs = keys(dvs)
        num_dvs = length(dvs)
        for i in 1:num_dvs
            @assert dvs[i] isa Int
            @assert dvs[i] == i
        end

        fntydvdeps = keys(dv_deps)
        @assert length(dv_deps) == num_dvs
        for i in 1:num_dvs
            dv_dep_i = dv_deps[i]
            @assert fntydvs[i] == fntydvdeps[i]
            @assert length(dv_dep_i) <= num_ivs
            fntydv_dep_i = keys(dv_dep_i)
            for j in 1:length(dv_dep_i)
                iv_symbol = fntydv_dep_i[j]
                iv_index = dv_dep_i[j]
                @assert iv_symbol isa Symbol
                @assert iv_symbol in fntyivs
                @assert iv_index == ivs[iv_symbol]
            end
        end

        if FType == VectorOfParameterizedMDFApplyFuncType
            f_types = typeof.(f)
            @assert length(f) == num_dvs
            all_fastchains = all(map(ft_i->ft_i <: FastChain, f_types))
            if all_fastchains
                initθ = DiffEqFlux.initial_params.(f)
                param_lengths = length.(initθ)
                param_lengths_cumulative = cumsum(vcat([0], param_lengths...)) 
                param_indices = [param_lengths_cumulative[i] + 1 : param_lengths_cumulative[i+1] for i in 1:num_dvs]
                extra = Dict(:param_indices => param_indices)
                fs_mdf = map(1:num_dvs) do i
                    dv_i = (fntydvs[i],)
                    iv_i = keys(dv_deps[i])
                    MultiDimensionalFunction(f[i], ParameterizedMatrixApplyFuncType(), iv_i, dv_i, (iv_i,))
                end
                f = fs_mdf
            else
                throw("currently only supports vector of FastChains")
            end

        else
            extra = Nothing
        end

        new{typeof(f), FType, typeof(extra), IVsNT, DVsNT, DVDepsNT}(f, ftype, extra, ivs, dvs, dv_deps)
    end

    function MultiDimensionalFunction(f, ftype::FType, ivs, dvs, dv_deps) where {FType <: AbstractApplyFuncType}

        for iv in ivs
            @assert iv isa Symbol
        end

        for dv in dvs
            @assert dv isa Symbol
        end

        @assert length(dv_deps) == length(dvs)
        for dv_dep in dv_deps
            @assert length(dv_dep) <= length(ivs)
            for iv in dv_dep
                @assert iv isa Symbol
                @assert iv in ivs
            end
        end

        ivs_nt = NamedTuple{ivs}(1:length(ivs))
        dvs_nt = NamedTuple{dvs}(1:length(dvs))
        dv_deps_tup = tuple(map(i->ivs_nt[dv_deps[i]], 1:length(dvs))...)
        dv_deps_nt = NamedTuple{dvs}(dv_deps_tup)

        MultiDimensionalFunction(f, ftype, ivs_nt, dvs_nt, dv_deps_nt)
    end
end

# Neural network

@generated function (f::MultiDimensionalFunction{Func, VectorOfParameterizedMDFApplyFuncType})(θ, x::Union{<:Real, AbstractVector{<:Real}}...; flat=Val{false}) where {Func}
    quote
        @assert length(x) == length(f.ivs) # assuming uniform local indexing
        num_dvs = length(f.dvs)
        param_indices = f.extra_data[:param_indices]
        f_evals = Vector{Array{eltype(θ)}}(undef, length(f.dvs))
        for i in 1:num_dvs
            ivs_i = f.dv_deps[i]
            # slice out the portion of θ that is appropriate and then call the sub MDF with that slice on the appropriate data
            x_i = x[collect(ivs_i)]
            param_indices_i = param_indices[i]
            f_evals[i] = f.f[i]((@view θ[param_indices_i]), x_i...; flat=flat)
        end
        return f_evals
    end
end

@generated function (f::MultiDimensionalFunction{Func, ParameterizedMatrixApplyFuncType})(θ, x::Union{<:Real, AbstractVector{<:Real}}...; flat=Val{false}) where {Func}
    quote
        # need to make a 2D array with the cartesian product of all the broadcasted 
        # and then apply the function to the point array
        # and reshape into the correct shape if flat is nothing
        if flat == Val{false}
            point_array, resize_info = cartesian_product(x...; flat=Val{true}, resize_info=Val{true}, debug=Val{false}) 
            transformed_array = f.f(point_array, θ)
            reshaped_transformed_array = reshape(transformed_array, :, resize_info...)
            return reshaped_transformed_array
        elseif flat == Val{true}
            point_array = cartesian_product(x...; flat=Val{true}, resize_info=Val{false}, debug=Val{false}) 
            transformed_array = f.f(point_array, θ)
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


