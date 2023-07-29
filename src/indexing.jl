
# This file is adapted from JuliaArrays/StaticArrays.jl License is MIT:
# https://github.com/JuliaArrays/StaticArrays.jl/blob/master/LICENSE.md

#######################################
## Multidimensional scalar indexing  ##
#######################################

# Note: all indexing behavior defaults to dense, linear indexing

@propagate_inbounds function Base.getindex(a::TupleVector{N}, ind::Int) where N
    @boundscheck checkbounds(a, ind)
    _getindex_scalar(Val(N), a, ind)
end

function _getindex_scalar(::Val{N}, a::TupleVector, ind::Int) where N
    Base.@_propagate_inbounds_meta
    return a[ind]
end

@propagate_inbounds function Base.setindex!(a::TupleVector{N}, value, ind::Int) where N
    @boundscheck checkbounds(a, ind)
    _setindex!_scalar(Val(N), a, value, ind)
    return a
end

@generated function _setindex!_scalar(::Val{N}, a::TupleVector, value, ind::Int) where N
    Base.@_propagate_inbounds_meta
    a[ind] = value
end

#########################
## Indexing utilities  ##
#########################

@inline index_size(::Val, ::Int) = Val(1)
@inline index_size(::Val, a::TupleVector{N}) where N = Val(N)
@inline index_size(s::Val, ::Colon) = s
@inline index_size(s::Val, a::SOneTo{n}) where n = Val(n)

@inline index_sizes(::S, ind) where {S<:Val} = index_size(S, ind)
@inline index_sizes(::S, inds...) where {S<:Val} = map(index_size, S, inds)

@inline index_sizes(::Int) = Val(1)
@inline index_sizes(a::TupleVector{N}) where N = Val(N)

@inline index_sizes2() = ()
@inline index_sizes2(::Int, inds...) = (Val(()), index_sizes2(inds...)...)
@inline index_sizes2(a::TupleVector{N}, inds...) where N = (Val(N), index_sizes2(inds...)...)

out_index_size(ind_size::Type{Val{N}}) where N = Val(N)
linear_index_size(ind_size::Type{Val{N}}) where N = N

out_index_size2(ind_sizes::Type{<:Val}...) = Val(_out_index_size2((), ind_sizes...))
@inline _out_index_size2(t::Tuple) = t
@inline _out_index_size2(t::Tuple, ::Type{Val{S}}, ind_sizes...) where {S} = _out_index_size2((t..., S...), ind_sizes...)

linear_index_size2(ind_sizes::Type{<:Val}...) = _linear_index_size2((), ind_sizes...)
@inline _linear_index_size2(t::Tuple) = t
@inline _linear_index_size2(t::Tuple, ::Type{Val{S}}, ind_sizes...) where {S} = _linear_index_size2((t..., prod(S)), ind_sizes...)

untuple(::Val{t}) where t = Val(t[1])

_ind(::Int, ::Type{Int}) = :ind
_ind(i::Int, ::Type{<:TupleVector}) = :(ind[$i])
_ind(j::Int, ::Type{Colon}) = j
_ind(j::Int, ::Type{<:SOneTo}) = j

_ind2(i::Int, ::Int, ::Type{Int}) = :(inds[$i])
_ind2(i::Int, j::Int, ::Type{<:TupleVector}) = :(inds[$i][$j])
_ind2(i::Int, j::Int, ::Type{Colon}) = j
_ind2(i::Int, j::Int, ::Type{<:SOneTo}) = j

################################
## Non-scalar linear indexing ##
################################

@inline function Base.getindex(a::TupleVector{N}, ::Colon) where N
    _getindex(a::TupleVector, Val(N), :)
end

@generated function _getindex(a::TupleVector, s::Val{N}, ::Colon) where N
    exprs = [:(a[$i]) for i = 1:N]
    return quote
        Base.@_inline_meta
        @inbounds return similar_type(a,s)(tuple($(exprs...)))
    end
end

@propagate_inbounds function Base.getindex(a::TupleVector, inds::TupleVector{N,Int}) where N
    _getindex(a, Val(N), inds)
end

@generated function _getindex(a::TupleVector, s::Val{N}, inds::TupleVector{N, Int}) where N
    exprs = [:(a[inds[$i]]) for i = 1:N]
    return quote
        Base.@_propagate_inbounds_meta
        similar_type(a, s)(tuple($(exprs...)))
    end
end

@inline function Base.setindex!(a::TupleVector{N}, v, ::Colon) where N
    _setindex!(a::TupleVector, v, Val(N), :)
    return a
end

@generated function _setindex!(a::TupleVector, v, ::Val{L}, ::Colon) where {L}
    exprs = [:(a[$i] = v) for i = 1:L]
    return quote
        Base.@_inline_meta
        @inbounds $(Expr(:block, exprs...))
    end
end

@generated function _setindex!(a::TupleVector, v::AbstractVector, ::Val{L}, ::Colon) where {L}
    exprs = [:(a[$i] = v[$i]) for i = 1:L]
    return quote
        Base.@_propagate_inbounds_meta
        @boundscheck if length(v) != L
            throw(DimensionMismatch("tried to assign $(length(v))-element array to length-$L destination"))
        end
        @inbounds $(Expr(:block, exprs...))
    end
end

@generated function _setindex!(a::TupleVector, v::TupleVector{M}, ::Val{L}, ::Colon) where {M,L}
    exprs = [:(a[$i] = v[$i]) for i = 1:L]
    return quote
        Base.@_propagate_inbounds_meta
        @boundscheck if M != L
            throw(DimensionMismatch("tried to assign $M-element array to length-$L destination"))
        end
        $(Expr(:block, exprs...))
    end
end

@propagate_inbounds function Base.setindex!(a::TupleVector, v, inds::TupleVector{N,Int}) where N
    _setindex!(a, v, Val(N), inds)
    return a
end

@generated function _setindex!(a::TupleVector, v, ::Val{N}, inds::TupleVector{N,Int}) where N
    exprs = [:(a[inds[$i]] = v) for i = 1:N]
    return quote
        Base.@_propagate_inbounds_meta
        similar_type(a, s)(tuple($(exprs...)))
    end
end

@generated function _setindex!(a::TupleVector, v::AbstractVector, ::Val{N}, inds::TupleVector{N,Int}) where N
    exprs = [:(a[inds[$i]] = v[$i]) for i = 1:N]
    return quote
        Base.@_propagate_inbounds_meta
        @boundscheck if length(v) != $N
            throw(DimensionMismatch("tried to assign $(length(v))-element array to length-$N destination"))
        end
        $(Expr(:block, exprs...))
    end
end

@generated function _setindex!(a::TupleVector, v::TupleVector{M}, ::Val{N}, inds::TupleVector{N,Int}) where {N,M}
    exprs = [:(a[inds[$i]] = v[$i]) for i = 1:N]
    return quote
        Base.@_propagate_inbounds_meta
        @boundscheck if M != N
            throw(DimensionMismatch("tried to assign $M-element array to length-$N destination"))
        end
        $(Expr(:block, exprs...))
    end
end

###########################################
## Multidimensional non-scalar indexing  ##
###########################################

# To intercept `A[i1, ...]` where all `i` indexes have static sizes,
# create a wrapper used to mark non-scalar indexing operations.
# We insert this at a point in the dispatch hierarchy where we can intercept any
# `typeof(A)` (specifically, including dynamic arrays) without triggering ambiguities.

struct TupleIndexing{I}
    ind::I
end
unwrap(i::TupleIndexing) = i.ind

function Base.to_indices(A, I::Tuple{Vararg{Union{Integer, CartesianIndex, TupleVector{N,Int} where N}}})
    inds = to_indices(A, axes(A), I)
    return map(TupleIndexing, inds)
end

# getindex

Base.@propagate_inbounds function Base.getindex(a::TupleVector{N}, ind::Union{Int, TupleVector{M, Int} where M, SOneTo, Colon}) where N
    _getindex(a, index_sizes(Val(N), ind), ind)
end

function Base._getindex(::IndexStyle, A::AbstractVector, ind::TupleIndexing)
    return StaticVectors._getindex(A, index_sizes(unwrap(ind)), unwrap(ind))
end

function Base._getindex(::IndexStyle, A::AbstractArray, i1::TupleIndexing, I::TupleIndexing...)
    inds = (unwrap(i1), map(unwrap, I)...)
    return StaticVectors._getindex(A, index_sizes2(inds...), inds)
end

@generated function _getindex(a::AbstractVector, ind_size::Val, ind)
    newsize = out_index_size(ind_size)
    linearsize = linear_index_size(ind_size)
    exprs = [:(getindex(a, $(_ind(i, ind)))) for i ∈ 1:linearsize]
    quote
        Base.@_propagate_inbounds_meta
        similar_type(a, $newsize)(tuple($(exprs...)))
    end
end

@generated function _getindex(a::AbstractArray, ind_sizes::Tuple{Vararg{Val}}, inds)
    newsize = untuple(out_index_size2(ind_sizes.parameters...))
    linearsizes = linear_index_size2(ind_sizes.parameters...)
    exprs = Array{Expr}(undef, linearsizes)

    # Iterate over input indices
    ind_types = inds.parameters
    current_ind = ones(Int,length(linearsizes))
    more = !isempty(exprs)
    while more
        exprs_tmp = [_ind2(i, current_ind[i], ind_types[i]) for i = 1:length(linearsizes)]
        exprs[current_ind...] = :(getindex(a, $(exprs_tmp...)))

        # increment current_ind
        current_ind[1] += 1
        for i ∈ 1:length(linearsizes)
            if current_ind[i] > linearsizes[i]
                if i == length(linearsizes)
                    more = false
                    break
                else
                    current_ind[i] = 1
                    current_ind[i+1] += 1
                end
            else
                break
            end
        end
    end

    quote
        Base.@_propagate_inbounds_meta
        similar_type(a, $newsize)(tuple($(exprs...)))
    end
end

# setindex!

@propagate_inbounds function Base.setindex!(a::TupleVector{N}, value, ind::Union{Int, TupleVector{M, Int} where M, Colon}) where N
    _setindex!(a, value, index_sizes(N, ind), ind)
end

function Base._setindex!(::IndexStyle, a::AbstractVector, value, ind::TupleIndexing)
    return StaticVectors._setindex!(a, value, index_sizes(ind), unwrap(ind))
end

function Base._setindex!(::IndexStyle, a::AbstractArray, value, i1::TupleIndexing, I::TupleIndexing...)
    inds = (unwrap(i1), map(unwrap, I)...)
    return StaticVectors._setindex!(a, value, index_sizes2(inds...), inds)
end

# setindex! from a scalar
@generated function _setindex!(a::AbstractVector, value, ind_sizes::Val, ind)
    linearsize = linear_index_size(ind_size)
    exprs = [:(setindex!(a, value, $(_ind(i, ind)))) for i ∈ 1:linearsize]
    quote
        Base.@_propagate_inbounds_meta
        $(exprs...)
        return a
    end
end

@generated function _setindex!(a::AbstractArray, value, ind_sizes::Tuple{Vararg{Val}}, inds)
    linearsizes = linear_index_size2(ind_sizes.parameters...)
    exprs = Array{Expr}(undef, linearsizes)

    # Iterate over input indices
    ind_types = inds.parameters
    current_ind = ones(Int,length(ind_types))
    more = !isempty(exprs)
    while more
        exprs_tmp = [_ind2(i, current_ind[i], ind_types[i]) for i = 1:length(ind_types)]
        exprs[current_ind...] = :(setindex!(a, value, $(exprs_tmp...)))

        # increment current_ind
        current_ind[1] += 1
        for i ∈ 1:length(linearsizes)
            if current_ind[i] > linearsizes[i]
                if i == length(linearsizes)
                    more = false
                    break
                else
                    current_ind[i] = 1
                    current_ind[i+1] += 1
                end
            else
                break
            end
        end
    end

    quote
        Base.@_propagate_inbounds_meta
        $(exprs...)
        return a
    end
end

# setindex! from an array
@generated function _setindex!(a::AbstractVector, v::AbstractVector, ind_size::Val, ind)
    linearsize = linear_index_size(ind_size)
    # Iterate over input indices
    exprs = [:(setindex!(a, v[$j], $(_ind(j, ind)))) for j ∈ 1:linearsize]
    quote
        Base.@_propagate_inbounds_meta
        if length(v) != $linearsize
            newsize = $linearsize
            throw(DimensionMismatch("tried to assign $(length(v))-element array to $newsize destination"))
        end
        $(exprs...)
        return a
    end
end

@generated function _setindex!(a::AbstractArray, v::AbstractArray, ind_sizes::Tuple{Vararg{Val}}, inds)
    linearsizes = linear_index_size2(ind_sizes.parameters...)
    exprs = Array{Expr}(undef, linearsizes)

    # Iterate over input indices
    ind_types = inds.parameters
    current_ind = ones(Int,length(ind_types))
    more = true
    j = 1
    while more
        exprs_tmp = [_ind2(i, current_ind[i], ind_types[i]) for i = 1:length(ind_types)]
        exprs[current_ind...] = :(setindex!(a, v[$j], $(exprs_tmp...)))

        # increment current_ind
        current_ind[1] += 1
        for i ∈ 1:length(linearsizes)
            if current_ind[i] > linearsizes[i]
                if i == length(linearsizes)
                    more = false
                    break
                else
                    current_ind[i] = 1
                    current_ind[i+1] += 1
                end
            else
                break
            end
        end
        j += 1
    end

    quote
        Base.@_propagate_inbounds_meta
        if length(v) != $(prod(linearsizes))
            newsize = $(linearsizes)
            throw(DimensionMismatch("tried to assign $(length(v))-element array to $newsize destination"))
        end
        $(exprs...)
        return a
    end
end

# checkindex

Base.checkindex(B::Type{Bool}, inds::AbstractUnitRange, i::TupleIndexing{T}) where T = Base.checkindex(B, inds, unwrap(i))

# unsafe_view

# unsafe_view need only deal with vargs of `TupleIndexing`, as wrapped by to_indices.
# i1 is explicitly specified to avoid ambiguities with Base
Base.unsafe_view(A::AbstractVector, ind::TupleIndexing) = Base.unsafe_view(A, unwrap(ind))
Base.unsafe_view(A::AbstractArray, i1::TupleIndexing, indices::TupleIndexing...) = Base.unsafe_view(A, unwrap(i1), map(unwrap, indices)...)

# Views of views need a new method for Base.SubArray because storing indices
# wrapped in TupleIndexing in field indices of SubArray causes all sorts of problems.
# Additionally, in some cases the SubArray constructor may be called directly
# instead of unsafe_view so we need this method too (Base._maybe_reindex
# is a good example)
# the tuple indices has to have at least one element to prevent infinite
# recursion when viewing a zero-dimensional array (see issue #705)
Base.SubArray(A::AbstractVector, indices::Tuple{TupleIndexing, Vararg{TupleIndexing}}) = Base.SubArray(A, map(unwrap, indices))
