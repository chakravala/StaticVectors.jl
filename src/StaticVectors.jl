
# this file is inspired from StaticArrays.jl
# https://github.com/JuliaArrays/StaticArrays.jl

module StaticVectors

@inline inv(z) = Base.inv(z)
@inline /(a,b) = Base.:/(a,b)
@inline -(a,b) = Base.:-(a,b)
-(x) = Base.:-(x)
-(x::Symbol) = :(-$x)
for (OP,op) ∈ ((:∏,:*),(:∑,:+))
    @eval @inline $OP(x...) = Base.$op(x...)
end

import LinearAlgebra

export TupleVector, Values, Variables, FixedVector

import Base: @propagate_inbounds, @_inline_meta, @pure

if haskey(ENV,"STATICJL")
    import StaticArrays: SVector, MVector, SizedVector, StaticVector, _diff
    const Values,Variables,FixedVector,TupleVector = SVector,MVector,SizedVector,StaticVector
else

abstract type TupleVector{N,T} <: AbstractVector{T} end

# Being a member of TupleMatrixLike or TupleVectorLike implies that Val(A)
# returns a static Val instance (none of the dimensions are Dynamic). The converse may not be true.
# These are akin to aliases like StridedArray and in similarly bad taste, but the current approach
# in Base necessitates their existence.
const TupleMatrixLike{n,T} = Union{
    LinearAlgebra.Transpose{T, <:TupleVector{n,T}},
    LinearAlgebra.Adjoint{T, <:TupleVector{n,T}},
    LinearAlgebra.Diagonal{T, <:TupleVector{n,T}},
}
const TupleVectorLike{n,T} = Union{TupleVector{n,T}, TupleMatrixLike{n,T}}

@pure Base.length(::T) where T<:TupleVector{N} where N = N
@pure Base.length(::Type{<:TupleVector{N}}) where N = N
@pure Base.lastindex(::T) where T<:TupleVector{N} where N = N
@pure Base.size(::T) where T<:TupleVector{N} where N = (N,)
@pure Base.size(::Type{<:TupleVector{N}}) where N = (N,)
@pure @inline Base.size(t::T,d::Int) where T<:TupleVector{N} where N = d > 1 ? 1 : length(t)
@pure @inline Base.size(t::Type{<:TupleVector},d::Int) = d > 1 ? 1 : length(t)

include("SOneTo.jl")
include("util.jl")
include("traits.jl")
include("Values.jl")
include("Variables.jl")
include("FixedVector.jl")
include("initializers.jl")
include("convert.jl")
include("abstractvector.jl")
include("indexing.jl")
include("broadcast.jl")
include("mapreduce.jl")
include("arraymath.jl")
include("linalg.jl")

const SVector,MVector,SizedVector = Values,Variables,FixedVector

end

"""
    count(a::Int,b::Int) -> Values

Generates an incremental `Int` list of `Values` from `a` to `b`.
"""
@pure Base.count(a::Int,b::Int) = countvalues(a,b)

"""
    countvalues(a::Int,b::Int) -> Values

Generates an incremental `Int` list of `Values` from `a` to `b`.
"""
@pure countvalues(a::Int,b::Int) = Values{max(0,b-a+1),Int}(a:b...)

"""
    evenvalues(a::Int,b::Int) -> Values

Generates an `Int` list of `Values` by skipping in twos from `a` to `b`.
"""
@pure evenvalues(a::Int,b::Int) = Values{((b-a)÷2)+1,Int}(a:2:b...)

"""
    evens(a::Int,b::Int) -> Values

Generates an `Int` list of `Values` by skipping in twos from `a` to `b`.
"""
const evens = evenvalues

end # module
