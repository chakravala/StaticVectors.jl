# StaticVectors.jl

*Statically sized tuple vectors for Julia*

Subtypes of `TupleVector` are provided as a light weight alternative to [StaticArrays.jl](https://github.com/JuliaArrays/StaticArrays.jl).
This package was created for use in [AbstractTensors.jl](https://github.com/chakravala/AbstractTensors.jl), [Grassmann.jl](https://github.com/chakravala/Grassmann.jl), [FieldAlgebra.jl](https://github.com/chakravala/FieldAlgebra.jl), [Similitude.jl](https://github.com/chakravala/Similitude.jl), [Geophysics.jl](https://github.com/chakravala/Geophysics.jl), and various other related repositories developed by [chakravala](https://github.com/chakravala).

**TupleVector** provides a framework for implementing statically sized tuple vectors
in Julia, using the abstract type `TupleVector{N,T} <: AbstractVector{T}`.
Subtypes of `TupleVector` will provide fast implementations of common array and
linear algebra operations. Note that here "statically sized" means that the
size can be determined from the *type*, and "static" does **not** necessarily
imply `immutable`.

The package also provides some concrete static vector types: `Values` which may be used as-is (or else embedded in your own type).
Mutable versions `Variables` are also exported, as well
as `FixedVector` for annotating standard `Vector`s with static size information.

### Quick start

Add *StaticVectors* from the [Pkg REPL](https://docs.julialang.org/en/latest/stdlib/Pkg/#Getting-Started-1), i.e., `pkg> add StaticVectors`. Then:
```julia
using StaticVectors

# Create Values using various forms, using constructors, functions or macros
v1 = Values(1, 2, 3)
v1.v === (1, 2, 3) # Values uses a tuple for internal storage
v2 = Values{3,Float64}(1, 2, 3) # length 3, eltype Float64
v5 = zeros(Values{3}) # defaults to Float64
v7 = Values{3}([1, 2, 3]) # Array conversions must specify size

# Can get size() from instance or type
size(v1) == (3,)
size(typeof(v1)) == (3,)

# Supports all the common operations of AbstractVector
v7 = v1 + v2
v8 = sin.(v2)

# Indexing can also be done using static vectors of integers
v1[1] === 1
v1[:] === v1
typeof(v1[[1,2,3]]) <: Vector # Can't determine size from the type of [1,2,3]
```

### Approach

The package provides a range of different useful built-in `TupleVector` types,
which include mutable and immutable vectors based upon tuples, vectors based upon
structs, and wrappers of `Vector`. There is a relatively simple interface for
creating your own, custom `TupleVector` types, too.

This package also provides methods for a wide range of `AbstractVector` functions,
specialized for (potentially immutable) `TupleVector`s. Many of Julia's
built-in method definitions inherently assume mutability, and further
performance optimizations may be made when the size of the vector is known to the
compiler. One example of this is by loop unrolling, which has a substantial
effect on small arrays and tends to automatically trigger LLVM's SIMD
optimizations. In combination with intelligent fallbacks to
the methods in Base, we seek to provide a comprehensive support for statically
sized vectors, large or small, that hopefully "just works".

`TupleVector` is directly inspired from [StaticArrays.jl](https://github.com/JuliaArrays/StaticArrays.jl).
