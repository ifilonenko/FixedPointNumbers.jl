# __precompile__()

module FixedPointNumbers

using Base: reducedim_initarray

import Base: ==, <, <=, -, +, *, /, ~, isapprox,
             convert, promote_rule, show, showcompact, isinteger, abs, decompose,
             isnan, isinf, isfinite,
             zero, oneunit, one, typemin, typemax, realmin, realmax, eps, sizeof, reinterpret,
             float, trunc, round, floor, ceil, bswap,
             div, fld, rem, mod, mod1, fld1, min, max, minmax,
             start, next, done, r_promote, reducedim_init, rand,
             dot
if isdefined(Base, :rem1)
    import Base: rem1
end
using Base: @pure

### Rounding Schemes
abstract type RoundingScheme end
abstract type Saturated <: RoundingScheme end
abstract type Exact <: RoundingScheme end
abstract type SatAndRandomized <: Saturated end
abstract type SatAndNearestNeighbor <: Saturated end
abstract type ExactAndRandomized <: Exact end
abstract type ExactAndNearestNeighbor <: Exact end

### FixedPoints
# T => BaseType
# f => Number of Bytes reserved for fractional part
abstract type FixedPoint{T <: Integer, f} <: Real end
# T => BaseType
# s => Scale factor that can be any Real number
abstract type ScaledFixedPoint{T <: Integer, f, s, r <: RoundingScheme} <: Real end
const GenericFixedPoint = Union{FixedPoint, ScaledFixedPoint}

export
    RoundingScheme,
    Saturated,
    Exact,
    SatAndRandomized,
    SatAndNearestNeighbor,
    ExactAndRandomized,
    ExactAndNearestNeighbor,
    GenericFixedPoint,
    FixedPoint,
    ScaledFixedPoint,
    Fixed,
    Normed,
    Scaled,
# "special" typealiases
    # Q and U typealiases are exported in separate source files
# literal constructor constants
    uf8,
    uf10,
    uf12,
    uf14,
    uf16,
# Functions
    scaledual,
    widen1,
    maxf,
    get_T,
    get_f,
    get_s,
    get_r,
# Constants
    _log2_10

reinterpret(x::GenericFixedPoint) = x.i
reinterpret(::Type{T}, x::FixedPoint{T,f}) where {T,f} = x.i
reinterpret(::Type{T}, x::ScaledFixedPoint{T,f,s,r}) where {T,f,s,r <: RoundingScheme} = x.i

# construction using the (approximate) intended value, i.e., N0f8
*(x::Real, ::Type{X}) where {X<:FixedPoint} = X(x)

# comparison
==(x::T, y::T) where {T <: GenericFixedPoint} = x.i == y.i
 <(x::T, y::T) where {T <: GenericFixedPoint} = x.i < y.i
<=(x::T, y::T) where {T <: GenericFixedPoint} = x.i <= y.i

"""
    isapprox(x::FixedPoint, y::FixedPoint; rtol=0, atol=max(eps(x), eps(y)))

For FixedPoint numbers, the default criterion is that `x` and `y` differ by no more than `eps`, the separation between adjacent fixed-point numbers.
"""
function isapprox(x::T, y::T; rtol=0, atol=max(eps(x), eps(y))) where {T <: FixedPoint}
    maxdiff = T(atol+rtol*max(abs(x), abs(y)))
    rx, ry, rd = reinterpret(x), reinterpret(y), reinterpret(maxdiff)
    abs(signed(widen1(rx))-signed(widen1(ry))) <= rd
end
function isapprox(x::FixedPoint, y::FixedPoint; rtol=0, atol=max(eps(x), eps(y)))
    isapprox(promote(x, y)...; rtol=rtol, atol=atol)
end
function isapprox(x::ScaledFixedPoint, y::Real,num_its=100000,rtol=0)
    diff = float(eps(x)) + rtol*float(eps(x))
    abs(mean([float(x) for a in 1:num_its]) - y) <= float(diff)
end

# predicates
isinteger(x::FixedPoint{T,f}) where {T,f} = (x.i&(1<<f-1)) == 0
isinteger(x::ScaledFixedPoint{T,f,s,r}) where {T,f,s,r} = (x.i&(1<<f-1)) == 0
isinteger(x::Real) = x-floor(x) == 0.0
isinteger(x::R,tolerance::R) where {R<:Real} = (x-floor(x)) <= tolerance


# traits
typemax(::Type{T}) where {T <: FixedPoint} = T(typemax(rawtype(T)), 0)
typemin(::Type{T}) where {T <: FixedPoint} = T(typemin(rawtype(T)), 0)
realmin(::Type{T}) where {T <: FixedPoint} = eps(T)
realmax(::Type{T}) where {T <: FixedPoint} = typemax(T)

widen1(::Type{Int8})   = Int16
widen1(::Type{UInt8})  = UInt16
widen1(::Type{Int16})  = Int32
widen1(::Type{UInt16}) = UInt32
widen1(::Type{Int32})  = Int64
widen1(::Type{UInt32}) = UInt64
widen1(::Type{Int64})  = Int128
widen1(::Type{UInt64}) = UInt128
widen1(::Type{UInt128}) = UInt128
widen1(x::Integer) = x % widen1(typeof(x))

# Max number of floating bits
maxf(::Type{Int8})=7
maxf(::Type{Int16})=15
maxf(::Type{Int32})=31
maxf(::Type{Int64})=63

const ShortInts = Union{Int8,UInt8,Int16,UInt16}

floattype(::Type{FixedPoint{T,f}}) where {T <: ShortInts,f} = Float32
floattype(::Type{ScaledFixedPoint{T,f,s,r}}) where {T <: ShortInts,f,s, r<: RoundingScheme} = Float32
floattype(::Type{FixedPoint{T,f}}) where {T <: Integer,f} = Float64
floattype(::Type{ScaledFixedPoint{T,f,s,r}}) where {T <: Integer,f,s, r<: RoundingScheme} = Float64
floattype(::Type{F}) where {F <: GenericFixedPoint} = floattype(supertype(F))
floattype(x::GenericFixedPoint) = floattype(typeof(x))

nbitsfrac(::Type{FixedPoint{T,f}}) where {T <: Integer,f} = f
nbitsfrac(::Type{ScaledFixedPoint{T,f,s,r}}) where {T <: Integer,f,s,r<: RoundingScheme} = f
nbitsfrac(::Type{F}) where {F <: GenericFixedPoint} = nbitsfrac(supertype(F))

rawtype(::Type{FixedPoint{T,f}}) where {T <: Integer,f} = T
rawtype(::Type{ScaledFixedPoint{T,f,s,r}}) where {T <: Integer,f,s,r<: RoundingScheme} = T
rawtype(::Type{F}) where {F <: GenericFixedPoint} = rawtype(supertype(F))
rawtype(x::GenericFixedPoint) = rawtype(typeof(x))

# This IOBuffer is used during module definition to generate typealias names
_iotypealias = IOBuffer()

# Printing. These are used to generate type-symbols, so we need them
# before we include any files.
function showtype(io::IO, ::Type{X}) where {X <: FixedPoint}
    print(io, typechar(X))
    f = nbitsfrac(X)
    m = sizeof(X)*8-f-signbits(X)
    print(io, m, 'f', f)
    io
end
function showtype(io::IO, ::Type{X}) where {X <: ScaledFixedPoint}
    print(io, typechar(X))
    f = nbitsfrac(X)
    m = sizeof(X)*8-f-signbits(X)
    print(io, m, 'f', f)
    io
end
function show(io::IO, x::ScaledFixedPoint{T,f,s,r}) where {T,f,s,r <: RoundingScheme}
    showcompact(io, x)
    showtype(io, typeof(x))
    print(io, "($s)")
    io
end
function show(io::IO, x::FixedPoint{T,f}) where {T,f}
    showcompact(io, x)
    showtype(io, typeof(x))
end

const _log2_10 = 3.321928094887362
showcompact(io::IO, x::ScaledFixedPoint{T,f,s,r}) where {T,f,s,r <: RoundingScheme} = show(io, round(convert(Float64,x),ceil(Int,f/_log2_10)))
showcompact(io::IO, x::FixedPoint{T,f}) where {T,f} = show(io, round((Float64,x), ceil(Int,f/_log2_10)))

if VERSION >= v"0.7.0-DEV.1790"
    function Base.showarg(io::IO, a::Array{T}, toplevel) where {T<:FixedPoint}
        toplevel || print(io, "::")
        print(io, "Array{")
        showtype(io, T)
        print(io, ",$(ndims(a))}")
        toplevel && print(io, " with eltype ", T)
    end
end

include("fixed.jl")
include("normed.jl")
include("scaled.jl")
include("deprecations.jl")

eps(::Type{T}) where {T <: GenericFixedPoint} = T(oneunit(rawtype(T)),0)
eps(::T) where {T <: GenericFixedPoint} = eps(T)
sizeof(::Type{T}) where {T <: GenericFixedPoint} = sizeof(rawtype(T))

# Promotions for reductions
const Treduce = Float64
r_promote(::typeof(+), x::FixedPoint{T}) where {T} = Treduce(x)
r_promote(::typeof(*), x::FixedPoint{T}) where {T} = Treduce(x)

reducedim_init(f::typeof(identity),
                              op::typeof(+),
                              A::AbstractArray{T}, region) where {T <: FixedPoint} =
    reducedim_initarray(A, region, zero(Treduce))
reducedim_init(f::typeof(identity),
                              op::typeof(*),
                              A::AbstractArray{T}, region) where {T <: FixedPoint} =
    reducedim_initarray(A, region, oneunit(Treduce))

for f in (:div, :fld, :fld1)
    @eval begin
        $f(x::T, y::T) where {T <: FixedPoint} = $f(reinterpret(x),reinterpret(y))
    end
end
for f in (:rem, :mod, :mod1, :rem1, :min, :max)
    if f === :rem1 && !isdefined(Base, :rem1)
        continue
    end
    @eval begin
        $f(x::T, y::T) where {T <: FixedPoint} = T($f(reinterpret(x),reinterpret(y)),0)
    end
end

# When multiplying by a float, reduce two multiplies to one.
# Particularly useful for arrays.
scaledual(Tdual::Type, x) = oneunit(Tdual), x
scaledual(b::Tdual, x) where {Tdual <: Number} = b, x
scaledual(Tdual::Type, x::Union{T,AbstractArray{T}}) where {T <: FixedPoint} =
    (Tdual, 1/oneunit(T)), reinterpret(rawtype(T), x)
scaledual(b::Tdual, x::Union{T,AbstractArray{T}}) where {Tdual <: Number,T <: FixedPoint} =
    (Tdual, b/oneunit(T)), reinterpret(rawtype(T), x)

@noinline function throw_error(::Type{T}, x) where {T <: FixedPoint}
    n = 2^(8*sizeof(T))
    bitstring = sizeof(T) == 1 ? "an 8-bit" : "a $(8*sizeof(T))-bit"
    io = IOBuffer()
    showcompact(io, typemin(T)); Tmin = String(take!(io))
    showcompact(io, typemax(T)); Tmax = String(take!(io))
    throw(ArgumentError("$T is $bitstring type representing $n values from $Tmin to $Tmax; cannot represent $x"))
end

rand(::Type{T}) where {T <: FixedPoint} = reinterpret(T, rand(rawtype(T)))
rand(::Type{T}, sz::Dims) where {T <: FixedPoint} = reinterpret(T, rand(rawtype(T), sz))

end # module
