module FixedPointNumbers

using Base: reducedim_initarray

import Base: ==, <, <=, -, +, *, /, ~, isapprox,
             convert, promote_rule, show, showcompact, isinteger, abs, decompose,
             isnan, isinf, isfinite,
             zero, oneunit, one, typemin, typemax, realmin, realmax, eps, sizeof, reinterpret,
             float, trunc, round, floor, ceil, bswap,
             div, fld, rem, mod, mod1, fld1, min, max, minmax,
             start, next, done, r_promote, reducedim_init, rand
if isdefined(Base, :rem1)
    import Base: rem1
end
using Base: @pure

# T => BaseType
# f => Number of Bytes reserved for fractional part
abstract type FixedPoint{T <: Integer, f} <: Real end

abstract type ScaledFixedPoint{T <: Integer, f, s} <: Real end
const GenericFixedPoint = Union{FixedPoint, ScaledFixedPoint}

export
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
    scaledual

reinterpret(x::GenericFixedPoint) = x.i
reinterpret(::Type{T}, x::FixedPoint{T,f}) where {T,f} = x.i
reinterpret(::Type{T}, x::ScaledFixedPoint{T,f,s}) where {T,f,s} = x.i

# construction using the (approximate) intended value, i.e., N0f8
*(x::Real, ::Type{X}) where {X<:GenericFixedPoint} = X(x)

# comparison
==(x::T, y::T) where {T <: GenericFixedPoint} = x.i == y.i
 <(x::T, y::T) where {T <: GenericFixedPoint} = x.i  < y.i
<=(x::T, y::T) where {T <: GenericFixedPoint} = x.i <= y.i
"""
    isapprox(x::FixedPoint, y::FixedPoint; rtol=0, atol=max(eps(x), eps(y)))

For FixedPoint numbers, the default criterion is that `x` and `y` differ by no more than `eps`, the separation between adjacent fixed-point numbers.
"""
function isapprox(x::T, y::T; rtol=0, atol=max(eps(x), eps(y))) where {T <: GenericFixedPoint}
    maxdiff = T(atol+rtol*max(abs(x), abs(y)))
    rx, ry, rd = reinterpret(x), reinterpret(y), reinterpret(maxdiff)
    abs(signed(widen1(rx))-signed(widen1(ry))) <= rd
end
function isapprox(x::GenericFixedPoint, y::GenericFixedPoint; rtol=0, atol=max(eps(x), eps(y)))
    isapprox(promote(x, y)...; rtol=rtol, atol=atol)
end

# predicates
isinteger(x::FixedPoint{T,f}) where {T,f} = (x.i&(1<<f-1)) == 0
isinteger(x::ScaledFixedPoint{T,f,s}) where {T,f,s} = (x.i&(1<<f-1)) == 0

# traits
typemax(::Type{T}) where {T <: GenericFixedPoint} = T(typemax(rawtype(T)), 0)
typemin(::Type{T}) where {T <: GenericFixedPoint} = T(typemin(rawtype(T)), 0)
realmin(::Type{T}) where {T <: GenericFixedPoint} = eps(T)
realmax(::Type{T}) where {T <: GenericFixedPoint} = typemax(T)

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

const ShortInts = Union{Int8,UInt8,Int16,UInt16}

floattype(::Type{FixedPoint{T,f}}) where {T <: ShortInts,f} = Float32
floattype(::Type{ScaledFixedPoint{T,f,s}}) where {T <: ShortInts,f,s} = Float32
floattype(::Type{FixedPoint{T,f}}) where {T <: Integer,f} = Float64
floattype(::Type{ScaledFixedPoint{T,f,s}}) where {T <: Integer,f,s} = Float64
floattype(::Type{F}) where {F <: GenericFixedPoint} = floattype(supertype(F))
floattype(x::GenericFixedPoint) = floattype(typeof(x))

nbitsfrac(::Type{FixedPoint{T,f}}) where {T <: Integer,f} = f
nbitsfrac(::Type{ScaledFixedPoint{T,f,s}}) where {T <: Integer,f,s} = f
nbitsfrac(::Type{F}) where {F <: GenericFixedPoint} = nbitsfrac(supertype(F))

rawtype(::Type{FixedPoint{T,f}}) where {T <: Integer,f} = T
rawtype(::Type{ScaledFixedPoint{T,f,s}}) where {T <: Integer,f,s} = T
rawtype(::Type{F}) where {F <: GenericFixedPoint} = rawtype(supertype(F))
rawtype(x::GenericFixedPoint) = rawtype(typeof(x))

# This IOBuffer is used during module definition to generate typealias names
_iotypealias = IOBuffer()

# Printing. These are used to generate type-symbols, so we need them
# before we include any files.
function showtype(io::IO, ::Type{X}) where {X <: GenericFixedPoint}
    print(io, typechar(X))
    f = nbitsfrac(X)
    m = sizeof(X)*8-f-signbits(X)
    print(io, m, 'f', f)
    io
end
function show(io::IO, x::FixedPoint{T,f}) where {T,f}
    showcompact(io, x)
    showtype(io, typeof(x))
end
function show(io::IO, x::ScaledFixedPoint{T,f,s}) where {T,f,s}
    showcompact(io, x)
    showtype(io, typeof(x))
end
const _log2_10 = 3.321928094887362
showcompact(io::IO, x::FixedPoint{T,f}) where {T,f} = show(io, round(convert(Float64,x), ceil(Int,f/_log2_10)))
showcompact(io::IO, x::ScaledFixedPoint{T,f,s}) where {T,f,s} = show(io, round(convert(Float64,x), ceil(Int,f/_log2_10)))

if VERSION >= v"0.7.0-DEV.1790"
    function Base.showarg(io::IO, a::Array{T}, toplevel) where {T<:GenericFixedPoint}
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
r_promote(::typeof(+), x::GenericFixedPoint) = Treduce(x)
r_promote(::typeof(*), x::GenericFixedPoint)= Treduce(x)

reducedim_init(f::typeof(identity),
                              op::typeof(+),
                              A::AbstractArray{T}, region) where {T <: GenericFixedPoint} =
    reducedim_initarray(A, region, zero(Treduce))
reducedim_init(f::typeof(identity),
                              op::typeof(*),
                              A::AbstractArray{T}, region) where {T <: GenericFixedPoint} =
    reducedim_initarray(A, region, oneunit(Treduce))

for f in (:div, :fld, :fld1)
    @eval begin
        $f(x::T, y::T) where {T <: GenericFixedPoint} = $f(reinterpret(x),reinterpret(y))
    end
end
for f in (:rem, :mod, :mod1, :rem1, :min, :max)
    if f === :rem1 && !isdefined(Base, :rem1)
        continue
    end
    @eval begin
        $f(x::T, y::T) where {T <: GenericFixedPoint} = T($f(reinterpret(x),reinterpret(y)),0)
    end
end

# When multiplying by a float, reduce two multiplies to one.
# Particularly useful for arrays.
scaledual(Tdual::Type, x) = oneunit(Tdual), x
scaledual(b::Tdual, x) where {Tdual <: Number} = b, x
scaledual(Tdual::Type, x::Union{T,AbstractArray{T}}) where {T <: GenericFixedPoint} =
    convert(Tdual, 1/oneunit(T)), reinterpret(rawtype(T), x)
scaledual(b::Tdual, x::Union{T,AbstractArray{T}}) where {Tdual <: Number,T <: GenericFixedPoint} =
    convert(Tdual, b/oneunit(T)), reinterpret(rawtype(T), x)

@noinline function throw_converterror(::Type{T}, x) where {T <: GenericFixedPoint}
    n = 2^(8*sizeof(T))
    bitstring = sizeof(T) == 1 ? "an 8-bit" : "a $(8*sizeof(T))-bit"
    io = IOBuffer()
    showcompact(io, typemin(T)); Tmin = String(take!(io))
    showcompact(io, typemax(T)); Tmax = String(take!(io))
    throw(ArgumentError("$T is $bitstring type representing $n values from $Tmin to $Tmax; cannot represent $x"))
end

rand(::Type{T}) where {T <: GenericFixedPoint} = reinterpret(T, rand(rawtype(T)))
rand(::Type{T}, sz::Dims) where {T <: GenericFixedPoint} = reinterpret(T, rand(rawtype(T), sz))

end # module
