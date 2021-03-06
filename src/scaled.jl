# Scaled{T,f,s,r} maps Ints with a scale factor of s s.t. the domain is
# [-s*2^(f-1,...,-s,0,s,...s*(2^(f-1)-1))] where f is number of bits
using Distributions
struct Scaled{T<:Signed,f,s,r <: RoundingScheme} <: ScaledFixedPoint{T,f,s,r}
    i::T
    # constructor for manipulating the representation;
    # selected by passing an extra dummy argument
    Scaled{T,f,s,r}(i::Integer,_) where {T <:Integer,f,s,r<: RoundingScheme} = new{T,f,s,r}(i%T)
    Scaled{T,f,s,r}(x::Scaled) where {T <:Integer,f,s,r<: RoundingScheme} =
        new{T,f,s,r}(float(x.i))
    Scaled{T,f,s,r}(x::A) where {T <:Integer,f,s,r<: RoundingScheme, A<:Array} =
        convert(Array{Scaled{T,f,s,r}},x)
    Scaled{T,f,s,r}(x) where {T,f,s,r <: RoundingScheme} =
        convert(Scaled{T,f,s,r}, x)
end
export
    eq,
    le,
    leq,
    up,
    sat_add,
    subdomain,
    get_T,
    get_f,
    get_s,
    get_r

reinterpret(::Type{Scaled{T,f,s,r}},x::T) where {T <: Signed ,f,s,r} = Scaled{T,f,s,r}(x,0)
zero(::Type{Scaled{T,f,s,r}}) where {T <: Signed ,f,s,r} = Scaled{T,f,s,r}(zero(T),0)
function oneunit(::Type{T}) where {T <: Scaled}
    T(typemax(rawtype(T)) >> (8*sizeof(T)-nbitsfrac(T)), 0)
end
one(::Type{T}) where {T <: Scaled} = oneunit(T)
zero(x::Scaled) = zero(typeof(x))
oneunit(x::Scaled) =  one(typeof(x))
one(x::Scaled) = oneunit(x)

typechar(::Type{X}) where {X <: Scaled} = 'S'
signbits(::Type{X}) where {X <: Scaled} = 1
# datatype comparison
eq(x::Type{T1}, y::Type{T2}) where {T1,T2 <: Signed} =   sizeof(x)==sizeof(y)
le(x::Type{T1}, y::Type{T2}) where {T1,T2 <: Signed} = sizeof(x) <sizeof(y)
leq(x::Type{T1}, y::Type{T2}) where {T1,T2 <: Signed} =  sizeof(x)<=sizeof(y)

# datatype upgrade when necessary
up(x::Type{T},b1::I,b2::I) where {T <: Signed,I <: Integer} = b1+b2 > maxf(x) ? widen1(x) : x

# basic operators
-(x::Scaled{T,f,s,r}) where {T,f,s,r} = Scaled{T,f,s,r}(-x.i,0) # negation
abs(x::Scaled{T,f,s,r}) where {T,f,s,r} = Scaled{T,f,s,r}(abs(x.i),0) # absolute value
# Saturated Add when T,f,s,r are the same
function sat_add(x::Scaled{T,f,s,r}, y::Scaled{T,f,s,r}) where {T <:Signed,f,s,r<: RoundingScheme}
    Scaled{T,f,s,r}(clamp(Int32(x.i)+Int32(y.i),-(1<<f),(1<<f-1)),0)
end
# Exact adds: only when f's are different
+(x::Scaled{T1,f1,s,r}, y::Scaled{T2,f2,s,r}) where {T1,T2 <:Signed,f1,f2,s,r<: RoundingScheme} =
    Scaled{leq(T1,T2) ? up(T2,f1,f2) : up(T1,f1,f2),max(f1,f2)+2,s,r}(Int64(x.i)+Int64(y.i),0)
-(x::Scaled{T1,f1,s,r}, y::Scaled{T2,f2,s,r}) where {T1,T2 <:Signed,f1,f2,s,r<: RoundingScheme} =
    Scaled{leq(T1,T2) ? up(T2,f1,f2) : up(T1,f1,f2),max(f1,f2)+2,s,r}(Int64(x.i)-Int64(y.i),0)


# Multiplication
*(x::Scaled{T1,f1,s1,r}, y::Scaled{T2,f2,s2,r}) where {T1,T2 <:Signed,f1,s1,f2,s2,r<: RoundingScheme} =
    float(Scaled{leq(T1,T2) ? up(T2,f1,f2) : up(T1,f1,f2),f1+f2+1,s1*s2,r}(Base.widemul(x.i,y.i),0))

# Type Multiplication
*(x::Type{Scaled{T1,f1,s1,r}}, y::Type{Scaled{T2,f2,s2,r}}) where {T1,T2 <:Signed,f1,s1,f2,s2,r<: RoundingScheme} =
    Scaled{leq(T1,T2) ? up(T2,f1,f2) : up(T1,f1,f2),f1+f2+1,s1*s2,r}

# extended functions
# Whether a Scaled is within the same subdomain as another Scaled
subdomain(x::Type{Scaled{T1,f1,s1,r}}, y::Type{Scaled{T2,f2,s2,r}},tol::Real) where {T1,T2 <:Signed,f1,s1,f2,s2,r<: RoundingScheme} =
    isinteger(s1/s2,tol) && log2(s1/s2) <= f2-f1+1

# conversions
convert(::Type{Scaled{T,f,s,r}}, x::Integer) where {T,f,s,r <: Exact} = Scaled{T,f,s,r}(T(x/s),0)
convert(::Type{Scaled{T,f,s,r}}, x::Integer) where {T,f,s,r <: Saturated} = Scaled{T,f,s,r}(T(clamp(x/s,-(1<<f),(1<<f-1))),0)

function convert(::Type{Scaled{T,f,s,r}}, x::AbstractFloat) where {T,f,s,r <: RoundingScheme}
    if r==ExactAndRandomized
        return Scaled{T,f,s,r}(convert(T,floor(Base.widemul(x,(1./s)) + rand(Uniform(0,1)))),0)
    elseif r==ExactAndNearestNeighbor
        return Scaled{T,f,s,r}(convert(T,floor(Base.widemul(x,(1./s)) + 0.5)),0)
    elseif r==SatAndRandomized
        return Scaled{T,f,s,r}(convert(T,clamp(floor(Base.widemul(x,(1./s)) + rand(Uniform(0,1))),-(1<<f),(1<<f-1))),0)
    elseif r==SatAndNearestNeighbor
        return Scaled{T,f,s,r}(convert(T,clamp(floor(Base.widemul(x,(1./s)) + 0.5),-(1<<f),(1<<f-1))),0)
    end
end
convert(::Type{TF}, x::Scaled{T,f,s,r}) where {TF <: AbstractFloat,T,f,s,r} = Float64(x.i)*s

# get methods
get_T(::Type{Scaled{T,f,s,r}}) where {T,f,s,r} = T
get_f(::Type{Scaled{T,f,s,r}}) where {T,f,s,r} = f
get_s(::Type{Scaled{T,f,s,r}}) where {T,f,s,r} = s
get_r(::Type{Scaled{T,f,s,r}}) where {T,f,s,r} = r

#promotions
promote_rule(ft::Type{Scaled{T,f,s,r}}, ::Type{Scaled}) where {T <: Integer,f,s,r <: RoundingScheme} = Scaled{T,f,s,r}
promote_rule(ft::Type{Scaled{T,f,s,r}}, ::Type{TI}) where {T <: Integer,f,s,r <: RoundingScheme,TI <: Integer} = Scaled{T,f,s,r}
promote_rule(::Type{Scaled{T,f,s,r}}, ::Type{TF}) where {T <: Integer,f,s,r <: RoundingScheme,TF <: AbstractFloat} = TF
promote_rule(::Type{Scaled{T,f,s,r}}, ::Type{Rational{TR}}) where {T <: Integer,f,s,r <: RoundingScheme,TR} = Rational{TR}
