# Scaled{T,s,r} maps Ints with a scale factor of s s.t. the domain is
# [-s*2^(b-1,...,-s,0,s,...s*(2^(b-1)-1))] where b is number of bits
using Distributions
struct Scaled{T<:Signed,s,r <: RoundingScheme} <: ScaledFixedPoint{T,s,r}
    i::T
    # constructor for manipulating the representation;
    # selected by passing an extra dummy argument
    Scaled{T,s,r}(i::Integer,_) where {T <:Integer,s,r<: RoundingScheme} = new{T,s,r}(i%T)
    Scaled{T,s,r}(x) where {T,s,r <: RoundingScheme} =
        convert(Scaled{T,s,r}, x)
end
typechar(::Type{X}) where {X <: Scaled} = 'S'
# datatype comparison
# TODO: Not needed until we start working with situations with different T values
# eq(x::Type{T}, y::Type{T}) where {T <: Signed} =   sizeof(x)==sizeof(y)
# less(x::Type{T}, y::Type{T}) where {T <: Signed} = sizeof(x) <sizeof(y)
# leq(x::Type{T}, y::Type{T}) where {T <: Signed} =  sizeof(x)<=sizeof(y)

# basic operators
-(x::Scaled{T,s,r}) where {T,s,r} = Scaled{T,s,r}(-x.i,true) # negation
abs(x::Scaled{T,s,r}) where {T,s,r} = Scaled{T,s,r}(abs(x.i),true) # absolute value
# TODO: Non-saturated adds require domain to be Scaled{max(b1,b2)+1,s,r}, +1 isn't possible
# Saturated Add when T and s are the same
+(x::Scaled{T,s,r}, y::Scaled{T,s,r}) where {T <:Signed,s,r<: RoundingScheme} =
    Scaled{T,s,r}(clamp(x.i+y.i,-(1<<maxf(T)),(1<<maxf(T))-1))
-(x::Scaled{T,s,r}, y::Scaled{T,s,r}) where {T <:Signed,s,r<: RoundingScheme} =
    Scaled{T,s,r}(clamp(x.i-y.i,-(1<<maxf(T)),(1<<maxf(T))-1))
# Multiplication with the assumption that they are have the same T
# TODO: Arbitrary T: this requires T to include Int24 and Int40 or arbitrary sizes
*(x::Scaled{T,s1,r}, y::Scaled{T,s2,r}) where {T <:Signed,s1,s2,r<: RoundingScheme} =
    Scaled{widen1(T),s1+s2,r}(Base.widemul(x.i,y.i))
# TODO: Dot Product
# TODO: Matrix Multiplication

# conversions
convert(::Type{Scaled{T,s,r}}, x::Integer) where {T,s,r} = Scaled{T,s,r}(x%T,0)
function convert(::Type{Scaled{T,s,r}}, x::AbstractFloat) where {T,s,r <: RoundingScheme}
    if r==Randomized
        return Scaled{T,s,r}(convert(T,floor(Base.widemul(x,(1./s)) + rand(Uniform(0,1)))),0)
    elseif r==NearestNeighbor
        return Scaled{T,s,r}(convert(T,floor(Base.widemul(x,(1./s)) + 0.5)),0)
    end
end
convert(::Type{TF}, x::Scaled{T,s}) where {TF <: AbstractFloat,T,s} = float(x.i)

#promotions
promote_rule(ft::Type{Scaled{T,s,r}}, ::Type{TI}) where {T <: Integer,s,r <: RoundingScheme,TI <: Integer} = Scaled{T,s,r}
promote_rule(::Type{Scaled{T,s,r}}, ::Type{TF}) where {T <: Integer,s,r <: RoundingScheme,TF <: AbstractFloat} = TF
promote_rule(::Type{Scaled{T,s,r}}, ::Type{Rational{TR}}) where {T <: Integer,s,r <: RoundingScheme,TR} = Rational{TR}
