module Symbolic

using ..Types: RationalComplex, ComplexOffsetVector, raw_complex
using ..Ellipticals: EllipticPraecursor, getindex, QSpecial
using ..SpecialWeierstrass: Weierstrass

using OffsetArrays

###############################################################################

struct BoundedVector{T <: Number}
    vector ::OffsetArray{T, 1, Array{T, 1}}
end

function BoundedVector{T}(indices::UnitRange{Int}) where T <: Number
    vector = OffsetVector{T}(undef, indices)
    BoundedVector{T}(vector)
end

function Base.getindex(bv::BoundedVector{T}, key::Int) where T
    if key in axes(bv.vector, 1)
        @inbounds return bv.vector[key]
    else
        return zero(T)
    end
end

function Base.setindex!(bv::BoundedVector{T}, val::T, key::Int) where T
    if key in axes(bv.vector, 1)
        @inbounds bv.vector[key] = val
    end
end

lastindex(bv::BoundedVector{T}) where T = Base.lastindex(bv.vector)

Base.fill!(bv::BoundedVector{T}, val::T) where T = Base.fill!(bv.vector, val)

include("EllipticalForm.jl")
include("PolynomialForm.jl")

end # module Symbolic