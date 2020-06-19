module Types

using OffsetArrays

const Float = Float64
export Float

# OffsetArray{ComplexF64}
const ComplexOffsetMatrix = OffsetArray{Complex{Float},2, Array{Complex{Float}, 2}}
export ComplexOffsetMatrix

# OffsetVector{ComplexF64}
const ComplexOffsetVector = OffsetArray{Complex{Float},1,Array{Complex{Float}, 1}}
export ComplexOffsetVector

# OffsetVector{Int}
const IntOffsetVector = OffsetArray{Int, 1, Array{Int, 1}}
export IntOffsetVector

const RationalComplex = Complex{Rational{Int}}
export RationalComplex

const RComplex2IntDict = Dict{RationalComplex, Int}
export RComplex2IntDict

const RComplex2OffsetMatrix = Dict{RationalComplex,ComplexOffsetMatrix}
export RComplex2OffsetMatrix

const Variable = Int
export Variable

const Coefficient = ComplexF64
export Coefficient

const Pole = RationalComplex
export Pole

const Power = Int
export Power

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

Base.lastindex(bv::BoundedVector) = Base.lastindex(bv.vector)

Base.fill!(bv::BoundedVector{T}, val::T) where T = Base.fill!(bv.vector, val)

###############################################################################

mutable struct CashedVector{N <: Number}
    array ::OffsetVector{N}
    last_cashed ::Int
end

function CashedVector{N}(range ::UnitRange{Int}) where {N <: Number}
    array = OffsetVector{N}(undef, range)
    last_cashed = firstindex(array) - 1
    CashedVector(array, last_cashed)
end

function Base.getindex(v::CashedVector{N} , index::Int) where {N}
    if index <= v.last_cashed
       return v.array[index] 
    else
        error("Access to undefined values")
    end
end

function Base.setindex!(v::CashedVector{N}, val::N, index::Int) where {N}
    if index == v.last_cashed+1
        v.array[index] = val
        v.last_cashed += 1
    elseif index <= v.last_cashed
        v.array[index] = val
    else
        error("Access to undefined values")
    end
end

last_cashed(v::CashedVector{N}) where N = v.last_cashed

not_cashed(v::CashedVector{N}, ind::Int) where N = (ind > v.last_cashed)

###############################################################################

end # module TypeSynonims