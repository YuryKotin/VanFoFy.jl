module TypeSynonyms

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

end # module TypeSynonims