module TypeSynonyms

using OffsetArrays

# OffsetArray{ComplexF64}
const ComplexOffsetMatrix = OffsetArray{Complex{Float64},2, Array{Complex{Float64}, 2}}
export ComplexOffsetMatrix

# OffsetVector{ComplexF64}
const ComplexOffsetVector = OffsetArray{Complex{Float64},1,Array{Complex{Float64}, 1}}
export ComplexOffsetVector

# OffsetVector{Int64}
const IntOffsetVector = OffsetArray{Int64, 1, Array{Int64, 1}}
export IntOffsetVector

const RationalComplex = Complex{Rational{Int64}}
export RationalComplex

const RComplex2IntDict = Dict{RationalComplex, Int64}
export RComplex2IntDict

const RComplex2OffsetMatrix = Dict{RationalComplex,ComplexOffsetMatrix}
export RComplex2OffsetMatrix

const Variable = Int64
export Variable

const Coefficient = ComplexF64
export Coefficient

end # module TypeSynonims