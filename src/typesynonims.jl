# OffsetArray{ComplexF64}
const ComplexOffsetMatrix = OffsetArray{Complex{Float64},2, Array{Complex{Float64}, 2}}

# OffsetVector{ComplexF64}
const ComplexOffsetVector = OffsetArray{Complex{Float64},1,Array{Complex{Float64}, 1}}

# OffsetVector{Int64}
const IntOffsetVector = OffsetArray{Int64, 1, Array{Int64, 1}}

const RationalComplex = Complex{Rational{Int64}}

const WeiDeltas = Dict{RationalComplex, Int64, Bool}

const RComplexToOffsetMatrix = Dict{RationalComplex,ComplexOffsetMatrix}
