module Symbolic

using ..Types: RationalComplex, ComplexOffsetVector
using ..Ellipticals: EllipticPraecursor, getindex, raw_complex, Weierstrass, QSpecial

using OffsetArrays

include("EllipticalForm.jl")
include("PolynomialForm.jl")

end # module Symbolic