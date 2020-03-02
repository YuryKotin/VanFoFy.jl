module LongitudinalShearCase

using ..TypeSynonyms
using ..Symbolic

include("Lshear-fiber.jl")
include("Lshear-cohesive.jl")

struct LongitudinalShear
    fibers     ::Dict{RationalComplex, Fiber}
    cohesive   ::Cohesive
    slae_left  ::Matrix{Float64}
    slae_right ::Vector{Float64}
end

end # module LongitudinalShear