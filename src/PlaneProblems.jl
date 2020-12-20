module PlaneProblems

using ..Types: RationalComplex
using ..VarPolyForms: VarPolyForm
using ..Input: FiberData, LayerData

include("PlaneFiber.jl")
include("PlaneCohesive.jl")

struct PlaneProblem
    inclusions               :: Vector{PlaneInclusion}
    cohesive                 :: PlaneCohesive
end

end # module PlaneProblems