module PlaneProblems

using ..Types: RationalComplex
using ..VarPolyForms: VarPolyForm, VarPolyFormBox
using ..Input: FiberData, LayerData, InclusionData

include("PlaneFiber.jl")
include("PlaneCohesive.jl")

struct PlaneProblem
    inclusions               :: Vector{PlaneInclusion}
    cohesive                 :: PlaneCohesive
end

end # module PlaneProblems