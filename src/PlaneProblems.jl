module PlaneProblems

include("PlaneFiber.jl")
include("PlaneCohesive.jl")

struct PlaneProblem
    doubly_periodic_function :: DoublyPeriodicFunction
    inclusions               :: Vector{PlaneInclusion}
    cohesive                 :: PlaneCohesive
end

end # module PlaneProblems