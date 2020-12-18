module General

using ..Types: RationalComplex

struct Contour
    center    ::RationalComplex
    radius    ::Float64
    max_power ::Int
end

abstract type ProblemType end

abstract type Fiber{P <: ProblemType} end

abstract type Cohesive{P <: ProblemType} end

end # module General