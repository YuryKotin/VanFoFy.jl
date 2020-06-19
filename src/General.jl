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
#=
abstract type SubMatrixType{P <: ProblemType} end
export SubMatrixType

struct SubMatrix_Zero{P <: ProblemType} <: SubMatrixType{P}
    contour ::Contour
end
export SubMatrix_Zero

struct SubMatrix_Fiber{P <: ProblemType} <: SubMatrixType{P}
    contour ::Contour
    fiber   ::Fiber{P}
end
export SubMatrix_Fiber

struct SubMatrix_Cohesive{P <: ProblemType} <: SubMatrixType{P}
    contour  ::Contour
    cohesive ::Cohesive{P}
end
export SubMatrix_Cohesive

abstract type SubMatrix_Averaging{P <: ProblemType} <: SubMatrixType{P} end
export SubMatrix_Averaging

struct ReferenceMatrix{P <: ProblemType}
    matrix ::Matrix{ SubMatrixType{P} }
end
export ReferenceMatrix

=#
end # module General