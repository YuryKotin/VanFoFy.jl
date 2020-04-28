module Design
#=
using OffsetArrays

const Float = Float64

#==============================================================================
Polynomials
==============================================================================#

# DONE
struct VarLinForm{T}
    form ::Dict{Int, T}
end

#DONE
function add!(dest<:VarLinForm, source<:VarLinForm, factor=1) end

#DONE
function add_conjugated!(dest<:VarLinForm, source<:VarLinForm, factor=1) end

#DONE
function mul!(form<:VarLinForm, factor) end

#######################################

struct PolynomialForm{N <: Number}
    terms ::OffsetVector{ VarLinForm{N} }
end

function add!(dest<:PolynomialForm, source<:PolynomialForm, factor=1) end

function mul!(poly<:PolynomialForm, factor) end

function mul_by_power!(poly<:PolynomialForm, factor) end

function conjugate(poly<:PolynomialForm) end

function z_conj_diff(poly<:PolynomialForm) end

function matrix_form(poly<:PolynomialForm) end

#######################################

struct NormedPolynomial{N <: Number, FL <: AbstractFloat}
    poly ::PolynomialForm{N}
    r_norm ::FL
end

function add!(dest<:NormedPolynomial, source<:NormedPolynomial, factor=1) end

function empty!(normed_poly<:NormedPolynomial) end

function conjugate(normed_poly<:NormedPolynomial, r_contour) end

function re_conjugate!(normed_poly<:NormedPolynomial, r_old, r_new) end

function z_conj_diff(normed_poly<:NormedPolynomial, r_contour) end

function matrix_form(normed_poly<:NormedPolynomial) end

#######################################

struct PlaneLayer
    E       ::Float
    ν       ::Float
    r_inner ::Float
    r_outer ::Float
    ϕ       ::NormedPolynomial{Complex{Float}, Float}
    z_bar_Φ ::NormedPolynomial{Complex{Float}, Float}
    bar_ψ   ::NormedPolynomial{Complex{Float}, Float}
end

function couple!(dest::PlaneLayer, source::PlaneLayer) end

function displacements(layer::PlaneLayer) end

function forces(layer::PlaneLayer) end

#######################################

struct PlaneFiber
    layers ::Vector{PlaneLayer{Float}}
    var_indices ::UnitRange{Int}
    max_power ::Int
end

displacements(fiber::PlaneFiber) = displacements(fiber.layers[end])

forces(fiber::PlaneFiber) = forces(fiber.layers[end])

#######################################

struct PlaneCohesive end

function displacements(cohesive::PlaneCohesive) end

function forces(cohesive::PlaneCohesive) end

#######################################

struct PlaneCell
    fibers ::Vector{ PlaneFiber{Float} }
    cohesive ::PlaneCohesive{Float}
    slae_A ::Matrix{Float}
end

function fill_slae_from_fibers!(cell::PlaneCell)

end

###############################################################################
###############################################################################
###############################################################################

struct Contour
    center    ::Complex{ Rational{Int} }
    radius    ::Float
    max_power ::Int
end

abstract type ProblemType end

abstract type Fiber{P <: ProblemType} end

abstract type Cohesive{P <: ProblemType} end

abstract type SubMatrixType{P <: ProblemType} end

struct SubMatrix_Zero{P <: ProblemType} <: SubMatrixType{P}
    contour ::Contour
end

struct SubMatrix_Fiber{P <: ProblemType} <: SubMatrixType{P}
    contour ::Contour
    fiber   ::Fiber{P}
end

struct SubMatrix_Cohesive{P <: ProblemType} <: SubMatrixType{P}
    contour  ::Contour
    cohesive ::Cohesive{P}
end

abstract type SubMatrix_Averaging{P <: ProblemType} <: SubMatrixType{P} end

struct ReferenceMatrix{P <: ProblemType}
    matrix ::Matrix{ SubMatrixType{P} }
end




#==============================================================================
Ellipticals
==============================================================================#

struct EllipticPraecursor{FL <: AbstractFloat}
    terms :: OffsetVector{ Complex{FL} }
end

=#
end # module Design