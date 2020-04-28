module PlaneProblems
#=
using ..Input: CellData
using ..TypeSynonyms
using ..Symbolic: NormedPolynomial, add!, z_conj_diff, re_conjugate!


using ..General

abstract type PlaneProblem <: ProblemType end

struct PlaneLayer
    E       ::Float
    ν       ::Float
    r_inner ::Float
    r_outer ::Float
    ϕ       ::NormedPolynomial{Complex{Float}, Float}
    z_bar_Φ ::NormedPolynomial{Complex{Float}, Float}
    bar_ψ   ::NormedPolynomial{Complex{Float}, Float}
end

function couple!(dest::PlaneLayer, source::PlaneLayer) 
    if dest.r_inner != source.r_outer
        error("Layers don't touch")
    end

    E1 = source.E
    ν1 = source.ν
    G1 = E1/(2*(1+ν1))
    κ1 = 3 - 4*ν1

    E2 = dest.E
    ν2 = dest.ν
    κ2 = 3 - 4*ν2
    G2 = E2/(2*(1+ν2))

    empty!(dest.ϕ)
    empty!(dest.z_bar_Φ)
    empty!(dest.bar_ψ)

    add!(dest.ϕ, source.ϕ,       (1+κ1*G2/G1)/(1+κ2))
    add!(dest.ϕ, source.z_bar_Φ, (1-G2/G1)   /(1+κ2))
    add!(dest.ϕ, source.bar_ψ,   (1-G2/G1)   /(1+κ2))

    dest.z_bar_Φ = z_conj_diff(dest.ϕ, dest.r_inner)

    add!(dest.bar_ψ, source.ϕ)
    add!(dest.bar_ψ, source.z_bar_Φ)
    add!(dest.bar_ψ, source.bar_ψ)
    add!(dest.bar_ψ, dest.ϕ, -1)
    add!(dest.bar_ψ, dest.z_bar_Φ, -1)

    re_conjugate!(dest.bar_ψ,   dest.r_inner, dest.r_outer)
    re_conjugate!(dest.z_bar_Φ, dest.r_inner, dest.r_outer)
end

function displacements(layer::PlaneLayer)
    E = layer.E
    ν = layer.ν
    G = E/(2*(1+ν))
    κ = 3 - 4*ν
    
    displ = NormedPolynomial{Complex{Float}, Float}()
    add!(displ, layer.ϕ,        κ/(2G))
    add!(displ, layer.z_bar_Φ, -1/(2G))
    add!(displ, layer.bar_ψ,   -1/(2G))

    return displ
end

function forces(layer::PlaneLayer)
    f = NormedPolynomial{Complex{Float}, Float}()

    add!(f, layer.ϕ)
    add!(f, layer.z_bar_Φ, -1)
    add!(f, layer.bar_ψ,   -1)

    return f
end

#######################################

struct PlaneFiber
    layers ::Vector{PlaneLayer}
end

struct PlaneFiber_Displacements <: Fiber{PlaneProblem}
    fiber ::PlaneFiber
end

struct PlaneFiber_Displacements <: Fiber{PlaneProblem}
    fiber ::PlaneFiber
end



#=
using ..Symbolic: Polynomial, max_abs_index, add!

struct Layer
    E::Float64
    ν::Float64

    ϕ         ::Polynomial
    z_bar_Phi ::Polynomial
    bar_ψ     ::Polynomial
end

function fill_slae!(layer::Layer, slae) 
    ϕ = layer.ϕ
    ψ = layer.ψ

    E = layer.E
    ν = layer.ν
    G = E / (2.0*(1.0+ν))
    κ = 3.0-4.0*ν
    
    max_index = maximum(max_abs_index(ϕ), max_abs_index(ψ))
    displacements = Polynomial(-max_index : max_index)
    stresses      = Polynomial(-max_index : max_index)
    
    add_plain!(displacements, ϕ, κ/G)
    add_z_conj_diff!(displacements, ϕ, -1.0/G)
    add_conj!(displacements, ψ, -1.0/G)

    add_plain!(stresses, ϕ)
    add_z_conj_diff!(stresses, ϕ)
    add_conj!(stresses, ψ)

    # TODO fill slae
end

struct Fiber
    layers ::Vector{Layer}
end

function fill_slae!(fiber::Fiber, slae) 
    fill_slae!(fiber.layers[end], slae) 
end

struct Cohesive

end

struct PlaneProblem
    fibers    ::Dict{Pole, Fiber}
    cohesive  ::Cohesive

    slae_rows ::Dict{Pole, UnitRange{Int64}}
    slae_left ::Matrix{Float64}
    slae_right::Vector{Float64}
end

function PlaneProblem(cell::CellData) 
    # TODO
end

function fill_slae!(cohesive::Cohesive, pole::Pole, slae) 
    # TODO
end

function fill_slae!(plane_problem::PlaneProblem)
    slae_left = plane_problem.slae_left
    rows = plane_problem.slae_rows

    for pole in keys(plane_problem.fibers)
        slae_view = view(slae_left[rows[pole]])
        fill_slae!(plane_problem.fibers[pole], slae_view)
        fill_slae!(plane_problem.cohesive, pole, slae_view)
    end

    # TODO averaging
end

=#
=#
end # module PlaneProblems