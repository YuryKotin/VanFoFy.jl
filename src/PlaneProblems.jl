module PlaneProblems

using ..Input: CellData
using ..TypeSynonyms

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

end # module PlaneProblems