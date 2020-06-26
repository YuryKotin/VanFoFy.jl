module PlaneProblems

using ..FunctionalTerms: PolynomialTerm
using ..Types: VarLinForm, IntOffsetVector

using OffsetArrays

struct PlaneLayer
    E       ::Float64
    ν       ::Float64
    r_inner ::Float64
    r_outer ::Float64
    ϕ       ::VarLinForm{PolynomialTerm}
    z_bar_Φ ::VarLinForm{PolynomialTerm}
    bar_ψ   ::VarLinForm{PolynomialTerm}
end

function PlaneLayer(E::Float64, 
                    ν::Float64, 
                    r_inner::Float64, 
                    r_outer::Float64, 
                    r_fiber::Float64, 
                    var_indices::IntOffsetVector) 

    n_vars = size(var_indices, 1)
    if n_vars % 2 == 1
        error("Odd number of variables complex parts")
    end

    max_power = n_vars - 1
    bottom_var = firstindex(var_indices)
    top_var    = lastindex(var_indices)

    ϕ = VarLinForm(
                    OffsetVector([PolynomialTerm(-max_power:max_power+2, r_fiber) 
                                for v in var_indices], 
                                bottom_var : top_var)
        )
    
    ψ = VarLinForm(
                    OffsetVector([PolynomialTerm(-max_power:max_power+2, r_fiber) 
                                for v in var_indices], 
                                bottom_var : top_var)
        )       

    for v in bottom_var : top_var
        p = (v - bottom_var) ÷ 2
        ϕ_poly = ϕ[v]
        ψ_poly = ψ[v]
        if (v - bottom_var) % 2 == 0
            ϕ_poly[p] = 1.0+0.0im
            ψ_poly[p] = 1.0+0.0im
        else
            ϕ_poly[p] = 1.0im
            ψ_poly[p] = 1.0im
        end
    end

    if r_inner != 0.0

    end
end

function PlaneLayer(E, ν, r_outer, prev_layer) end

#=
using ..Input: CellData
using ..Types
using ..FunctionalTerms: NormedPolynomial, add!, z_conj_diff
using ..General: ProblemType, Fiber

abstract type PlaneProblem <: ProblemType end

struct PlaneLayer
    E       ::Float64
    ν       ::Float64
    r_inner ::Float64
    r_outer ::Float64
    ϕ       ::NormedPolynomial
    z_bar_Φ ::NormedPolynomial
    bar_ψ   ::NormedPolynomial
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

    dest.z_bar_Φ = z_conj_diff(dest.ϕ, dest.r_inner, 1.0)

    add!(dest.bar_ψ, source.ϕ,       +1.0)
    add!(dest.bar_ψ, source.z_bar_Φ, +1.0)
    add!(dest.bar_ψ, source.bar_ψ,   +1.0)
    add!(dest.bar_ψ, dest.ϕ,         -1.0)
    add!(dest.bar_ψ, dest.z_bar_Φ,   -1.0)
end

function displacements(layer::PlaneLayer)
    E = layer.E
    ν = layer.ν
    G = E/(2*(1+ν))
    κ = 3 - 4*ν
    
    displ = similar(layer.ϕ)
    add!(displ, layer.ϕ,        κ/(2G))
    add!(displ, layer.z_bar_Φ, -1/(2G))
    add!(displ, layer.bar_ψ,   -1/(2G))

    return displ
end

function forces(layer::PlaneLayer)
    force = similar(layer.ϕ)

    add!(force, layer.ϕ,       +1.0)
    add!(force, layer.z_bar_Φ, -1.0)
    add!(force, layer.bar_ψ,   -1.0)

    return force
end

#######################################

struct PlaneFiber
    layers ::Vector{PlaneLayer}
end

struct PlaneFiber_Displacements <: Fiber{PlaneProblem}
    fiber ::PlaneFiber
end

struct PlaneFiber_Forces <: Fiber{PlaneProblem}
    fiber ::PlaneFiber
end

=#

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