module Design

using OffsetArrays

#==============================================================================
Polynomials
==============================================================================#

# DONE
struct VarLinForm{T}
    form ::Dict{Int, T}
end

#DONE
function add!(dest, source, factor=1) end

#DONE
function add_conjugated!(dest, source, factor=1) end

#DONE
function mul!(form, factor) end

#######################################

struct PolynomialForm{N <: Number}
    terms ::OffsetVector{ VarLinForm{N} }
end

function add!(dest, source, factor=1) end

function mul!(poly, factor) end

function mul_by_power!(poly, factor) end

function conjugate(poly) end

function z_conj_diff(poly) end

function matrix_form(poly) end

#######################################

struct NormedPolynomial{N <: Number, FL <: AbstractFloat}
    poly ::PolynomialForm{N}
    r_norm ::FL
end

function add!(dest::NormedPolynomial, source::NormedPolynomial, factor=1) end

function empty!(normed_poly::NormedPolynomial) end

function conjugate(normed_poly::NormedPolynomial, r_contour) end

function re_conjugate!(normed_poly::NormedPolynomial, r_old, r_new) end

function z_conj_diff(normed_poly::NormedPolynomial, r_contour) end

function matrix_form(normed_poly::NormedPolynomial) end

#######################################

struct PlaneLayer{FL <: AbstractFloat}
    E       ::FL
    ν       ::FL
    r_inner ::FL
    r_outer ::FL
    ϕ       ::NormedPolynomial{Complex{FL}, FL}
    z_bar_Φ ::NormedPolynomial{Complex{FL}, FL}
    bar_ψ   ::NormedPolynomial{Complex{FL}, FL}
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

    add!(dest.bar_ψ, source.ϕ,       1/(1+κ2))
    add!(dest.bar_ψ, source.z_bar_Φ, 1/(1+κ2))
    add!(dest.bar_ψ, source.bar_ψ,   1/(1+κ2))
    add!(dest.bar_ψ, dest.ϕ,        -1/(1+κ2))
    add!(dest.bar_ψ, dest.z_bar_Φ,      -1/(1+κ2))

    re_conjugate!(dest.bar_ψ,   dest.r_inner, dest.r_outer)
    re_conjugate!(dest.z_bar_Φ, dest.r_inner, dest.r_outer)
end

#==============================================================================
Ellipticals
==============================================================================#

struct EllipticPraecursor{FL <: AbstractFloat}
    terms :: OffsetVector{ Complex{FL} }
end


end # module Design