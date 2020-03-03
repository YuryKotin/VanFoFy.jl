module Symbolic

using OffsetArrays

using ..TypeSynonyms

abstract type GeneralFunction end

differentiate!(f::GeneralFunction) = 
    error("Differentiation is not defined in the concrete type")

###############################################################################

struct VarLinForm
    form ::Dict{Variable, Coefficient}
end

VarLinForm() = VarLinForm(Dict{Variable, Coefficient}())

Base.getindex(vlf::VarLinForm, key::Variable) = get(vlf.dict, key, 0.0im)

Base.setindex!(vlf::VarLinForm, val::Coefficient, key::Variable) = setindex!(vlf.dict, val, key)

###############################################################################

struct Polynomial <: GeneralFunction
    terms ::OffsetVector{VarLinForm}
end

function Polynomial(ind_range) 
    Polynomial(
        OffsetVector{VarLinForm}([VarLinForm() for i in ind_range], ind_range)
        )
end

function max_abs_index(poly::Polynomial)
    maximum(abs(firstindex(poly.terms)), abs(lastindex(poly.terms)))
end

function add_plain!(dest::Polynomial, source::Polynomial, factor::Float64=1.0)
    for i in eachindex(source.terms)
        for v in keys(source[i])
            dest.terms[i][v] += source.terms[i][v] * factor
        end
    end
end

function add_conj!(dest::Polynomial, source::Polynomial, factor::Float64=1.0)
    for i in eachindex(source.terms)
        for v in keys(source[i])
            dest.terms[-i][v] += source.terms[i][v] * factor
        end
    end
end

function add_z_conj_diff!(dest::Polynomial, source::Polynomial, factor::Float64=1.0)
    for i in eachindex(source.terms)
        for v in keys(source[i])
            dest.terms[2-i][v] += i * conj(source.terms[i][v]) * factor
        end
    end
end

###############################################################################

"""
  (z/r)^n

power -  n,

Norming factor ''r'' is stored elsewhere.
"""
mutable struct ZPowN <: GeneralFunction
    power  ::Int64
end

function differentiate!(f::ZPowN, r_norm::Float64) 
    # [(z-a)^n/r^n]' = (n/r) z^{n-1}/r^{n-1}
    multiplier = f.power / r_norm
    f.power -= 1
    return multiplier
end

###############################################################################

struct Polynomial <: GeneralFunction
    terms ::Dict{Int64, ComplexF64}
end

###############################################################################

"""
  r^{n+2} ℘^{(n)}(z)/(n+1)!

deriv -  n.

Norming factor ''r'' is stored elsewhere.
"""
mutable struct WpDerivative <: GeneralFunction
    deriv::Int64
end

function differentiate!(f::WpDerivative, r_norm::Float64) 
    # [r_k^{n+2} ℘^{(n)}(z-ζ_k)/(n+1)!]' =
    # [(n+2) / r_k] *
    # r_k^{n+3} ℘^{(n+1)}(z-ζ_k)/(n+2)!
    multiplier = (f.deriv+2) / r_norm
    f.deriv += 1
    return multiplier
end

###############################################################################

"""
  r^{n+2} Q^{(n)}(z)/(n+1)!

deriv -  n.

Norming factor ''r'' is stored elsewhere.
"""
mutable struct QSpecial <: GeneralFunction
    deriv::Int64
end

function differentiate!(f::QSpecial, r_norm::Float64) 
    # [r_k^{n+2} Q^{(n)}(z-ζ_k)/(n+1)!]' =
    # [(n+2) / r_k] *
    # r_k^{n+3} Q^{(n+1)}(z-ζ_k)/(n+2)!
    multiplier = (f.deriv+2) / r_norm
    f.deriv += 1
    return multiplier
end

###############################################################################

mutable struct Variable
    number      ::Int64
    is_computed ::Bool
    value       ::Float64
end


mutable struct Term{F <: GeneralFunction}
    variable    ::Variable
    coefficient ::ComplexF64
    functional  ::F
end

function differentiate!(term::Term)
    multiplier = differentiate!(term.functional)
    term.coefficient *= multiplier
    nothing
end

###############################################################################

abstract type FunctionalForm end

struct LinearForm{F<:GeneralFunction} <: FunctionalForm
    terms          ::Vector{Term{F}}
    shift_rational ::RationalComplex
    shift_raw      ::ComplexF64
    norm_factor    ::Float64
    conjugated     ::Bool
    mul_by_z       ::Bool
end

function conjugate!(form::LinearForm)
    for term in form.terms
        term.coefficient = conj(term.coefficient)
    end
    form.conjugated = !form.conjugated
    nothing
end

function mul!(form::LinearForm, factor::Float64)
    for term in form.terms
        term.coefficient *= factor
    end
    nothing
end

struct SolutionForm <:FunctionalForm
    terms ::Vector{LinearForm}
end

function differentiate!(form <:FunctionalForm)
    map(differentiate!, form.terms)
    nothing
end

function conjugate!(form <:FunctionalForm)
    map(conjugate!, form.terms)
    nothing
end

function mul!(form <:FunctionalForm, factor::Float64)
    f!(t) = mul!(t, factor)
    map(f!, form.terms)
    nothing
end

end # module Symbolic