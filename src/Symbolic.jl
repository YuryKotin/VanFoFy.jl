module Symbolic

using OffsetArrays

using ..TypeSynonyms

abstract type GeneralFunction end

differentiate!(f::GeneralFunction) = 
error("Differentiation is not defined in the concrete type")

abstract type GeneralPolynomial{T <: Number} end

function add!(dest::GP, source::GP, factor::T) where GP <: GeneralPolynomial{T} where T
    if matched_powers(dest, source)
        for p in powers(dest)
            dest[p] += source[p] * factor
        end
    else
        error("Polynomials powers don't match")
    end
end

###############################################################################

mutable struct Binomial{T} <: GeneralPolynomial{T}
    termA ::T
    termB ::T
    powA  ::Int
    powB  ::Int
end

function Base.getindex(b::Binomial{T}, key::Int) where T
    if key == b.powA 
        return b.termA 
    end
    if key == b.powB 
        return b.termB 
    end
    return zero(T)
end

function Base.setindex!(b::Binomial{T}, value::T, key::Int) where T
    if key == b.powA 
        b.termA = value
        return 
    end
    if key == b.powB 
        b.termB = value 
        return
    end
    error("Key doesn't match Binomial powers")
end

powers(b::Binomial) = (b.powA, b.powB)

matched_powers(b1::Binomial, b2::Binomial) = 
(b1.powA == b2.powA) && (b1.powB == b2.powB)

###############################################################################

struct Polynomial{T <: Number}
    terms ::Dict{Int, T}
end

function Polynomial{T}() where T <: Number 
    Polynomial(Dict{Int, T}())
end

Base.getindex(poly::Polynomial{T}, key::Int) where T <: Number = 
get(poly.terms, key, zero(T))

Base.setindex!(poly::Polynomial{T}, value::T, key::Int) where T <: Number =
Base.setindex!(poly.terms, value, key)

###

eachpower(p::Polynomial) = keys(p.terms)

function max_abs_index(poly::Polynomial)
    powers = eachpower(poly)
    max_power = maximum(powers)
    min_power = minimum(powers)
    maximum(abs(max_power), abs(min_power))
end

function add!(dest::Polynomial{T}, source::Polynomial{T}, factor::T) where T <: Number
    for power in eachpower(source)
        dest[power] += source[power] * factor
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

struct VarLinForm
    form ::Dict{Variable, Coefficient}
end

VarLinForm() = VarLinForm(Dict{Variable, Coefficient}())

Base.getindex(vlf::VarLinForm, key::Variable) = get(vlf.dict, key, 0.0im)

Base.setindex!(vlf::VarLinForm, val::Coefficient, key::Variable) = setindex!(vlf.dict, val, key)

variables(vlf::VarLinForm) = keys(vlf.form)

function add!(dest::VarLinForm, source::VarLinForm, factor=1)
    for var in variables(source)
        dest[var] += source[var] * factor
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