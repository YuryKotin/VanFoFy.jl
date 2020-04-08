module Symbolic

using OffsetArrays

using ..TypeSynonyms


struct VarLinForm
    form ::Dict{Variable, Coefficient}
end

VarLinForm() = VarLinForm(Dict{Variable, Coefficient}())

Base.getindex(vlf::VarLinForm, key::Variable) = get(vlf.form, key, 0.0im)

Base.setindex!(vlf::VarLinForm, val::Coefficient, key::Variable) = setindex!(vlf.form, val, key)

variables(vlf::VarLinForm) = keys(vlf.form)

function add!(dest::VarLinForm, source::VarLinForm, factor=1)
    for var in variables(source)
        dest[var] += source[var] * factor
    end
end

function add_conjugated!(dest::VarLinForm, source::VarLinForm, factor=1)
    for var in variables(source)
        dest[var] += conj(source[var] * factor)
    end
end


function mul!(form::VarLinForm, factor)
    for var in variables(form)
        form[var] *= factor
    end
end

###############################################################################

struct PolynomialForm
    terms ::OffsetVector{VarLinForm}
end

function PolynomialForm(index_range::UnitRange) 
    PolynomialForm(
        OffsetVector{VarLinForm}([VarLinForm() for i in index_range], index_range)
    )
end

Base.getindex(poly::PolynomialForm, key::Int) = poly.terms[key]

firstindex(poly::PolynomialForm) = firstindex(poly.terms)
lastindex(poly::PolynomialForm) = lastindex(poly.terms)
eachpower(poly::PolynomialForm) = eachindex(poly.terms)

function add!(dest::PolynomialForm, source::PolynomialForm, factor=1)
    for power in eachpower(source)
        add!(dest[power], source[power],  factor)
    end
end

function mul!(poly::PolynomialForm, factor)
    for power in eachpower(poly)
        mul!(poly[power], factor)
    end
end

"""
    mul_by_power!(poly::PolynomialForm, factor)

Transform z^n to (z^n * factor^n) for each polynomial term.
Powers indices must encompass 0.
"""
function mul_by_power!(poly::PolynomialForm, factor)
    bottom = firstindex(poly)
    top  = lastindex(poly)
    if (bottom <= 0) && (top >= 0)
        min_bound = minimum(-bottom, top)
        factor_pow = factor

        # Transform central symmetrical part of polynomial
        for power in 1 : min_bound
            mul!(poly[power], factor_pow)
            mul!(poly[-power], 1/factor_pow)
            factor_pow *= factor
        end

        # Transform margins if any
        for power in min_bound+1 : top
            mul!(poly[power], factor_pow)
            factor_pow *= factor
        end
        for power in min_bound+1 : -bottom
            mul!(poly[-power], 1/factor_pow)
            factor_pow *= factor
        end
    else
        error("Indices don't encompass 0. This situation doesn't treated here")
    end
end

"""
    conjugate(poly::PolynomialForm)

Create conjugated form of a given polynomial form.
"""
function conjugate(poly::PolynomialForm)
    res = PolynomialForm(-eachpower(poly))
    for i in eachpower(res)
        add_conjugated!(res[i], poly[-i])
    end
    return res
end

"""
    z_conj_diff(poly::PolynomialForm)

For given polynomial form ϕ(z) creates z bar Φ(z).
"""
function z_conj_diff(poly::PolynomialForm)
    bottom = firstindex(poly)
    top    = lastindex(poly)
    res = PolynomialForm(-top+2 : -bottom+2)

    for n in eachpower(poly)
        add_conjugated!(res[-n+2], poly[n], n)
    end
    return res
end

function matrix_form(poly::PolynomialForm)
    bottom = firstindex(poly)
    top    = lastindex(poly)
    min_var_n, max_var_n = extrema(variables(poly[bottom]))
    for n in bottom+1 : top
        min_n, max_n = extrema(variables(poly[n]))
        (min_n < min_var_n) && (min_var_n = min_n)
        (max_n < max_var_n) && (max_var_n = max_n)
    end

    matrix = OffsetArray{ComplexF64, 2}(undef, bottom:top, min_var_n:max_var_n)
    fill!(matrix, 0.0im)

    for pow in eachpower(poly)
        for var in variables(poly[pow])
            matrix[pow, var] = poly[pow][var]
        end
    end

    return matrix
end















#= 

###############################################################################

abstract type GeneralFunction end

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
 =#
end # module Symbolic