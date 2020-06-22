struct VarLinForm
    form ::ComplexOffsetVector
end

struct PolynomialForm
    terms ::OffsetVector{ VarLinForm }
end

struct NormedPolynomial
    poly ::PolynomialForm
    r_norm ::Float64
end

###############################################################################
## VarLinForm
###############################################################################

VarLinForm(var_range::UnitRange) = VarLinForm(OffsetVector([0.0im for v in var_range], var_range))

Base.getindex(vlf::VarLinForm, key::Int) = vlf.form[key]

Base.setindex!(vlf::VarLinForm, val::ComplexF64, key::Int) = setindex!(vlf.form, val, key)

variables(vlf::VarLinForm) = eachindex(vlf.form)

function add!(dest::VarLinForm, source::VarLinForm, factor)
    @__dot__ dest += source * factor
end

function add_conjugated!(dest::VarLinForm, source::VarLinForm, factor)
    @__dot__ dest += conj(source) * factor
end

function mul!(vlf::VarLinForm, factor)
    vlf .*= factor
end

empty!(vlf::VarLinForm) = fill!(vlf.form, 0.0im)

###############################################################################
## PolynomialForm
###############################################################################

function PolynomialForm(index_range::UnitRange)
    PolynomialForm(
        OffsetVector{VarLinForm}([VarLinForm() for i in index_range], index_range)
    )
end

Base.getindex(poly::PolynomialForm, key::Int) = poly.terms[key]

firstindex(poly::PolynomialForm) = firstindex(poly.terms)
lastindex(poly::PolynomialForm) = lastindex(poly.terms)
eachpower(poly::PolynomialForm) = eachindex(poly.terms)

function empty!(form::PolynomialForm)
    for p in eachpower(form)
        empty!(form.terms[p])
    end
end

function similar(poly::PolynomialForm)
    index_range = firstindex(poly) : lastindex(poly)
    PolynomialForm(index_range)
end

function similar_plane_conj(poly::PolynomialForm)
    index_range = -lastindex(poly)+2 : -firstindex(poly)+2
    PolynomialForm(index_range)
end

function plane_conj_invarianted(poly::PolynomialForm)
    bottom = firstindex(poly)
    top = lastindex(poly)
    top == -bottom + 2
end

function add!(dest::PolynomialForm, source::PolynomialForm, factor)
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
"""
function mul_by_power!(poly::PolynomialForm, factor)
    bottom = firstindex(poly)
    top  = lastindex(poly)

    for power in 1 : top
        factor_pow = factor
        mul!(poly[power], factor_pow)
        factor_pow *= factor
    end

    for power in -1 : -1 : bottom
        factor_pow = 1.0/factor
        mul!(poly[power], factor_pow)
        factor_pow /= factor
    end
end

"""
    conjugate(poly::PolynomialForm)

Create conjugated form of a given polynomial form.
"""
function conjugate(poly::PolynomialForm)
    conj_range = -lastindex(poly) : -firstindex(poly)
    res = PolynomialForm(conj_range)
    for i in eachpower(res)
        add_conjugated!(res[i], poly[-i], 1.0)
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

###############################################################################
## NormedPolynomial
###############################################################################

NormedPolynomial(index_range::UnitRange, r_norm::Float64) = NormedPolynomial(PolynomialForm(index_range), r_norm)
similar(npoly::NormedPolynomial) = NormedPolynomial(similar(npoly.poly), npoly.r_norm)

empty!(npoly::NormedPolynomial) = empty!(npoly.poly)
firstindex(npoly::NormedPolynomial) = firstindex(npoly.poly)
lastindex(npoly::NormedPolynomial) = lastindex(npoly.poly)

Base.getindex(npoly::NormedPolynomial, key::Int) = npoly.poly[key]

plane_conj_invarianted(npoly::NormedPolynomial) = plane_conj_invarianted(npoly.poly)

function add!(dest::NormedPolynomial, source::NormedPolynomial, factor) 
    if source.r_norm == dest.r_norm
        add!(dest.poly, source.poly, factor)
    else
        error("Norming radii don't match")
    end
end

function conjugate(npoly::NormedPolynomial, r_contour) 
    c_npoly = NormedPolynomial(conjugate(npoly.poly), npoly.r_norm)
    r = npoly.r_norm
    R = r_contour
    mul_by_power!(c_npoly, (R/r)^2)
    return c_npoly
end

function z_conj_diff(normed_poly::NormedPolynomial, r_contour) end

#=

###############################################################################

abstract type GeneralFunction end

struct FunctionalTerm{F <: GeneralFunction}
    func   ::F
    variable ::Int
end



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