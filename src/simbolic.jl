abstract type GeneralFunction end

differentiate!(f::GeneralFunction) = 
    error("Differentiation is not defined in the concrete type")

###############################################################################

"""
  (z-ζ_k)^n/r_k^n

power -  n,
center - ζ_k,
norm -   r_k.
"""
mutable struct ZPowN <: GeneralFunction
    power  ::Int64
    center ::ComplexF64
    norm   ::Float64
end

function differentiate!(f::ZPowN) 
    # [(z-a)^n/r^n]' = (n/r) z^{n-1}/r^{n-1}
    multiplier = f.power / f.norm
    f.power -= 1
    return multiplier
end

###############################################################################

"""
  r_k^{n+2} ℘^{(n)}(z-ζ_k)/(n+1)!


deriv -  n,
center - ζ_k,
norm -   r_k.
"""
mutable struct WpDerivative <: GeneralFunction
    deriv::Int64
    center::ComplexF64
    norm::Float64
end

function differentiate!(f::WpDerivative) 
    # [r_k^{n+2} ℘^{(n)}(z-ζ_k)/(n+1)!]' =
    # [(n+2) / r_k] *
    # r_k^{n+3} ℘^{(n+1)}(z-ζ_k)/(n+2)!
    multiplier = (f.deriv+2) / f.norm
    f.deriv += 1
    return multiplier
end

###############################################################################

"""
  r_k^{n+2} Q^{(n)}(z-ζ_k)/(n+1)!


deriv -  n,
center - ζ_k,
norm -   r_k.
"""
mutable struct QSpecial <: GeneralFunction
    deriv::Int64
    center::ComplexF64
    norm::Float64
end

function differentiate!(f::QSpecial) 
    # [r_k^{n+2} Q^{(n)}(z-ζ_k)/(n+1)!]' =
    # [(n+2) / r_k] *
    # r_k^{n+3} Q^{(n+1)}(z-ζ_k)/(n+2)!
    multiplier = (f.deriv+2) / f.norm
    f.deriv += 1
    return multiplier
end

###############################################################################

mutable struct Term{F <: GeneralFunction}
    variable    ::Int64
    coefficient ::ComplexF64
    functional  ::F
    conjugated  ::Bool
    mul_by_z    ::Bool
end

function differentiate!(term::Term)
    multiplier = differentiate!(term.functional)
    term.coefficient *= multiplier
    nothing
end

function conjugate!(term::Term)
    term.coefficient = conj(term.coefficient)
    term.conjugated = !term.conjugated
    nothing
end

function mul!(term::Term, factor::Float64)
    term.coefficient *= factor
    nothing
end

###############################################################################

struct LinearForm
    terms ::Vector{Term}
end

function differentiate!(form::LinearForm)
    foreach(t->differentiate!(t), form.terms)
    nothing
end

function conjugate!(form::LinearForm)
    foreach(t->conjugate!(t), form.terms)
    nothing
end

function mul!(form::LinearForm, factor::Float64)
    foreach(t->mul!(t, factor), form.terms)
    nothing
end