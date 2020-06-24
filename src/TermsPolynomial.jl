"
Полином ``a_1 (z/r)^{k_1} + ... + a_n (z/r)^{k_n}``
"
struct PolynomialTerm <: FunctionalTerm
    "Коэффициенты полинома"
    coeffs :: ComplexOffsetVector
    "Нормирующий радиус"
    norm_r :: Float64
    "Простой конструктор"
    PolynomialTerm(coeffs::ComplexOffsetVector, norm_r::Float64) = new(coeffs, norm_r)
    "Конструктор"
    function PolynomialTerm(ind::UnitRange, norm_r::Float64)
        coeffs = OffsetVector{ComplexF64}(undef, ind)
        fill!(coeffs, 0.0im)
        new(coeffs, norm_r)
    end
end

Base.getindex(term::PolynomialTerm, i::Int) = getindex(term.coeffs, i)
Base.setindex!(term::PolynomialTerm, val::ComplexF64, i::Int) = setindex!(term.coeffs, val, i)
Base.firstindex(term::PolynomialTerm) = firstindex(term.coeffs)
Base.lastindex(term::PolynomialTerm) = lastindex(term.coeffs)
Base.eachindex(term::PolynomialTerm) = eachindex(term.coeffs)

"""
    differentiate(term :: PolynomialTerm) ::PolynomialTerm

Дифферециирование полинома. Возвращает новый полином.
"""
function differentiate(term :: PolynomialTerm) ::PolynomialTerm
    bottom = firstindex(term)
    top = lastindex(term)
    # Индексы дифференциированного полинома
    if top == 0 && bottom == 0
        # Производная константы остается константой
        d_bottom = 0
        d_top    = 0
        d_coeffs = OffsetVector([0.0im,], d_bottom : d_top)
        return PolynomialTerm(d_coeffs, term.norm_r)
    elseif top == 0
        # Производная от многочлена, заканчивающегося константой,
        # становится многочленом с максимум -2-й степенью
        d_bottom = bottom - 1
        d_top    = -2
    else
        # В обычном случае просто уменьшаем степени на единицу,
        # если только нижняя степень не нулевая
        d_bottom = bottom != 0 ? bottom - 1 : 0
        d_top = top - 1
    end

    d_coeffs = OffsetVector{ComplexF64}(undef, d_bottom : d_top)
    d_term = PolynomialTerm(d_coeffs, term.norm_r)
    
    for k in d_bottom : d_top
        d_term[k] = term[k+1] * (k+1) / term.norm_r
    end
    
    return d_term
end

"""
    function conjugate(term::PolynomialTerm, conj_r::Float64) ::PolynomialTerm

Сопряжение полинома по контуру радиусом conj_r. Возвращает новый полином.
"""
function conjugate(term::PolynomialTerm, conj_r::Float64) ::PolynomialTerm
    bottom = firstindex(term)
    top = lastindex(term)
    # Индексы сопряженного полинома
    c_bottom = -top
    c_top = -bottom

    conj_term = PolynomialTerm(c_bottom:c_top, term.norm_r)

    R = conj_r
    r = term.norm_r
    # z = r e^{i θ}
    # bar (z/r) = (1/r) (R^2 / z) = (R/r)^2 (r/z)
    factor = (R / r)^2
    
     # bar (z/r)^n = (R^2/r^2)^n (r/z)^n
    f_pow = factor
    for n in 1 : top
        conj_term[-n] = f_pow * conj(term[n])
        f_pow *= factor
    end

    f_pow = 1.0 / factor
    for n in -1 : -1 : bottom
        conj_term[-n] = f_pow * conj(term[n])
        f_pow /= factor
    end

    conj_term[0] = conj(term[0])

    return conj_term
end

"""
    function z_conj_diff(term::PolynomialTerm, conj_r::Float64) ::PolynomialTerm

Для многочлена term(z) вычисляет многочлен z bar term'(z)
"""
function z_conj_diff(term::PolynomialTerm, conj_r::Float64) ::PolynomialTerm
    bottom = firstindex(term)
    top = lastindex(term)

    # Индексы дифференциированного полинома
    if top == 0 && bottom == 0
        # Производная константы становится нулем
        # Ноль на что ни умножай, все равно ноль
        d_bottom = 0
        d_top    = 0
        d_coeffs = OffsetVector([0.0im,], 0 : 0)
        return PolynomialTerm(d_coeffs, term.norm_r)
    elseif top == 0
        # Производная от многочлена, заканчивающегося константой,
        # становится многочленом с максимум -2-й степенью
        d_bottom = bottom - 1
        d_top    = -2
    else
        # В обычном случае просто уменьшаем степени на единицу,
        # если только нижняя степень не нулевая
        d_bottom = bottom != 0 ? bottom - 1 : 0
        d_top = top - 1
    end

    # Индексы сопряженного полинома
    c_bottom = -d_top
    c_top = -d_bottom
    
    # Индексы умноженного на z полинома
    z_bottom = c_bottom + 1
    z_top = c_top +1

    z_term = PolynomialTerm(z_bottom:z_top, term.norm_r)
    
    R = conj_r
    r = term.norm_r
    # z = r e^{i θ}
    # bar (z/r) = (1/r) (R^2 / z) = (R/r)^2 (r/z)
    factor = (R / r)^2
    # bar (z/r)^n = (R^2/r^2)^n (r/z)^n

    # a_n z^n (d)-> n a_n z^{n-1} (c)-> n bar a_n z^{1-n} (z)-> n bar a_n z^{2-n}
    f_pow = factor
    for n in 1 : d_top+1
        z_term[2-n] = f_pow * n * conj(term[n])
        f_pow *= factor
    end

    f_pow = 1.0 / factor
    for n in -1 : -1 : bottom
        z_term[2-n] = f_pow * n * conj(term[n])
        f_pow /= factor
    end

    if bottom <= 1 <= top
        z_term[1] = conj(term[1])
    end

    return z_term
end

#=
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

Base.firstindex(poly::PolynomialForm) = firstindex(poly.terms)
Base.lastindex(poly::PolynomialForm) = lastindex(poly.terms)
eachpower(poly::PolynomialForm) = eachindex(poly.terms)

function Base.empty!(form::PolynomialForm)
    for p in eachpower(form)
        empty!(form.terms[p])
    end
end

function Base.similar(poly::PolynomialForm)
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
Base.similar(npoly::NormedPolynomial) = NormedPolynomial(similar(npoly.poly), npoly.r_norm)

Base.empty!(npoly::NormedPolynomial) = empty!(npoly.poly)
Base.firstindex(npoly::NormedPolynomial) = firstindex(npoly.poly)
Base.lastindex(npoly::NormedPolynomial) = lastindex(npoly.poly)

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
=#

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