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
        n_ind = size(ind, 1)
        coeffs = OffsetVector( zeros(ComplexF64, n_ind), ind)
        new(coeffs, norm_r)
    end
end

Base.getindex(term::PolynomialTerm, i::Int) = term.coeffs[i]
Base.setindex!(term::PolynomialTerm, val::ComplexF64, i::Int) = setindex!(term.coeffs, val, i)
Base.firstindex(term::PolynomialTerm) = firstindex(term.coeffs)
Base.lastindex(term::PolynomialTerm) = lastindex(term.coeffs)
Base.eachindex(term::PolynomialTerm) = firstindex(term.coeffs) : lastindex(term.coeffs)
Base.ndims(term::PolynomialTerm) = 1
Base.axes(term::PolynomialTerm) = axes(term.coeffs)
function Base.similar(term::PolynomialTerm) 
    PolynomialTerm(eachindex(term), term.norm_r)
end

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
    reconjugate(term::PolynomialTerm, old_conj_r::Float64, new_conj_r::Float64)

Пересопряжение полинома по другому контуру
"""
function reconjugate(term::PolynomialTerm, old_conj_r::Float64, new_conj_r::Float64)
    rc_term = similar(term)
    
    # bar (z/r)^n = (R1/r)^{2n} (r/z)^n =
    #             = (R2/r)^{2n} (r/z)^n =
    #             = (R2/R1)^{2n} (R1/r)^{2n} (r/z)^n
    R1 = old_conj_r
    R2 = new_conj_r

    factor = (R1 / R2)^2
    f_pow = factor
    
    bottom = firstindex(term)
    top = lastindex(term)
    for i in 1 : top
        rc_term[i] = term[i] * f_pow
        f_pow *= factor
    end

    rc_term[0] = term[0]

    f_pow = 1/factor
    for i in -1 : -1 : bottom
        rc_term[i] = term[i] * f_pow
        f_pow /= factor
    end
    
    return rc_term
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

    # После дифференцииования все степени сдвигаются на -1,
    # поэтому множитель при сопряжении должен быть не (R/r)^{2n}, а (R/r)^{2n-2}
    f_pow = 1.0 # (R/r)^0
    
    # a_n z^n (d)-> n a_n z^{n-1} (c)-> n bar a_n z^{1-n} (z)-> n bar a_n z^{2-n}
    for n in 1 : d_top+1
        z_term[2-n] = f_pow * n * conj(term[n])
        f_pow *= factor
    end

    f_pow = 1.0 / factor^2 # (R/r)^{-4}
    for n in -1 : -1 : bottom
        z_term[2-n] = f_pow * n * conj(term[n])
        f_pow /= factor
    end

    if bottom <= 1 <= top
        z_term[1] = conj(term[1])
    end

    return z_term
end

"""
    add!(dest::PolynomialTerm, source::PolynomialTerm, factor)

Сложение двух нормированных полиномов. dest += source * factor.

Выход за границы массивов не проверяется.
"""
function add!(dest::PolynomialTerm, source::PolynomialTerm, factor)
    r = dest.norm_r
    R = source.norm_r
    # z/r = (z/R) * (R/r)
    if r == R
        # На всякий случай, чтобы избежать ошибок округления
        norm_factor = 1.0
    else
        norm_factor = R / r
    end

    # (z/r)^n = (z/R)^n * (R/r)^n
    term_factor = 1.0 # (R/r)^n
    for n in 0 :lastindex(source)
        dest[n] += source[n] * factor * term_factor
        term_factor *= norm_factor
    end

    term_factor = 1.0 / norm_factor
    for n in -1 : -1 :firstindex(source)
        dest[n] += source[n] * factor * term_factor
        term_factor /= norm_factor
    end
end

"""
    add_term_series!(output, term::PolynomialTerm; factor::ComplexF64)
Прибавление полинома к контейнеру output
"""
function add_term_series!(  output, 
                            term  ::PolynomialTerm; 
                            factor::ComplexF64)
    for i in eachindex(term)
        output[i] += term[i] * factor
    end
end

###############################################################################

"
Линейно-полиномиальная форма вида
``v_1 P_1(z) + ... +  v_N P_N(z)``,
где 
- ``v_m`` - переменные линейной формы
- ``P_m(z)`` - полином
"
const PolynomialSolution = VarLinForm{PolynomialTerm}

"""
    z_conj_diff(sol::PolynomialSolution, conj_r::Float64)

Для линейно-полиномиальной формы ϕ(z) вычисляет многочлен z bar ϕ'(z)
"""
function z_conj_diff(sol::PolynomialSolution, conj_r::Float64)
    VarLinForm{PolynomialTerm}(
        OffsetVector(
            [
                z_conj_diff(sol[v], conj_r) for v in eachindex(sol)
            ],
            axes(sol, 1)
        )
    )
end

"""
    function conjugate(sol::PolynomialSolution, conj_r::Float64)

Сопряжение линейно-полиномиальной формы по контуру радиусом conj_r. Возвращает новое решение.
"""
function conjugate(sol::PolynomialSolution, conj_r::Float64)
    VarLinForm{PolynomialTerm}(
        OffsetVector(
            [
                conjugate(sol[v], conj_r) for v in eachindex(sol)
            ],
            axes(sol, 1)
        )
    )
end

"""
    reconjugate(sol::PolynomialSolution, old_conj_r::Float64, new_conj_r::Float64)

Пересопряжение линейно-полиномиальной формы по другому контуру
"""
function reconjugate(sol::PolynomialSolution, old_conj_r::Float64, new_conj_r::Float64)
    VarLinForm{PolynomialTerm}(
        OffsetVector(
            [
                reconjugate(sol[v], old_conj_r, new_conj_r) for v in eachindex(sol)
            ],
            axes(sol, 1)
        )
    )
end

"""
    add!(dest::PolynomialSolution, source::PolynomialSolution, factor)

Сложение двух линейно-полиномиальных форм. dest += source * factor.

Выход за границы массивов не проверяется.
"""
function add!(dest::PolynomialSolution, source::PolynomialSolution, factor)
    for v in eachindex(source)
        add!(dest[v], source[v], factor)
    end
end

"""
    similar(source::VarLinForm{PolynomialTerm})

Создание аналогичной линейно-полиномиальной формы
"""
function Base.similar(source::VarLinForm{PolynomialTerm})
    form = VarLinForm{PolynomialTerm}(
        OffsetVector(
            [
                PolynomialTerm(
                    firstindex(source[v]) : lastindex(source[v]),
                    source[v].norm_r
                )
                for v in eachindex(source)
            ],
            firstindex(source) : lastindex(source)
        )

    )
    return form
end

eachpower(sol::PolynomialSolution) = eachindex(first(sol))