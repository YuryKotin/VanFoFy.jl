#==============================================================================
Специальная функция Q(z)

# Export

- struct QSpecial:
        Экземпляр специальной функции Q(z)

- QSpecial(℘::Weierstrass):
        Конструктор

- Base.getindex(Q::QSpecial, rz::RationalComplex, n_deriv::Int):
        Кэшируемое вычисление n-й производной функции Q(z) в заданной точке

# Description

α Q(z) = 0.5℘'(z) + ζ(z)℘(z) - β z ℘(z) + β ζ(z) - g_2 z / 12

Нужна для придания двояко-периодичности выражению z bar Φ(z) + Ψ(z).
Все вычисления кэшируются так же, как и у эллиптической функции Вейерштрасса.
===============================================================================#

"
Специальная функция Q(z)
"
module SpecialQ

using ..Types: differentiate!, raw_complex, Lattice
using ..Types: RationalComplex, CashedVector, ComplexOffsetMatrix
using ..Types: last_cashed, not_cashed
using ..SpecialWeierstrass: Weierstrass

using OffsetArrays

"
Экземпляр специальной функции Q(z)

# Конструктор

    QSpecial(℘::Weierstrass)

Экземпляр специальной функции Q(z) строится на основании экземпляра 
эллиптической функции Вейерштрасса ℘(z).
"
struct QSpecial #<: EllipticFunction
    "Экземпляр эллиптической функции Вейерштрасса"
    ℘              ::Weierstrass
    "Решетка периодов"
    lattice        ::Lattice
    # Константы
    α              ::ComplexF64
    β              ::ComplexF64
    γ1             ::ComplexF64
    γ3             ::ComplexF64
    "Кэшированные ряды производных в различных точках ячейки"
    cash           ::Dict{RationalComplex, CashedVector{ComplexF64}}
    "Разложения производных в ряд Тейлора около полюса"
    derivs_series  ::ComplexOffsetMatrix
    # Конструктор
    function QSpecial(℘::Weierstrass)
        # Вычисление многочисленных констант
        lattice = ℘.lattice
        
        ω1 = lattice.ω1
        ω3 = lattice.ω3
        
        w1 = 2ω1
        w3 = 2ω3
        
        η1 = ℘.η1
        η3 = ℘.η3
    
        α = (2im*π)/(conj(w1)*w3-conj(w3)*w1)
        β = 2(η3*conj(w1)-η1*conj(w3))/(conj(w1)*w3-conj(w3)*w1)
    
        g2 = ℘.g2
        # γ1 = (1/α) (2 β η1 - (1/6) g2 ω1)
        # γ3 = (1/α) (2 β η3 - (1/6) g2 ω3)
        γ1 = (2*β*η1-g2*ω1/6)/α
        γ3 = (2*β*η3-g2*ω3/6)/α
        
        # Кэш для вычислений производных
        cash = Dict{RationalComplex, CashedVector{ComplexF64}}()
    
        # Разложение функции Q(z) в ряд Тейлора около полюса
        # Используется разложение в ряд Лорана функции ℘(z) около полюса
        Q_series = OffsetVector{ComplexF64}(undef, 0:2*(℘.max_derivative+2))
        fill!(Q_series, 0.0im)
        ℘_series = ℘.℘_series
        K = lastindex(℘_series)

        # По следующей формуле собирается разложение в ряд функции Q(z)
        # α Q(z) = 0.5℘'(z) + ζ(z)℘(z) - β z ℘(z) + β ζ(z) - g_2 z / 12
        for k in 0 : K-1
        # 0.5 ℘'(z)
            Q_series[k] += 0.5 * (k+1) * ℘_series[k+1]
        # +ζ(z) ℘(z)
        ## +(1/z) * ℘(z)
            Q_series[k] += ℘_series[k+1]
        ## +(1/z^2) ζ(z)
            Q_series[k] -= ℘_series[k+1]/(k+2)
        ##
            for m in 0 : K-k-1
                Q_series[m+k+1] -= ℘_series[k] * ℘_series[m] / (m+1)
            end
        # -β z ℘(z)
            Q_series[k+1] -= β * ℘_series[k]
        # +β ζ(z)
            Q_series[k+1] -= β * ℘_series[k] / (k+1)
        end
        # -g_2 z / 12
        Q_series[1] -= ℘.g2 / 12.0
    
        Q_series[:] ./= α
    
        # Разложение производных ``Q^{(n)}(z)/(n+1)!`` в ряд Тейлора около полюса
        derivs_series = OffsetArray{ComplexF64, 2}(undef, 0:℘.max_derivative, 0:℘.max_derivative+2)
        for k in 0:℘.max_derivative
            Q_series[:] ./= k+1
            for n in 0 : ℘.max_derivative+2
                derivs_series[k,n] = Q_series[n]
            end
            differentiate!(Q_series)
        end
    
        ######################################################################
    
        new(℘, lattice, α, β, γ1, γ3, cash, derivs_series)
    end
end

"""
    cash_elder_derivs!(output::CashedVector{ComplexF64}, Q::QSpecial, rz::RationalComplex, n_deriv::Int)

Кэширование недостающих старших производных
"""
function cash_elder_derivs!(output::CashedVector{ComplexF64}, Q::QSpecial, 
                              rz::RationalComplex, n_deriv::Int)
    ℘ = Q.℘
    α = Q.α
    β = Q.β

    z = raw_complex(Q.lattice, rz)
    az = abs(z)
    
    last_cashed_deriv = last_cashed(output)
    for n in (last_cashed_deriv+1)+1 : (n_deriv+1)
        qn = 0.0im
        qn += 0.5 * ℘[rz, n]
        qn += -℘[rz, -1] * ℘[rz, n-1] / n
        for k in 1 : n-1
            qn += -0.5 * ℘[rz, k-1] * ℘[rz, n-k-1] / (k*(n-k))
        end
        qn += -az * z * β * ℘[rz, n-1] / n
        qn += -az^2 * β * ℘[rz, n-2] / (n-1)

        output[n-1] = qn / α
    end
end

"""
    cash_first_derivs!(output::CashedVector{ComplexF64}, Q::QSpecial, rz::RationalComplex)

Начальное заполнение кэша нормированными значениями Q(z) и Q'(z).
"""
function cash_first_derivs!(output::CashedVector{ComplexF64}, Q::QSpecial, rz::RationalComplex)
    ℘ = Q.℘
    α = Q.α
    β = Q.β
    g2 = ℘.g2

    z = raw_complex(Q.lattice, rz)
    az = abs(z)

    Q0 = 0.0im
    Q0 += 0.5*℘[rz, 1]
    Q0 += -℘[rz, -1] * ℘[rz, 0]
    Q0 += -az * z * β * ℘[rz, 0]
    Q0 += -az^2 * β * ℘[rz, -1]
    Q0 += -az^3 * g2 * z / 12.0

    Q1 = 0.0im
    Q1 += ℘[rz, 2] / 3.0
    Q1 += -0.5 * ℘[rz, -1] * ℘[rz, 1]
    Q1 += -0.5 * az * z * β * ℘[rz, 1]
    Q1 += -az^2 * β * ℘[rz, 0]
    Q1 += -az^4 * g2 / 12.0

    output[0] = Q0 / α
    output[1] = Q1 / α
end

"""
    Base.getindex(Q::QSpecial, rz::RationalComplex, n_deriv::Int)

Кэшируемое вычисление n-й производной функции Q(z) в заданной точке

# Output

``|z|^(n+3)/(n+1)! Q^{(n)}(z)``
"""
function Base.getindex(Q::QSpecial, rz::RationalComplex, n_deriv::Int)
    if rz == 0
        error("Can't compute Q special functions at zero point")
    end

    if rz in keys(Q.cash)
        derivs_array = Q.cash[rz]
        sign = 1.0
    elseif -rz in keys(Q.cash)
        derivs_array = Q.cash[-rz]
        sign = (n_deriv % 2 == 0 ? -1.0 : 1.0)
    else
        derivs_array = CashedVector{ComplexF64}(0:2*Q.℘.max_derivative+1)
        Q.cash[rz] = derivs_array
        cash_first_derivs!(derivs_array, Q, rz)
        sign = 1.0
    end

    if not_cashed(derivs_array, n_deriv)
        cash_elder_derivs!(derivs_array, Q, rz, n_deriv)
    end
    return derivs_array[n_deriv] * sign
end

end # module SpecialQ