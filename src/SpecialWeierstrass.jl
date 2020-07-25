#==============================================================================
Эллиптические функции Вейерштрасса (ЭФВ)

# Export

- `struct Weierstrass`: 
        Экземпляр эллиптической функции Вейерштрасса с предвычисленными параметрами и 
        кэшированными значениями

 - `WeierstrassData(l::Lattice, max_derivative ::Int)`: 
        Конструктор

 - `Base.getindex(w::Weierstrass, rz::RationalComplex, n_deriv::Int)`: 
        Кэшируемое вычисление n-й производной ЭФВ в заданной точке

# Description

ЭФВ нужны для построения двояко-периодических аналитических функций

Здесь реализован кэш в виде словаря с ключами в виде комплексных рациональных чисел.
В методе Ван Фо Фы требуются многочисленные вычисления ЭФВ и их производных на
смещениях между центрами волокон.
Если волокна расположены в виде регулярной решетки, то кэшированые вычислений дает
значительный выигрыш в скорости.
Однако если вычислять смещения через числа с плавающей точкой, то из-за ошибок округления
математически одинаковые смещения получат разные ключи, что приведет к лишним вычислениям
в практически совпадающих точках.
Поэтому центры волокон задаются рационально-комплексными числами, смещения тогда получаются 
рационально-комплексными и только такие числа принимаются в качестве аргументов.
==============================================================================#

"
Эллиптические функции Вейерштрасса
"
module SpecialWeierstrass

using ..SpecialTheta: Theta, theta, theta12
using ..Types: RationalComplex, ComplexOffsetVector, ComplexOffsetMatrix
using ..Types: CashedVector, last_cashed, not_cashed
using ..Types: Lattice, raw_complex, differentiate!

using OffsetArrays

"""
Экземпляр эллиптической функции Вейерштрасса с предвычисленными параметрами и кэшированными значениями

# Конструктор

    WeierstrassData(l::Lattice, max_derivative ::Int)

## Arguments
 - `l::Lattice`: решетка периодов
 - `max_derivative::Int`: максимальная требуемая производная
"""
struct Weierstrass #<: EllipticFunction
    "Решетка периодов"
    lattice  ::  Lattice
    "Корень характеристического уравнения"
    e1  ::  ComplexF64
    "Первый инвариант"
    g2  ::  ComplexF64
    "Второй инвариант"
    g3  ::  ComplexF64
    "Приращение ζ-функции на полупериоде ω1"
    η1  ::  ComplexF64
    "Приращение ζ-функции на полупериоде ω3"
    η3  ::  ComplexF64
    "Множитель в формуле для ℘-функции = frac{π θ_3(0) θ_4(0)}{2ω_1})^2"
    ℘_factor  ::  ComplexF64
    "Множитель в формуле для σ-функции"
    σ_factor  ::  ComplexF64
    "Экземпляр тэта-функции"
    θ  ::  Theta
    "Кэшированные ряды производных от ЭФВ в различных точках ячейки"
    cash  ::  Dict{RationalComplex, CashedVector{ComplexF64}}
    "Разложение ЭФВ в ряд Лорана около полюса"
    ℘_series  ::  ComplexOffsetVector
    "Разложения производных ЭФВ ``℘^{(n)}(z)/(n+1)!`` в ряд Лорана около полюса"
    derivs_series  ::  ComplexOffsetMatrix
    "Номер самой старшей производной"
    max_derivative  ::  Int
    "Конструктор"
    function Weierstrass(l::Lattice, max_derivative ::Int)
        # Вычисление многочисленных констант
        ϵ = eps()
        θ = Theta(l, ϵ)

        ω1 = l.ω1
        ω3 = l.ω3
        tolerance = eps(imag(ω3))

        th2 = theta(θ, th_k=2, d_n=0, z=0.0im)
        th4 = theta(θ, th_k=4, d_n=0, z=0.0im)
        th24 = th2^4
        th44 = th4^4
        e_factor = π^2 / (12 * ω1^2)
        e1 = e_factor * (th24 + 2th44)
        e2 = e_factor * (th24 - th44)
        e3 = -e_factor * (2th24 + th44)

        g2 = -4 * (e2*e3 + e3*e1 + e1*e2)
        g3 = 4 * e1 * e2 * e3

        d3th1 = theta(θ, th_k=1, d_n=3, z=0.0im)
        d1th1 = theta(θ, th_k=1, d_n=1, z=0.0im)

        η1 = -π^2 / (12ω1) * d3th1 / d1th1
        η3 = (η1 * ω3 - 1im*π/2) / ω1

        th3 = theta(θ, th_k=3, d_n=0, z=0.0im)
        ℘_factor = (π * th3 * th4 / (2ω1))^2

        σ_factor = 2ω1 / (π*d1th1)

        # Кэш для вычислений производных
        cash = Dict{RationalComplex, CashedVector{ComplexF64}}()

        # Ненулевые коэффициенты при положительных степенях разложения ЭФВ в ряд Лорана около полюса
        c = OffsetVector{ComplexF64}(undef, 2:max_derivative+2)
        c[2] = g2/20.0
        c[3] = g3/28.0
        for n in 4:lastindex(c)
            c[n] = 0.0im
            for m in 2:n-2
                c[n] += c[m] * c[n-m]
            end
            c[n] *= 3/((2n+1)*(n-3))
        end

        # Часть по неотрицательным степеням разложения ЭФВ в ряд Лорана около полюса
        laurent_praecursor = OffsetVector{ComplexF64}(undef, 0:2*(max_derivative+2))
        fill!(laurent_praecursor, 0.0im)
        for n in 2 : max_derivative+2
            laurent_praecursor[2n-2] = c[n]
        end
        # Копия разложения, поскольку `laurent_praecursor` дальше будет дифференциироваться
        ℘_series = deepcopy(laurent_praecursor)

        # Разложения производных ЭФВ ``℘^{(n)}(z)/(n+1)!`` в ряд Лорана около полюса
        # Получаются путем последовательного дифференциирования разложения ℘-функции
        derivs_series = OffsetArray{ComplexF64, 2}(undef, -1:max_derivative, 0:max_derivative+2)
        # факториал (k+1)!
        fact = 1.0
        for k in 0:max_derivative
            for n in 0 : max_derivative+2
                derivs_series[k,n] = laurent_praecursor[n] / fact
            end
            differentiate!(laurent_praecursor)
            fact *= k+2
        end
        # Разложение ζ-функции Вейерштрасса
        # Получаются путем интегрирования разложения ℘-функции
        derivs_series[-1,0] = 0.0im
        for n in 1 : max_derivative+2
            derivs_series[-1,n] = derivs_series[0,n-1]/n
        end

        ##########
        new(l, e1, g2, g3, η1, η3, ℘_factor, σ_factor, θ, cash, 
            ℘_series, derivs_series, max_derivative)
    end
end

"""
    cash_elder_derivs!(output::CashedVector{ComplexF64}, n_deriv::Int)

Кэширование недостающих старших производных
"""
function cash_elder_derivs!(output::CashedVector{ComplexF64}, n_deriv::Int)
    # wp_{n+2}(z) = frac{6}{(n+1)(n+2)} sum_{k=0}^n wp_k(z) wp_{n-k}(z),
    # n >= 1.
    last_cashed_deriv = last_cashed(output)
    for n in last_cashed_deriv-1 : n_deriv-2
        res = 0.0im
        for k in 0 : n
            res += output[k] * output[n-k]
        end
        res *= 6.0 / ((n+1) * (n+2))
        output[n+2] = res
    end
end

"""
    cash_first_derivs!(w::Weierstrass, rz::RationalComplex, output::CashedVector{ComplexF64})

Начальное заполнение кэша нормированными производными от -1-й до 2-й включительно в заданной точке
"""
function cash_first_derivs!(w::Weierstrass, rz::RationalComplex, output::CashedVector{ComplexF64})
    ω1 = w.lattice.ω1
    e1 = w.e1
    η1 = w.η1
    g2 = w.g2
    ℘_factor = w.℘_factor
    σ_factor = w.σ_factor

    z = raw_complex(w.lattice,rz)
    # Нормирующий множитель, равный модулю z
    # С таким множителем все производные ЭФВ по модулю примерно равны 1
    norm_factor = abs(z)

    u = (π / (2ω1)) * z
    th1, th2, dth1, dth2 = theta12(w.θ, u)
    
    # Трюк для улучшения точности вычислений при малых |z|
    # |th1_norm| ≈ 1
    th1_norm = th1 / norm_factor

    ζ = (π/(2ω1)) * (dth1/th1_norm)
    ζ += norm_factor * (η1 / ω1) * z

    ℘ = ℘_factor * (th2 / th1_norm)^2
    ℘ += norm_factor^2 * e1

    ℘_prime = th1 * th2 * dth2 - dth1 * th2 * th2
    ℘_prime /= th1_norm^3
    ℘_prime *= (π / ω1) * ℘_factor

    ℘_pp = 6*℘*℘ - 0.5 * g2 * norm_factor^4

    output[-1] = -ζ
    output[0] = ℘             # /0!
    output[1] = ℘_prime       # /1!
    output[2] = ℘_pp / 2.0    # /2!
end

"""
    Base.getindex(w::Weierstrass, rz::RationalComplex, n_deriv::Int)

Кэшируемое вычисление n-й производной ЭФВ в заданной точке

# Output

``|z|^(n+2)/(n+1)! ℘^{(n)}(z)``
"""
function Base.getindex(w::Weierstrass, rz::RationalComplex, n_deriv::Int)
    if rz == 0
        error("Can't compute Weierstrass elliptic functions at zero point")
    end

    if rz in keys(w.cash)
        # Если для точки уже есть кэш, то берем значения оттуда
        vector = w.cash[rz]
        sign = 1.0
    elseif -rz in keys(w.cash)
        # Если кэш есть для противоположной точки, то берем значения оттуда
        # с поправкой на четность/нечетность производной
        vector = w.cash[-rz]
        sign = (n_deriv % 2 == 0) ? 1.0 : -1.0
    else
        # Если кэша не ни для самой точки, на для обратной,
        # то создаем новый кэш
        vector = CashedVector{ComplexF64}(-1:2*w.max_derivative+2)
        w.cash[rz] = vector
        # И заполняем его первоначальными производными
        cash_first_derivs!(w, rz, vector)
        sign = 1.0
    end

    if not_cashed(vector, n_deriv)
        # Если нужной производной в кэше нет, то кэшируем все недостающие
        cash_elder_derivs!(vector, n_deriv)
    end

    return vector[n_deriv]*sign
end

end # module SpecialWeierstrass