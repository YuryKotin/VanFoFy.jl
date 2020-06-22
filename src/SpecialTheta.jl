#==============================================================================
Тэта-функции

theta_1(z,q) = 2 sum_{n=0}^infty (-1)^n q^{(n+0.5)^2} sin((2n+1)z)
theta_2(z,q) = 2 sum_{n=0}^infty q^{(n+0.5)^2} cos((2n+1)z)
theta_3(z,q) = 1 + 2 sum_{n=1}^infty q^{n^2} cos(2nz)
theta_4(z,q) = 1 + 2 sum_{n=1}^infty (-1)^n q^{n^2} cos(2nz)

Тэта-функции нужны для вычисления функций Вейерштрасса.
Наиболее часто требуются θ1, θ2, θ1', θ2'.
Для ускорения счета заранее вычисляются последовательности q^{(n+0.5)^2} и (2n+1)q^{(n+0.5)^2}.

В этом файле:
 - struct Theta: предвычисленные параметры тэта-функций;
 - theta12(th, z): быстрое вычисление θ1, θ2, θ1', θ2';
 - theta(θ; th_k, d_n, z): медленное вычисление производных тэта-функций.
==============================================================================#

"Тэта-функции"
module SpecialTheta

using OffsetArrays

using ..Types: ComplexOffsetMatrix, Lattice

"""
Предвычисленные параметры тэта-функций

# Конструктор

    Theta(l::Lattice, ϵ::Float64)

## Arguments
- `l::Lattice`: решетка полупериодов
- `ϵ::Float64`: абсолютная погрешность вычисления функций
"""
struct Theta
    "Ном"
    nome ::ComplexF64
    "Предвычисленные степени нома для быстрого вычисления тета-функций"
    nome_powers ::ComplexOffsetMatrix
    "Количество предвычисленных степеней нома"
    powers_num ::Int
    "Точность вычисления тэта-функций"
    ϵ ::Float64
    "Конструктор"
    function Theta(l::Lattice, ϵ ::Float64)        
        ω3 = l.ω3 
        ω1 = l.ω1
        ###########################################################################
        # Вычисление нома
        τ = ω3 / ω1
        q = exp(1im * π * τ)

        ###########################################################################
        # Вычисление степеней нома
        #
        # В дальнейшем нам нужно часто вычислять значения функций θ1(z), θ2(z), θ1'(z), θ2'(z).
        # Эти вычисления можно значительно ускорить, если заранее вычислить необходимое количество
        # множителей вида q^{(n+0.5)^2} и (2n+1)q^{(n+0.5)^2}.

        # Для начала определяем количество таких степеней
        # Зададимся точкой, в которой будем контролировать погрешность вычислений
        # Возьмем угловую точку ячейки периодичности
        z0 = ω1+ω3
        
        # Установим максимальное число членов ряда, нужных для достижения такой погрешности
        max_n_terms = 10
        
        # Начинаем вычислять члены рядов Фурье для функций θ1'(z) и θ2'(z), 
        # пока оба члены не станут меньше требуемой погрешности или
        # пока их количество не превзойдет установленный максимум
        n = 0
        while true
            # theta_1(z,q) = 2 sum_{n=0}^infty (-1)^n q^{(n+0.5)^2} sin((2n+1)z)
            dt1 = 2*(2n+1)*cos((2n+1)*z0)*q^((n+0.5)^2)
            # theta_2(z,q) = 2 sum_{n=0}^infty q^{(n+0.5)^2} cos((2n+1)z)
            dt2 = 2*(2n+1)*sin((2n+1)*z0)*q^((n+0.5)^2)
            if (abs(dt1) < ϵ) && (abs(dt1) < ϵ)
                break
            end
            n += 1
            if n > max_n_terms
                error("Series don't converge")
            end
        end

        # Теперь вычисляем сами степени
        nome_powers = OffsetArray{ComplexF64}(undef, 1:2, 0:n)
        for k in 0 : n
            nome_powers[1, k] = q^((k+0.5)^2)
            nome_powers[2, k] = nome_powers[1, k] * (2k+1)
        end

        new(q, nome_powers, n, ϵ)
    end
end


"""
    theta12(th::Theta, z::ComplexF64)

Быстрое вычисление значений θ_1(z), θ_2(z), θ_1'(z), θ_2'(z) в заданной точке.
"""
function theta12(th ::Theta, z  ::ComplexF64)
    #####
    t1  = 0.0im
    t2  = 0.0im
    dt1 = 0.0im
    dt2 = 0.0im
    sign = 1.0

    N = th.powers_num
    for n in 0 : N
        z_n = (2n+1)*z
        sin_z = sin(z_n)
        cos_z = cos(z_n)
        # theta_1(z,q) = 2 sum_{n=0}^N (-1)^n q^{(n+0.5)^2} sin((2n+1)z)
        t1 += 2 * sign * th.nome_powers[1, n] * sin_z
        # theta_2(z,q) = 2 sum_{n=0}^N q^{(n+0.5)^2} cos((2n+1)z)
        t2 += 2 * th.nome_powers[1, n] * cos_z
        # theta_1'(z,q) = 2 sum_{n=0}^N (-1)^n (2n+1) q^{(n+0.5)^2} cos((2n+1)z))
        dt1 += 2 * sign * th.nome_powers[2, n] * cos_z
        # theta_2'(z,q) = -2 sum_{n=0}^N (2n+1) q^{(n+0.5)^2} sin((2n+1)z)
        dt2 += -2 * th.nome_powers[2, n] * sin_z

        sign = -sign
    end

    return t1, t2, dt1, dt2
end


""""
  theta(th::Theta; th_k::Int, d_n::Int, z::ComplexF64) ::ComplexF64

Вычисление значений производной тэта-функции.

# Arguments
- `th::Theta`: экземпляр параметров тэта-функций
- `th_k::Int`: номер тэта-функции (1, 2, 3, 4)
- `d_n::Int`: номер производной, >=0
- `z::ComplexF64`: точка
"""
function theta(th::Theta; th_k::Int, d_n::Int, z::ComplexF64) ::ComplexF64
    #####
    q = th.nome
    ϵ = th.ϵ

    if (th_k < 1) || (th_k > 4)
        error("Invalid theta function number")
    end

    # Пара sin z, cos z возвращается в саму себя после четырехкратного дифференциирования
    rem_d = d_n % 4
    # Знак-множитель после дифференциирования:
    # (sin z)' = +1.0 * cos z; (sin z)'' = -1.0 * sin z и т.д.
    d_sin_sign = (rem_d == 0) || (rem_d == 1) ? 1.0 : -1.0
    d_cos_sign = (rem_d == 0) || (rem_d == 3) ? 1.0 : -1.0

    # n - index of summation, sum - сумма ряда
    # В зависимости от номера тета-функции и её производной они имеют разные начальные значения
    if (th_k == 1) || (th_k == 2)
        sum = 0.0im
        n = 0
    else
        n = 1
        if d_n == 0
            sum = 1.0+0.0im
        else
            sum = 0.0im
        end
    end

    sign = (-1.0)^n
    while true
        # theta_1(z,q) = 2 sum_{n=0}^infty (-1)^n q^{(n+0.5)^2} sin((2n+1)z)
        if th_k == 1
            # term - очередной член ряда
            term = 2.0 * sign * (2n+1)^d_n * q^((n+0.5)^2) * d_sin_sign
            if d_n % 2 == 0
                term *= sin((2n+1)*z)
            else
                term *= cos((2n+1)*z)
            end
        end
        # theta_2(z,q) = 2 sum_{n=0}^infty q^{(n+0.5)^2} cos((2n+1)z)
        if th_k == 2
            term = 2.0 * (2n+1)^d_n * q^((n+0.5)^2) * d_cos_sign
            if d_n % 2 == 0
                term *= cos((2n+1)*z)
            else
                term *= sin((2n+1)*z)
            end
        end
        # theta_3(z,q) = 1 + 2 sum_{n=1}^infty q^{n^2} cos(2nz)
        if th_k == 3
            term = 2.0 * (2n)^d_n * q^(n^2) * d_cos_sign
            if d_n % 2 == 0
                term *= cos(2n*z)
            else
                term *= sin(2n*z)
            end
        end
        # theta_4(z,q) = 1 + 2 sum_{n=1}^infty (-1)^n q^{n^2} cos(2nz)
        if th_k == 4
            term = 2.0 * sign * (2n)^d_n * q^(n^2) * d_cos_sign
            if d_n % 2 == 0
                term *= cos(2n*z)
            else
                term *= sin(2n*z)
            end
        end

        sum += term
        n += 1
        sign = -sign

        if abs(term) < ϵ
            break
        end
    end

    return sum
end

end # module SpecialTheta