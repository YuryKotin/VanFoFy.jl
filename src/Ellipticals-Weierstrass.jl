###############################################################################
##### Weierstrass elliptic functions
###############################################################################


"""
Set of Weierstrass elliptic functions parameters for computational issues
"""
struct Weierstrass <: EllipticFunction
    ω1            ::  ComplexF64
    ω3            ::  ComplexF64
    e1            ::  ComplexF64
    g2            ::  ComplexF64
    g3            ::  ComplexF64
    η1            ::  ComplexF64
    # (\\frac{π θ_3(0) θ_4(0)}{2ω_1})^2
    # coefficient in ℘-function evaluation
    ℘_factor      ::  ComplexF64
    σ_factor      ::  ComplexF64
    # θ-function data
    θ             ::  Theta
    cash          ::  Dict{RationalComplex, CashedVector{ComplexF64}}
    max_derivative ::Int
end

"""
  WeierstrassData(; ω1::ComplexF64, ω3::ComplexF64)

Construct WeierstrassData object from lattice generators ω1 and ω3
"""
function Weierstrass(
        ω1::ComplexF64,
        ω3::ComplexF64,
        max_derivative ::Int
        )
    #####
    ϵ = eps()
    θ = Theta(ω1, ω3, ϵ)

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

    cash = Dict{RationalComplex, CashedVector{ComplexF64}}()

    return Weierstrass(ω1, ω3, e1, g2, g3, η1, ℘_factor, σ_factor, θ, cash, max_derivative)
end

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

raw_complex(w::Weierstrass, rz::RationalComplex) = real(rz)*w.ω1 + imag(rz)*w.ω3

function cash_first_derivs!(w::Weierstrass, rz::RationalComplex, output::CashedVector{ComplexF64})
    ω1 = w.ω1
    e1 = w.e1
    η1 = w.η1
    g2 = w.g2
    ℘_factor = w.℘_factor
    σ_factor = w.σ_factor

    z = raw_complex(w,rz)
    norm_factor = abs(z)

    u = (π / (2ω1)) * z
    th1, th2, dth1, dth2 = theta(w.θ, u)
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


function Base.getindex(w::Weierstrass, rz::RationalComplex, n_deriv::Int)
    if rz == 0
        error("Can't compute Weierstrass elliptic functions at zero point")
    end

    if rz in keys(w.cash)
        vector = w.cash[rz]
        sign = 1.0
    elseif -rz in keys(w.cash)
        vector = w.cash[-rz]
        sign = (n_deriv % 2 == 0) ? 1.0 : -1.0
    else
        vector = CashedVector{ComplexF64}(-1:w.max_derivative)
        cash_first_derivs!(w, rz, vector)
        sign = 1.0
    end

    if not_cashed(vector, n_deriv)
        cash_elder_derivs!(vector, n_deriv)
    end
    return vector[n_deriv]*sign
end

#=
"""
weierstrass_normalized!(; z::ComplexF64, norm_factor::Float64, w::WeierstrassData,
        upper_term::Int64, output::ComplexOffsetVector)

Compute sequence of normalized weierstrass elliptic functions derivatives

# Return
output[-2] = ``\\wp^{(-2)}(z)``,
output[-1] = ``\\wp^{(-1)}(z)``,
output[n]  = ``\\frac{r^{n+2}}{(n+1)!} \\wp^{(n)}(z)``, n ≥ 0.
"""
function weierstrass_normalized!(;
        z           ::ComplexF64,
        norm_factor ::Float64,
        w           ::Weierstrass,
        upper_term  ::Int64,
        output      ::ComplexOffsetVector
        )
    #####
    ω1 = w.ω1
    e1 = w.e1
    η1 = w.η1
    g2 = w.g2
    ℘_factor = w.℘_factor
    σ_factor = w.σ_factor

    u = (π / (2ω1)) * z
    th1, th2, dth1, dth2 = theta(w.θ, u)
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

    # wp_{n+2}(z) = frac{6}{(n+1)(n+2)} sum_{k=0}^n wp_k(z) wp_{n-k}(z),
    # n >= 1.
    for n in 1 : upper_term-2
        output[n+2] = 0.0im
        for k in 0 : n
            output[n+2] += output[k] * output[n-k]
        end
        output[n+2] *= 6.0 / ((n+1) * (n+2))
    end

    # for n in 1 : upper_term
    #     output[n] /= n+1
    # end

    output[-2] = -(η1 * (z^2) / (2ω1) + log(σ_factor * th1))

end


function expand_at_pole!(;
        output           ::ComplexOffsetVector,
        weierstrass_data ::Weierstrass,
        upper_term       ::Int64
        )
    #####
    if upper_term < 2
        error("Too few terms of expansion")
    end

    g2 = weierstrass_data.g2
    g3 = weierstrass_data.g3
    
    n_c = div(upper_term+2, 2)
    c = OffsetVector{ComplexF64}(undef, 2:n_c)

    c[2] = g2 / 20.0

    if n_c >= 3
        c[3] = g3 / 28.0
    end

    for n in 4 : n_c
        c[n] = 0
        for m in 2 : n-2
            c[n] += c[m] * c[n-m]
        end
        c[n] *= 3 / ((2n+1)*(n-3))
    end

    for n in 0 : upper_term
        output[n] = 0
    end
    for n in 2 : n_c
        output[2n-2] = c[n]
    end

    output[-1] = 0
    output[-2] = 0
end


function term_expansion_normalized!(;
        output         ::ComplexOffsetVector,
        dict_of_series ::RComplex2OffsetMatrix,
        delta          ::RationalComplex,
        n_derivative   ::Int64,
        upper_term     ::Int64,
        r_source       ::Float64,
        r_dest         ::Float64
        )

    # При решении задачи нужно вычислять разложения функций Вейерштрасса
    # в противоположных точках. Чтобы не делать двойную работу, значения вычисляются
    # только для z с положительной действительной частью, а для противоположного -
    # нечетные члены рядов берутся с противоположным знаком.
    if real(delta) < 0
        r = -r_dest  # (-1)^n r^n = (-r)^n
        delta = -delta
    else
        r = r_dest
        # delta = delta
    end

    # Массив разложений в ряд Тейлора производных ℘-функции Вейерштрасса
    # в точке delta.
    α = dict_of_series[delta]

    if delta != 0
        factor = (r_source / abs(delta))^(n_derivative+2)
    else
        factor = r_dest^(n_derivative+2)
    end

    for n in 0 : upper_term
        output[n] = factor * α[n, n_derivative]
        factor *= r
    end
end

=#