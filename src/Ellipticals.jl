module Ellipticals

using OffsetArrays

using ..TypeSynonyms: ComplexOffsetMatrix, ComplexOffsetVector, RComplex2IntDict
using ..TypeSynonyms: RationalComplex, RComplex2OffsetMatrix
using ..Input: FiberData, CellData

#####
##### Theta functions
#####

@doc """
  nome(; ω1::ComplexF64, ω3::ComplexF64)

Return the nome for the giving pair of lattice generators.
"""
function nome(;
        ω1::ComplexF64,
        ω3::ComplexF64
        )
    #####
    τ = ω3 / ω1
    q = exp(1im * π * τ)
end


@doc """"
  theta(; th_k::Integer, d_n::Integer,z::ComplexF64, q::ComplexF64, ϵ::Float64)

Compute ``\\frac{d^n}{dz^n} \\theta_k(z, q)`` with tolerance `ϵ`.
"""
function theta(;
        # theta-function number
        th_k ::Integer,
        # derivative number
        d_n  ::Integer,
        # point for computation
        z    ::ComplexF64,
        # nome
        q    ::ComplexF64,
        # tolerance
        ϵ    ::Float64)
   #####
   if (th_k < 1) || (th_k > 4)
       error("Invalid theta function number")
   end

   rem_d = d_n % 4
   d_sin_sign = (rem_d == 0) || (rem_d == 1) ? 1 : -1
   d_cos_sign = (rem_d == 0) || (rem_d == 3) ? 1 : -1

   res :: ComplexF64 = 0
   if (th_k == 1) || (th_k == 2)
       res = 0
       n = 0
   else
       n = 1
       if d_n == 0
           res = 1
       else
           res = 0
       end
   end

   sign = (-1)^n
   while true
       # theta_1(z,q) = 2 sum_{n=0}^infty (-1)^n q^{(n+0.5)^2} sin((2n+1)z)
       if th_k == 1
           term = 2 * sign * (2n+1)^d_n * q^((n+0.5)^2) * d_sin_sign
           if d_n % 2 == 0
               term *= sin((2n+1)*z)
           else
               term *= cos((2n+1)*z)
           end
       end
       # theta_2(z,q) = 2 sum_{n=0}^infty q^{(n+0.5)^2} cos((2n+1)z)
       if th_k == 2
           term = 2 * (2n+1)^d_n * q^((n+0.5)^2) * d_cos_sign
           if d_n % 2 == 0
               term *= cos((2n+1)*z)
           else
               term *= sin((2n+1)*z)
           end
       end
       # theta_3(z,q) = 1 + 2 sum_{n=1}^infty q^{n^2} cos(2nz)
       if th_k == 3
           term = 2 * (2n)^d_n * q^(n^2) * d_cos_sign
           if d_n % 2 == 0
               term *= cos(2n*z)
           else
               term *= sin(2n*z)
           end
       end
       # theta_4(z,q) = 1 + 2 sum_{n=1}^infty (-1)^n q^{n^2} cos(2nz)
       if th_k == 4
           term = 2 * sign * (2n)^d_n * q^(n^2) * d_cos_sign
           if d_n % 2 == 0
               term *= cos(2n*z)
           else
               term *= sin(2n*z)
           end
       end

       res += term
       n += 1
       sign = -sign

       if abs(term) < ϵ
           break
       end
   end

   return res
end

@doc """
  precompute_nome_powers(; q::ComplexF64, ϵ::Float64, z0::ComplexF64, max_n_terms ::Integer = 10)

Create array of precomputed nome powers for evaluating theta functions with given tolerance.

# RETURN:
- OffsetArray{ComplexF64}(undef, 1:2, 0:n)
"""
function precompute_nome_powers(;
        # nome
        q           ::ComplexF64,
        # tolerance
        ϵ           ::Float64,
        # point for tolerance checking
        z0          ::ComplexF64,
        # maximum number of terms
        max_n_terms ::Integer = 10
        )
    #####
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
    nome_powers = OffsetArray{ComplexF64}(undef, 1:2, 0:n)
    for k in 0 : n
        nome_powers[1, k] = q^((k+0.5)^2)
        nome_powers[2, k] = nome_powers[1, k] * (2k+1)
    end
    return nome_powers
end

@doc """
  theta_1d2d(; z::ComplexF64, nome_powers ::ComplexOffsetMatrix)

Compute θ_1(z), θ_2(z), θ_1'(z), θ_2'(z) with given precomputed array of nome powers.
"""
function theta_1d2d(;
        z           ::ComplexF64,
        # precomputed array of nome powers
        nome_powers ::ComplexOffsetMatrix
        )
    #####
    t1  :: ComplexF64 = 0
    t2  :: ComplexF64 = 0
    dt1 :: ComplexF64 = 0
    dt2 :: ComplexF64 = 0
    sign = 1.0

    N = lastindex(nome_powers, 2)
    for n in 0 : N
        z_n = (2n+1)*z
        sin_z = sin(z_n)
        cos_z = cos(z_n)
        # theta_1(z,q) = 2 sum_{n=0}^N (-1)^n q^{(n+0.5)^2} sin((2n+1)z)
        t1 += 2 * sign * nome_powers[1, n] * sin_z
        # theta_2(z,q) = 2 sum_{n=0}^N q^{(n+0.5)^2} cos((2n+1)z)
        t2 += 2 * nome_powers[1, n] * cos_z
        # theta_1'(z,q) = 2 sum_{n=0}^N (-1)^n (2n+1) q^{(n+0.5)^2} cos((2n+1)z))
        dt1 += 2 * sign * nome_powers[2, n] * cos_z
        # theta_2'(z,q) = -2 sum_{n=0}^N (2n+1) q^{(n+0.5)^2} sin((2n+1)z)
        dt2 += -2 * nome_powers[2, n] * sin_z

        sign = -sign
    end

    return t1, t2, dt1, dt2
end

#####
##### Weierstrass elliptic functions
#####


@doc """
Set of Weierstrass elliptic functions parameters for computational issues
"""
struct WeierstrassData
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
    # precomputed array of nome powers
    nome_powers   ::  OffsetArray{Complex{Float64},2,Array{Complex{Float64},2}}
end

@doc """
  WeierstrassData(; ω1::ComplexF64, ω3::ComplexF64)

Construct WeierstrassData object from lattice generators ω1 and ω3
"""
function WeierstrassData(;
        ω1::ComplexF64,
        ω3::ComplexF64
        )
    #####
    q = nome(ω1=ω1, ω3=ω3)
    tolerance = eps(imag(ω3))
    z_theta = complex(0.0, imag(ω3))
    nome_powers = precompute_nome_powers(q=q, ϵ=tolerance, z0=z_theta)

    th2 = theta(th_k=2, d_n=0, z=0.0im, q=q, ϵ=tolerance)
    th4 = theta(th_k=4, d_n=0, z=0.0im, q=q, ϵ=tolerance)
    th24 = th2^4
    th44 = th4^4
    e_factor = π^2 / (12 * ω1^2)
    e1 = e_factor * (th24 + 2th44)
    e2 = e_factor * (th24 - th44)
    e3 = -e_factor * (2th24 + th44)

    g2 = -4 * (e2*e3 + e3*e1 + e1*e2)
    g3 = 4 * e1 * e2 * e3

    d3th1 = theta(th_k=1, d_n=3, z=0.0im, q=q, ϵ=tolerance)
    d1th1 = theta(th_k=1, d_n=1, z=0.0im, q=q, ϵ=tolerance)

    η1 = -π^2 / (12ω1) * d3th1 / d1th1
    η3 = (η1 * ω3 - 1im*π/2) / ω1

    th3 = theta(th_k=3, d_n=0, z=0.0im, q=q, ϵ=tolerance)
    ℘_factor = (π * th3 * th4 / (2ω1))^2

    σ_factor = 2ω1 / (π*d1th1)

    return WeierstrassData(ω1, ω3, e1, g2, g3, η1, ℘_factor, σ_factor, nome_powers)
end

@doc """
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
        w           ::WeierstrassData,
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
    th1, th2, dth1, dth2 = theta_1d2d(z=u, nome_powers=w.nome_powers)
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

@doc """
  derivative_series_from_initial!(; initial_series::ComplexOffsetVector, 
    derivatives_series ::ComplexOffsetMatrix)

Transform initial series into array of derivatives' and integrals' series

Given long Taylor series
``f(z) = \\sum_{n=0}^N a_n z^n``
generates array of series
``\\frac{1}{k!}f^{(k)}(z) = \\sum_{n=0}^M b_{k,n} z^n, 0 < z ≤ K_2``,
``f^{(k)}(z) = \\sum_{n=0}^M b_{k,n} z^n, K_1 ≤ z < 0``,
where K1 and K2 are first and last indices if *derivative_series*.

Assuming  N ≥ M + (K_2 - K_1 + 1).
"""
function derivative_series_from_initial!(;
        initial_series     ::ComplexOffsetVector,
        derivatives_series ::ComplexOffsetMatrix
        )
    #####
    upper_term       = lastindex(derivatives_series, 1)
    upper_derivative = lastindex(derivatives_series, 2)
    lower_derivative = firstindex(derivatives_series, 2)

    #derivatives_series[0 : upper_term, 0] .= initial_series[0 : upper_term]
    for n in 0 : upper_term
        derivatives_series[n, 0] = initial_series[n]
    end

    for k in 1 : upper_derivative
        for n in 1 : (upper_term + upper_derivative - k + 1)
            initial_series[n-1] = initial_series[n] * n / (k+1)
        end
        initial_series[upper_term + upper_derivative - k + 1] = .0im

        #@views derivatives_series[:, k] .= initial_series[0 : upper_term]
        for n in 0 : upper_term
            derivatives_series[n, k] = initial_series[n]
        end
    end

    for k in -1 : -1 : lower_derivative
        derivatives_series[0, k] = initial_series[k]
        for n in 1 : upper_term
            derivatives_series[n, k] = derivatives_series[n-1, k+1] / n
        end
    end
end

function expand_at_pole!(;
        output           ::ComplexOffsetVector,
        weierstrass_data ::WeierstrassData,
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

function construct_dict_of_wei_series(;
        delta_to_nterms ::RComplex2IntDict,
        wei             ::WeierstrassData
        )
    #####
    max_n_terms = 0
    for v in values(delta_to_nterms)
        if v[1] > max_n_terms
            max_n_terms = v[1]
        end
    end

    wei_sequence = OffsetArray{ComplexF64}(undef, -2 : 2*max_n_terms)

    wei_series_dict = Dict{RationalComplex, ComplexOffsetMatrix}()

    for delta in keys(delta_to_nterms)
        upper_term = delta_to_nterms[delta]

        # Для вычисления разложений в ряд производных нужно разложить исходную
        # функцию с запасом членов.
        # Всего нужно n производных по n членов в каждой.
        diff_upper_term = 2 * upper_term

        if delta != 0
            z = real(delta)*wei.ω1 + imag(delta)*wei.ω3
            norm_factor = abs(z)

            weierstrass_normalized!(z=z,
                                 norm_factor=norm_factor,
                                 w=wei,
                                 upper_term=diff_upper_term,
                                 output=wei_sequence)
        else
            expand_at_pole!(output=wei_sequence,
                            weierstrass_data=wei,
                            upper_term=diff_upper_term)
        end

        series_array = OffsetArray{ComplexF64}(undef, 0:upper_term, -1:upper_term-2)

        derivative_series_from_initial!(initial_series=wei_sequence,
                                        derivatives_series=series_array)

        wei_series_dict[delta] = series_array
    end

    return wei_series_dict
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

function construct_wei_deltas(fibers::Vector{FiberData})
    delta_to_nterms = RComplex2IntDict()
    delta_to_nterms[0] = 0

    n_fibers = length(fibers)
    for i in 1 : n_fibers
        z1 = fibers[i].center
        n1 = fibers[i].n_terms
        
        for j in i+1 : n_fibers
            z2= fibers[j].center
            delta = z2 - z1
            if real(delta) < 0
                delta = -delta
            end

            n_terms = get(delta_to_nterms, delta, -1)
            n2 = fibers[j].n_terms
            n = max(n1, n2)
            if n_terms < n
                delta_to_nterms[delta] = n
            end
        end

        if delta_to_nterms[0] < n1
            delta_to_nterms[0] = n1
        end
    end

    return delta_to_nterms
end

struct Weierstrass
    data :: WeierstrassData
    expansions :: RComplex2OffsetMatrix
end

function Weierstrass(cell::CellData)
    ω1 = complex(cell.l1 / 2.0)
    ω3 = (cell.l2/2.0) * exp(1.0im * cell.γ)
    data = WeierstrassData(ω1, ω3)

    wei_deltas = construct_wei_deltas(cell.fibers)

    expansions = construct_dict_of_wei_series(wei_deltas, data)

    return Weierstrass(data, expansions)
end

struct Q_special
    wei :: Weierstrass
    expansions :: RComplex2OffsetMatrix
end

function Q_special(wei::Weierstrass)

end

end # module Ellipticals