###############################################################################
###  Q special function
###############################################################################

struct QSpecial <: EllipticFunction
    ℘              ::Weierstrass
    α              ::ComplexF64
    β              ::ComplexF64
    γ1             ::ComplexF64
    γ3             ::ComplexF64
    cash           ::Dict{RationalComplex, CashedVector{ComplexF64}}
    max_derivative ::Int
end

function QSpecial(℘::Weierstrass, max_derivative::Int)
    ω1 = ℘.ω1
    ω3 = ℘.ω3
    
    w1 = 2ω1
    w3 = 2ω3
    
    η1 = ℘.η1
    η3 = (η1*ω3 - 0.5im*π)/ω1

    α = (2im*π)/(conj(w1)*w3-conj(w3)*w1)
    β = 2(η3*conj(w1)-η1*conj(w3))/(conj(w1)*w3-conj(w3)*w1)

    g2 = ℘.g2
    γ1 = (2*β*η1-g2*w1/12)/α
    γ3 = (2*β*η3-g2*w3/12)/α

    cash = Dict{RationalComplex, CashedVector{ComplexF64}}()

    QSpecial(℘, α, β, γ1, γ3, cash, max_derivative)
end

function cash_Q_elder_derivs!(output::CashedVector{ComplexF64}, Q::QSpecial, 
                              rz::RationalComplex, n_deriv::Int)
    ℘ = Q.℘
    α = Q.α
    β = Q.β

    z = raw_complex(℘, rz)
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

function cash_Q_first_derivs!(output::CashedVector{ComplexF64}, Q::QSpecial, rz::RationalComplex)
    ℘ = Q.℘
    α = Q.α
    β = Q.β
    g2 = ℘.g2

    z = raw_complex(℘, rz)
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


function Base.getindex(Q::QSpecial, rz::RationalComplex, n_deriv::Int)
    if rz == 0
        error("Can't compute Q special functions at zero point")
    end

    if rz in keys(Q.cash)
        derivs_array = Q.cash[rz]
    else
        derivs_array = CashedVector{ComplexF64}(0:Q.max_derivative)
        cash_Q_first_derivs!(derivs_array, Q, rz)
    end

    if not_cashed(derivs_array, n_deriv)
        cash_Q_elder_derivs!(derivs_array, Q, rz, n_deriv)
    end
    return derivs_array[n_deriv]
end