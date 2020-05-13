###############################################################################
##### Theta functions
###############################################################################


struct Theta
    nome        ::ComplexF64
    nome_powers ::Array{ComplexF64, 2}
    powers_num  ::Int64
    ϵ           ::Float64
end

function Theta(ω1::ComplexF64,
               ω3::ComplexF64,
               ϵ ::Float64)
    q = nome(ω1, ω3)
    nome_powers = precompute_nome_powers(q, ϵ, ω1+ω3)
    powers_num = size(nome_powers, 2)

    Theta(q, nome_powers, powers_num, ϵ)
end

"""
  nome(ω1::ComplexF64, ω3::ComplexF64)

Return the nome for the giving pair of lattice generators.
"""
function nome(
        ω1::ComplexF64,
        ω3::ComplexF64
        )
    #####
    τ = ω3 / ω1
    q = exp(1im * π * τ)
end


"""
  precompute_nome_powers(; q::ComplexF64, ϵ::Float64, z0::ComplexF64, max_n_terms ::Integer = 10)

Create array of precomputed nome powers for evaluating theta functions with given tolerance.

# RETURN:
- Array{ComplexF64, 2}(2, n)
"""
function precompute_nome_powers(
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
    nome_powers = Array{ComplexF64, 2}(undef, 2, n+1)
    for k in 0 : n
        nome_powers[1, k+1] = q^((k+0.5)^2)
        nome_powers[2, k+1] = nome_powers[1, k+1] * (2k+1)
    end

    return nome_powers
end

@doc """
  theta12(th::Theta, z::ComplexF64)

Compute θ_1(z), θ_2(z), θ_1'(z), θ_2'(z) with given precomputed array of nome powers.
"""
function theta12(
        th ::Theta,
        z  ::ComplexF64
        )
    #####
    t1  :: ComplexF64 = 0
    t2  :: ComplexF64 = 0
    dt1 :: ComplexF64 = 0
    dt2 :: ComplexF64 = 0
    sign = 1.0

    N = th.powers_num-1
    for n in 0 : N
        z_n = (2n+1)*z
        sin_z = sin(z_n)
        cos_z = cos(z_n)
        # theta_1(z,q) = 2 sum_{n=0}^N (-1)^n q^{(n+0.5)^2} sin((2n+1)z)
        t1 += 2 * sign * th.nome_powers[1, n+1] * sin_z
        # theta_2(z,q) = 2 sum_{n=0}^N q^{(n+0.5)^2} cos((2n+1)z)
        t2 += 2 * th.nome_powers[1, n+1] * cos_z
        # theta_1'(z,q) = 2 sum_{n=0}^N (-1)^n (2n+1) q^{(n+0.5)^2} cos((2n+1)z))
        dt1 += 2 * sign * th.nome_powers[2, n+1] * cos_z
        # theta_2'(z,q) = -2 sum_{n=0}^N (2n+1) q^{(n+0.5)^2} sin((2n+1)z)
        dt2 += -2 * th.nome_powers[2, n+1] * sin_z

        sign = -sign
    end

    return t1, t2, dt1, dt2
end


""""
  theta(θ::Theta; th_k::Int, d_n::Int,z::ComplexF64)

Compute ``\\frac{d^n}{dz^n} \\theta_k(z, q)`` with tolerance `ϵ`.
"""
function theta(
        # Theta data struct
        th ::Theta;
        # theta-function number
        th_k ::Int,
        # derivative number
        d_n  ::Int,
        # point for computation
        z    ::ComplexF64,
        )
   #####
    q = th.nome
    ϵ = th.ϵ

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
