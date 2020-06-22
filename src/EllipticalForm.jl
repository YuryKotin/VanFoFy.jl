abstract type EllipticalTerm end

struct WeierstrassTerm <: EllipticalTerm
    deriv      ::Int
    pole       ::RationalComplex
    num_factor ::ComplexF64
    norm_r     ::Float64
end

struct QSpecialTerm <: EllipticalTerm
    deriv      ::Int
    pole       ::RationalComplex
    num_factor ::ComplexF64
    norm_r     ::Float64
end

struct ZTerm <: EllipticalTerm
    num_factor ::ComplexF64
end

struct ConstTerm <: EllipticalTerm
    num_factor ::ComplexF64
end

###############################################################################

function differentiate(term::EllipticalTerm) end

function differentiate(term::WeierstrassTerm)
    # (r^{n+2}/(n+1)! ℘^{(n)}(z))' = [(n+2)/r] [r^(n+3)/(n+2)! ℘^{(n+1)}(z)]
    deriv = term.deriv + 1
    pole = term.pole
    num_factor = term.num_factor * (term.deriv + 2) / term.norm_r
    norm_r = term.norm_r
    WeierstrassTerm(deriv, pole, num_factor, norm_r)
end

function differentiate(term::QSpecialTerm)
    # (r^{n+3}/(n+1)! Q^{(n)}(z))' = [(n+2)/r] [r^(n+4)/(n+2)! Q^{(n+1)}(z)]
    deriv = term.deriv + 1
    pole = term.pole
    num_factor = term.num_factor * (term.deriv + 2) / term.norm_r
    norm_r = term.norm_r
    QSpecialTerm(deriv, pole, num_factor, norm_r)
end

differentiate(term::ZTerm) = ConstTerm(term.num_factor)

differentiate(term::ConstTerm) = ConstTerm(0.0im)

###############################################################################

const EllipticalLinForm = OffsetArray{EllipticalTerm, 1, Array{EllipticalTerm, 1}}

function differentiate(form::EllipticalLinForm)
    d = OffsetVector{EllipticalTerm}(undef, axes(form)[1])
    for i in eachindex(form)
        d[i] = differentiate(form[i])
    end
    return d
end

###############################################################################

function add_term_series!(output       ::BoundedVector{ComplexF64}, 
                            term       ::WeierstrassTerm; 
                            point      ::RationalComplex,
                            factor     ::ComplexF64 = 1.0+0.0im, 
                            norm_r     ::Float64,
                            power_shift::Int = 0, 
                            conjugated ::Bool = false,
                            praecursor ::EllipticPraecursor)
    
    ℘ = praecursor.℘
    Δ = point - term.pole
    if Δ == 0

        r = term.norm_r
        R = norm_r
        n = term.deriv
        
        addend = n % 2 == 0 ? 1.0+0.0im : -1.0+0.0im
        addend *= (r / R)^(n+2) * term.num_factor
        if !conjugated
            output[-(n+2)+power_shift] += addend * factor
        else
            output[n+2+power_shift] += conj(addend) * factor
        end
        
        rn = r^(n+2)
        Rk = 1.0
        K = lastindex(output)
        for k in 0 : K
            addend = term.num_factor * rn * Rk * ℘.derivs_series[n,k]
            if !conjugated
                output[k+power_shift] += addend * factor
            else
                output[-k+power_shift] += conj(addend) * factor
            end
            Rk *= R
        end

    else        
        z = raw_complex(℘.lattice, Δ)
        rΔ = term.norm_r / abs(z)
        RΔ = norm_r / abs(z)
        K = lastindex(output)
        n = term.deriv
        
        rΔ_n = (n == -1 ? rΔ : rΔ^(n+2) / (n+1))
        RΔ_k = 1.0
        binom = 1.0
        
        for k in 0 : K
            addend = term.num_factor * rΔ_n * RΔ_k * ℘[Δ, n+k]
            if n == -1
                addend /= (k == 0 ? 1 : k)
            else
                addend *= binom
            end

            if !conjugated
                output[k+power_shift] += addend * factor
            else
                output[-k+power_shift] += conj(addend) * factor
            end
            
            binom *= (n+k+1)/(k+1)
            RΔ_k *= RΔ
        end
    end
end

#######################################

function add_term_series!(output::BoundedVector{ComplexF64}, 
                        term       ::QSpecialTerm; 
                        point      ::RationalComplex, 
                        factor     ::ComplexF64 = 1.0+0.0im, 
                        norm_r     ::Float64,
                        power_shift::Int = 0, 
                        conjugated ::Bool = false,
                        praecursor ::EllipticPraecursor)
    Q = praecursor.Q
    ℘ = praecursor.℘
    Δ = point - term.pole

    r = term.norm_r
    R = norm_r
    K = lastindex(output)
    n = term.deriv
    if Δ == 0
        rn = r^(n+3)
        Rk = 1.0
        for k in 0 : K
            addend = term.num_factor * rn * Rk * Q.derivs_series[n,k]
            if !conjugated
                output[k+power_shift] += addend * factor
            else
                output[-k+power_shift] += conj(addend) * factor
            end
            Rk *= R
        end
    else
        z = raw_complex(Q.lattice, Δ)
        rΔ = r / abs(z)
        RΔ = R / abs(z)

        rΔ_n = rΔ^(n+3)
        RΔ_k = 1.0
        binom = 1.0
        
        for k in 0 : K
            addend = term.num_factor * rΔ_n * RΔ_k * Q[Δ, n+k]
            addend *= binom

            if !conjugated
                output[k+power_shift] += addend * factor
            else
                output[-k+power_shift] += conj(addend) * factor
            end
            
            binom *= (n+k+2)/(k+1)
            RΔ_k *= RΔ
        end
    end
end

#######################################

function add_term_series!(output::BoundedVector{ComplexF64}, 
                        term       ::ZTerm; 
                        point      ::RationalComplex, 
                        factor     ::ComplexF64 = 1.0+0.0im, 
                        norm_r     ::Float64,
                        power_shift::Int = 0, 
                        conjugated ::Bool = false,
                        praecursor ::EllipticPraecursor)
    z = raw_complex(praecursor.l, point)
    if conjugated
        output[0+power_shift] += conj(z * term.num_factor) * factor
        output[1+power_shift] += (abs(z)^2 / norm_r) * conj(term.num_factor) * factor
    else
        output[0+power_shift] += z * term.num_factor * factor
        output[1+power_shift] += norm_r * term.num_factor * factor
    end
end

#######################################

function add_term_series!(output::BoundedVector{ComplexF64}, 
                        term       ::ConstTerm; 
                        point      ::RationalComplex, 
                        factor     ::ComplexF64 = 1.0+0.0im, 
                        norm_r     ::Float64,
                        power_shift::Int = 0, 
                        conjugated ::Bool = false,
                        praecursor ::EllipticPraecursor)
    if conjugated
        output[0+power_shift] += conj(term.num_factor) * factor
    else
        output[0+power_shift] += term.num_factor * factor
    end
end