
abstract type FunctionForm end

abstract type PeriodicForm <: FunctionForm end

struct ZForm <: FunctionForm
    coeffs ::ComplexOffsetVector
end

struct WeierstrassForm <: PeriodicForm
    coeffs ::ComplexOffsetVector
    derivs ::IntOffsetVector
    pole   ::RationalComplex
    r      ::Float64
end

struct QSpecialForm <: PeriodicForm
    coeffs ::ComplexOffsetVector
    derivs ::IntOffsetVector
    pole   ::RationalComplex
end

#=============================================================================#

function ZForm(var_range::UnitRange)
    coeffs = OffsetVector{ComplexF64}(undef, var_range)
    fill!(coeffs, 0.0im)
    ZForm(coeffs)
end

function WeierstrassForm(var_range, pole, r)
    coeffs = OffsetVector{ComplexF64}(undef, var_range)
    derivs = OffsetVector{ComplexF64}(undef, var_range)
    fill!(coeffs, 0.0im)
    fill!(derivs, 0)
    
    WeierstrassForm(coeffs, derivs, pole, r)
end

function QSpecialForm(var_range, pole)
    coeffs = OffsetVector{ComplexF64}(undef, var_range)
    derivs = OffsetVector{ComplexF64}(undef, var_range)
    fill!(coeffs, 0.0im)
    fill!(derivs, 0)
    
    QSpecialForm(coeffs, derivs, pole)
end

#=============================================================================#

@forward ZForm.coeffs           axes, getindex, setindex!
@forward WeierstrassForm.coeffs axes, getindex, setindex!
@forward QSpecialForm.coeffs    axes, getindex, setindex!

#-----------------------------------------------------------------------------#

variables(form<:FunctionForm) = UnitRange(axes(form, 1))

#=============================================================================#

conjQ(z, Q) = Q ? conj(z) : z

#=============================================================================#

function series_on_circle!(input::ZForm, pole::RationalComplex, r::Float64, 
                    output::PolynomForm, instance::PeriodicInstance,
                    is_conj::Bool=false)
    #---------------------------------------------------
    ps = is_conj ? -1 : 1

    ξ = raw_complex(pole, instance)
    vars = variables(input)
    for v in vars
        output[ps * 0, v] = conjQ(input[v] * ξ, is_conj)
        output[ps * 1, v] = conjQ(input[v] * r, is_conj)
    end
end

#=============================================================================#

function taylor(periodic<:PeriodicForm, Δ, instance) end

#=============================================================================#

function series_on_circle_pos!(input<:PeriodicForm,
                            pole::RationalComplex, r::Float64, 
                            output::PolynomForm, instance::PeriodicInstance,
                            is_conj::Bool=false)
    #------------------------------------------------------
    ps = is_conj ? -1 : 1
    M  = is_conj ? -first(powers(output)) : last(powers(output))
    
    Δ = pole - input.pole
    rho = taylor(input, Δ, instance)
    if Δ == 0
        for v in variables(input)
            n = input.derivs[v]
            
            rn = r^(n+2)
            for m in 0 : M
                val = input[v] * rn * rho[n,m]
                output[ps * m, v] = conjQ(val, is_conj)
                rn *= r
            end
        end
    else # Δ ≠ 0
        absΔ = abs(raw_complex(Δ, input.mother))
        for v in variables(input)
            n = input.derivs[v]
            Rn = (input.r / absΔ)^(n+2)
            rΔ = r / absΔ
            rm = 1.0
            for m in 0 : M
                val = input[v] * Rn * rm * rho[n,m]
                output[ps * m, v] = conjQ(val, is_conj)
                rm *= rΔ
            end
        end
    end
end

function series_on_circle!(input::WeierstrassForm, pole::RationalComplex, r::Float64, 
                    output::PolynomForm, special::PeriodicInstance; is_conj::Bool=false)
    #------------------------------------------------------
    series_on_circle_pos!(input, pole, r, output, special, is_conj)
    
    ps = is_conj ? -1 : 1
    
    if pole == input.pole
        for v in variables(input)
            n = input.derivs[v]
            val = input[v] * (-1)^(n % 2)
            output[ps * (-n-2), v] = conjQ(val, is_conj)
        end
    end
end

#-----------------------------------------------------------------------------#

function series_on_circle!(input::QSpecialForm,
                    pole::RationalComplex, r::Float64, 
                    output::PolynomForm, special::PeriodicInstance,
                    is_conj::Bool=false)
    #------------------------------------------------------
    series_on_circle_pos!(input, pole, r, output, special, is_conj)
end
