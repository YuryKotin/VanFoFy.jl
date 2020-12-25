
abstract type VarFunctionForm end

abstract type VarDoublyPeriodicForm <: VarFunctionForm end

struct VarZForm <: VarFunctionForm
    coeffs ::ComplexOffsetVector
end

struct VarWeierstrassForm <: VarDoublyPeriodicForm
    coeffs ::ComplexOffsetVector
    derivs ::IntOffsetVector
    pole   ::RationalComplex
    r      ::Float64
end

struct VarQSpecialForm <: VarDoublyPeriodicForm
    coeffs ::ComplexOffsetVector
    derivs ::IntOffsetVector
    pole   ::RationalComplex
end
#-----------------------------------------------------------------------------#
#=============================================================================#

function VarZForm(var_range::UnitRange)
    coeffs = OffsetVector{ComplexF64}(undef, var_range)
    fill!(coeffs, 0.0im)
    VarZForm(coeffs)
end

function VarWeierstrassForm(var_range, pole, r)
    coeffs = OffsetVector{ComplexF64}(undef, var_range)
    derivs = OffsetVector{ComplexF64}(undef, var_range)
    fill!(coeffs, 0.0im)
    fill!(derivs, 0)
    
    VarWeierstrassForm(coeffs, derivs, pole, r)
end

function VarQSpecialForm(var_range, pole)
    coeffs = OffsetVector{ComplexF64}(undef, var_range)
    derivs = OffsetVector{ComplexF64}(undef, var_range)
    fill!(coeffs, 0.0im)
    fill!(derivs, 0)
    
    VarQSpecialForm(coeffs, derivs, pole)
end

#=============================================================================#

@forward VarZForm.coeffs           axes, getindex, setindex!
@forward VarWeierstrassForm.coeffs axes, getindex, setindex!
@forward VarQSpecialForm.coeffs    axes, getindex, setindex!

#-----------------------------------------------------------------------------#

variables(form<:VarFunctionForm) = UnitRange(axes(form, 1))

#=============================================================================#

conjQ(z, Q) = Q ? conj(z) : z

#=============================================================================#

function on_circle!(input::VarZForm, pole::RationalComplex, r::Float64, 
                    output::VarPolyForm, instance::PeriodicInstance,
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

#-----------------------------------------------------------------------------#
#=============================================================================#

function taylor(periodic<:VarDoublyPeriodicForm, Δ, instance) end

#=============================================================================#

function on_circle_pos_pow!(input<:VarDoublyPeriodicForm,
                            pole::RationalComplex, r::Float64, 
                            output::VarPolyForm, instance::PeriodicInstance,
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
                output[ps * m, v] = conjQ(
                    input.coeffs[v] * rn * 
                    rho[m,n],
                    is_conj)
                rn *= r
            end
        end
    else # Δ ≠ 0
        absΔ = abs(raw_complex(Δ, input.mother))
        for v in variables
            n = input.derivs[v]
            Rn = (input.r / absΔ)^(n+2)
            rΔ = r / absΔ
            rm = 1.0
            for m in 0 : M
                output[ps * m, v] = conjQ(
                    input.coeffs[v] * Rn * rm * 
                    rho[m,n],
                    is_conj)
                rm *= rΔ
            end
        end
    end
end

function on_circle!(input::VarWeierstrassForm, pole::RationalComplex, r::Float64, 
                    output::VarPolyForm, instance::PeriodicInstance; is_conj::Bool=false)
    #------------------------------------------------------
    on_circle_pos_pow!(input, pole, r, output, instance, is_conj)
    
    ps = is_conj ? -1 : 1
    
    Δ = pole - input.pole
    if Δ == 0
        for v in variables(input)
            n = input.derivs[v]
            
            output[ps * (-n-2), v] = conjQ(
                input.coeffs[v] * (-1)^(n % 2), 
                is_conj)
        end
    end
end

#-----------------------------------------------------------------------------#

function on_circle!(input::VarQSpecialForm; 
                    pole::RationalComplex, r::Float64, 
                    output::VarPolyForm, instance::PeriodicInstance,
                    is_conj::Bool=false)
    #------------------------------------------------------
    on_circle_pos_pow!(input, pole, r, output, instance, is_conj)
end
