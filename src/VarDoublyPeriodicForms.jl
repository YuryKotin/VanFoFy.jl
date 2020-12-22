
abstract type VarDoublyPeriodicForm end

struct VarZForm <: VarDoublyPeriodicForm
    coeffs ::ComplexOffsetVector
    mother ::DoublyPeriodicFunction
end

struct VarWeierstrassForm <: VarDoublyPeriodicForm
    coeffs ::ComplexOffsetVector
    derivs ::IntOffsetVector
    pole   ::RationalComplex
    r      ::Float64
    mother ::DoublyPeriodicFunction
end

struct VarQSpecialForm <: VarDoublyPeriodicForm
    coeffs ::ComplexOffsetVector
    derivs ::IntOffsetVector
    pole   ::RationalComplex
    mother ::DoublyPeriodicFunction
end
#-----------------------------------------------------------------------------#
#=============================================================================#

function VarZForm(var_range::UnitRange, mother)
    coeffs = OffsetVector{ComplexF64}(undef, var_range)
    fill!(coeffs, 0.0im)
    VarZForm(coeffs, mother)
end

function VarWeierstrassForm(var_range, pole, r, mother)
    coeffs = OffsetVector{ComplexF64}(undef, var_range)
    derivs = OffsetVector{ComplexF64}(undef, var_range)
    fill!(coeffs, 0.0im)
    fill!(derivs, 0)
    
    VarWeierstrassForm(coeffs, derivs, pole, r, mother)
end

function VarQSpecialForm(var_range, pole, mother)
    coeffs = OffsetVector{ComplexF64}(undef, var_range)
    derivs = OffsetVector{ComplexF64}(undef, var_range)
    fill!(coeffs, 0.0im)
    fill!(derivs, 0)
    
    VarQSpecialForm(coeffs, derivs, pole, mother)
end

#=============================================================================#

@forward VarZForm.coeffs           axes, getindex, setindex!
@forward VarWeierstrassForm.coeffs axes, getindex, setindex!
@forward VarQSpecialForm.coeffs    axes, getindex, setindex!

#-----------------------------------------------------------------------------#

variables(form<:VarDoublyPeriodicForm) = UnitRange(axes(form, 1))

#=============================================================================#

conjQ(z, Q) = Q ? conj(z) : z

#=============================================================================#

function on_circle!(input::VarZForm, pole::RationalComplex, r::Float64, 
                    output::VarPolyForm; is_conj::Bool=false)
    #---------------------------------------------------
    ps = is_conj ? -1 : 1

    ξ = raw_complex(pole, input.mother)
    vars = variables(input)
    for v in vars
        output[ps * 0, v] = conjQ(input[v] * ξ, is_conj)
        output[ps * 1, v] = conjQ(input[v] * r, is_conj)
    end
end

#-----------------------------------------------------------------------------#
#=============================================================================#

function ρ(mother, Δ, n, m) end
function ϱ(mother, Δ, n, m) end

#=============================================================================#

function on_circle!(input::VarWeierstrassForm, pole::RationalComplex, r::Float64, 
                    output::VarPolyForm; is_conj::Bool=false)
    #------------------------------------------------------
    ps = is_conj ? -1 : 1
    
    M = last(powers(output))
    
    Δ = pole - input.pole
    if Δ == 0
        for v in variables(input)
            n = input.derivs[v]
            
            output[ps * (-n-2), v] = conjQ(
            input.coeffs[v] * (-1)^(n % 2), 
            is_conj)
            
            rn = r^(n+2)
            for m in 0 : M
                output[ps * m, v] = conjQ(
                    input.coeffs[v] * rn * 
                    ρ(input.mother, Δ, n, m),
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
                    ρ(input.mother, Δ, n, m),
                    is_conj)
                rm *= rΔ
            end
        end
    end
end

#-----------------------------------------------------------------------------#

function on_circle!(
    in::VarQSpecialForm, 
    pole::RationalComplex, 
    r::Float64, 
    out::VarPolyForm; 
    conjQ::Bool=false)
end
