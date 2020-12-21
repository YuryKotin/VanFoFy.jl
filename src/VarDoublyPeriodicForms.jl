
abstract type VarDoublyPeriodicForm end

struct VarZForm <: VarDoublyPeriodicForm
    coeffs ::ComplexOffsetVector
    mother ::DoublyPeriodicFunction
end

struct VarWeierstrassForm <: VarDoublyPeriodicForm
    coeffs ::ComplexOffsetVector
    derivs ::IntOffsetVector
    pole   ::RationalComplex
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

function VarWeierstrassForm(var_range, pole, mother)
    coeffs = OffsetVector{ComplexF64}(undef, var_range)
    derivs = OffsetVector{ComplexF64}(undef, var_range)
    fill!(coeffs, 0.0im)
    fill!(derivs, 0)
    
    VarWeierstrassForm(coeffs, derivs, pole, mother)
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

function on_circle!(input::VarZForm, pole::RationalComplex, r::Float64, 
                    output::VarPolyForm; conjQ::Bool=false)
    #---------------------------------------------------
    ξ = raw_complex(pole, input.mother)
    
    vars = variables(input)
    @views begin
        if conjQ
            @. output[vars, 0] = conj(input[vars] * ξ)
            @. output[vars,-1] = conj(input[vars] * r)
        else
            @. output[vars, 0] =      input[vars] * ξ
            @. output[vars, 1] =      input[vars] * r
        end
    end
end

#-----------------------------------------------------------------------------#
#=============================================================================#

function ρ(mother, Δ, n, m) end

#=============================================================================#


function on_circle!(input::VarWeierstrassForm, pole::RationalComplex, r::Float64, 
                    output::VarPolyForm; conjQ::Bool=false)
    #------------------------------------------------------
    Δ = pole = input.pole
    
end

#-----------------------------------------------------------------------------#

function on_circle!(
    in::VarQSpecialForm, 
    pole::RationalComplex, 
    r::Float64, 
    out::VarPolyForm; 
    conjQ::Bool=false)
end
