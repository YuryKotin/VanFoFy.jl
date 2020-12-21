struct VarPolyForm
    data ::ComplexOffsetMatrix
end

#-----------------------------------------------------------------------------#

function VarPolyForm(var_range, pow_range)
    data = OffsetArray{ComplexF64, 2}(undef, var_range, pow_range)
    fill!(data, 0.0im)
    VarPolyForm(data)
end

#-----------------------------------------------------------------------------#

@forward VarPolyForm.data getindex, setindex!, fill!, axes

#-----------------------------------------------------------------------------#

Base.similar(form::VarPolyForm) = VarPolyForm(similar(form.data))

#-----------------------------------------------------------------------------#

variables(form::VarPolyForm) = axes(form.data, 1)
powers(form::VarPolyForm) =    axes(form.data, 2)

#=============================================================================#

struct VarPolyFormBox
    # f = ϕ(z)/. z-> rτ
    f ::VarPolyForm 
    # zbF = z \bar Φ /.{z -> rτ, \bar z -> r/τ}
    zbF ::VarPolyForm 
    # by = \bar ψ(z)/. \bar z -> r/τ
    by ::VarPolyForm 
    # displ = (κ f - zbF - by)/(2G)
    displ ::VarPolyForm
    # force = f + zbF + by
    force ::VarPolyForm
end

#-----------------------------------------------------------------------------#

function VarPolyFormBox(form::VarPolyForm)
    f = similar(form)
    zbF = similar(form)
    by = similar(form)
    displ = similar(form)
    force = similar(form)
    VarPolyFormBox(f, zbF, by, displ, force)
end

#=============================================================================#

function on_circle!(input::VarPolyForm, r::Float64, output::VarPolyForm)
    axes(input) == axes(output) || error("In and out containers have different shapes")
    @views for p in powers(input)
        factor = r^p
        @. output[:,p] = input[:,p] * factor
    end
end

#-----------------------------------------------------------------------------#

function on_circle_conj!(input::VarPolyForm, r::Float64, output::VarPolyForm)
    variables(input) == variables(output) || error("In and out containers have different variables")
    
    first_pow_in  = first(powers(input))
    last_pow_in   = last( powers(input))
    first_pow_out = first(powers(output))
    last_pow_out  = last( powers(output))
    first_pow_in == -last_pow_out  || error("In and out containers have different powers")
    last_pow_in  == -first_pow_out || error("In and out containers have different powers")
    
    @views for p in powers(input)
        factor = r^p
        @. output[:,-p] = conj(input[:,p]) * factor
    end
end

#-----------------------------------------------------------------------------#

function from_circle(input::VarPolyForm, r::Float64)
    output = similar(input)
    @views for p in powers(input)
        factor = r^p
        @. output[:,p] = input[:,p] / factor
    end
    return output
end

#-----------------------------------------------------------------------------#

function from_circle_conj(input::VarPolyForm, r::Float64)
    conj_powers = -last( powers(input)): -first(powers(input))
    output = VarPolyForm(variables(input), conj_powers)

    @views for p in powers(input)
        factor = r^p
        @. output[:,-p] = conj(input[:,p]) / factor
    end
    return output
end
