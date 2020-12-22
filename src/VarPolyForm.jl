struct VarPolyForm
    data ::ComplexOffsetMatrix
end

#-----------------------------------------------------------------------------#

function VarPolyForm(var_range, pow_range)
    data = OffsetArray{ComplexF64, 2}(undef, pow_range, var_range)
    fill!(data, 0.0im)
    VarPolyForm(data)
end

#-----------------------------------------------------------------------------#

@forward VarPolyForm.data getindex, setindex!, fill!, axes

#-----------------------------------------------------------------------------#

Base.similar(form::VarPolyForm) = VarPolyForm(similar(form.data))

#-----------------------------------------------------------------------------#

powers(form::VarPolyForm) =    axes(form.data, 1)
variables(form::VarPolyForm) = axes(form.data, 2)

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
    axes(input) == axes(output) || 
        error("In and out containers have different shapes")
    
    pows = powers(input)
    vars = variables(input)
    rp = r^first(pows)
    for v in vars
        factor = rp
        for p in pows
            output[v,p] = input[v,p] * factor
            factor *= r
        end
    end
end

#-----------------------------------------------------------------------------#

function on_circle_conj!(input::VarPolyForm, r::Float64, output::VarPolyForm)
    variables(input) == variables(output) || 
        error("In and out containers have different variables")
    first(powers(input)) == -last(powers(output))  || 
        error("In and out containers have different powers")
    last(powers(input))  == -first(powers(output)) || 
        error("In and out containers have different powers")
    
    pows = powers(input)
    vars = variables(input)
    rp = r^first(pows)
    for v in vars
        factor = rp
        for p in pows
            output[v,-p] = conj(input[v,p]) * factor
            factor *= r
        end
    end
end

#-----------------------------------------------------------------------------#

function from_circle(input::VarPolyForm, r::Float64)
    output = similar(input)

    pows = powers(input)
    vars = variables(input)
    rp = r^first(pows)
    for v in vars
        factor = rp
        for p in pows
            output[v,p] = input[v,p] / factor
            factor *= r
        end
    end

    return output
end

#-----------------------------------------------------------------------------#

function from_circle_conj(input::VarPolyForm, r::Float64)
    conj_powers = -last( powers(input)): -first(powers(input))
    output = VarPolyForm(variables(input), conj_powers)

    pows = powers(input)
    vars = variables(input)
    rp = r^first(pows)
    for v in vars
        factor = rp
        for p in pows
            output[v,-p] = conj(input[v,p]) / factor
            factor *= r
        end
    end

    return output
end
