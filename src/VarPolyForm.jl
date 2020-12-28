struct PolynomForm
    coeffs ::ComplexOffsetMatrix
end

#-----------------------------------------------------------------------------#

function PolynomForm(var_range, pow_range)
    data = OffsetArray{ComplexF64, 2}(undef, pow_range, var_range)
    fill!(data, 0.0im)
    PolynomForm(data)
end

#-----------------------------------------------------------------------------#

@forward PolynomForm.coeffs getindex, setindex!, fill!, axes

#-----------------------------------------------------------------------------#

Base.similar(form::PolynomForm) = PolynomForm(similar(form.coeffs))

#-----------------------------------------------------------------------------#

powers(form::PolynomForm) =    axes(form.coeffs, 1)
variables(form::PolynomForm) = axes(form.coeffs, 2)

#=============================================================================#

struct VarPolyFormBox
    # f = ϕ(z)/. z-> rτ
    f ::PolynomForm 
    # zbF = z \bar Φ /.{z -> rτ, \bar z -> r/τ}
    zbF ::PolynomForm 
    # by = \bar ψ(z)/. \bar z -> r/τ
    by ::PolynomForm 
    # displ = (κ f - zbF - by)/(2G)
    displ ::PolynomForm
    # force = f + zbF + by
    force ::PolynomForm
end

#-----------------------------------------------------------------------------#

function VarPolyFormBox(form::PolynomForm)
    f = similar(form)
    zbF = similar(form)
    by = similar(form)
    displ = similar(form)
    force = similar(form)
    VarPolyFormBox(f, zbF, by, displ, force)
end

#=============================================================================#

conjQ(z, Q) = Q ? conj(z) : z

#=============================================================================#

function series_on_circle!(input::PolynomForm, r::Float64, output::PolynomForm, is_conj::Bool=false)
    variables(input) == variables(output) || 
        error("In and out containers have different variables")
    if is_conj
        powers(input) == conj_range(powers(output))  || 
            error("In and out containers have different powers")
    else
        powers(input) == powers(output)  || 
            error("In and out containers have different powers")
    end    
    
    ps = is_conj ? -1 : 1
    
    pows = powers(input)
    vars = variables(input)
    rp = r^first(pows)
    for v in vars
        factor = rp
        for p in pows
            val = input[v,p] * factor
            output[ps * p, v] = conjQ(val, is_conj)
            factor *= r
        end
    end
end

#-----------------------------------------------------------------------------#

function on_circle_conj!(input::PolynomForm, r::Float64, output::PolynomForm)
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

function from_circle(input::PolynomForm, r::Float64)
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

function from_circle_conj(input::PolynomForm, r::Float64)
    conj_powers = -last( powers(input)): -first(powers(input))
    output = PolynomForm(variables(input), conj_powers)

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
