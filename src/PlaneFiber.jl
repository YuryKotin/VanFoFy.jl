struct PlaneLayer
    E ::Float64
    ν ::Float64
    r ::Float64
    # polynomial form of solution
    ϕ ::VarPolyForm
    χ ::VarPolyForm
end

function PlaneLayer(E, ν, r, r0, first_index, max_power)
    pow_range_ϕ = -max_power     : max_power + 2
    pow_range_χ = -max_power - 2 : max_power
    
    ar = OffsetVector{Int}(undef, 0:max_power+2)
    br = OffsetVector{Int}(undef, 0:max_power)
    ai = similar(ar)
    bi = similar(br)

    ind = first_index
    for i in axes(ar,1)
        ar[i] = ind; ind += 1
        ai[i] = ind; ind += 1
    end
    for i in axes(br,1)
        br[i] = ind; ind += 1
        bi[i] = ind; ind += 1
    end 
    
    var_range = first(ar) : last(bi)
    
    ϕ = VarPolyForm(var_range, pow_range_ϕ)
    χ = VarPolyForm(var_range, pow_range_χ)
    
    for i in axes(ar,1)
        ϕ[ar[i],i] = 1.0 + 0.0im
        ϕ[ai[i],i] = 0.0 + 1.0im
    end
    for i in axes(br,1)
        χ[br[i],i] = 1.0 + 0.0im
        χ[bi[i],i] = 0.0 + 1.0im
    end
    
    if r0 > 0
        # TODO Hollow layer
    end

    PlaneLayer(E, ν, r, ϕ, χ)    
end

function displacements_forces(layer::PlaneLayer, buffer::VarPolyFormBox)
    for i in 1 : 5
        axes(buffer[i]) == axes(layer.ϕ) || 
            error("Temp container has wrong shape")
    end

    E = layer.ϕ
    ν = layer.ν
    r = layer.r
    ϕ = layer.ϕ
    χ = layer.χ
    
    # f = ϕ(z)/. z-> rτ
    f = buffer[1]
    @views for p in powers(ϕ)
        factor = r^p
        @. f[:,p] = ϕ[:,p] * factor
    end
    
    # zbF = z \bar Φ /.{z -> rτ, \bar z -> r/τ}
    zbF = buffer[2]
    @views for p in powers(ϕ)
        factor = p * r^p
        @. zbF[:,2-p] = conj(ϕ[:,p]) * factor
    end
    
    # by = \bar χ(z)/. \bar z -> r/τ
    bc = buffer[3]
    @views for p in powers(χ)
        factor = r^p
        @. bc[:,-p] = χ[:,p] * factor
    end

    G = E / (2(1+ν))
    κ = 3 - 4ν

    displ = buffer[4]
    force = buffer[5]
    @views for p in powers(f)
        @. displ[:,p] = (κ*f[:,p] - zbF[:,p] - bc[:,p]) / (2G)
        @. force[:,p] =    f[:,p] + zbF[:,p] + bc[:,p]
    end
    
    return
end

function PlaneLayer(E, ν, r, prev_layer::PlaneLayer, buffer::VarPolyFormBox)
    displacements_forces(prev_layer, buffer)
    displ = buffer[4]
    force = buffer[5]
    
    G = E / (2(1+ν))
    κ = 3 - 4ν
    r1 = prev_layer.r
    
    # f = ϕ(z)/. z-> rτ
    f = buffer[1]
    @views for p in powers(ϕ)
        @. f[:,p] = (1 / (1+κ))*(2G*displ[:,p] + force[:,p])
    end
    
    # zbF = z \bar Φ /.{z -> rτ, \bar z -> r/τ}
    zbF = buffer[2]
    @views for p in powers(ϕ)
        @. zbF[:,2-p] = conj(f[:,p]) * p
    end
    
    # by = \bar χ(z)/. \bar z -> r/τ
    bc = buffer[3]
    @views for p in powers(ϕ)
        @. bc[:,p] = force[:,p] - f[:,p] - zbF[:,p]
    end
        
    ϕ = similar(prev_layer.ϕ)
    χ = similar(prev_layer.χ)
    @views for p in powers(ϕ)
        factor = r1^p
        @. ϕ[:,p] = f[:,p] / factor
    end
    @views for p in powers(χ)
        factor = r1^p
        @. χ[:,p] = bc[:,-p] / factor
    end

    PlaneLayer(E, ν, r, ϕ, χ)    
end

struct PlaneFiber
    layers ::Vector{PlaneLayer}
    r0     ::Float64
    buffer ::VarPolyFormBox
end

function PlaneFiber(data::FiberData, first_index, max_power)
    r0 = data.r0
    
    n_layers = size(data.layers, 1)
    layers = Vector{PlaneLayer}(undef, n_layers)
    
    ld = data.layers[1]
    layers[1] = PlaneLayer(
        ld.E, ld.ν, ld.r,
        r0, first_index, max_power)
    
    buffer = [similar(layers[1].ϕ) for _ in 1:5]
    
    if n_layers > 1
        for n in 2 : n_layers
            ld = data.layers[n]
            layers[n] = PlaneLayer(ld.E, ld.ν, ld.r, layers[n-1], buffer)
        end
    end

    displacements_forces(layers[end], buffer)

    PlaneFiber(layers, r0, buffer)
end

struct PlaneInclusion
    fiber  ::PlaneFiber
    center ::ComplexF64
end

function PlaneInclusion(data::InclusionData, first_index, center)
    fiber = PlaneFiber(data.fiber, first_index, data.max_power)
    PlaneInclusion(fiber, center)
end